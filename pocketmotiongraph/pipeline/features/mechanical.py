"""Mechanical / physico‑chemical features (stubs)."""
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd
from Bio.PDB import PDBParser
import os
import tempfile
import subprocess
import mdtraj as md
KD_HYDROPHOBICITY = {
    "A": 1.8,  "C": 2.5,  "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5,  "K": -3.9, "L": 3.8,
    "M": 1.9,  "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2,  "W": -0.9, "Y": -1.3
}

def compute_hydrophobicity(sequence: str) -> list[float]:
    """
    Map sequence to Kyte-Doolittle hydrophobicity values.
    """
    return [KD_HYDROPHOBICITY.get(res, 0.0) for res in sequence]

def compute_bfactor(pdb_path: str) -> list[float]:
    """
    Extract B-factor (temperature factor) for each CA atom from a PDB file.
    Returns list[float] for all CA atoms in order.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    b_factors = []
    for atom in structure.get_atoms():
        if atom.name == "CA":
            b_factors.append(atom.bfactor)
    return b_factors

import numpy as np
from scipy.spatial.distance import cdist
from Bio.PDB import PDBParser
from scipy.spatial import ConvexHull

def calculate_pocket_volume(
        pdb_path: str,
        center: tuple,
        radius: float = 5.5,
        grid_spacing: float = 1.0,
        atom_radius: float = 1.65  # Van der Waals radius for carbon
) -> float:
    """
    Calculate pocket volume using core POVME algorithm with BioPython.

    Args:
        pdb_path: Path to PDB file
        center: (x, y, z) center of inclusion sphere
        radius: Radius of inclusion sphere in Å
        grid_spacing: Grid resolution in Å
        atom_radius: Minimum distance from protein atoms

    Returns:
        Pocket volume in Å³
    """
    # 1. Generate inclusion sphere points
    x = np.arange(center[0]-radius, center[0]+radius+grid_spacing, grid_spacing)
    y = np.arange(center[1]-radius, center[1]+radius+grid_spacing, grid_spacing)
    z = np.arange(center[2]-radius, center[2]+radius+grid_spacing, grid_spacing)

    # Create 3D grid
    xx, yy, zz = np.meshgrid(x, y, z)
    grid_points = np.vstack([xx.ravel(), yy.ravel(), zz.ravel()]).T

    # Filter to points within sphere
    dist_from_center = np.linalg.norm(grid_points - np.array(center), axis=1)
    inclusion_points = grid_points[dist_from_center <= radius]

    # 2. Load protein structure with BioPython
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_path)
    protein_coords = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    protein_coords.append(atom.get_coord())

    protein_coords = np.array(protein_coords)

    # 3. Remove points near protein atoms
    if len(protein_coords) > 0:
        # Calculate distances between grid points and protein atoms
        dist_matrix = cdist(inclusion_points, protein_coords)
        min_dists = np.min(dist_matrix, axis=1)

        # Filter points that are too close to protein atoms
        pocket_points = inclusion_points[min_dists > atom_radius]
    else:
        pocket_points = inclusion_points

    # 4. Calculate volume
    if len(pocket_points) > 4:
        hull = ConvexHull(pocket_points)
        return hull.volume
    elif len(pocket_points) > 0:
        # Approximate volume for small pockets
        return len(pocket_points) * (grid_spacing ** 3)
    else:
        return 0.0

def measure_drift_from_pdbs(pdb_initial,
                                pdb_final,
                                ligand_resnames=["UNL", "MOL", "LIG", "UNK", "VER", "MIL", "DIS", "TAC", "BEA"]):
        traj_initial = md.load(pdb_initial)
        traj_final = md.load(pdb_final)

        print(f"measuring drift {pdb_initial} ---- {pdb_final}")

        if traj_initial.n_atoms != traj_final.n_atoms:
            raise ValueError("Atom count mismatch between initial and final frame.")

        # ---- 2. ALIGN final → initial on protein Cα --------------------------
        prot_ca = traj_initial.topology.select("protein and name CA")
        if prot_ca.size == 0:
            raise ValueError("No protein Cα atoms found for alignment.")

        traj_final.superpose(traj_initial,
                             atom_indices=prot_ca,    # fit atoms    (mobile)
                             ref_atom_indices=prot_ca) # reference atoms (target)

        # ---- 3. pick ligand atoms -------------------------------------------
        ligand_atoms = None
        for resname in ligand_resnames:
            atoms = traj_initial.topology.select(f"resname {resname}")
            if len(atoms) > 0:
                ligand_atoms = atoms
                break

        if ligand_atoms is None:
            raise ValueError("No ligand atoms found.")

        # ---- 4. COMs and drift ----------------------------------------------
        com_initial = traj_initial.xyz[0][ligand_atoms].mean(axis=0)
        com_final = traj_final.xyz[0][ligand_atoms].mean(axis=0)
        drift_vector = com_final - com_initial
        drift_distance = np.linalg.norm(drift_vector)

        return drift_distance * 10, drift_vector, len(ligand_atoms)


CONTACTS = [
    "HydrogenBond",
    "Hydrophobic",
    "SaltBridge",
    "PiStacking",
    "PiCation",
    "Halogen",
    "WaterBridge",
]

def parse_plip_interactions(csv_path: Path, tag: str):
    """
    Parses the PLIP summary file for a given tag ('initial' or 'final').

    Returns:
    - `residue_to_counts`: dict[int, np.ndarray] of shape [len(CONTACTS)]
    - `global_counts`: dict[str, float], key is f"{Type}_{tag}"
    """
    df = pd.read_csv(csv_path)
    df_tag = df[df["Complex"].str.contains(tag)]

    residue_to_counts = defaultdict(lambda: np.zeros(len(CONTACTS)))
    global_counts = defaultdict(float)

    for _, row in df_tag.iterrows():
        match = re.match(r"(\d+)([A-Z]+)?", str(row["Residue"]))
        if match:
            resid = int(match.group(1))
            resname = match.group(2)  # e.g., "THR"
        else:
            continue  # or log warning

        typ = row["Type"]
        if typ in CONTACTS:
            idx = CONTACTS.index(typ)
            residue_to_counts[resid][idx] += 1
            global_counts[f"{typ}_{tag}"] += 1.0

    return dict(residue_to_counts), dict(global_counts)

def residue_fingerprint(csv_path: Path, tag: str) -> dict[int, np.ndarray]:
    """Returns {resid: fingerprint_vector}"""
    fp, _ = parse_plip_interactions(csv_path, tag)
    return fp

def plip_features_from_summary(csv_path: Path) -> dict[str, float]:
    """Returns global counts like {'Hbond_final': 5.0}"""
    _, counts_final = parse_plip_interactions(csv_path, "final")
    _, counts_init  = parse_plip_interactions(csv_path, "initial")
    return {**counts_init, **counts_final}

# Usage example
if __name__ == "__main__":
    # pdb_path = "/home/zhenli/git/valleyfevermutation/mutation_pipeline/md/AFR1/simulation_explicit/pocket14/replica_0/AFR1_prepared_Milbemycin_pocket14_complex_clean.pdb"
    # center = (29.7694, 12.4363, 6.8341)
    # # the result 286 vs 297 is close to prankweb
    # try:
    #     volume = calculate_pocket_volume(pdb_path, center, radius=5.0)
    #     print(f"Calculated pocket volume: {volume:.2f} Å³")
    # except Exception as e:
    #     print(f"Error calculating volume: {e}")

    # pdb_init = "/home/zhenli/git/valleyfevermutation/mutation_pipeline/md/AFR1/simulation_explicit/pocket1/replica_1/AFR1_prepared_Milbemycin_pocket1_complex_recombined_complex_explicit_initial_frame.pdb"
    # pdb_final = "/home/zhenli/git/valleyfevermutation/mutation_pipeline/md/AFR1/simulation_explicit/pocket1/replica_1/AFR1_prepared_Milbemycin_pocket1_complex_recombined_complex_explicit_final_frame.pdb"
    #
    # drift, _, _ = measure_drift_from_pdbs(pdb_init, pdb_final)
    # print(drift)


    replica_dir = Path("/git/valleyfevermutation/mutation_pipeline/md/AFR1/simulation_explicit/pocket1/replica_1")
    csv_path    = Path("/home/zhenli/git/valleyfevermutation/mutation_pipeline/md/AFR1/simulation_explicit/pocket1/replica_1/plip_results/all_plip_interactions_summary.csv")

    row = plip_features_from_summary(csv_path)
    row["protein_id"] = "AFR1"
    row["pocket_id"]  = "pocket1"
    row["stable"]     = 1   # <-- placeholder label
    for r in row:
        print(f"{r}")
    #
    # df = pd.DataFrame([row])
    # df.to_csv("plip_feature_table.csv", index=False)
    # print(df.head(1).T)
