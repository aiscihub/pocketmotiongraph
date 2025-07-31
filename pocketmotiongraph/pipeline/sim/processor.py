"""Run MD-trajectory analysis (stub)."""
from __future__ import annotations
import logging
import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.coordinates.PDB import PDBWriter
import os
from pathlib import Path
import tempfile
from pipeline.config import PROJECT_DIR
from pipeline.features.mechanical import calculate_pocket_volume, measure_drift_from_pdbs
from pipeline.io.prankweb import (
    load_prankweb_pocket_residues,
    load_pocket_center,
    load_pocket_score,
)
from pipeline.io.gff import parse_deeptmhmm_gff
from pipeline.sim.domain_cfg import DOMAIN_RANGES


@dataclass
class SimulationProcessor:
    replica_dir: Path
    protein_id: str
    pocket_id: str
    ligand_name: str = "Milbemycin"
    stability_threshold: float = 0.982
    logger: logging.Logger = None
    pocket_residues: list[int] = field(default_factory=list, init=False)
    universe: mda.Universe = field(default=None, init=False)
    target_n : int = 20
    tm_ranges: list = None

    def __post_init__(self):
        """Initializes logger and loads the Universe after the dataclass is created."""
        if self.logger is None:
            self.logger = logging.getLogger(__name__)
        try:
            top_path = (
                self.replica_dir
                / f"{self.protein_id}_prepared_{self.ligand_name}_{self.pocket_id}_complex_recombined_complex_explicit_initial_frame.pdb"
            )
            traj_path = (
                self.replica_dir
                / f"{self.protein_id}_prepared_{self.ligand_name}_{self.pocket_id}_complex_recombined_complex_explicit_trajectory.dcd"
            )
            self.universe = mda.Universe(str(top_path), str(traj_path))
            gff_path = (
                    PROJECT_DIR
                    / "mutation_pipeline"
                    / "dataset"
                    / "proteins"
                    / "gff"
                    / f"{self.protein_id}.gff3"
            )

            if gff_path.exists():
                tm_dict = parse_deeptmhmm_gff(gff_path)
                tm_ranges = tm_dict.get(self.protein_id, [])
                self.tm_ranges = tm_ranges
        except Exception as e:
            self.logger.error(f"Failed to load Universe for {self.replica_dir}: {e}")
            self.universe = None

    def run(self) -> dict[str, Any] | None:
        """Load MD trajectory, calculate all metrics, and return them in a dictionary."""
        if self.universe is None or len(self.universe.trajectory) < 2:
            self.logger.error(
                f"Universe not loaded or trajectory too short for {self.replica_dir}. Skipping."
            )
            return None

        prank_csv = (
            PROJECT_DIR
            / "mutation_pipeline"
            / "outputs"
            / "prank"
            / f"{self.protein_id}_relaxed.pdb_predictions.csv"
        )
        self.pocket_residues = load_prankweb_pocket_residues(
            prank_csv, pocket_name=self.pocket_id
        )
        if not self.pocket_residues:
            self.logger.error(
                f"Could not load pocket residues for {self.pocket_id}. Skipping."
            )
            return None

        domain_cfg = DOMAIN_RANGES.get(self.protein_id, {})

        # --- 1. Calculate Per-Frame Metrics (as lists) ---
        pocket_rmsd_per_frame = self._compute_pocket_rmsd_per_frame()
        if pocket_rmsd_per_frame is None:
            return None

        ca_coords = [
            self.universe.select_atoms("protein and name CA").positions.copy()
            for ts in self.universe.trajectory
        ]
        nbd_dists, gate_apertures = self._compute_dynamic_metrics(domain_cfg)

        # --- 2. Calculate Global Metrics (as single values) ---
        ligand_drift = self._compute_ligand_drift_pdb()

        global_rmsd_drift = self._compute_rmsd_drift()
        global_energy_variance = self._compute_energy_variance()

        pocket_center = load_pocket_center(prank_csv, pocket_name=self.pocket_id)
        pdb_for_volume = (
            self.replica_dir
            / f"{self.protein_id}_prepared_{self.ligand_name}_{self.pocket_id}_complex_recombined_complex.pdb"
        )
        cavity_volume_val = calculate_pocket_volume(str(pdb_for_volume), pocket_center)

         # --- NEW --- TM-helix overlap (static, per-pocket)
        tmh_overlap = self._compute_tmh_overlap(tm_ranges=self.tm_ranges)  # 0-1

# --- 3. Assemble Final Features Dictionary ---
        features = {
            # Per-frame data
            "pocket_rmsd_per_frame": pocket_rmsd_per_frame,
            "ca_coords": ca_coords,
            # Identity
            "protein_id": self.protein_id,
            "pocket_id": self.pocket_id,

            # Global data (single values)
            "tmh_overlap":         tmh_overlap,
            "rmsd_drift": global_rmsd_drift,
            "energy_variance": global_energy_variance,
            "cavity_volume": cavity_volume_val,
            "gate_aperture_median": np.median(gate_apertures)  if gate_apertures else 0.0,
            "nbd_distance_median": np.median(nbd_dists) if nbd_dists else 0.0,
            "ligand_drift":        ligand_drift,
        }

        self._dump_frames()
        return features

    def _compute_dynamic_metrics(
        self, domain_cfg: dict
    ) -> tuple[list[float], list[float]]:
        nbd_dists, gate_apertures = [], []
        if self.tm_ranges is None:
            gff_path = (
                PROJECT_DIR
                / "mutation_pipeline"
                / "dataset"
                / "proteins"
                / "gff"
                / f"{self.protein_id}.gff3"
            )
            tm_ranges = []
            if gff_path.exists():
                tm_dict = parse_deeptmhmm_gff(gff_path)
                tm_ranges = tm_dict.get(self.protein_id, [])
                self.tm_ranges = tm_ranges
        else:
            tm_ranges = self.tm_ranges
        for ts in self.universe.trajectory:
            try:
                nbd1_sel = self.universe.select_atoms(
                    f"resid {domain_cfg['NBD1'][0]}-{domain_cfg['NBD1'][1]} and name CA"
                )
                nbd2_sel = self.universe.select_atoms(
                    f"resid {domain_cfg['NBD2'][0]}-{domain_cfg['NBD2'][1]} and name CA"
                )
                if nbd1_sel.n_atoms > 0 and nbd2_sel.n_atoms > 0:
                    nbd_dists.append(
                        np.linalg.norm(
                            nbd1_sel.center_of_mass() - nbd2_sel.center_of_mass()
                        )
                    )
                else:
                    nbd_dists.append(0.0)
            except (KeyError, IndexError):
                nbd_dists.append(0.0)
            try:
                tm_idx1, tm_idx2 = domain_cfg.get("gate_tm_pair", (None, None))
                if (
                    tm_idx1 is not None
                    and tm_idx1 < len(tm_ranges)
                    and tm_idx2 < len(tm_ranges)
                ):
                    gate1_sel = self.universe.select_atoms(
                        f"resid {tm_ranges[tm_idx1][0]}-{tm_ranges[tm_idx1][1]} and name CA"
                    )
                    gate2_sel = self.universe.select_atoms(
                        f"resid {tm_ranges[tm_idx2][0]}-{tm_ranges[tm_idx2][1]} and name CA"
                    )
                    if gate1_sel.n_atoms > 0 and gate2_sel.n_atoms > 0:
                        gate_apertures.append(
                            np.linalg.norm(
                                gate1_sel.center_of_mass() - gate2_sel.center_of_mass()
                            )
                        )
                    else:
                        gate_apertures.append(0.0)
                else:
                    gate_apertures.append(0.0)
            except (KeyError, IndexError):
                gate_apertures.append(0.0)
        return nbd_dists, gate_apertures

    def _compute_pocket_rmsd_per_frame(self) -> np.ndarray | None:
        try:
            pocket_sel_str = "protein and name CA and resid " + " ".join(
                map(str, self.pocket_residues)
            )
            rmsd_calc = RMSD(
                self.universe, self.universe, select=pocket_sel_str, ref_frame=0
            )
            rmsd_calc.run()
            return rmsd_calc.results.rmsd[:, 2]
        except Exception as e:
            self.logger.error(f"RMSD calculation failed: {e}")
            return None

    def _compute_rmsd_drift(self) -> float | None:
        try:
            ref = self.universe.select_atoms("protein or resname MIL")
            if ref.n_atoms == 0:
                return 0.0
            rmsd_calc = RMSD(self.universe, ref, select="protein or resname MIL")
            rmsd_calc.run()
            rmsd_values = rmsd_calc.results.rmsd[:, 2]
            drift = np.abs(rmsd_values - rmsd_values[0]).mean()
            return float(drift)
        except Exception:
            return 0.0

    def _compute_energy_variance(self) -> float | None:
        try:
            log_path = (
                self.replica_dir
                / f"{self.protein_id}_prepared_{self.ligand_name}_{self.pocket_id}_complex_recombined_complex_interaction_energies.txt"
            )
            if not log_path.exists():
                return 0.0
            energies = np.loadtxt(log_path, skiprows=1, usecols=1, ndmin=1)
            #return np.var(energies) if energies.size > 0 else 0.0
            #Changet to standard
            #quick answer
            #Switch to option B (standard deviation) and keep your current per-protein z-scaling.
            #Physical meaning. σ answers “how much does the interaction energy fluctuate, on average?”—easy to discuss in a paper.
            #Dynamic range. σ (hundreds of kJ) already reduces the spread by ≈ √10³ without destroying linear relationships.
            #Stability. After z-scaling every protein still contributes numbers of order 1, but outliers in the raw trace can’t pull the variance estimate into 10⁵–10⁶ territory.
            return float(np.std(energies, ddof=0)) if energies.size > 0 else 0.0
        except Exception:
            return 0.0


    def _compute_ligand_drift_pdb(self) -> float:
        pdb_initial = os.path.join(self.replica_dir, f"{self.protein_id}_prepared_{self.ligand_name}_{self.pocket_id}_complex_recombined_complex_explicit_stripped_initial_frame.pdb")
        pdb_final = os.path.join(self.replica_dir, f"{self.protein_id}_prepared_{self.ligand_name}_{self.pocket_id}_complex_recombined_complex_explicit_stripped_final_frame.pdb")
        if not Path(pdb_final).is_file():
            self.logger.warning(f"[ligand_drift] Missing final frame: {pdb_final}")
            return None
        try:
            drift, _, _ = measure_drift_from_pdbs(pdb_initial, pdb_final)
            print(f"[info] Ligand Drift Distance  : {drift:.3f} ")
            return drift
        except Exception as e:
            self.logger.warning(f"[ligand_drift] Error measuring drift: {e}")
            return None

    def _compute_tmh_overlap(self, tm_ranges: list[tuple[int, int]]) -> float:
        """
        Fraction of pocket residues that fall inside any DeepTMHMM helix.
        Returns 0.0 when either list is empty.
        """
        if not tm_ranges or not self.pocket_residues:
            return 0.0

        tm_set = {resid for start, end in tm_ranges for resid in range(start, end + 1)}
        pocket_set = set(self.pocket_residues)

        overlap = len(pocket_set & tm_set) / len(pocket_set)
        return float(overlap)

    def _dump_frames(self):
        frames_dir = self.replica_dir / "frames"
        frames_dir.mkdir(exist_ok=True)
        with PDBWriter(str(frames_dir / "dummy.pdb"), multiframe=False) as W:
            pass  # create writer once to get topology; we'll reopen below

        stride = max(1, len(self.universe.trajectory) // self.target_n)
        for i, ts in enumerate(self.universe.trajectory[::stride]):
            out = frames_dir / f"frame_{i:04d}.pdb"
            if len(list(frames_dir.glob("frame_*.pdb"))) >= self.target_n:
                self.logger.info(f"[SKIP] Already found {self.target_n} frames in {frames_dir}")
                return
            with PDBWriter(str(out), multiframe=False) as W:
                W.write(self.universe.atoms)
            self.logger.info(f"saving {out.stem}")

    # For testing only
    def force_load(self, topology: str, trajectory: str):
        """Forcefully load a custom MD trajectory and topology for testing."""
        try:
            self.universe = mda.Universe(topology, trajectory)
            self.logger.info(f"[TEST] Successfully force-loaded trajectory from {trajectory}")
        except Exception as e:
            self.logger.error(f"[TEST] Failed to force load test trajectory: {e}")
            self.universe = None