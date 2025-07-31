"""ESM-2 embedding and feature extraction helpers."""
import numpy as np
import torch
from pathlib import Path

import esm
from Bio.PDB import PDBParser, is_aa
from MDAnalysis import Universe

from pipeline.features.mechanical import compute_hydrophobicity, compute_bfactor, residue_fingerprint

restype_3to1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z', 'UNK': 'X'
}

# --- Helper Functions for Feature Generation ---

def compute_sasa(seq: str) -> list[float]:
    """Placeholder: estimate SASA based on hydrophobicity."""
    hydrophobic_residues = {'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W'}
    return [0.2 if aa in hydrophobic_residues else 0.8 for aa in seq]

def compute_secondary_structure(seq: str) -> torch.Tensor:
    """Placeholder: one-hot encode secondary structure based on residue type."""
    helix_formers = {'A', 'E', 'L'}
    sheet_formers = {'V', 'I', 'Y'}
    result = []
    for aa in seq:
        if aa in helix_formers: result.append([1, 0, 0])
        elif aa in sheet_formers: result.append([0, 1, 0])
        else: result.append([0, 0, 1])
    return torch.tensor(result, dtype=torch.float32)

def get_sidechain_volume(seq: str) -> list[float]:
    """Get approximate side-chain volume for each residue."""
    volume_table = {'A': 67, 'R': 148, 'N': 96, 'D': 91, 'C': 86, 'Q': 114, 'E': 109, 'G': 48, 'H': 118, 'I': 124, 'L': 124, 'K': 135, 'M': 124, 'F': 135, 'P': 90, 'S': 73, 'T': 93, 'W': 163, 'Y': 141, 'V': 105}
    return [volume_table.get(aa, 100) for aa in seq]

# --- Core Feature and Embedding Functions ---

def get_sequence_and_resids_from_pdb(pdb_path: Path) -> tuple[str, list[int]]:
    """Extracts sequence and matching residue IDs from a PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))
    sequence = ""
    resid_list = []
    for residue in structure.get_residues():
        if is_aa(residue):
            sequence += restype_3to1.get(residue.get_resname(), 'X')
            resid_list.append(residue.get_id()[1])
    return sequence, resid_list

def get_ca_resids_from_md(pdb_path: Path) -> list[int]:
    """Extracts residue numbers for all C-alpha atoms."""
    u = Universe(str(pdb_path))
    return [atom.resid for atom in u.select_atoms("protein and name CA")]

def build_esm_resids_from_pdb(pdb_path: Path, esm_model, batch_converter) -> torch.Tensor | None:
    """Align ESM embeddings to MD structure residues."""
    seq, esm_resids = get_sequence_and_resids_from_pdb(pdb_path)
    ca_resids = get_ca_resids_from_md(pdb_path)
    if not seq: return None

    data = [(pdb_path.stem, seq)]
    _, _, tokens = batch_converter(data)
    tokens = tokens.to(next(esm_model.parameters()).device)

    with torch.no_grad():
        results = esm_model(tokens, repr_layers=[33], return_contacts=False)
    full_emb = results["representations"][33][0, 1 : len(seq) + 1]

    esm_resid_to_index = {resid: i for i, resid in enumerate(esm_resids)}
    matched_indices = [esm_resid_to_index.get(r, -1) for r in ca_resids]

    # Handle cases where a CA resid might not be in the ESM sequence
    final_indices = []
    for idx in matched_indices:
        if idx != -1:
            final_indices.append(idx)
        else:
            # Fallback: if a resid is missing, just append the last valid index
            # This is a simple strategy to maintain tensor shape.
            if final_indices:
                final_indices.append(final_indices[-1])
            else: # If the very first resid is missing
                final_indices.append(0)

    return full_emb[final_indices]

def normalize_feature(feat: torch.Tensor) -> torch.Tensor:
    feat = feat.float()
    return (feat - feat.mean(0)) / feat.std(0).clamp_min(1e-6)

def extract_features(pdb_path: Path, esm_model, batch_converter, graph_path: Path, override = False) -> torch.Tensor:
    """
    Consolidated function to build the complete node feature vector (x).
    """
    # Save feature file next to graph
    x_path = graph_path.with_name(graph_path.stem.replace("_graph", "_x") + ".pt")

    # Always delete stale feature to force regeneration

    # 1. ESM embedding
    esm_features = build_esm_resids_from_pdb(pdb_path, esm_model, batch_converter)
    if esm_features is None:
        raise ValueError(f"Failed to extract ESM embeddings from {pdb_path}")

    seq, _ = get_sequence_and_resids_from_pdb(pdb_path)
    device = esm_features.device

    try:
        # 2. Other physicochemical features (normalized)
        hydro_feat  = normalize_feature(torch.tensor(compute_hydrophobicity(seq), device=device).unsqueeze(-1))
        bfac_feat   = normalize_feature(torch.tensor(compute_bfactor(str(pdb_path)), device=device).unsqueeze(-1))
        sasa_feat   = normalize_feature(torch.tensor(compute_sasa(seq), device=device).unsqueeze(-1))
        ss_feat     = normalize_feature(compute_secondary_structure(seq).to(device))  # already a tensor
        volume_feat = normalize_feature(torch.tensor(get_sidechain_volume(seq), device=device).unsqueeze(-1))
        # Concatenate all features
        # Normalize ESM embeddings across residues
        esm_features = (esm_features - esm_features.mean(0)) / esm_features.std(0).clamp_min(1e-6)

        # print("Feature stats:",
        #       "hydro:", hydro_feat.mean().item(), hydro_feat.std().item(),
        #       "bfacor:", bfac_feat.mean().item(), bfac_feat.std().item())

        x = torch.cat([esm_features, hydro_feat, bfac_feat, sasa_feat, ss_feat, volume_feat], dim=-1)
        x = x.cpu()
    except Exception as e:
        # This error can happen if a feature list doesn't match the number of atoms
        raise RuntimeError(f"Feature dimension mismatch in {pdb_path}: {e}")
    if x_path.exists():
        if override:
            x_path.unlink()
        else:
            return x
    torch.save(x, x_path)
    return x
