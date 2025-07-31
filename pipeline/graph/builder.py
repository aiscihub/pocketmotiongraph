"""Convert feature dict to torch_geometric.Data."""
from pathlib import Path
# builder.py  – add near top
from pathlib import Path
import json

import pandas as pd
from torch_geometric.data import Data
from torch_cluster import radius_graph
import torch, numpy as np

from pipeline.features.mechanical import residue_fingerprint, CONTACTS, plip_features_from_summary
from pipeline.graph.domain_map_gff import cached_domain_map
from pipeline.graph.g_config import MAX_PROTEINS, _pid_lookup, PROTEIN_LIST, EXPECTED_U_RAW_DIM, \
    U_RAW_FEATURES, THRESHOLDS, NODE_FEATURE_DIM, NEIGHBOR_CUTOFF, NODE_FEATURE_DIM_STATIC, CONTACT_CUTOFF, FLAG_DIM, \
    DOMAIN_FLAGS

# ------------------------------------------------------------
# Normalise raw contact counts -> 0‒1  (log1p  then /log1p(5))
# ------------------------------------------------------------
_LOG_CONST = np.log1p(5.0)        # assume max ≈ 5 contacts / residue / type

def _norm_counts(vec: np.ndarray) -> np.ndarray:
    """vec shape = [7] raw integers 0,1,2…  → float 0–1"""
    return np.log1p(vec) / _LOG_CONST

# --- NEW quick contact counter (ligand ↔ protein heavy atoms ≤ 4 Å) --------
def _contact_count(xyz: torch.Tensor, lig_atom_mask: torch.Tensor) -> float:
    """
    Return the number of heavy-atom pairs (ligand i, protein j)
    whose Euclidean distance ≤ 4.0 Å in *this frame*.
        xyz            – (N_atoms, 3) float tensor, Å
        lig_atom_mask  – bool mask for ligand atoms (same length as xyz)
    """
    lig_xyz  = xyz[lig_atom_mask]
    prot_xyz = xyz[~lig_atom_mask]
    if lig_xyz.numel() == 0 or prot_xyz.numel() == 0:
        return 0.0
    d2 = torch.cdist(lig_xyz, prot_xyz, p=2)        # (N_lig, N_prot)
    return float((d2 <= 4.0).sum())

# --------------------------------------------------------------------------- #
#   BUILDER with PLIP fingerprints merged into existing node / global feats   #
# --------------------------------------------------------------------------- #

def _augment_with_plip(features: dict, node_feats: torch.Tensor) -> torch.Tensor:
    """Concatenate PLIP fingerprints [init, final, delta] to each node feature.
    The residue ID order is taken from the provided `frame_path` (a stripped PDB).
    """
    ref_pdb: Path = features.get("plip_ref_pdb") or features.get("frame_path")
    if ref_pdb is None or not Path(ref_pdb).exists():
        return node_feats
    csv_path = features.get("plip_csv_path") or (ref_pdb.parent / "plip_results" / "all_plip_interactions_summary.csv")
    if not csv_path.exists():
        return node_feats

    fp_init  = residue_fingerprint(csv_path, "initial")
    fp_final = residue_fingerprint(csv_path, "final")
    vec_len  = len(CONTACTS)

    # Use Bio.PDB to extract residue ids (numbers) in CA order from *ref_pdb*
    from Bio.PDB import PDBParser
    res_order = []
    keys = []
    for res in PDBParser(QUIET=True).get_structure("ref", str(ref_pdb)).get_residues():
        if any(atom.name == "CA" for atom in res):
            res_order.append(res.id[1])  # integer resid
            keys.append((res.get_parent().id, res.id[1], res.id[2]))

    if len(res_order) != node_feats.shape[0]:
        # Cannot safely align – skip augmentation
        return node_feats

    plip_vecs = []
    for resid in res_order:
        v_i = _norm_counts(fp_init.get(resid, np.zeros(vec_len)))
        v_f = _norm_counts(fp_final.get(resid, np.zeros(vec_len)))
        plip_vecs.append(np.concatenate([v_i, v_f, v_f - v_i]))


    plip_tensor = torch.tensor(plip_vecs, dtype=torch.float32)
    if plip_tensor.shape[0] != node_feats.shape[0]:
        # atom mapping mismatch – skip augmentation
        return node_feats

    # print(f"[debug] node_feats.shape: {node_feats.shape}")
    # print(f"[debug] plip_tensor.shape: {plip_tensor.shape}")
    #hits = (plip_tensor.sum(dim=1) > 0).sum().item()
    #print(f"[PLIP] matched {hits}/{len(keys)} residues with ≥1 contact")
    features["_cached_node_keys"] = keys         # pass to build_pyg_graph
    return torch.cat([node_feats, plip_tensor], dim=1)

def _assemble_global_u(features: dict) -> torch.Tensor:
    """Concatenate original u_raw scalars with PLIP global counts."""
    # -- new first two slots -------------------------------------------------
    rmsd_drift_this_frame = float(features.get("pocket_rmsd_drift",
                                               features.get("rmsd_drift", 0.0)))
    print(f"rmsd_drift_this_frame --- {rmsd_drift_this_frame}")
    base_vals = [
        rmsd_drift_this_frame,
        float(features.get("energy_variance", 0.0)),
        float(features.get("gate_aperture_median", 0.0)),
        float(features.get("nbd_distance_median", 0.0)),
        float(features.get("cavity_volume", 0.0)),
        float(features.get("tmh_overlap", 0.0)),
        float(features.get("ligand_drift", 0.0)),
    ]

    # ------------------ NEW per-frame contact_count -----------------------
    try:
        # Optional PLIP CSV (initial/final snapshot)
        csv_path = features.get("plip_csv_path")
        g_counts = {}
        contact_cnt = 0
        if csv_path and csv_path.exists():
            g_counts = plip_features_from_summary(csv_path)

            # Contact count — try direct coords, else fallback to ramp from PLIP init/final
            contact_types = sorted(set(k.split("_")[0] for k in g_counts if "_" in k))
            init = sum(g_counts.get(f"{t}_initial", 0.0) for t in contact_types)
            final = sum(g_counts.get(f"{t}_final", 0.0) for t in contact_types)
            i = int(features.get("frame_idx", 0))
            N = int(features.get("frames_per_replica", 1))
            t = i / max(N - 1, 1)
            #print(f"{init} {final} -- {contact_cnt}")
            contact_cnt = init + t * (final - init)
        else:
            print(f"didn't find file {csv_path}")
            contact_cnt = 0

        if contact_cnt > 0 :
            print(f"contact_cnt = {contact_cnt}")

        base_vals.append(contact_cnt)  # slot 7 → contact count
        base_u = torch.tensor(base_vals, dtype=torch.float32)

        plip_u = torch.tensor([g_counts.get(f"{t}_{phase}", 0.0)
                                   for phase in ("initial", "final")
                                   for t in CONTACTS], dtype=torch.float32)
        # FIXED: Domain flag calculation (multi-hot encoding)
        one_hot = np.zeros(len(DOMAIN_FLAGS), dtype=np.float32)
        pocket_residues = features.get("pocket_residue_ids") or []

        # Get residue keys if pocket_residues not provided
        if not pocket_residues and "_cached_node_keys" in features:
            pocket_residues = [resid for (_, resid, _) in features["_cached_node_keys"]]

        # Only process if we have pocket residues
        if pocket_residues:
            protein_id = features.get("protein_id")
            if protein_id:
                # Get domain map for this protein
                domain_map = cached_domain_map(protein_id)

                for resid in pocket_residues:
                    dom = domain_map.get(resid)
                    if dom in DOMAIN_FLAGS:
                        idx = DOMAIN_FLAGS.index(dom)
                        one_hot[idx] = 1.0  # Set flag if domain exists in pocket
                       # print(f"{resid} ---- {dom} = 1")
        return torch.cat([base_u, plip_u, torch.tensor(one_hot, dtype=torch.float32)])
    except Exception as e:
        raise e
# ------------------------------------------------------------------ #
#  Protein-ID ➜ integer (0 … MAX_PROTEINS-1) ➜ 1-element LongTensor  #
#  This feeds the nn.Embedding inside PocketStabilityEGNN.           #
# ------------------------------------------------------------------ #
# module-level singleton {pdb_id: idx}

def _pid_index(pdb_id: str) -> int:
    """Return a stable, <=MAX_PROTEINS-1 index for each unique protein id."""
    if pdb_id not in _pid_lookup:
        new_idx = len(_pid_lookup)
        if new_idx >= MAX_PROTEINS:
            raise ValueError(
                f"MAX_PROTEINS={MAX_PROTEINS} exhausted – add capacity before "
                f"adding more proteins (offending id = {pdb_id})"
            )
        _pid_lookup[pdb_id] = new_idx
    return _pid_lookup[pdb_id]


def build_pyg_graph(features: dict) -> Data:
    """Assemble a PyG Data object with PLIP and existing features."""

    # --- node & edge tensors -------------------------------------------------
    x_esm = features["x"]
    if not torch.is_tensor(x_esm):
        x_esm = torch.tensor(x_esm, dtype=torch.float32)

    # augment with PLIP residue fingerprints
    x = _augment_with_plip(features, x_esm)

    pos = torch.tensor(np.asarray(features["ca_coords"][0]), dtype=torch.float32)
    edge_index = radius_graph(pos, r=NEIGHBOR_CUTOFF, loop=False)
    row, col   = edge_index
    edge_attr  = torch.norm(pos[row] - pos[col], dim=1).unsqueeze(1)

    # --- global feature vector ---------------------------------------------
    u_raw = _assemble_global_u(features)

    # label (regression – drift per frame)
    y_reg = torch.tensor([float(features["y_data"][0])]).view(1, 1)

    data = Data(
        x=x,
        edge_index=edge_index,
        edge_attr=edge_attr.float(),
        pos=pos,
        u_raw=u_raw,
        y=y_reg,
    )

    data.meta = dict(
        protein_id = features["protein_id"],
        pocket_id  = features.get("pocket_id",  "UNK"),
        replica    = features.get("replica_id", "UNK"),
        frame_idx  = int(features.get("frame_idx", -1)),
    )

    if "_cached_node_keys" in features:
       data.node_keys = features["_cached_node_keys"]      # [(chain,resid,icode)]
    else:
        if isinstance(x, torch.Tensor) and x.shape[0] > 0:
            # Re-run residue key extraction for reference
            ref_pdb: Path = features.get("plip_ref_pdb") or features.get("frame_path")
            if ref_pdb and Path(ref_pdb).exists():
                from Bio.PDB import PDBParser
                keys = []
                for res in PDBParser(QUIET=True).get_structure("ref", str(ref_pdb)).get_residues():
                    if any(atom.name == "CA" for atom in res):
                        keys.append((res.get_parent().id, res.id[1], res.id[2]))  # (chain, resid, icode)
                if len(keys) == x.shape[0]:
                    data.node_keys = keys
    # Attach numeric residue ID for attention visualization
    if hasattr(data, "node_keys"):
        # node_keys = [(chain, resid, icode)]  →  use resid as integer ID
        data.residue_id = torch.tensor(
            [resid for (_, resid, _) in data.node_keys],
            dtype=torch.long
        )

    # Attach a string identifier for saving .csv (e.g., 'CIMG_00533')
    data.pdb_id = features.get("protein_id", "UNK")
    data.pid = torch.tensor([ _pid_index(data.pdb_id) ], dtype=torch.long)

    validate_features(data)

    return data

def load_contact_residues(csv_path: Path) -> set[int]:
    """
    Returns the set of residue IDs that PLIP flagged as contacting the ligand
    in ANY frame of this replica.
    """
    if not csv_path.exists():
        return set()

    df = pd.read_csv(csv_path, usecols=["Residue"])
    resid_set: set[int] = set()
    for res in df["Residue"].dropna().astype(str):
        digits = "".join(filter(str.isdigit, res))
        if digits:
            resid_set.add(int(digits))
    return resid_set

def write_contact_sidecar(data, out_graph_path, plip_csv_path, cutoff=CONTACT_CUTOFF):
    """
    Builds edge-level is_contact tensor based on PLIP summary residue hits.
    """
    sidecar = Path(out_graph_path).with_suffix(".contact.pt")
    if sidecar.exists():
        return

    plip_residues = load_contact_residues(plip_csv_path)
    if not plip_residues:
        return

    row, col = data.edge_index
    dists = data.edge_attr[:, 0]
    mask = (dists <= cutoff)

    # residue_id is already attached during graph building
    res_ids = data.residue_id.tolist()

    is_contact = torch.zeros(row.size(0), dtype=torch.bool)
    for e in mask.nonzero(as_tuple=True)[0]:
        i_node = row[e].item()
        j_node = col[e].item()
        if res_ids[i_node] in plip_residues or res_ids[j_node] in plip_residues:
            is_contact[e] = True

    n_positive = is_contact.sum().item()
    n_total    = is_contact.numel()

    if n_positive == 0:
        print(f"[warn] No contacts found for {sidecar.stem} — check PLIP coverage or residue IDs.")
    elif n_positive > 0.5 * n_total:
        print(f"[warn] >50% of edges marked contact? ({n_positive}/{n_total}) → possible misalignment")
    else:
        print(f"[ok] Contact labels: {n_positive} / {n_total} edges = {100.0 * n_positive / n_total:.2f}%")

    torch.save(is_contact, sidecar)
    print(f"[contact.pt] saved: {sidecar}")

# -----------------------------------------------------------

def validate_features(data):
    print(f"Validating features for graph")
    u = data.u_raw
    #print(f" data.u_raw.shape[0] = {data.u_raw.shape[0]}")
    assert data.u_raw.shape[0] == EXPECTED_U_RAW_DIM

    if data.pdb_id not in PROTEIN_LIST:
        raise ValueError(f"[Graph Build] Unrecognized protein ID: '{data.pdb_id}'. Must be one of: {sorted(PROTEIN_LIST)}")

    if data.x.shape[1] != NODE_FEATURE_DIM_STATIC:
        raise ValueError(f"[Graph Build] Expected node feature dim = {NODE_FEATURE_DIM}, got {data.x.shape[1]}")

    for name, idx in U_RAW_FEATURES.items():
        try:
            val = float(u[idx])
        except Exception:
            raise RuntimeError(f"{name} missing or not float")

        if not torch.isfinite(torch.tensor(val)):
            raise RuntimeError(f"{name} is NaN/Inf: {val}")
        else:
            lo, hi = THRESHOLDS.get(name, (-float('inf'), float('inf')))
            if not (lo <= val <= hi):
                raise RuntimeError(f"{name} out of expected range: {val:.2f} not in [{lo}, {hi}]")
