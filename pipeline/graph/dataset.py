from pathlib import Path

from torch.utils.data import Dataset
import torch

from pipeline.graph.g_config import Y_SCALE


class DriftDataset(Dataset):
    """
    Optional feature–masking: pass a list of names to `drop`.
    Allowed names  = ['rmsd_drift', 'energy_variance', 'gate_aperture_median',
                      'nbd_distance_median', 'cavity_volume', 'tmh_overlap',
                      'ligand_drift', 'all', 'ligand_only']
    """
    _feat_index = {                 # index inside u_raw
        "rmsd_drift":            0,
        "energy_variance":       1,
        "gate_aperture_median":  2,
        "nbd_distance_median":   3,
        "cavity_volume":         4,
        "tmh_overlap":           5,
        "ligand_drift":          6,
    }


    def __init__(self, graph_dir, *, keep_x=None,
                 scaler=None, mask_drift=True, drop=None):
        if drop is None:
            drop = []
        self.keep_x = (None if keep_x is None else
                       torch.as_tensor(keep_x, dtype=torch.long))
        self.paths = []
        self.meta = []
        for f in sorted(Path(graph_dir).rglob("frame_*_graph.pt")):
            try:
                g = torch.load(f, map_location='cpu')
                drift = float(g.u_raw[self._feat_index["ligand_drift"]])
                if drift > 10.0:
                    continue  # skip bad graphs
                self.paths.append(f)
                self.meta.append({
                    "protein_id": g.meta["protein_id"],
                    "pocket_id":  g.meta["pocket_id"],
                })
            except Exception as e:
                print(f"[WARN] Skipping {f} due to load error: {e}")
                continue

        self.scaler = scaler
        self.mask   = mask_drift
        self.drop   = set(drop or [])

    def __len__(self):
        return len(self.paths)

    def __getitem__(self, idx):
        g = torch.load(self.paths[idx], map_location='cpu')
        g.y_scaled = g.y * Y_SCALE

        u = g.u_raw.clone()

        # optional z-scaling (μ,σ) fitted on *train* only
        if self.scaler is not None:
            mu, std = self.scaler[g.meta["protein_id"]]
            u   = (u - mu) / std     # drift & envar only

        # --- on-the-fly feature masking ----------------------------------
        if self.keep_x is not None:
            g.x = g.x[:, self.keep_x]

        # attach to graph (no file write)
        g.u = u
        if len(self.keep_x) > 0:
            z_energy = u[1:2].view(1, 1).expand(g.x.size(0), 1).to(g.x.device)
            g.x = torch.cat([g.x, z_energy], dim=-1)

        return self._attach_contacts(g, self.paths[idx])

    """
        Mixin: if a sibling *.contact.pt* file exists, load it and
        attach `g.is_contact` (Bool[E]) aligned with edge_index.
        """
    def _attach_contacts(self, g, graph_path):
        contact_path = graph_path.with_suffix(".contact.pt")
        if contact_path.exists():
            g.is_contact = torch.load(contact_path, map_location='cpu')
        else:
            # Fallback – create an all-false tensor so code won't break
            g.is_contact = torch.zeros(g.edge_index.size(1), dtype=torch.bool)
        return g

    def set_phys_stats_by_protein(self, mu_sigma_dict, fallback=False):
        """
        Set per-protein normalization for u_raw vector.
        Modifies each graph in-place by adding a normalized `u` attribute.
        """
        self.scaler = {}

        for meta in self.meta:
            pid = meta["protein_id"]
            if pid in mu_sigma_dict:
                self.scaler[pid] = (
                    mu_sigma_dict[pid]["mu"],
                    mu_sigma_dict[pid]["std"].clamp_min(1e-6),
                )
            elif fallback:
                self.scaler[pid] = (mu_sigma_dict["_GLOBAL_"]["mu"],
                                    mu_sigma_dict["_GLOBAL_"]["std"])
            else:
                raise KeyError(f"[DriftDataset] No mu/std for protein {pid}")
