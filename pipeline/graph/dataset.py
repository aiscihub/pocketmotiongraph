# dataset.py
from pathlib import Path
from typing import Optional
import warnings
import torch
from torch.utils.data import Dataset
from torch_geometric.nn import radius_graph

from pipeline.graph.g_config import Y_SCALE


def _ensure_angstroms(pos: torch.Tensor, sample: int = 512):
    """Detect nm vs Å via NN-distance; return (pos_in_Å, flag)."""
    assert pos.ndim == 2 and pos.size(1) == 3, f"pos must be [N,3], got {tuple(pos.shape)}"
    N = pos.size(0)
    k = min(N, sample)
    idx = torch.arange(k, device=pos.device)  # deterministic for reproducibility
    with torch.no_grad():
        d = torch.cdist(pos.index_select(0, idx), pos.index_select(0, idx), p=2)
        d[d == 0] = float("inf")
        nn_med = float(d.min(dim=1).values.median())
    if 0.25 <= nn_med <= 0.70:     # ~0.38 nm CA–CA
        return pos * 10.0, "nm->A"
    if 2.5 <= nn_med <= 5.5:       # ~3.8 Å CA–CA
        return pos, "A"
    return pos, "unknown"


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

    def __init__(
            self,
            graph_dir,
            *,
            keep_x=None,
            scaler=None,
            mask_drift=True,
            drop=None,
            # --- transform knobs ---
            fix_units: bool = False,
            cutoff_A: float = 10.0,
            max_neighbors: Optional[int] = 64,
            sanity_checks: bool = True,
    ):
        if drop is None:
            drop = []
        self.keep_x = (None if keep_x is None else torch.as_tensor(keep_x, dtype=torch.long))
        self.paths = []
        self.meta = []
        self._cache = {}          # path -> Data (transformed)
        self._cache_size = 20000  # set to len(dataset) if you have RAM; or a smaller LRU

        for f in sorted(Path(graph_dir).rglob("frame_*_graph.pt")):
            try:
                g = torch.load(f, map_location="cpu")
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

        # transform knobs
        self.fix_units = fix_units
        self.cutoff_A = float(cutoff_A)
        self.max_neighbors = max_neighbors
        self.sanity_checks = sanity_checks

    def __len__(self):
        return len(self.paths)

    def __getitem__(self, idx):

        p = self.paths[idx]
        g_cached = self._cache_get(p)
        if g_cached is not None:
            return g_cached.clone()

        g = torch.load(self.paths[idx], map_location="cpu")

        # --------- INLINE TRANSFORM (units -> Å, rebuild edges once) ----------
        if self.fix_units:
            already_ok = getattr(g, "_units_A", False) and getattr(g, "_cutoff_A", None) == self.cutoff_A
            if not already_ok:
                # 1) Scale positions to Å if needed
                g.pos, flag = _ensure_angstroms(g.pos)

                # 2) If we converted nm->Å, also convert label y and u_raw['ligand_drift']
                if flag == "nm->A":
                    if hasattr(g, "y") and g.y is not None:
                        g.y = g.y * 10.0
                    try:
                        li = self._feat_index["ligand_drift"]
                        # ensure tensor
                        if not torch.is_tensor(g.u_raw):
                            g.u_raw = torch.as_tensor(g.u_raw, dtype=torch.float32)
                        g.u_raw[li] = g.u_raw[li] * 10.0
                    except Exception:
                        pass  # don't break if u_raw missing or sized differently

                # 3) Rebuild radius graph with Å cutoff
                old_E = int(g.edge_index.size(1)) if hasattr(g, "edge_index") else -1
                edge_index = radius_graph(
                    g.pos, r=self.cutoff_A, loop=False, max_num_neighbors=self.max_neighbors
                )
                row, col   = edge_index
                edge_attr  = (g.pos[row] - g.pos[col]).norm(dim=1, keepdim=True).to(torch.float32)

                g.edge_index   = edge_index
                g.edge_attr    = edge_attr
                g._units_A     = True
                g._cutoff_A    = self.cutoff_A
                g._unit_flag   = flag
                g._edges_changed = (old_E != -1) and (int(edge_index.size(1)) != old_E)

                # 4) Light sanity checks (non-fatal warnings)
                if self.sanity_checks:
                    with torch.no_grad():
                        # NN median should be ~3.8 Å
                        N = g.pos.size(0)
                        k = min(N, 512)
                        idx_k = torch.arange(k)
                        d = torch.cdist(g.pos.index_select(0, idx_k), g.pos.index_select(0, idx_k))
                        d[d == 0] = float("inf")
                        nn_med = float(d.min(dim=1).values.median())
                        if not (3.0 < nn_med < 5.5):
                            warnings.warn(f"[units] NN median {nn_med:.2f} Å (flag={flag}) looks off")
                        # No edge beyond cutoff
                        dmax = float(edge_attr.max())
                        if dmax > self.cutoff_A + 1e-4:
                            warnings.warn(f"[units] Found edge {dmax:.2f} Å > cutoff {self.cutoff_A} Å")

        # ------------- build u (after any possible unit conversion) ------------
        u = g.u_raw.clone() if torch.is_tensor(g.u_raw) else torch.as_tensor(g.u_raw, dtype=torch.float32)

        # optional z-scaling (μ,σ) fitted on *train* only
        if self.scaler is not None:
            mu, std = self.scaler[g.meta["protein_id"]]
            u = (u - mu) / std  # assumes mu/std vector matches u_raw layout

        # --- on-the-fly feature masking ----------------------------------
        if self.keep_x is not None:
            g.x = g.x[:, self.keep_x]

        # attach normalized global vector; keep raw in g.u_raw
        g.u = u

        # append z_energy as an extra node channel if keep_x provided
        if (self.keep_x is not None) and (self.keep_x.numel() > 0):
            # energy_variance index is 1 in _feat_index
            z_energy = u[1:2].view(1, 1).expand(g.x.size(0), 1).to(g.x.device)
            g.x = torch.cat([g.x, z_energy], dim=-1)

        # ---- attach contacts (safe: skip if edges changed) -------------------
        g = self._attach_contacts_safe(g, self.paths[idx])

        # compute y_scaled *after* any unit conversion
        if hasattr(g, "y") and g.y is not None:
            g.y_scaled = g.y * Y_SCALE

        # cache the transformed sample in-memory
        self._cache_put(p, g)

        return g

    def _attach_contacts_safe(self, g, graph_path):
        """
        If a sibling *.contact.pt* exists, attach g.is_contact aligned with edge_index.
        If edges were rebuilt by our inline transform, skip to avoid misalignment.
        """
        contact_path = graph_path.with_suffix(".contact.pt")
        if contact_path.exists() and not getattr(g, "_edges_changed", False):
            g.is_contact = torch.load(contact_path, map_location="cpu")
        else:
            # Fallback – create an all-false tensor so code won't break
            E = g.edge_index.size(1)
            g.is_contact = torch.zeros(E, dtype=torch.bool)
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
                self.scaler[pid] = (
                    mu_sigma_dict["_GLOBAL_"]["mu"],
                    mu_sigma_dict["_GLOBAL_"]["std"],
                )
            else:
                raise KeyError(f"[DriftDataset] No mu/std for protein {pid}")

    # simple in-memory cache (FIFO)
    def _cache_put(self, key, g):
        if len(self._cache) >= self._cache_size:
            # simple FIFO eviction; replace with LRU if you want
            self._cache.pop(next(iter(self._cache)))
        self._cache[key] = g

    def _cache_get(self, key):
        return self._cache.get(key, None)
