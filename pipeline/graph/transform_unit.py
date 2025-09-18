# transforms_units.py
from __future__ import annotations
import warnings
import torch
from torch_geometric.data import Data
from torch_geometric.nn import radius_graph


def ensure_angstroms(pos: torch.Tensor, sample: int = 512):
    """
    Detect nm vs Å via nearest-neighbor distance.
    Returns (pos_in_Å, flag: 'A' | 'nm->A' | 'unknown').
    """
    assert pos.ndim == 2 and pos.size(1) == 3, f"pos must be [N,3], got {tuple(pos.shape)}"
    N = pos.size(0)
    # deterministic sample (first K) to keep transforms reproducible
    k = min(N, sample)
    idx = torch.arange(k, device=pos.device)
    with torch.no_grad():
        d = torch.cdist(pos.index_select(0, idx), pos.index_select(0, idx), p=2)
        d[d == 0] = float("inf")
        nn_med = float(d.min(dim=1).values.median())
    if 0.25 <= nn_med <= 0.70:     # ~0.38 nm CA–CA
        return pos * 10.0, "nm->A"
    if 2.5 <= nn_med <= 5.5:       # ~3.8 Å CA–CA
        return pos, "A"
    return pos, "unknown"


class ToAngstromsRadiusGraph:
    """
    PyG transform: ensure coordinates are in Å, (re)build radius graph with a given cutoff,
    and refresh edge_attr (= edge length in Å). Runs once per Data (cached).
    """
    def __init__(self, cutoff_A: float = 10.0, max_neighbors: int | None = 64,
                 check_nn: bool = True):
        self.cutoff_A = float(cutoff_A)
        self.max_neighbors = max_neighbors
        self.check_nn = check_nn

    def __call__(self, data: Data) -> Data:
        # Skip if already processed for this cutoff
        if getattr(data, "_units_A", False) and getattr(data, "_cutoff_A", None) == self.cutoff_A:
            return data

        # 1) Units → Å
        posA, flag = ensure_angstroms(data.pos)
        data.pos = posA

        # 2) Build radius graph (single-graph Data: no batch arg)
        edge_index = radius_graph(
            data.pos, r=self.cutoff_A, loop=False,
            max_num_neighbors=self.max_neighbors
        )
        row, col = edge_index
        edge_attr = (data.pos[row] - data.pos[col]).norm(dim=1, keepdim=True).to(torch.float32)

        data.edge_index = edge_index
        data.edge_attr  = edge_attr
        data._units_A   = True
        data._cutoff_A  = self.cutoff_A
        data._unit_flag = flag

        # 3) Light sanity checks (non-fatal warnings)
        with torch.no_grad():
            # Nearest-neighbor (should be ~3.8 Å)
            if self.check_nn:
                N = data.pos.size(0)
                k = min(N, 512)
                idx = torch.arange(k)
                d = torch.cdist(data.pos.index_select(0, idx), data.pos.index_select(0, idx))
                d[d == 0] = float("inf")
                nn_med = float(d.min(dim=1).values.median())
                if not (3.0 < nn_med < 5.5):
                    warnings.warn(f"[units] NN median {nn_med:.2f} Å (flag={flag}) looks off")

            # No edge should exceed cutoff (tiny epsilon for FP)
            dmax = float(edge_attr.max())
            if dmax > self.cutoff_A + 1e-4:
                warnings.warn(f"[units] Found edge {dmax:.2f} Å > cutoff {self.cutoff_A} Å")

        return data
