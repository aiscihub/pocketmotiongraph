# g_config_patched.py – dynamic feature‑group generator

from __future__ import annotations

from typing import Dict, List
import torch

__all__ = ["make_feature_groups"]


# -----------------------------------------------------------------------------
# Constants that *never* change in the current pipeline
# -----------------------------------------------------------------------------
PHYS_DIM   = 9    # hydropathy, B‑factor, SASA, volume, 3‑way SS one‑hot
PLIP_DIM   = 21   # 7 contact types × (initial, final, Δ)


def make_feature_groups(
        g_example: "torch_geometric.data.Data", *,
        len_domain_flags: int = 4,
        include_energy_slot: bool = True,
) -> Dict[str, List[int]]:
    """Return a dict mapping *feature‑group name*  → list of column indices.

    Parameters
    ----------
    g_example : a loaded PyG graph whose ``x`` tensor has the *full* set of
        node features produced by ``builder.py`` (before any ablation).
    len_domain_flags : number of domain‑flag columns present in ``g.u``;
        default = 4, adjust if you change the flag list.
    include_energy_slot : if ``True`` reserve a placeholder column index
        *at the end* that the Dataset class may append for energy variance.
    """
    x_dim = g_example.x.shape[1]

    # Determine presence and size of PLIP block (appended by builder)
    remaining_after_phys = x_dim - PHYS_DIM
    has_plip = remaining_after_phys >= PLIP_DIM

    # Compute ESM dimension by removing phys and PLIP dims
    dim_esm = max(
        remaining_after_phys - (PLIP_DIM if has_plip else 0),
        0
    )

    groups: Dict[str, List[int]] = {}
    idx = 0

    # 1) ESM embeddings (optional)
    if dim_esm > 0:
        groups["esm_embeddings"] = list(range(idx, idx + dim_esm))
        idx += dim_esm

    # 2) Physicochemical + SS (always 9)
    groups["physicochemical"] = list(range(idx, idx + PHYS_DIM))
    idx += PHYS_DIM

    # 3) PLIP node fingerprint (optional)
    if has_plip:
        groups["plip_fingerprint"] = list(range(idx, idx + PLIP_DIM))
        idx += PLIP_DIM

    # 4) Energy variance slot (appended later by Dataset)
    if include_energy_slot:
        groups["energy_variance"] = [idx]
        idx += 1

    return groups
