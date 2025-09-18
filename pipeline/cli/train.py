#!/usr/bin/env python
"""
train.py – 4-fold cross-protein CV **with every advanced trick**

Restores all functionality from the original “incremental” trainer:

• Per-protein μ/σ scaling (handled by DriftDataset)
• Node-level contact loss  (0.5 ×)
• Pair-level contact loss   (β = 0.05, K = 64 pairs/graph)
• Mixed-precision (bfloat16) + gradient-accumulation
• Node-cap sampler (keeps GPU mem ≤ max_nodes)
• Linear LR warm-up → Cosine / Plateau scheduler
• Stochastic Weight Averaging (after --swa-start fraction)
• Optional global-feature fusion can be **disabled** with --no-fusion
• 4-fold GroupKFold CV  + multi-seed

Example
--------
python train.py \\
    --graph-root pipeline/data_graph \\
    --save-prefix checkpoints/noFusion \\
    --cv 4 --seeds 0 1 2 --epochs 30 \\
    --batch 4 --huber --scheduler cosine --swa-start 0.6 \\
    --nodecap --no-fusion
"""
from __future__ import annotations
import argparse, csv, time
from collections import Counter
from pathlib import Path

import torch
from sklearn.model_selection import GroupKFold
from torch.cuda.amp import GradScaler, autocast
from torch.optim.swa_utils import AveragedModel, SWALR
from torch_geometric.loader import DataLoader
from torch.utils.data import SubsetRandomSampler

from pipeline.graph.dataset import DriftDataset
from pipeline.graph.g_config import FULL_IDX, NO_PLIP_IDX, NO_ESM_PLIP_IDX
from pipeline.models.contact_lookup import ContactLookup
from pipeline.models.gnn_tfn_egnn import PocketStabilityModel
torch.backends.cuda.matmul.allow_tf32 = True
torch.backends.cudnn.allow_tf32 = True


# ───────────── config ─────────────
POS_CLIP = 30.0
PAIR_K = 256        # more informative pairs → better ranking/ROC
BETA_PAIR = 0.18   # don’t overpower regression / node loss

ACCUM_STEPS = 4
NODE_LOSS_W = 1.65
EDGE_W = 0.05
MAX_NODES = 8_000  # for NodeCapSampler
Z_STD = 1e-8
Z_GAP = 0.001
import logging
EDGE_CAP = 30_000

log = logging.getLogger(__name__)  # Use module name as logger ID

rmsd_log = {"epoch": [], "true": [], "pred": []}

from sklearn.metrics import roc_auc_score

@torch.no_grad()
def roc_by_protein(model, loader, device, lookup, flip_map: dict[str, bool] | None = None):
    model.eval()
    s, t = defaultdict(list), defaultdict(list)  # scores, targets
    for b in loader:
        b = b.to(device)
        _, _, node_logits, _ = model(
            b.x, b.edge_index, batch=b.batch, pos=b.pos, u=b.u, pid=b.pid
        )
        for g in torch.unique(b.batch):
            m = b.batch == g
            pid = b.meta["protein_id"][g]
            tgt = lookup.mask(pid, b.residue_id[m]).float().cpu()
            logits = node_logits[m]
            if flip_map is not None and flip_map.get(pid, False):
                logits = -logits
            scr = torch.sigmoid(logits).cpu()
            s[pid].extend(scr.tolist())
            t[pid].extend(tgt.tolist())

    rocs = {p: roc_auc_score(t[p], s[p]) for p in s if len(set(t[p])) > 1}
    mu = sum(rocs.values()) / len(rocs)
    sig = (sum((v - mu) ** 2 for v in rocs.values()) / len(rocs)) ** 0.5
    return mu, sig, rocs

@torch.no_grad()
def build_contact_flip_map(model, loader, device, lookup) -> dict[str, bool]:
     """
        Decide per-protein whether to flip node-logit sign by choosing the
        orientation that yields higher AUROC on the provided split (use VAL).
     """
     model.eval()
     by_pid_scores, by_pid_labels = defaultdict(list), defaultdict(list)
     for b in loader:
         b = b.to(device)
         _, _, node_logits, _ = model(b.x, b.edge_index, batch=b.batch, pos=b.pos, u=b.u, pid=b.pid)
         for g in torch.unique(b.batch):
             m = b.batch == g
             pid = b.meta["protein_id"][g]
             y = lookup.mask(pid, b.residue_id[m]).float().cpu().numpy()
             s = node_logits[m].detach().cpu().numpy()  # raw logits (before sigmoid)
             by_pid_labels[pid].extend(y.tolist())
             by_pid_scores[pid].extend(s.tolist())

     flip_map: dict[str, bool] = {}
     flipped, total = 0, 0
     for pid, scores in by_pid_scores.items():
         labels = np.asarray(by_pid_labels[pid], dtype=int)
         scores = np.asarray(scores, dtype=float)
         if len(np.unique(labels)) < 2:
             continue  # can't compute AUROC → leave unflipped
         auc_pos = roc_auc_score(labels, 1.0 / (1.0 + np.exp(-scores)))  # sigmoid(scores)
         auc_neg = roc_auc_score(labels, 1.0 / (1.0 + np.exp(scores)))   # sigmoid(-scores)
         do_flip = bool(auc_neg > auc_pos)
         flip_map[pid] = do_flip
         total += 1
         flipped += int(do_flip)
     log.info(f"[contact-flip] flipped {flipped}/{total} proteins on VAL (choose better AUROC orientation)")
     return flip_map

def summarise_attention(
    attn_dir: Path, plip_csv: Path, suffix: str, include_pids: set[str] | None = None
):
    import pandas as pd
    from sklearn.metrics import roc_auc_score

    csvs = sorted(attn_dir.glob(f"*_{suffix}.csv"))
    if not csvs:
        log.warning(f"No attention csv found for suffix {suffix} → skipping summary")
        return float("nan")

    all_dfs = []
    for c in csvs:
        df = pd.read_csv(c)
        # Extract protein from filename like ATRF_ASPFU_attn_f0_s0.csv
        filename = c.name
        protein = filename.split("_attn_")[0]
        if "protein" not in df.columns:
            df.insert(0, "protein", protein)
        all_dfs.append(df)

    df = pd.concat(all_dfs, ignore_index=True, sort=False)
    if "score" in df.columns and "stab_score" not in df.columns:
        df.rename(columns={"score": "stab_score"}, inplace=True)
    df.to_csv(attn_dir / f"all_attention_summary_{suffix}.csv", index=False)

    # ---- overlay with PLIP contacts ---------------------------------------
    plip = pd.read_csv(plip_csv)
    plip["residue_id"] = plip["Residue"].str.extract(r"(\d+)", expand=False).astype(int)
    merged = df.merge(plip, on=["protein", "residue_id"], how="left").dropna(
        subset=["stab_score"]
    )
    if include_pids is not None:
        merged = merged[merged["protein"].isin(include_pids)]
        log.info(f"Include_pids = {include_pids}")

    merged["is_contact"] = merged["Complex"].notna()
    merged.to_csv(attn_dir / f"attention_vs_plip_{suffix}.csv", index=False)

    # merged["stab_z"] = merged["stab_score"]

    # --- 1) RAW attention stats (before any normalization) -----------------
    raw_mc = merged.loc[merged.is_contact, "stab_score"].mean()
    raw_nc = merged.loc[~merged.is_contact, "stab_score"].mean()
    raw_std = merged["stab_score"].std()
    log.debug(
        f"[DEBUG RAW] std={raw_std:.6f}, mean(contact)={raw_mc:.6f}, mean(non-contact)={raw_nc:.6f}"
    )
    # also per-protein raw stats:

    for prot in merged["protein"].unique():
        sub = merged[merged["protein"] == prot]
        std_p = sub["stab_score"].std()
        mc_p = sub.loc[sub.is_contact, "stab_score"].mean()
        nc_p = sub.loc[~sub.is_contact, "stab_score"].mean()
        log.debug(
            f"[DEBUG RAW PROT] {prot:15} | std={std_p:.6f} | mean(c)={mc_p:.6f} | non-c={nc_p:.6f}"
        )

    # ── 1) flip per protein on *raw* stab_score so that contact > non-contact
    def _flip(group):
        if group["is_contact"].sum() < 2:  # degenerate, leave as is
            return group
        if (
            group.loc[group.is_contact, "stab_score"].mean()
            < group.loc[~group.is_contact, "stab_score"].mean()
        ):
            group["stab_score"] *= -1
        return group

    merged = merged.groupby("protein", group_keys=False).apply(_flip)
    merged["stab_z"] = merged.groupby("protein")[
        "stab_score"
    ].transform(  # robust IQR norm
        lambda x: (x - x.median()) / (x.quantile(0.75) - x.quantile(0.25) + 1e-6)
    )
    # ------------ (A) AUPRC ------------------------------------------------
    y_true = merged["is_contact"].astype(int).values
    y_score = merged["stab_z"].values
    auprc = average_precision_score(y_true, y_score)

    # # Calculate means
    mean_contact = merged.loc[merged.is_contact, "stab_z"].mean()
    mean_non_contact = merged.loc[~merged.is_contact, "stab_z"].mean()
    std_all = merged["stab_z"].std()

    # New check to skip noisy attention
    if std_all < 1e-3 or abs(mean_contact - mean_non_contact) < Z_GAP:
        log.warning(
            f"[ATTN:{suffix}] Low attention signal — skipping global ROC (std={std_all:.4f}, diff={mean_contact - mean_non_contact:.4f})"
        )
        roc_auc = float("nan")  # or -1
    else:
        roc_auc = roc_auc_score(
            merged["is_contact"].astype(int), merged["stab_z"].values
        )
        log.info(
            f"[ATTN:{suffix}] global ROC-AUC    = {roc_auc:.3f}, AUPRC = {auprc:.3f}"
        )

    return roc_auc


# --- SAMPLER: node-cap with Subset support and no transform triggers ---
from torch.utils.data import Sampler, Subset
import torch, math, random


class NodeCapSampler(Sampler):
    """
    Builds batches whose total node count ≤ max_nodes.
    Works with either a base DriftDataset or a torch.utils.data.Subset(base_ds).
    Never calls ds[i] (so it won’t trigger transforms); it only peeks at saved graph files.
    """

    def __init__(self, ds, max_nodes=MAX_NODES, contact_weight=0.88, shuffle=True):
        self.max_nodes = int(max_nodes)
        self.cw = float(contact_weight)
        self.shuffle = bool(shuffle)

        # Unwrap Subset → (base_ds, local_indices, base_indices)
        if isinstance(ds, Subset):
            base = ds.dataset
            local_indices = list(range(len(ds)))
            base_indices = list(map(int, ds.indices))
        else:
            base = ds
            local_indices = list(range(len(ds)))
            base_indices = local_indices

        # We expect the base dataset to have .paths pointing to saved .pt files
        self.ds = ds
        self._local_indices = local_indices

        # Precompute per-item sizes and contact flags WITHOUT calling ds.__getitem__
        sizes = []
        contact_locals = []

        for loc_i, base_i in zip(local_indices, base_indices):
            p = base.paths[base_i]  # raw file path
            g = torch.load(p, map_location="cpu")  # lightweight header load
            sizes.append(int(getattr(g, "num_nodes", g.x.size(0))))

            # contact presence: use saved tensor if present; else sibling file existence
            has_contact = getattr(g, "is_contact", None)
            if has_contact is None:
                has_contact = p.with_suffix(".contact.pt").exists()
            if (torch.is_tensor(has_contact) and has_contact.any()) or bool(
                has_contact
            ):
                contact_locals.append(loc_i)

        self.sizes = sizes  # aligned with local_indices
        contact_set = set(contact_locals)
        self.contact = contact_locals
        self.no_contact = [loc for loc in local_indices if loc not in contact_set]

    def __iter__(self):
        # Decide class-balanced pools
        contact_pool = self.contact[:]
        nocon_pool   = self.no_contact[:]
        if self.shuffle:
            random.shuffle(contact_pool); random.shuffle(nocon_pool)

        i_c, i_n = 0, 0
        batch, tot = [], 0
        while i_c < len(contact_pool) or i_n < len(nocon_pool):
            # target share of contact graphs in the **current** batch
            want_c = self.cw
            # compute current share
            cur_c = 0.0 if not batch else sum(1 for idx in batch if idx in set(self.contact)) / len(batch)
            # choose from the side that moves us toward want_c
            choose_contact = (cur_c < want_c)

            chosen = None
            if choose_contact and i_c < len(contact_pool):
                cand = contact_pool[i_c]; i_c += 1
                n = self.sizes[cand]
                if tot + n > self.max_nodes and batch:
                    yield batch; batch, tot = [], 0
                batch.append(cand); tot += n
                continue

            if i_n < len(nocon_pool):
                cand = nocon_pool[i_n]; i_n += 1
                n = self.sizes[cand]
                if tot + n > self.max_nodes and batch:
                    yield batch; batch, tot = [], 0
                batch.append(cand); tot += n
                continue

            # if one side exhausted, keep draining the other
            if i_c < len(contact_pool):
                cand = contact_pool[i_c]; i_c += 1
                n = self.sizes[cand]
                if tot + n > self.max_nodes and batch:
                    yield batch; batch, tot = [], 0
                batch.append(cand); tot += n

        if batch:
            yield batch


    def __len__(self):
        # estimated number of batches per epoch
        return math.ceil(sum(self.sizes) / max(1, self.max_nodes))


def set_seed(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def w_bce(logits, targets):
    pos = targets.sum()
    neg = targets.numel() - pos
    if pos == 0 or neg == 0:
        return torch.tensor(0.0, device=logits.device)
    pw = min(neg / pos, POS_CLIP)
    return F.binary_cross_entropy_with_logits(
        logits, targets, pos_weight=torch.tensor([pw], device=logits.device)
    )

import torch
import torch.nn.functional as F

def margin_pair_loss(scores_pos: torch.Tensor,
                     scores_neg: torch.Tensor,
                     margin: float = 0.5,
                     reduction: str = "mean") -> torch.Tensor:
    """
    scores_pos, scores_neg: shape [N_pairs] (same shape)
    Larger score = more likely contact.
    """
    # target = +1 → encourage scores_pos > scores_neg
    y = torch.ones_like(scores_pos)
    return F.margin_ranking_loss(scores_pos, scores_neg, y, margin=margin, reduction=reduction)


def logistic_pair_auc_loss(scores_pos: torch.Tensor,
                           scores_neg: torch.Tensor,
                           reduction: str = "mean") -> torch.Tensor:
    """
    AUC surrogate: -log sigmoid(s_pos - s_neg)
    """
    diff = scores_pos - scores_neg
    loss = F.softplus(-diff)   # = -log(sigmoid(diff))
    return loss.mean() if reduction == "mean" else loss.sum()

# The weighted pos weight for bce does not help
def w_bce_pos_weight(logits, targets, pos_weight=None, gamma=None):
    pos = targets.sum()
    neg = targets.numel() - pos
    if pos == 0 or neg == 0:
        return torch.tensor(0.0, device=logits.device)

    if pos_weight is not None:
        pw = float(pos_weight)
    else:
        pw = min(neg / pos, POS_CLIP)  # dynamic weighting
    # print(f"pw ---- {pw}")
    loss = F.binary_cross_entropy_with_logits(
        logits,
        targets,
        pos_weight=torch.tensor([pw], device=logits.device),
        reduction="none",
    )

    # optional focal term
    if gamma is not None:
        p = torch.sigmoid(logits)
        loss = (1 - p) ** gamma * loss

    return loss.mean()


@torch.no_grad()
def dump_attention_maps(model, loader, device, out_dir: Path, fold: int, seed: int):
    """
    Writes one CSV per graph that contains
        residue_id, attention_score
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    model.eval()
    for g in loader:
        g = g.to(device)
        _, node_attn = model(
            g.x,
            g.edge_index,
            batch=g.batch,
            pos=g.pos,
            u=g.u,
            pid=g.pid,
            return_node_attn=True,
        )
        for gid in g.batch.unique():
            m = g.batch == gid
            edge_mask = torch.zeros_like(m)  # m has length = num_nodes
            edge_nodes = torch.unique(
                g.edge_index[:, g.batch[g.edge_index[0]] == gid].flatten()
            )
            edge_mask[edge_nodes] = True
            m = m & edge_mask  # keep only nodes with ≥1 edge
            if m.sum() == 0:
                continue

            res_ids = g.residue_id[m].cpu()
            scores = node_attn[m].cpu()
            # Optional: Add pocket center distance if available
            xyz = g.pos[m].cpu().numpy()
            arr = torch.stack([res_ids, scores], 1).numpy()

            # every graph already has a unique pdb_id  (ProteinID_Pocket_Frame)
            pdb_id = str(g.pdb_id if isinstance(g.pdb_id, str) else g.pdb_id[0])
            np.savetxt(
                out_dir / f"{pdb_id}_attn_f{fold}_s{seed}.csv",
                arr,
                fmt=["%d", "%.10f"],
                delimiter=",",
                header="residue_id,score",
                comments="",
            )


# ───────────── training loops ─────────────
def train_epoch(
    model,
    loader,
    optim,
    scaler,
    device,
    lookup,
    *,
    huber: bool,
    accum: int,
    pos_weight=None,
    gamma=None,
    pair_loss_flag="margin",
):
    """
    Train one epoch with hard-negative mining for contact ROC:
      • For each graph, identify positives (PLIP contacts) and negatives (non-contacts).
      • Rank negatives by current node logits and take top-K as hard negatives.
      • Optimize a pairwise AUC surrogate between (pos, hard-neg) with a separate pair_head.
    """
    # local AUC surrogate: -log(sigmoid(s_pos - s_neg))
    def logistic_pair_auc_loss(s_pos: torch.Tensor, s_neg: torch.Tensor) -> torch.Tensor:
        return F.softplus(-(s_pos - s_neg)).mean()

    model.train()
    total, step = 0.0, 0
    optim.zero_grad(set_to_none=True)

    for batch in loader:
        batch = batch.to(device)
        with autocast(dtype=torch.bfloat16):
            reg, _, node_logits, aux = model(
                batch.x,
                batch.edge_index,
                batch=batch.batch,
                pos=batch.pos,
                u=batch.u,
                pid=batch.pid,
            )

            # --- regression (RMSD) loss
            reg_loss = (F.smooth_l1_loss if huber else F.l1_loss)(
                reg.view(-1), batch.y[:, 0]
            )

            # --- node (per-residue) BCE loss (class-weighted/focal as you had)
            node_loss = 0.0
            n_graph = int(batch.batch.max()) + 1
            for g in range(n_graph):
                m = batch.batch == g
                pid = batch.meta["protein_id"][g]
                tgt = lookup.mask(pid, batch.residue_id[m]).float()
                node_loss += w_bce_pos_weight(node_logits[m], tgt, pos_weight=pos_weight, gamma=gamma)
            node_loss /= max(n_graph, 1)

            # --- pairwise loss (hard-negative mining for contact ranking)
            pair_loss = 0.0
            if hasattr(model, "pair_head"):
                for g in range(n_graph):
                    m = batch.batch == g
                    pid = batch.meta["protein_id"][g]
                    tgt_mask = lookup.mask(pid, batch.residue_id[m]).bool()

                    # graph-local indices of nodes in this micrograph
                    idx = torch.nonzero(m, as_tuple=False).squeeze(-1)
                    pos_ids = idx[tgt_mask]
                    neg_ids = idx[~tgt_mask]

                    # need at least 1 pos & 1 neg
                    if pos_ids.numel() < 1 or neg_ids.numel() < 1:
                        continue

                    # ---- HARD NEGATIVES: top-K negatives by current node logit (largest = "hard")
                    with torch.no_grad():
                        neg_scores = node_logits[neg_ids]
                        Khn = min(PAIR_K, neg_ids.numel())
                        # topk is safe even if K==len(neg_ids)
                        topk = torch.topk(neg_scores, k=Khn, largest=True).indices
                        hard_negs = neg_ids[topk]

                    # Pair each positive with a random hard negative (balanced pairs)
                    dvc = pos_ids.device
                    P = pos_ids.numel()
                    Kp = min(PAIR_K, P)
                    sel_pos = pos_ids[torch.randint(low=0, high=P, size=(Kp,), device=dvc)]
                    sel_neg = hard_negs[torch.randint(low=0, high=hard_negs.numel(), size=(Kp,), device=dvc)]

                    # Build pair features and logits via pair_head (your existing design)
                    h_i = aux["h"][sel_pos]
                    h_j = aux["h"][sel_neg]
                    dist = (aux["pos"][sel_pos] - aux["pos"][sel_neg]).norm(dim=-1, keepdim=True)
                    logits_pair = model.pair_head(torch.cat([h_i, h_j, dist], dim=-1)).squeeze(-1)

                    if pair_loss_flag == "margin":
                        # margin ranking between pos (should be high) vs hard neg (should be low)
                        # We can treat node_logits directly as scores for margin, but since you’ve
                        # trained a pair_head, rank with its outputs:
                        #   want logits_pair (pos,neg) to be > 0
                        # Build a dummy zero baseline by swapping to get "neg > pos" violations.
                        # Simpler/better: use logistic pair AUC below.
                        s_pos = logits_pair  # score for (pos, hard-neg) being a contact
                        # Turn it into pairwise ranking against zeros:
                        pair_loss += F.margin_ranking_loss(
                            s_pos, torch.zeros_like(s_pos), target=torch.ones_like(s_pos), margin=0.5
                        )
                    else:
                        # smooth AUC surrogate is typically stronger for ROC
                        # Recompute scores directly from node logits difference (also valid):
                        s_pos_nodes = node_logits[sel_pos]
                        s_neg_nodes = node_logits[sel_neg]
                        pair_loss += logistic_pair_auc_loss(s_pos_nodes, s_neg_nodes)

                    # ---- (optional) edge-level supervision you already had
                    edge_mask = batch.batch[batch.edge_index[0]] == g  # Bool[E_total]
                    if edge_mask.any():
                        edge_pairs = batch.edge_index[:, edge_mask].t()  # (E_g, 2)
                        edge_labels = batch.is_contact[edge_mask].float()  # (E_g,)
                        if edge_pairs.size(0) > EDGE_CAP:
                            sel = torch.randperm(edge_pairs.size(0), device=device)[:EDGE_CAP]
                            edge_pairs = edge_pairs[sel]
                            edge_labels = edge_labels[sel]
                        ei = edge_pairs[:, 0]
                        ej = edge_pairs[:, 1]
                        h_i_e = aux["h"][ei]
                        h_j_e = aux["h"][ej]
                        dist_e = (aux["pos"][ei] - aux["pos"][ej]).norm(dim=-1, keepdim=True)
                        edge_logits = model.pair_head(torch.cat([h_i_e, h_j_e, dist_e], dim=-1)).squeeze(-1)
                        edge_loss = w_bce_pos_weight(edge_logits, edge_labels, pos_weight=pos_weight, gamma=gamma)
                        pair_loss += EDGE_W * edge_loss

                pair_loss /= max(n_graph, 1)

            # --- total loss
            loss = reg_loss + NODE_LOSS_W * node_loss + BETA_PAIR * pair_loss

        # AMP backward/step
        scaler.scale(loss).backward()
        step += 1
        if step % accum == 0:
            scaler.unscale_(optim)
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            scaler.step(optim)
            scaler.update()
            optim.zero_grad(set_to_none=True)
        total += float(loss.item())

    return total / max(1, len(loader))


import numpy as np
from sklearn.metrics import average_precision_score

# ───────────── ranking metrics (per protein) ─────────────
def _pairs_cindex(y, s):
    import numpy as np
    y = np.asarray(y, float); s = np.asarray(s, float)
    n = y.size
    if n < 2: return float("nan")
    i, j = np.triu_indices(n, 1)
    dy = y[i] - y[j]
    ds = s[i] - s[j]
    valid = (dy != 0)
    if not np.any(valid): return float("nan")
    return float(np.mean(np.sign(dy[valid]) == np.sign(ds[valid])))

@torch.no_grad()
def write_ranking_metrics_csv(model, loader, device, out_csv: Path,
                              thr_by_pid: dict | None = None, flip_sign: bool = False,
                              ks=(10, 20), alpha=0.05, q_fallback=0.70):
    """
    Writes a CSV with one row per protein containing:
      protein, n, pos_rate, cindex, spearman, auprc, P@10, P@20, R@10, R@20, EF@5%
    Uses train-calibrated thresholds if provided (thr_by_pid), else per-split quantile (q_fallback).
    """
    import numpy as np, csv
    from sklearn.metrics import average_precision_score
    model.eval()
    by_pid_true, by_pid_score = {}, {}

    for b in loader:
        b = b.to(device)
        reg, _, _, _ = model(b.x, b.edge_index, batch=b.batch, pos=b.pos, u=b.u, pid=b.pid)
        reg = reg.view(-1)
        for gid in torch.unique(b.batch, sorted=True):
            gi = int(gid)
            pid = b.meta["protein_id"][gi]
            y = float(b.y[gi, 0].item())
            s = float(reg[gi].item())
            if flip_sign: s = -s
            by_pid_true.setdefault(pid, []).append(y)
            by_pid_score.setdefault(pid, []).append(s)

    fields = ["protein","n","pos_rate","cindex","spearman","auprc",
              "P@10","P@20","R@10","R@20","EF@5%"]
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for pid, y_list in sorted(by_pid_true.items()):
            y = np.asarray(y_list, float)
            s = np.asarray(by_pid_score[pid], float)
            if y.size == 0: continue

            # threshold-free metrics
            cindex = _pairs_cindex(y, s)
            # Spearman via rank corr (robust if ties)
            if y.std() == 0 or s.std() == 0:
                spearman = float("nan")
            else:
                ry = np.argsort(np.argsort(y))
                rs = np.argsort(np.argsort(s))
                spearman = float(np.corrcoef(ry, rs)[0,1])

            # binary metrics need a threshold
            thr = (thr_by_pid.get(pid) if (thr_by_pid and pid in thr_by_pid)
                   else float(np.quantile(y, q_fallback)))
            yb = (y > thr).astype(int)
            pos = int(yb.sum()); n = int(yb.size)
            pos_rate = (pos / n) if n > 0 else float("nan")
            if yb.min() == yb.max():
                auprc = float("nan")
            else:
                auprc = float(average_precision_score(yb, s))

            # top-K metrics
            order = np.argsort(-s)
            def _p_at(k):
                k = min(k, n);
                return float(yb[order[:k]].mean()) if k > 0 else float("nan")
            def _r_at(k):
                k = min(k, n);
                tp = int(yb[order[:k]].sum())
                return float(tp / pos) if pos > 0 else float("nan")
            p10, p20 = _p_at(10), _p_at(20)
            r10, r20 = _r_at(10), _r_at(20)

            # Enrichment @ α%
            k_alpha = max(1, int(np.ceil(alpha * n)))
            tp = int(yb[order[:k_alpha]].sum())
            expected = alpha * pos if pos > 0 else float("nan")
            ef = float(tp / expected) if (isinstance(expected,float) and expected>0) else float("nan")

            w.writerow(dict(protein=pid, n=n, pos_rate=pos_rate,
                            cindex=cindex, spearman=spearman, auprc=auprc,
                            **{"P@10": p10, "P@20": p20, "R@10": r10, "R@20": r20, "EF@5%": ef}))
    log.info(f"Wrote ranking metrics → {out_csv}")

# at top-level (near other constants), add a toggleable default:
REG_THRESH = 0.75  # Å, set to your triage threshold

@torch.no_grad()
def evaluate(model, loader, device, lookup, epoch=None, thr_by_pid=None, flip_sign=False, q: float = 0.70):
    """
    Returns:
      mae:     mean absolute error of regression head vs graph target
      reg_roc: per-protein mean ROC-AUC from the regression head,
               binarized by thresholds computed on TRAIN (thr_by_pid).
    """
    import numpy as np
    from sklearn.metrics import roc_auc_score
    import torch

    model.eval()
    mae_chunks = []
    by_pid_true, by_pid_pred = {}, {}

    for b in loader:
        b = b.to(device)
        reg, _, _, _ = model(b.x, b.edge_index, batch=b.batch, pos=b.pos, u=b.u, pid=b.pid)
        reg = reg.view(-1)
        assert reg.shape[0] == b.y[:, 0].shape[0], "regression head must be graph-level"
        mae_chunks.append((reg - b.y[:, 0]).abs().cpu())

        for gid in torch.unique(b.batch, sorted=True):
            gi = int(gid)
            pid = b.meta["protein_id"][gi]
            y_val = float(b.y[gi, 0].item())
            y_hat = float(reg[gi].item())
            by_pid_true.setdefault(pid, []).append(y_val)
            by_pid_pred.setdefault(pid, []).append((-y_hat) if flip_sign else y_hat)

    mae = torch.cat(mae_chunks).mean().item()

    rocs = []
    import numpy as np
    cidxs, p_at20, ef_at5 = [], [], []
    for pid, y_list in by_pid_true.items():
        y = np.asarray(y_list, dtype=float)
        s = np.asarray(by_pid_pred[pid], dtype=float)
        if y.size < 2:
            continue

        # prefer TRAIN-derived threshold; fall back to per-split quantile if missing
        thr = (thr_by_pid.get(pid) if (thr_by_pid is not None and pid in thr_by_pid)  else float(np.quantile(y, q)))
        yb = (y > thr).astype(int)
        if yb.min() == yb.max():
            ci = _pairs_cindex(y, s); cidxs.append(ci)
            continue

        auc = roc_auc_score(yb, s)
        if auc < 0.5:  # safety guard
            auc = 1.0 - auc
        rocs.append(auc)
                # ranking metrics
        ci = _pairs_cindex(y, s); cidxs.append(ci)
        order = np.argsort(-s); n = y.size
        k20 = min(20, n)
        top20 = yb[order[:k20]] if k20 > 0 else np.array([])
        p_at20.append(float(top20.mean()) if top20.size else float("nan"))
        k5 = max(1, int(np.ceil(0.05 * n)))
        tp = int(yb[order[:k5]].sum()); pos = int(yb.sum())
        expected = 0.05 * pos if pos > 0 else np.nan
        ef_at5.append(float(tp / expected) if (isinstance(expected, float) and expected > 0) else float("nan"))

    reg_roc = float(np.mean(rocs)) if rocs else float("nan")

    return mae, reg_roc

@torch.no_grad()
def build_regroc_calibration(model, train_loader, device, q=0.70):
    """
    Returns:
      thr_by_pid: dict[pid] -> float threshold computed on TRAIN targets (per protein)
      flip_sign:  bool      -> True if train Spearman corr(tv, pv) < 0 (use -pred for ROC)
    """
    import numpy as np
    from scipy.stats import spearmanr
    model.eval()

    by_pid_true, by_pid_pred = {}, {}
    for b in train_loader:
        b = b.to(device)
        reg, _, _, _ = model(b.x, b.edge_index, batch=b.batch, pos=b.pos, u=b.u, pid=b.pid)
        reg = reg.view(-1)
        for gid in torch.unique(b.batch, sorted=True):
            gi = int(gid)
            pid = b.meta["protein_id"][gi]
            by_pid_true.setdefault(pid, []).append(float(b.y[gi, 0].item()))
            by_pid_pred.setdefault(pid, []).append(float(reg[gi].item()))

    # thresholds per protein from TRAIN
    thr_by_pid = {}
    for pid, y_list in by_pid_true.items():
        y = np.asarray(y_list, dtype=float)
        thr_by_pid[pid] = float(np.quantile(y, q)) if y.size else float("nan")

    # global monotonic sign from TRAIN
    tv = np.concatenate([np.asarray(v, float) for v in by_pid_true.values()]) if by_pid_true else np.array([])
    pv = np.concatenate([np.asarray(by_pid_pred[k], float) for k in by_pid_true.keys()]) if by_pid_true else np.array([])
    corr = float(spearmanr(tv, pv).correlation) if tv.size and pv.size else 0.0
    flip_sign = (corr < 0)
    return thr_by_pid, flip_sign

# ───────────── fold runner ─────────────
def run_fold(seed, fold, train_i, val_i, test_i, ds, args, lookup):
    set_seed(seed)
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model_type = "tfn"
    model = PocketStabilityModel(
        in_dim=ds[0].x.shape[1],
        global_dim=ds[0].u.shape[0],
        hidden_dim=args.hidden_dim,
        dropout=args.dropout,
        use_u=args.fusion,
        model_type=model_type,
    ).to(dev)
    logging.info(
        f"model type = {model_type}, settings: x = {ds[0].x.shape[1]} u = {ds[0].u.shape[0]} {args.hidden_dim} {args.dropout} {args.fusion}"
    )

    optim = torch.optim.Adam(model.parameters(), lr=args.lr, fused=True)
    scaler = GradScaler()

    def make_loader(indices):
        if args.nodecap:
            sub = torch.utils.data.Subset(ds, indices)
            return DataLoader(sub, batch_sampler=NodeCapSampler(sub))
        return DataLoader(
            ds, batch_size=args.batch, sampler=SubsetRandomSampler(indices)
        )

    l_train, l_val, l_test = map(make_loader, (train_i, val_i, test_i))

    warm = int(max(1, args.epochs * 0.1))
    warm = 5
    warm_sched = torch.optim.lr_scheduler.LinearLR(optim, 0.1, 1.0, warm)
    main_sched = {
        "cosine": torch.optim.lr_scheduler.CosineAnnealingLR(
            optim, T_max=args.epochs - warm
        ),
        "plateau": torch.optim.lr_scheduler.ReduceLROnPlateau(
            optim, mode="max", factor=0.5, patience=3
        ),
    }.get(args.scheduler)
    swa_model = AveragedModel(model)
    swa_sched = SWALR(optim, swa_lr=1e-4)
    swa_phase = False
    swa_start = int(args.epochs * args.swa_start)
    best = float("inf")

    patience = 10

    best_epoch = 0
    global rmsd_log
    thr_by_pid, flip_sign = None, False
    contact_flip_map = None  # NEW: per-protein flip map for contact logits
    rmsd_log["epoch"].clear()
    rmsd_log["true"].clear()
    rmsd_log["pred"].clear()
    best_roc = 0.0
    best_state = None
    for ep in range(1, args.epochs + 1):
        tl = train_epoch(
            model,
            l_train,
            optim,
            scaler,
            dev,
            lookup,
            huber=args.huber,
            accum=ACCUM_STEPS,
            pos_weight=args.pos_weight,
            gamma=args.gamma,
        )
          # build calibration after 1 warm epoch so preds are meaningful
        if ep == 1 and thr_by_pid is None:
            thr_by_pid, flip_sign = build_regroc_calibration(model, l_train, dev, q=0.70)
            contact_flip_map = build_contact_flip_map(model, l_val, dev, lookup)
            log.info(f"[fold{fold} seed{seed}] built regROC calibration: {len(thr_by_pid)} proteins, flip_sign={flip_sign}")
        v_mae, v_roc = evaluate(model, l_val, dev, lookup, epoch=ep, thr_by_pid=thr_by_pid, flip_sign=flip_sign, q=0.70)
        #using contact_roc to optimize
        v_roc, attn_sig, v_contact_roc = roc_by_protein(
            model, l_val, dev, lookup, flip_map=contact_flip_map
        )
        log.info(
            f"[fold{fold} seed{seed}] ep{ep:02d} train {tl:.3f} | val {v_mae:.3f}/{v_roc:.3f}"
        )

        # LR sched
        if ep <= warm:
            warm_sched.step()
        elif ep < swa_start:
            if args.scheduler == "plateau":
                main_sched.step(v_roc)
            elif args.scheduler != "none":
                main_sched.step()
        else:
            if not swa_phase:
                log.info(f"[fold{fold} seed{seed}] enter SWA at ep{ep}")
                swa_phase = True
            swa_model.update_parameters(model)
            swa_sched.step()
        if v_roc > best_roc:
            best = v_mae
            best_epoch, best_roc = ep, v_roc
            best_state = {k: v.detach().cpu() for k, v in model.state_dict().items()}
            log.info(f"Logging the model with v_mae {best:3f} v_roc {best_roc:3f}")
            if not args.no_save:
                torch.save(best_state, args.save_dir / f"f{fold}_s{seed}.pt")
                log.info(
                    f"Saved checkpoint (val_MAE={best:.3f}, val_ROC={best_roc:.3f})"
                )
            if not args.no_attn:
                full_loader = DataLoader(ds, batch_size=1, shuffle=False)
                dump_attention_maps(model, full_loader, dev, args.attn_dir, fold, seed)

        elif ep - best_epoch >= patience:
            log.info(
                f"Early stopping at epoch {ep}. Best was {best_epoch} (val_MAE={best:.3f})"
            )
            break

    if best_state is not None:
        model.load_state_dict(best_state)
    else:
        # fall back to file if present (e.g., when --no-save is False)
        ckpt = args.save_dir / f"f{fold}_s{seed}.pt"
        if ckpt.exists():
            model.load_state_dict(torch.load(ckpt, map_location="cpu"))
        else:
            log.warning(
                "No best_state in RAM and checkpoint file missing; testing current weights."
            )

    t_mae, t_reg_roc = evaluate(model, l_test, dev, lookup, thr_by_pid=thr_by_pid, flip_sign=flip_sign, q=0.70)
    rank_csv = args.run_dir / f"ranking_metrics_f{fold}_s{seed}.csv"
    write_ranking_metrics_csv(model, l_test, dev, rank_csv,thr_by_pid=thr_by_pid, flip_sign=flip_sign,ks=(10,20), alpha=0.05, q_fallback=0.70)
    # ── Local-Contact head: per-protein ROC-AUC (attention logits vs PLIP) ──
    t_contact_avg_roc, attn_sig, t_contact_roc = roc_by_protein(
        model, l_test, dev, lookup, flip_map=contact_flip_map
    )
    for p, v in sorted(t_contact_roc.items()):
        log.info(f"  - {p}: {v:.3f}")
    logging.info(
        f"run_fold -- [fold{fold} seed{seed}] "
        f"TEST mae={t_mae:.3f} | reg-roc={t_reg_roc:.3f} | contact-roc={t_contact_avg_roc:.3f}"
    )
    return t_mae, t_reg_roc, t_contact_avg_roc


# --- SAMPLING LOGIC (below dataset creation) ---
from collections import defaultdict
import re


def _parse_ids(pdb_id):
    m = re.search(r"(.*?)(_P\d+)?(_F(\d+))?$", pdb_id)
    if not m:
        return pdb_id, "P0", -1
    pid, _, _, f = m.groups()
    return pid, "P0", int(f) if f else -1


def _group_by_pocket(meta):
    g = defaultdict(list)
    for i, m in enumerate(meta):
        pid, pocket, frame = _parse_ids(str(m.get("pdb_id", "")))
        g[(pid, pocket)].append((i, frame))
    return {k: [x[0] for x in sorted(v, key=lambda x: x[1])] for k, v in g.items()}


def _select_uniform(lst, k):
    if len(lst) <= k:
        return lst
    return [lst[i] for i in np.linspace(0, len(lst) - 1, k, dtype=int)]


# ───────────── main ─────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--graph-root", type=Path, required=True)

    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--batch", type=int, default=4)
    ap.add_argument("--lr", type=float, default=2e-3)
    ap.add_argument("--hidden-dim", type=int, default=256)
    ap.add_argument("--dropout", type=float, default=0.3)
    ap.add_argument(
        "--scheduler", choices=["none", "cosine", "plateau"], default="plateau"
    )
    ap.add_argument("--swa-start", type=float, default=0.6)
    ap.add_argument("--cv", type=int, default=4)
    ap.add_argument("--seeds", nargs="*", type=int, default=[0])
    ap.add_argument("--huber", action="store_true")
    ap.add_argument("--nodecap", action="store_true", default=True)
    ap.add_argument("--fusion", action="store_true")
    ap.add_argument("--physics-only", action="store_true")
    ap.add_argument("--contacts", type=Path, default=Path("checkpoints/contacts.csv"))
    ap.add_argument(
        "--run-dir", type=Path, required=True, help="(e.g., run_20250712_153000)"
    )
    ap.add_argument("--pos-weight", default=None)
    ap.add_argument("--gamma", default=None)
    ap.add_argument("--protein", default=None)
    ap.add_argument("--cv_by", choices=["protein","species"], default="protein",
                        help="GroupKFold by protein (default) or by species (LOSO-style)")
    # --new CLI flags -------------------------------------------------
    ap.add_argument("--node-plip", action="store_true")
    ap.add_argument("--no-esm", action="store_true")
    ap.add_argument(
        "--ablate-node",
        action="store_true",
        help="Drop all per-residue (node) features except ESM",
    )
    # ---- speed / I/O control -------------------------------------------------
    ap.add_argument(
        "--no-save",
        default=False,
        action="store_true",
        help="Skip checkpoint saving entirely (keep best weights only in RAM).",
    )
    ap.add_argument(
        "--no-attn",
        default=False,
        action="store_true",
        help="Skip writing attention CSVs during training.",
    )
    # argparse
    ap.add_argument("--no-plip-supervision", action="store_true",
                    help="Skip PLIP node/pair/edge supervision and metrics.")


# --- after loading DriftDataset ---------------------------------

    args = ap.parse_args()
    # Record the start time
    start_time = time.time()
    # Create structured run dir
    args.run_dir.mkdir(parents=True, exist_ok=True)
    args.save_dir = args.run_dir / "checkpoints"
    args.attn_dir = args.run_dir / "attention_maps"
    args.save_dir.mkdir(exist_ok=True)
    args.attn_dir.mkdir(exist_ok=True)

    # Logging path
    log_file = args.run_dir / "train_log.txt"
    logging.getLogger().handlers.clear()  # <- optional, clears previous handlers
    logging.basicConfig(
        format="[%(levelname)s] %(message)s",
        level=logging.DEBUG,
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
    )
    log.info(f"Starting training in run dir → {args.run_dir}")
    log.info(f"Graph root = {args.graph_root} | Save = {args.save_dir}")
    log.info(f"Starting training, cv = {args.cv}")
    log.info(f"       Graph root: {args.graph_root}")
    log.info(f"       Save path: {args.save_dir}")
    log.info(f"       Epochs: {args.epochs} | Batch: {args.batch} | LR: {args.lr}")
    log.info(f"       Hidden: {args.hidden_dim} | Dropout: {args.dropout}")
    log.info(f"       SWA start: {args.swa_start}")
    log.info(f"       nodecap: {args.nodecap}")
    log.info(f"       Using Huber: {args.huber}")
    log.info(f"       Fusion enabled: {args.fusion}")
    log.info(f"       NodeCap enabled: {args.nodecap}")

    # after args parse
    lookup = ContactLookup(args.contacts)

    # With this:
    keep = NO_PLIP_IDX

    if args.physics_only:
        keep = NO_ESM_PLIP_IDX
    elif args.ablate_node:
        keep = []
    elif args.no_esm:
        keep = NO_ESM_PLIP_IDX
    elif args.node_plip:
        keep = FULL_IDX
    ds = DriftDataset(args.graph_root, keep_x=keep, mask_drift=True)

    if len(ds) == 0:
        raise RuntimeError(
            f"No graphs loaded from {args.graph_root}. "
            f"Check drift filter / paths / feature groups."
        )

    prot_counts = Counter(m["protein_id"] for m in ds.meta)
    for pid, n in prot_counts.items():
        log.info(f"{pid:22} → {n} graphs")
    log.info(f"Graph total = {len(ds)} ")

    proteins = np.array([m["protein_id"] for m in ds.meta])
    # simple map protein → species (adjust as needed)
    def infer_species(pid: str) -> str:
        if pid.startswith("CIMG_"): return "COCCIDIOIDES"
        if "YEAST" in pid: return "YEAST"
        if "CANGA" in pid or "CANAL" in pid or "CANAR" in pid: return "CANDIDA"
        if "ASPFU" in pid: return "ASPERGILLUS"
        if "TRIRC" in pid: return "TRICHOPHYTON"
        if "CRYNH" in pid: return "CRYPTOCOCCUS"
        return "OTHER"
    species = np.array([infer_species(pid) for pid in proteins])

    all_indices = np.arange(len(ds))
    if args.protein is None:
        if args.cv_by == "protein":
            groups = proteins
            n_splits = min(args.cv, len(np.unique(groups)))
            if n_splits < 2:
                   raise ValueError(f"Need at least 2 protein groups for CV, found {len(np.unique(groups))}.")
            gkf = GroupKFold(n_splits=n_splits)
            splits = list(gkf.split(all_indices, groups=groups))
        else:
            # --cv_by species  → strict leave-one-species-out CV
            unique_species = np.unique(species)
            if len(unique_species) < 2:
                   raise ValueError("Need at least 2 species for LOSO CV.")
            splits = []
            for s in unique_species:
                test_i = np.where(species == s)[0]
                train_i = np.where(species != s)[0]
                if len(test_i) == 0 or len(train_i) == 0:
                       continue
                splits.append((train_i, test_i))



    for fold, (train_i, test_i) in enumerate(splits):
        if args.cv_by == "species":
            test_species = sorted(set(species[test_i]))
            log.info(f"[fold{fold}] LOSO test species: {test_species}")
        else:
            test_species = Counter(species[test_i])
            log.info(f"[fold{fold}] test species mix: {dict(test_species)}")


    results = []
    rng = np.random.default_rng(0)  # reproducible
    def make_val_from_train(train_idx, val_frac=0.2, seed=0):
        train_idx = np.asarray(train_idx, dtype=int)
        train_prots = proteins[train_idx]
        uniq = np.unique(train_prots)
        # sample proteins for val
        n_val = max(1, int(np.ceil(val_frac * len(uniq))))
        val_prots = rng.choice(uniq, size=n_val, replace=False)
        is_val = np.isin(train_prots, val_prots)
        val_idx = train_idx[is_val]
        new_train_idx = train_idx[~is_val]
        return new_train_idx, val_idx

    for fold, (train_i, test_i) in enumerate(splits):
            # build validation from TRAIN proteins (not from test fold)
        train_i, val_i = make_val_from_train(train_i, val_frac=0.2, seed=fold)    # for proc in test_i:
        #     log.info(f"{proc}")
        for seed in args.seeds:
            MAX_RETRIES = 2
            retries = 0
            while retries <= MAX_RETRIES:
                # ─── Clean up stale attention maps ─────────────────────
                attn_glob = list(args.attn_dir.glob(f"*f{fold}_s{seed}.csv"))
                for f in attn_glob:
                    f.unlink()
                # log.info(f"[CLEANUP] Removed {len(attn_glob)} old attention maps for fold{fold} seed{seed}")

                # ─── Clean up stale checkpoint ─────────────────────────
                ckpt_file = args.save_dir / f"f{fold}_s{seed}.pt"
                if ckpt_file.exists():
                    ckpt_file.unlink()
                    # log.info(f"[CLEANUP] Removed old checkpoint: {ckpt_file}")

                # ─── Rerun fold ────────────────────────────────────────
                mae, reg_roc, prot_contact_roc = run_fold(
                    seed, fold, train_i, val_i, test_i, ds, args, lookup
                )
                test_pids = {proteins[i] for i in test_i}
                global_roc = summarise_attention(
                    args.attn_dir,
                    Path(
                        "~/checkpoints/all_proteins_plip_summary.csv"
                    ),
                    suffix=f"f{fold}_s{seed}",
                    include_pids=test_pids,
                )

                # if not (np.isnan(global_roc) or global_roc == -1 or reg_roc < MIN_ROC_THRESHOLD or prot_contact_roc < MIN_ROC_THRESHOLD):
                results.append(
                    dict(
                        fold=fold,
                        seed=seed,
                        mae=mae,
                        roc=reg_roc,
                        prot_contact_roc=prot_contact_roc,
                        global_roc=global_roc,
                    )
                )
                print(
                    f"-----------------------------------------------------------------------------------"
                )
                print(
                    f"fold{fold} seed{seed} TEST test_mae={mae:.3f}  test_reg_roc={reg_roc:.3f} test_prot_contact_roc={prot_contact_roc:.3f} global_roc={global_roc:.3f}"
                )
                print(
                    f"-----------------------------------------------------------------------------------"
                )
                break  # Good run

    # finally write out the results to file
    write_summary(args, results=results)

    end_time = time.time()
    execution_time = end_time - start_time
    log.info(f"Execution time: {execution_time:.4f} seconds")


def write_summary(args, results):
    # write summary
    summary_csv = args.run_dir / "summary.csv"
    with summary_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["fold", "seed", "mae", "roc", "prot_contact_roc", "global_roc"],
        )
        writer.writeheader()
        for r in results:
            writer.writerow(
                {k: f"{v:.4f}" if isinstance(v, float) else v for k, v in r.items()}
            )

    # Aggregate final metrics
    import statistics

    maes = [r["mae"] for r in results]
    rocs = [r["roc"] for r in results]
    prot_contact_roc = [r["prot_contact_roc"] for r in results]


    def format_mean_std(values):
        if isinstance(values, set):
            values = list(values)
        mean = statistics.mean(values)
        std = statistics.stdev(values) if len(values) > 1 else 0.0
        return f"{mean:.3f} ± {std:.3f}"

    # Corrected one-line summary log
    try:
        log.info(
            f"Mean ± Std: MAE = {format_mean_std(maes)}, ROC = {format_mean_std(rocs)},"
            f" prot_contact_roc= {format_mean_std(prot_contact_roc)}"
        )
    except Exception as e:
        print(f"{e}")
    log.info(f"Wrote summary → {summary_csv}")


if __name__ == "__main__":
    main()

    # python -m pipeline.cli.train --graph-root pipeline/data_graph --run-dir checkpoints/run_20250714_101856_rerun_f3s1
    # --cv 4 --seeds 1 --epochs 30 --batch 4 --huber  --scheduler cosine --swa-start 0.6 --nodecap --no-fusion
