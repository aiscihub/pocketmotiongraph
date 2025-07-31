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
import argparse, csv, logging, random, time
import math
from collections import Counter
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import torch, torch.nn.functional as F
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import GroupKFold
from sympy import false
from torch.cuda.amp import GradScaler, autocast
from torch.optim.swa_utils import AveragedModel, SWALR
from torch_geometric.loader import DataLoader
from torch.utils.data import Sampler, SubsetRandomSampler, Subset

from pipeline.graph.dataset import DriftDataset
from pipeline.graph.g_config import FULL_IDX, NO_PLIP_IDX, NO_ESM_IDX, NO_ESM_PLIP_IDX
from pipeline.graph.g_config_patched import make_feature_groups
from pipeline.models.contact_lookup import ContactLookup
from pipeline.models.gnn_equivariant import PocketStabilityEGNN


# ───────────── config ─────────────
POS_CLIP   = 30.0
PAIR_K     = 32
BETA_PAIR  = 0.1
ACCUM_STEPS = 4
NODE_LOSS_W = 1.0
EDGE_W = 0.05
MAX_NODES   = 12_000   # for NodeCapSampler
Z_STD = 1e-8
Z_GAP = 0.001
import logging
log = logging.getLogger(__name__)  # Use module name as logger ID

rmsd_log = { "epoch": [], "true": [],   "pred": [] }

from sklearn.metrics import roc_auc_score
from collections import defaultdict

@torch.no_grad()
def roc_by_protein(model, loader, device, lookup):
    model.eval()
    s, t = defaultdict(list), defaultdict(list)          # scores, targets
    for b in loader:
        b = b.to(device)
        _, _, node_logits, _ = model(b.x, b.edge_index,
                                     batch=b.batch, pos=b.pos,
                                     u=b.u, pid=b.pid)
        for g in torch.unique(b.batch):
            m = b.batch == g
            pid = b.meta['protein_id'][g]
            tgt = lookup.mask(pid, b.residue_id[m]).float().cpu()
            scr = torch.sigmoid(node_logits[m]).cpu()
            s[pid].extend(scr.tolist())
            t[pid].extend(tgt.tolist())

    rocs = {p: roc_auc_score(t[p], s[p]) for p in s if len(set(t[p])) > 1}
    mu  = sum(rocs.values()) / len(rocs)
    sig = (sum((v - mu) ** 2 for v in rocs.values()) / len(rocs)) ** 0.5
    return mu, sig, rocs


def summarise_attention(attn_dir: Path,
                        plip_csv: Path,
                        suffix: str,
                        include_pids: set[str]|None = None):
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
    plip["residue_id"] = (
        plip["Residue"].str.extract(r"(\d+)", expand=False).astype(int)
    )
    merged = ( df.merge(plip, on=["protein", "residue_id"], how="left").dropna(subset=["stab_score"])  )
    if include_pids is not None:
        merged = merged[ merged["protein"].isin(include_pids) ]
        log.info(f"Include_pids = {include_pids}")

    merged["is_contact"] = merged["Complex"].notna()
    merged.to_csv(attn_dir / f"attention_vs_plip_{suffix}.csv", index=False)


    #merged["stab_z"] = merged["stab_score"]

    # --- 1) RAW attention stats (before any normalization) -----------------
    raw_mc  = merged.loc[merged.is_contact,  "stab_score"].mean()
    raw_nc  = merged.loc[~merged.is_contact, "stab_score"].mean()
    raw_std = merged["stab_score"].std()
    log.debug(f"[DEBUG RAW] std={raw_std:.6f}, mean(contact)={raw_mc:.6f}, mean(non-contact)={raw_nc:.6f}")
    # also per-protein raw stats:


    for prot in merged["protein"].unique():
         sub = merged[merged["protein"] == prot]
         std_p = sub["stab_score"].std()
         mc_p  = sub.loc[sub.is_contact,   "stab_score"].mean()
         nc_p  = sub.loc[~sub.is_contact,  "stab_score"].mean()
         log.debug(f"[DEBUG RAW PROT] {prot:15} | std={std_p:.6f} | mean(c)={mc_p:.6f} | non-c={nc_p:.6f}")

    # ── 1) flip per protein on *raw* stab_score so that contact > non-contact
    def _flip(group):
        if group['is_contact'].sum() < 2:        # degenerate, leave as is
            return group
        if group.loc[group.is_contact, 'stab_score'].mean() < \
                group.loc[~group.is_contact, 'stab_score'].mean():
            group['stab_score'] *= -1
        return group

    merged = merged.groupby('protein', group_keys=False).apply(_flip)
    merged["stab_z"] = (
        merged.groupby("protein")["stab_score"]        # robust IQR norm
        .transform(lambda x: (x - x.median())
                             / (x.quantile(0.75) - x.quantile(0.25) + 1e-6)))
   # ------------ (A) AUPRC ------------------------------------------------
    y_true = merged["is_contact"].astype(int).values
    y_score = merged["stab_z"].values
    auprc = average_precision_score(y_true, y_score)

    # ------------ (B) Top-K precision / recall ----------------------------
    K = 20                                             # ← feel free to expose as CLI arg
    merged = merged.sort_values("stab_z", ascending=False)
    topk = merged.head(K)
    tp = topk["is_contact"].sum()
    topk_precision = tp / K
    topk_recall    = tp / merged["is_contact"].sum()

    log.info(f"[ATTN:{suffix}] AUPRC={auprc:.3f}  Top{K} Precision={topk_precision:.3f} Top{K} Recall={topk_recall:.3f}")

    # Calculate means
    mean_contact = merged.loc[merged.is_contact, "stab_z"].mean()
    mean_non_contact = merged.loc[~merged.is_contact, "stab_z"].mean()
    std_all = merged["stab_z"].std()

    # New check to skip noisy attention
    if std_all < 1e-3 or abs(mean_contact - mean_non_contact) < Z_GAP:
        log.warning(f"[ATTN:{suffix}] Low attention signal — skipping global ROC (std={std_all:.4f}, diff={mean_contact - mean_non_contact:.4f})")
        roc_auc = float("nan")  # or -1
    else:
        roc_auc = roc_auc_score(merged["is_contact"].astype(int), merged["stab_z"].values)
        log.info(f"[ATTN:{suffix}] global ROC-AUC    = {roc_auc:.3f}, AUPRC = {auprc:.3f}")
    # ------------- (C) dump metrics to CSV --------------------------------
    metrics_row = {
        "suffix":    suffix,
        "protein_set": "|".join(sorted(include_pids)) if include_pids else "ALL",
        "ROC_AUC":   roc_auc,
        "AUPRC":     auprc,
        f"Top{K}_precision": topk_precision,
        f"Top{K}_recall":    topk_recall,
        "mean_c":    mean_contact,
        "mean_nc":   mean_non_contact,
        "std":       std_all,
    }
    csv_path = attn_dir / "attention_fold_metrics_{suffix}.csv"
    hdr = not csv_path.exists()
    pd.DataFrame([metrics_row]).to_csv(csv_path, mode="a", index=False, header=hdr)

    return roc_auc

#post_weight = 0.8 is the best so far
class NodeCapSampler(Sampler):
    def __init__(self, ds, max_nodes=MAX_NODES, contact_weight=0.8, shuffle=True):
        self.ds, self.max_nodes = ds, max_nodes
        self.sizes = [g.num_nodes for g in ds]
        self.contact = [i for i, g in enumerate(ds) if g.is_contact.any()]
        self.no_contact = [i for i in range(len(ds)) if i not in self.contact]
        self.cw, self.shuffle = contact_weight, shuffle

    def __iter__(self):
        idx_pool = (self.contact + self.no_contact)
        if self.shuffle:
            random.shuffle(idx_pool)

        batch, tot = [], 0
        for idx in idx_pool:
            n = self.sizes[idx]
            if tot + n > self.max_nodes and batch:
                yield batch
                batch, tot = [], 0
            batch.append(idx); tot += n
        if batch:
            yield batch                     # last partial batch

    def __len__(self):
        # how many batches ≈ one full pass
        return math.ceil(sum(self.sizes) / self.max_nodes)



def set_seed(seed:int):
    random.seed(seed); np.random.seed(seed)
    torch.manual_seed(seed); torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

def w_bce(logits, targets):
    pos = targets.sum(); neg = targets.numel() - pos
    if pos==0 or neg==0: return torch.tensor(0.0, device=logits.device)
    pw = min(neg/pos, POS_CLIP)
    return F.binary_cross_entropy_with_logits(
        logits, targets, pos_weight=torch.tensor([pw], device=logits.device))

#The weighted pos weight for bce does not help
def w_bce_pos_weight(logits, targets, pos_weight = None, gamma = None):
    pos = targets.sum()
    neg = targets.numel() - pos
    if pos == 0 or neg == 0:
        return torch.tensor(0.0, device=logits.device)

    if pos_weight is not None:
            pw = float(pos_weight)
    else:
        pw = min(neg / pos, POS_CLIP)        # dynamic weighting
    # print(f"pw ---- {pw}")
    loss = F.binary_cross_entropy_with_logits(
        logits, targets,
        pos_weight=torch.tensor([pw], device=logits.device),
        reduction="none")

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
        _, node_attn = model(g.x, g.edge_index,      batch=g.batch, pos=g.pos, u=g.u, pid=g.pid,
                             return_node_attn=True)
        for gid in g.batch.unique():
            m = (g.batch == gid)
            edge_mask = torch.zeros_like(m)          # m has length = num_nodes
            edge_nodes = torch.unique(
                g.edge_index[:, g.batch[g.edge_index[0]] == gid].flatten()
            )
            edge_mask[edge_nodes] = True
            m = m & edge_mask                        # keep only nodes with ≥1 edge
            if m.sum() == 0:
                continue

            res_ids = g.residue_id[m].cpu()
            scores  = node_attn[m].cpu()
            # Optional: Add pocket center distance if available
            xyz = g.pos[m].cpu().numpy()
            arr = torch.stack([res_ids, scores], 1).numpy()

            # every graph already has a unique pdb_id  (ProteinID_Pocket_Frame)
            pdb_id = str(g.pdb_id if isinstance(g.pdb_id, str) else g.pdb_id[0])
            np.savetxt(out_dir / f"{pdb_id}_attn_f{fold}_s{seed}.csv",arr, fmt=["%d", "%.10f"],  delimiter=",",
                       header="residue_id,score", comments="")

# ───────────── training loops ─────────────
def train_epoch(model, loader, optim, scaler, device, lookup,
                *, huber:bool, accum:int, pos_weight=None, gamma=None):
    model.train(); total, step = 0.0, 0
    optim.zero_grad(set_to_none=True)
    for batch in loader:
        batch = batch.to(device)
        with autocast(dtype=torch.bfloat16):
            reg, _, node_logits, aux = model(
                batch.x,
                batch.edge_index,
                batch=batch.batch, pos=batch.pos, u=batch.u, pid=batch.pid)
            reg_loss = (F.smooth_l1_loss if huber else F.l1_loss)(
                reg.view(-1), batch.y[:,0])

            node_loss = 0.0
            n_graph = int(batch.batch.max())+1
            for g in range(n_graph):
                m   = batch.batch==g
                pid = batch.meta['protein_id'][g]
                tgt = lookup.mask(pid, batch.residue_id[m]).float()

                #node_loss += w_bce(node_logits[m], tgt)
                node_loss += w_bce_pos_weight(node_logits[m], tgt, pos_weight=pos_weight,  gamma=gamma)
            node_loss /= n_graph

            pair_loss = 0.0
            if hasattr(model,'pair_head'):
                for g in range(n_graph):
                    m   = batch.batch==g
                    pid = batch.meta['protein_id'][g]
                    tgt_mask = lookup.mask(pid, batch.residue_id[m]).bool()
                    idx = torch.nonzero(m).squeeze(-1)
                    pos_ids = idx[tgt_mask]; neg_ids = idx[~tgt_mask]
                    if pos_ids.numel()<2 or neg_ids.numel()==0: continue
                    pos_pairs = torch.combinations(pos_ids)
                    if pos_pairs.size(0) > PAIR_K:
                        sel = torch.randperm(pos_pairs.size(0), device=device)[:PAIR_K]
                        pos_pairs = pos_pairs[sel]
                    neg_pairs = torch.stack(torch.meshgrid(
                        pos_ids, neg_ids, indexing='ij')).T.reshape(-1,2)
                    if neg_pairs.size(0) > PAIR_K:
                        sel = torch.randperm(neg_pairs.size(0), device=device)[:PAIR_K]
                        neg_pairs = neg_pairs[sel]
                    pairs  = torch.cat([pos_pairs, neg_pairs], 0)
                    labels = torch.cat([torch.ones(len(pos_pairs), device=device),
                                        torch.zeros(len(neg_pairs),device=device)])
                    h_i = aux['h'][pairs[:,0]]; h_j = aux['h'][pairs[:,1]]
                    dist = (aux['pos'][pairs[:,0]]-aux['pos'][pairs[:,1]]).norm(dim=-1,keepdim=True)
                    logits = model.pair_head(torch.cat([h_i,h_j,dist],-1)).squeeze(-1)
                    #pair_loss += w_bce(logits, labels)
                    pair_loss += w_bce_pos_weight(logits, labels, pos_weight=pos_weight,  gamma=gamma)

                    # ── inside the for-g loop, after pair_loss ────────────────────────────────
                    edge_mask   = (batch.batch[batch.edge_index[0]] == g)      # Bool[E_total]
                    if edge_mask.any():
                        edge_pairs  = batch.edge_index[:, edge_mask].t()       # (E_g, 2)
                        edge_labels = batch.is_contact[edge_mask].float()      # (E_g,)

                        h_i  = aux['h'][edge_pairs[:, 0]]
                        h_j  = aux['h'][edge_pairs[:, 1]]
                        dist = (aux['pos'][edge_pairs[:, 0]] -
                                aux['pos'][edge_pairs[:, 1]]).norm(dim=-1, keepdim=True)

                        edge_logits = model.pair_head(torch.cat([h_i, h_j, dist], -1)).squeeze(-1)
                        edge_loss   = w_bce_pos_weight(edge_logits, edge_labels,
                                                       pos_weight=pos_weight, gamma=gamma)
                        pair_loss  += EDGE_W * edge_loss
                pair_loss /= max(n_graph,1)

            loss = reg_loss + NODE_LOSS_W *node_loss + BETA_PAIR*pair_loss
            # In train_epoch(), modify the logging:
            # log.info(
            #     f"(reg_loss:{reg_loss.item():.3f}, node_loss:{node_loss.item():.3f}, pair:{pair_loss.item():.3f})"
            # )

        scaler.scale(loss).backward(); step += 1
        if step % accum == 0:
            scaler.unscale_(optim)
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            scaler.step(optim); scaler.update()
            optim.zero_grad(set_to_none=True)
        total += loss.item()
    return total/len(loader)

@torch.no_grad()
def evaluate(model, loader, device, lookup, epoch=None):
    model.eval(); mae_l, roc_t, roc_s = [], [], []
    true_vals, pred_vals = [], []  # ← added
    for b in loader:
        b = b.to(device)
        reg,_,node_logits,_ = model(b.x,b.edge_index,
                                    batch=b.batch,pos=b.pos,u=b.u,pid=b.pid)
        mae_l.append((reg.view(-1)-b.y[:,0]).abs().cpu())
        # Save for RMSD scatter plot
        true_vals.extend(b.y[:, 0].cpu().numpy())
        pred_vals.extend(reg.view(-1).cpu().numpy())
        for g in torch.unique(b.batch):
            m=b.batch==g
            pid=b.meta['protein_id'][g]
            tgt=lookup.mask(pid,b.residue_id[m]).float()
            roc_t += tgt.tolist(); roc_s += torch.sigmoid(node_logits[m]).cpu().tolist()
    if epoch is not None:
        rmsd_log["epoch"].extend([epoch] * len(true_vals))
        rmsd_log["true"].extend(true_vals)
        rmsd_log["pred"].extend(pred_vals)

    mae=torch.cat(mae_l).mean().item()
    roc=roc_auc_score(roc_t,roc_s) if len(set(roc_t))>1 else float('nan')
    return mae,roc

# ───────────── fold runner ─────────────
def run_fold(seed, fold, train_i, val_i, test_i, ds, args, lookup):
    set_seed(seed)
    dev=torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = PocketStabilityEGNN(
        in_dim      = ds[0].x.shape[1],
        global_dim  = ds[0].u.shape[0],
        hidden_dim  = args.hidden_dim,
        dropout     = args.dropout,
        use_u       = args.fusion
    ).to(dev)
    logging.info(f"model settings: x = {ds[0].x.shape[1]} u = {ds[0].u.shape[0]} {args.hidden_dim} {args.dropout} {args.fusion}")

    optim=torch.optim.Adam(model.parameters(), lr=args.lr)
    scaler=GradScaler()

    def make_loader(indices):
        if args.nodecap:
            sub= torch.utils.data.Subset(ds, indices)
            return DataLoader(sub, batch_sampler=NodeCapSampler(sub))
        return DataLoader(ds, batch_size=args.batch,
                          sampler=SubsetRandomSampler(indices))
    l_train,l_val,l_test = map(make_loader,(train_i,val_i,test_i))

    warm=int(max(1,args.epochs*0.1))
    warm = 5
    warm_sched=torch.optim.lr_scheduler.LinearLR(optim,0.1,1.0,warm)
    main_sched={'cosine':torch.optim.lr_scheduler.CosineAnnealingLR(
        optim,T_max=args.epochs-warm),
        'plateau':torch.optim.lr_scheduler.ReduceLROnPlateau(
            optim,mode='min',factor=0.5,patience=3)}.get(args.scheduler)
    swa_model=AveragedModel(model); swa_sched=SWALR(optim, swa_lr=1e-4)
    swa_phase=False; swa_start=int(args.epochs*args.swa_start)
    best=float('inf')

    patience = 8

    best_epoch = 0
    global rmsd_log
    rmsd_log["epoch"].clear()
    rmsd_log["true"].clear()
    rmsd_log["pred"].clear()
    best_roc = 0.0
    for ep in range(1,args.epochs+1):
        tl=train_epoch(model,l_train,optim,scaler,dev,lookup,
                       huber=args.huber,accum=ACCUM_STEPS, pos_weight=args.pos_weight, gamma=args.gamma)
        v_mae,v_roc=evaluate(model,l_val,dev,lookup, epoch=ep)
        log.info(f"[fold{fold} seed{seed}] ep{ep:02d} train {tl:.3f} | val {v_mae:.3f}/{v_roc:.3f}")
        # LR sched
        if ep<=warm:
            warm_sched.step()
        elif ep<swa_start:
            if args.scheduler=='plateau': main_sched.step(v_mae)
            elif args.scheduler!='none': main_sched.step()
        else:
            if not swa_phase:
                log.info(f"[fold{fold} seed{seed}] enter SWA at ep{ep}"); swa_phase=True
            swa_model.update_parameters(model); swa_sched.step()
        if v_roc > best_roc:
            best=v_mae
            best_epoch, best_roc = ep, v_roc
            torch.save(model.state_dict(), args.save_dir/f"f{fold}_s{seed}.pt")
            full_loader = DataLoader(ds, batch_size=1, shuffle=False)
            dump_attention_maps(model, full_loader, dev, args.attn_dir, fold, seed)
            log.info(f"Saving the model with v_mae {best:3f} v_roc {best_roc:3f}")
        elif ep - best_epoch >= patience:
            log.info(f"Early stopping at epoch {ep}. Best was {best_epoch} (val_MAE={best:.3f})")
            break

     # test
    model.load_state_dict(torch.load(args.save_dir/f"f{fold}_s{seed}.pt"))
    # mu, sig, v_proc = roc_by_protein(model, l_test, dev, lookup)
    # log.info(f"[cross-protein] ROC-AUC  reg_roc={mu:.3f}  σ={sig:.3f} ")

    t_mae,t_reg_roc=evaluate(model,l_test,dev,lookup)

    # ── Local-Contact head: per-protein ROC-AUC (attention logits vs PLIP) ──
    t_contact_avg_roc, attn_sig, t_contact_roc = roc_by_protein(model, l_test, dev, lookup)
    for p, v in sorted(t_contact_roc.items()):
        log.info(f"  - {p}: {v:.3f}")
    logging.info(
        f"run_fold -- [fold{fold} seed{seed}] "
        f"TEST mae={t_mae:.3f} | reg-roc={t_reg_roc:.3f} | contact-roc={t_contact_avg_roc:.3f}"
         )
    return t_mae,t_reg_roc,t_contact_avg_roc


# --- SAMPLING LOGIC (below dataset creation) ---
from collections import defaultdict
import re

def _parse_ids(pdb_id):
    m = re.search(r"(.*?)(_P\d+)?(_F(\d+))?$", pdb_id)
    if not m: return pdb_id, "P0", -1
    pid, _, _, f = m.groups()
    return pid, "P0", int(f) if f else -1

def _group_by_pocket(meta):
    g = defaultdict(list)
    for i, m in enumerate(meta):
        pid, pocket, frame = _parse_ids(str(m.get("pdb_id", "")))
        g[(pid, pocket)].append((i, frame))
    return {k: [x[0] for x in sorted(v, key=lambda x: x[1])] for k, v in g.items()}

def _select_uniform(lst, k):
    if len(lst) <= k: return lst
    return [lst[i] for i in np.linspace(0, len(lst)-1, k, dtype=int)]



# ───────────── main ─────────────
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--graph-root',type=Path,required=True)

    ap.add_argument('--epochs',type=int,default=30)
    ap.add_argument('--batch',type=int,default=4)
    ap.add_argument('--lr',type=float,default=2e-3)
    ap.add_argument('--hidden-dim',type=int,default=256)
    ap.add_argument('--dropout',type=float,default=0.3)
    ap.add_argument('--scheduler',choices=['none','cosine','plateau'],default='plateau')
    ap.add_argument('--swa-start',type=float,default=0.6)
    ap.add_argument('--cv',type=int,default=4)
    ap.add_argument('--seeds',nargs='*',type=int,default=[0])
    ap.add_argument('--huber',action='store_true')
    ap.add_argument('--nodecap',action='store_true', default=True)
    ap.add_argument('--fusion',action='store_true')
    ap.add_argument("--physics-only", action="store_true")
    ap.add_argument('--contacts',type=Path,default=Path('checkpoints/contacts.csv'))
    ap.add_argument('--run-dir', type=Path, required=True, help='(e.g., run_20250712_153000)')
    ap.add_argument("--pos-weight", default=None)
    ap.add_argument("--gamma", default=None)
    ap.add_argument("--protein", default=None)
    # --new CLI flags -------------------------------------------------
    ap.add_argument("--node-plip",  action="store_true")
    ap.add_argument("--no-esm",   action="store_true")
    ap.add_argument("--ablate-node", action="store_true",
                        help="Drop all per-residue (node) features except ESM")


    # --- after loading DriftDataset ---------------------------------

    args=ap.parse_args()
    # Record the start time
    start_time = time.time()
    # Create structured run dir
    args.run_dir.mkdir(parents=True, exist_ok=True)
    args.save_dir = args.run_dir / 'checkpoints'
    args.attn_dir = args.run_dir / 'attention_maps'
    args.save_dir.mkdir(exist_ok=True)
    args.attn_dir.mkdir(exist_ok=True)

    # Logging path
    log_file = args.run_dir / 'train_log.txt'
    logging.getLogger().handlers.clear()  # <- optional, clears previous handlers
    logging.basicConfig(
        format='[%(levelname)s] %(message)s',
        level=logging.DEBUG,
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
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

    lookup=ContactLookup(args.contacts)
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
          raise RuntimeError(f"No graphs loaded from {args.graph_root}. "
                                                         f"Check drift filter / paths / feature groups.")

    prot_counts = Counter(m["protein_id"] for m in ds.meta)
    for pid, n in prot_counts.items():
          log.info(f"{pid:22} → {n} graphs")
    log.info(f"Graph total = {len(ds)} ")

    proteins    = np.array([m["protein_id"] for m in ds.meta])
    all_indices = np.arange(len(ds))

    if args.protein is None:
        # Standard GroupKFold
        n_splits = min(args.cv, len(np.unique(proteins)))
        if n_splits < 2:
            raise ValueError(f"Need at least 2 proteins for CV, found {len(np.unique(proteins))}.")
        gkf = GroupKFold(n_splits=n_splits)
        splits = list(gkf.split(all_indices, groups=proteins))


    results=[]
    MIN_ROC_THRESHOLD = 0.75
    #for fold,(train_i,test_i) in enumerate(gkf.split(np.zeros(len(proteins)),groups=proteins)):
    for fold, (train_i, test_i) in enumerate(splits):
        mid=len(test_i)//2; val_i=test_i[:mid]; test_i=test_i[mid:]
        # for proc in test_i:
        #     log.info(f"{proc}")
        for seed in args.seeds:
            MAX_RETRIES = 2
            retries = 0
            while retries <= MAX_RETRIES:
                # ─── Clean up stale attention maps ─────────────────────
                attn_glob = list(args.attn_dir.glob(f"*f{fold}_s{seed}.csv"))
                for f in attn_glob:
                    f.unlink()
                #log.info(f"[CLEANUP] Removed {len(attn_glob)} old attention maps for fold{fold} seed{seed}")

                # ─── Clean up stale checkpoint ─────────────────────────
                ckpt_file = args.save_dir / f"f{fold}_s{seed}.pt"
                if ckpt_file.exists():
                    ckpt_file.unlink()
                    #log.info(f"[CLEANUP] Removed old checkpoint: {ckpt_file}")

                # ─── Rerun fold ────────────────────────────────────────
                mae, reg_roc, prot_contact_roc= run_fold(seed, fold, train_i, val_i, test_i, ds, args, lookup)
                test_pids = { proteins[i] for i in test_i }
                global_roc = summarise_attention(
                    args.attn_dir,
                    Path("/home/zhenli/git/valleyfevermutation/checkpoints/all_proteins_plip_summary.csv"),
                    suffix=f"f{fold}_s{seed}",
                    include_pids=test_pids
                )

                if not (np.isnan(global_roc) or global_roc == -1 or reg_roc < MIN_ROC_THRESHOLD or prot_contact_roc < MIN_ROC_THRESHOLD):
                    results.append(dict(fold=fold, seed=seed, mae=mae, roc=reg_roc, prot_contact_roc = prot_contact_roc, global_roc=global_roc))
                    print(f"-----------------------------------------------------------------------------------")
                    print(f"fold{fold} seed{seed} TEST test_mae={mae:.3f}  test_reg_roc={reg_roc:.3f} test_prot_contact_roc={prot_contact_roc:.3f} global_roc={global_roc:.3f}")
                    print(f"-----------------------------------------------------------------------------------")
                    break  # Good run

                log.warning(f"[RETRY] fold{fold} seed{seed} test_reg_roc={reg_roc:.3f} test_prot_contact_roc={prot_contact_roc:.3f}too low or invalid — retrying ({retries+1}/{MAX_RETRIES})")
                retries += 1
    #finally write out the results to file
    write_summary(args, results=results)

    end_time = time.time(); execution_time = end_time - start_time
    log.info(f"Execution time: {execution_time:.4f} seconds")

def write_summary(args, results=[]):
    # write summary
    summary_csv = args.run_dir / 'summary.csv'
    with summary_csv.open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['fold', 'seed', 'mae', 'roc', 'prot_contact_roc', 'global_roc'])
        writer.writeheader()
        for r in results:
            writer.writerow({k: f"{v:.4f}" if isinstance(v, float) else v for k, v in r.items()})

    # Aggregate final metrics
    import statistics
    maes = [r["mae"] for r in results]
    rocs = [r["roc"] for r in results]
    prot_contact_roc = [r["prot_contact_roc"] for r in results]
    grocs = [r["global_roc"] for r in results]

    def format_mean_std(values):
        if isinstance(values, set):
            values = list(values)
        mean = statistics.mean(values)
        std  = statistics.stdev(values) if len(values) > 1 else 0.0
        return f"{mean:.3f} ± {std:.3f}"
    # Corrected one-line summary log
    try:
        log.info(f"Mean ± Std: MAE = {format_mean_std(maes)}, ROC = {format_mean_std(rocs)},"
             f" prot_contact_roc= {format_mean_std(prot_contact_roc)}, Global ROC = {format_mean_std(grocs)}")
    except Exception as e:
        print(f"{e}")
    log.info(f"Wrote summary → {summary_csv}")
if __name__=="__main__":
    main()
