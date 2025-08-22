#!/usr/bin/env python
"""
Portable trainer: cross-protein CV with advanced tricks, no local paths.
"""
from __future__ import annotations
import argparse, csv, logging, random, time, math, re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import GroupKFold
from torch.cuda.amp import GradScaler, autocast
from torch.optim.swa_utils import AveragedModel, SWALR
from torch_geometric.loader import DataLoader
from torch.utils.data import Sampler, SubsetRandomSampler, Subset

# ——— Your internal modules (kept; no absolute dirs) ———
from pipeline.graph.dataset import DriftDataset
from pipeline.graph.g_config import FULL_IDX, NO_PLIP_IDX, NO_ESM_IDX, NO_ESM_PLIP_IDX
from pipeline.graph.g_config_patched import make_feature_groups
from pipeline.models import PocketStabilityModel
from pipeline.models.contact_lookup import ContactLookup

# ───────────── config ─────────────
POS_CLIP      = 30.0
PAIR_K        = 32
BETA_PAIR     = 0.1
ACCUM_STEPS   = 4
NODE_LOSS_W   = 1.0
EDGE_W        = 0.05
MAX_NODES     = 12_000
Z_GAP         = 0.001

log = logging.getLogger(__name__)
rmsd_log = {"epoch": [], "true": [], "pred": []}

# ───────────── utils ─────────────
def set_seed(seed:int):
    random.seed(seed); np.random.seed(seed)
    torch.manual_seed(seed); torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

def w_bce_pos_weight(logits, targets, pos_weight=None, gamma=None):
    pos = targets.sum()
    neg = targets.numel() - pos
    if pos == 0 or neg == 0:
        return torch.tensor(0.0, device=logits.device)
    pw = float(pos_weight) if pos_weight is not None else min(neg / pos, POS_CLIP)
    loss = F.binary_cross_entropy_with_logits(
        logits, targets, pos_weight=torch.tensor([pw], device=logits.device), reduction="none"
    )
    if gamma is not None:
        p = torch.sigmoid(logits)
        loss = (1 - p) ** gamma * loss
    return loss.mean()

class NodeCapSampler(Sampler):
    def __init__(self, ds, max_nodes=MAX_NODES, shuffle=True):
        self.ds, self.max_nodes, self.shuffle = ds, max_nodes, shuffle
        self.sizes = [g.num_nodes for g in ds]

    def __iter__(self):
        idx_pool = list(range(len(self.ds)))
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
            yield batch

    def __len__(self):
        return math.ceil(sum(self.sizes) / self.max_nodes)

@torch.no_grad()
def evaluate(model, loader, device, lookup, epoch=None):
    model.eval(); mae_l, roc_t, roc_s = [], [], []
    true_vals, pred_vals = [], []
    for b in loader:
        b = b.to(device)
        reg, _, node_logits, _ = model(b.x, b.edge_index, batch=b.batch, pos=b.pos, u=b.u, pid=b.pid)
        mae_l.append((reg.view(-1) - b.y[:, 0]).abs().cpu())
        true_vals.extend(b.y[:, 0].cpu().numpy())
        pred_vals.extend(reg.view(-1).cpu().numpy())
        for g in torch.unique(b.batch):
            m   = b.batch == g
            pid = b.meta['protein_id'][g]
            tgt = lookup.mask(pid, b.residue_id[m]).float()
            roc_t += tgt.tolist()
            roc_s += torch.sigmoid(node_logits[m]).cpu().tolist()
    if epoch is not None:
        rmsd_log["epoch"].extend([epoch] * len(true_vals))
        rmsd_log["true"].extend(true_vals)
        rmsd_log["pred"].extend(pred_vals)
    mae = torch.cat(mae_l).mean().item()
    roc = roc_auc_score(roc_t, roc_s) if len(set(roc_t)) > 1 else float('nan')
    return mae, roc

@torch.no_grad()
def roc_by_protein(model, loader, device, lookup):
    model.eval()
    s, t = defaultdict(list), defaultdict(list)
    for b in loader:
        b = b.to(device)
        _, _, node_logits, _ = model(b.x, b.edge_index, batch=b.batch, pos=b.pos, u=b.u, pid=b.pid)
        for g in torch.unique(b.batch):
            m = b.batch == g
            pid = b.meta['protein_id'][g]
            tgt = lookup.mask(pid, b.residue_id[m]).float().cpu()
            scr = torch.sigmoid(node_logits[m]).cpu()
            s[pid].extend(scr.tolist()); t[pid].extend(tgt.tolist())
    rocs = {p: roc_auc_score(t[p], s[p]) for p in s if len(set(t[p])) > 1}
    mu  = sum(rocs.values()) / len(rocs) if rocs else float("nan")
    sig = (sum((v - mu) ** 2 for v in rocs.values()) / len(rocs)) ** 0.5 if rocs else float("nan")
    return mu, sig, rocs

def summarise_attention(attn_dir: Path, plip_csv: Path|None, suffix: str, include_pids: set[str]|None=None):
    """If plip_csv is None or missing, skip gracefully."""
    csvs = sorted(attn_dir.glob(f"*_{suffix}.csv"))
    if not csvs or not plip_csv or not plip_csv.exists():
        log.warning(f"[ATTN] Skipping attention summary (missing CSVs or PLIP csv).")
        return float("nan")
    dfs = []
    for c in csvs:
        df = pd.read_csv(c)
        protein = c.name.split("_attn_")[0]
        if "protein" not in df.columns:
            df.insert(0, "protein", protein)
        if "score" in df.columns and "stab_score" not in df.columns:
            df.rename(columns={"score": "stab_score"}, inplace=True)
        dfs.append(df)
    df = pd.concat(dfs, ignore_index=True, sort=False)
    plip = pd.read_csv(plip_csv)
    plip["residue_id"] = plip["Residue"].str.extract(r"(\d+)", expand=False).astype(int)
    merged = df.merge(plip, on=["protein", "residue_id"], how="left").dropna(subset=["stab_score"])
    if include_pids is not None:
        merged = merged[merged["protein"].isin(include_pids)]
    merged["is_contact"] = merged["Complex"].notna()

    # flip per protein on raw score so contact > non-contact
    def _flip(g):
        if g["is_contact"].sum() < 2: return g
        if g.loc[g.is_contact, "stab_score"].mean() < g.loc[~g.is_contact, "stab_score"].mean():
            g["stab_score"] *= -1
        return g
    merged = merged.groupby("protein", group_keys=False).apply(_flip)
    merged["stab_z"] = merged.groupby("protein")["stab_score"] \
        .transform(lambda x: (x - x.median()) / (x.quantile(0.75) - x.quantile(0.25) + 1e-6))

    y_true = merged["is_contact"].astype(int).values
    y_score = merged["stab_z"].values
    if merged["stab_z"].std() < 1e-3 or abs(
            merged.loc[merged.is_contact, "stab_z"].mean()
            - merged.loc[~merged.is_contact, "stab_z"].mean()) < Z_GAP:
        log.warning(f"[ATTN:{suffix}] Low attention signal — skipping global ROC")
        return float("nan")
    return roc_auc_score(y_true, y_score)

def dump_attention_maps(model, loader, device, out_dir: Path, fold: int, seed: int):
    out_dir.mkdir(parents=True, exist_ok=True)
    model.eval()
    for g in loader:
        g = g.to(device)
        _, node_attn = model(
            g.x, g.edge_index, batch=g.batch, pos=g.pos, u=g.u, pid=g.pid,
            return_node_attn=True
        )
        for gid in g.batch.unique():
            m = (g.batch == gid)
            # keep only nodes with ≥1 edge
            edge_nodes = torch.unique(g.edge_index[:, g.batch[g.edge_index[0]] == gid].flatten())
            m = m & torch.isin(torch.arange(len(g.batch), device=g.batch.device), edge_nodes)
            if m.sum() == 0: continue
            res_ids = g.residue_id[m].cpu()
            scores  = node_attn[m].cpu()
            pdb_id = str(g.pdb_id if isinstance(g.pdb_id, str) else g.pdb_id[0])
            np.savetxt(out_dir / f"{pdb_id}_attn_f{fold}_s{seed}.csv",
                       torch.stack([res_ids, scores], 1).numpy(),
                       fmt=["%d", "%.10f"], delimiter=",",
                       header="residue_id,score", comments="")

def train_epoch(model, loader, optim, scaler, device, lookup, *, huber: bool, accum: int, pos_weight=None, gamma=None):
    model.train(); total, step = 0.0, 0
    optim.zero_grad(set_to_none=True)
    for batch in loader:
        batch = batch.to(device)
        with autocast(dtype=torch.bfloat16):
            reg, _, node_logits, aux = model(batch.x, batch.edge_index, batch=batch.batch, pos=batch.pos, u=batch.u, pid=batch.pid)
            reg_loss = (F.smooth_l1_loss if huber else F.l1_loss)(reg.view(-1), batch.y[:, 0])
            # per-graph node contact loss
            node_loss = 0.0
            n_graph = int(batch.batch.max()) + 1
            for g in range(n_graph):
                m   = batch.batch == g
                pid = batch.meta['protein_id'][g]
                tgt = lookup.mask(pid, batch.residue_id[m]).float()
                node_loss += w_bce_pos_weight(node_logits[m], tgt, pos_weight=pos_weight, gamma=gamma)
            node_loss /= max(n_graph, 1)

            # pair/edge losses (as in your original)
            pair_loss = 0.0
            if hasattr(model, 'pair_head'):
                for g in range(n_graph):
                    m   = batch.batch == g
                    pid = batch.meta['protein_id'][g]
                    tgt_mask = lookup.mask(pid, batch.residue_id[m]).bool()
                    idx = torch.nonzero(m).squeeze(-1)
                    pos_ids = idx[tgt_mask]; neg_ids = idx[~tgt_mask]
                    if pos_ids.numel() < 2 or neg_ids.numel() == 0: continue
                    pos_pairs = torch.combinations(pos_ids)
                    if pos_pairs.size(0) > PAIR_K:
                        pos_pairs = pos_pairs[torch.randperm(pos_pairs.size(0), device=device)[:PAIR_K]]
                    neg_pairs = torch.stack(torch.meshgrid(pos_ids, neg_ids, indexing='ij')).T.reshape(-1, 2)
                    if neg_pairs.size(0) > PAIR_K:
                        neg_pairs = neg_pairs[torch.randperm(neg_pairs.size(0), device=device)[:PAIR_K]]
                    pairs  = torch.cat([pos_pairs, neg_pairs], 0)
                    labels = torch.cat([torch.ones(len(pos_pairs), device=device),
                                        torch.zeros(len(neg_pairs), device=device)])
                    h_i = aux['h'][pairs[:, 0]]; h_j = aux['h'][pairs[:, 1]]
                    dist = (aux['pos'][pairs[:, 0]] - aux['pos'][pairs[:, 1]]).norm(dim=-1, keepdim=True)
                    logits = model.pair_head(torch.cat([h_i, h_j, dist], -1)).squeeze(-1)
                    pair_loss += w_bce_pos_weight(logits, labels, pos_weight=pos_weight, gamma=gamma)

                    edge_mask = (batch.batch[batch.edge_index[0]] == g)
                    if edge_mask.any():
                        edge_pairs  = batch.edge_index[:, edge_mask].t()
                        edge_labels = batch.is_contact[edge_mask].float()
                        h_i  = aux['h'][edge_pairs[:, 0]]
                        h_j  = aux['h'][edge_pairs[:, 1]]
                        dist = (aux['pos'][edge_pairs[:, 0]] - aux['pos'][edge_pairs[:, 1]]).norm(dim=-1, keepdim=True)
                        edge_logits = model.pair_head(torch.cat([h_i, h_j, dist], -1)).squeeze(-1)
                        edge_loss   = w_bce_pos_weight(edge_logits, edge_labels, pos_weight=pos_weight, gamma=gamma)
                        pair_loss  += EDGE_W * edge_loss
                pair_loss /= max(n_graph, 1)

            loss = reg_loss + NODE_LOSS_W * node_loss + BETA_PAIR * pair_loss

        scaler.scale(loss).backward(); step += 1
        if step % accum == 0:
            scaler.unscale_(optim)
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            scaler.step(optim); scaler.update()
            optim.zero_grad(set_to_none=True)
        total += loss.item()
    return total / len(loader)

def run_fold(seed, fold, train_i, val_i, test_i, ds, args, lookup):
    set_seed(seed)
    dev = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    model = PocketStabilityModel(
        in_dim     = ds[0].x.shape[1],
        global_dim = ds[0].u.shape[0],
        hidden_dim = args.hidden_dim,
        dropout    = args.dropout,
        use_u      = args.fusion,
        model_type = args.model
    ).to(dev)
    log.info(f"[fold{fold} seed{seed}] model={args.model} dims: x={ds[0].x.shape[1]} u={ds[0].u.shape[0]} H={args.hidden_dim} drop={args.dropout} fusion={args.fusion}")

    optim  = torch.optim.Adam(model.parameters(), lr=args.lr)
    scaler = GradScaler()

    def make_loader(indices):
        if args.nodecap:
            sub = Subset(ds, indices)
            return DataLoader(sub, batch_sampler=NodeCapSampler(sub))
        return DataLoader(ds, batch_size=args.batch, sampler=SubsetRandomSampler(indices))

    l_train, l_val, l_test = map(make_loader, (train_i, val_i, test_i))

    warm = max(1, int(args.epochs * 0.1))
    warm_sched = torch.optim.lr_scheduler.LinearLR(optim, 0.1, 1.0, warm)
    main_sched = {
        'cosine' : torch.optim.lr_scheduler.CosineAnnealingLR(optim, T_max=args.epochs - warm),
        'plateau': torch.optim.lr_scheduler.ReduceLROnPlateau(optim, mode='min', factor=0.5, patience=3)
    }.get(args.scheduler)
    swa_model = AveragedModel(model)
    swa_sched = SWALR(optim, swa_lr=1e-4)
    swa_phase = False
    swa_start = int(args.epochs * args.swa_start)

    best_epoch, best_mae, best_roc = 0, float("inf"), 0.0
    patience = 8

    rmsd_log["epoch"].clear(); rmsd_log["true"].clear(); rmsd_log["pred"].clear()

    for ep in range(1, args.epochs + 1):
        tl = train_epoch(model, l_train, optim, scaler, dev, lookup, huber=args.huber, accum=ACCUM_STEPS,
                         pos_weight=args.pos_weight, gamma=args.gamma)
        v_mae, v_roc = evaluate(model, l_val, dev, lookup, epoch=ep)
        log.info(f"[fold{fold} seed{seed}] ep{ep:02d} train {tl:.3f} | val {v_mae:.3f}/{v_roc:.3f}")

        if ep <= warm:
            warm_sched.step()
        elif ep < swa_start:
            if args.scheduler == 'plateau': main_sched.step(v_mae)
            elif args.scheduler != 'none': main_sched.step()
        else:
            if not swa_phase:
                log.info(f"[fold{fold} seed{seed}] enter SWA at ep{ep}"); swa_phase = True
            swa_model.update_parameters(model); swa_sched.step()

        if v_roc > best_roc:
            best_epoch, best_mae, best_roc = ep, v_mae, v_roc
            torch.save(model.state_dict(), args.save_dir / f"f{fold}_s{seed}.pt")
            full_loader = DataLoader(ds, batch_size=1, shuffle=False)
            dump_attention_maps(model, full_loader, dev, args.attn_dir, fold, seed)
            log.info(f"Saved checkpoint @ ep{ep} (val_mae {best_mae:.3f} val_roc {best_roc:.3f})")
        elif ep - best_epoch >= patience:
            log.info(f"Early stopping @ ep{ep}. Best ep={best_epoch} (val_MAE={best_mae:.3f})")
            break

    model.load_state_dict(torch.load(args.save_dir / f"f{fold}_s{seed}.pt"))
    t_mae, t_reg_roc = evaluate(model, l_test, dev, lookup)

    t_contact_avg_roc, _, t_contact_roc = roc_by_protein(model, l_test, dev, lookup)
    for p, v in sorted(t_contact_roc.items()):
        log.info(f"  - {p}: {v:.3f}")

    log.info(f"[fold{fold} seed{seed}] TEST mae={t_mae:.3f} | reg-roc={t_reg_roc:.3f} | contact-roc={t_contact_avg_roc:.3f}")
    return t_mae, t_reg_roc, t_contact_avg_roc

def write_summary(args, results):
    summary_csv = args.run_dir / 'summary.csv'
    with summary_csv.open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['fold', 'seed', 'mae', 'roc', 'prot_contact_roc', 'global_roc'])
        writer.writeheader()
        for r in results:
            writer.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})
    import statistics as st
    def ms(vs):
        if not vs: return "nan ± nan"
        m = st.mean(vs); s = st.stdev(vs) if len(vs) > 1 else 0.0
        return f"{m:.3f} ± {s:.3f}"
    log.info(f"Mean ± Std: MAE={ms([r['mae'] for r in results])}, "
             f"ROC={ms([r['roc'] for r in results])}, "
             f"prot_contact_roc={ms([r['prot_contact_roc'] for r in results])}, "
             f"Global ROC={ms([r['global_roc'] for r in results])}")
    log.info(f"Wrote summary → {summary_csv}")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--graph-root', type=Path, required=True)
    ap.add_argument('--run-dir', type=Path, default=None, help='Output directory (default: runs/<timestamp>)')

    # model + train
    ap.add_argument('--model', choices=['egnn', 'egnn_rbf', 'egnn_rbf_att', 'tfn', 'tfn_full'], default='egnn_rbf_att')
    ap.add_argument('--epochs', type=int, default=30)
    ap.add_argument('--batch', type=int, default=4)
    ap.add_argument('--lr', type=float, default=2e-3)
    ap.add_argument('--hidden-dim', type=int, default=256)
    ap.add_argument('--dropout', type=float, default=0.3)
    ap.add_argument('--scheduler', choices=['none', 'cosine', 'plateau'], default='plateau')
    ap.add_argument('--swa-start', type=float, default=0.6)
    ap.add_argument('--cv', type=int, default=4)
    ap.add_argument('--seeds', nargs='*', type=int, default=[0])
    ap.add_argument('--huber', action='store_true')
    ap.add_argument('--nodecap', action='store_true')
    ap.add_argument('--fusion', action='store_true')
    ap.add_argument('--physics-only', action='store_true')
    ap.add_argument('--pos-weight', type=float, default=None)
    ap.add_argument('--gamma', type=float, default=None)
    ap.add_argument('--protein', default=None)

    # features / contacts / attention overlay
    ap.add_argument('--contacts', type=Path, default=None, help="PLIP contacts CSV (default: <graph-root>/contacts.csv)")
    ap.add_argument('--plip-csv', type=Path, default=None, help="Global PLIP summary for attention overlay; optional")
    args = ap.parse_args()

    # run dir
    if args.run_dir is None:
        ts = time.strftime("run_%Y%m%d_%H%M%S")
        args.run_dir = Path("runs") / ts
    args.run_dir.mkdir(parents=True, exist_ok=True)
    args.save_dir = args.run_dir / 'checkpoints'; args.save_dir.mkdir(exist_ok=True)
    args.attn_dir = args.run_dir / 'attention_maps'; args.attn_dir.mkdir(exist_ok=True)

    # logging
    log_file = args.run_dir / 'train_log.txt'
    logging.getLogger().handlers.clear()
    logging.basicConfig(
        format='[%(levelname)s] %(message)s',
        level=logging.INFO,
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )
    log.info(f"Graph root = {args.graph_root}")
    log.info(f"Save dir   = {args.save_dir}")

    # contacts default
    if args.contacts is None:
        candidate = args.graph_root / "contacts.csv"
        if candidate.exists():
            args.contacts = candidate
    if args.contacts is None or not args.contacts.exists():
        raise FileNotFoundError(
            f"--contacts CSV not found (looked for {args.contacts or '<graph-root>/contacts.csv'}). "
            "Provide --contacts to enable node/edge supervision."
        )

    lookup = ContactLookup(args.contacts)

    # feature groups
    if args.physics_only:
        keep = NO_ESM_PLIP_IDX
    else:
        keep = NO_PLIP_IDX
    ds = DriftDataset(args.graph_root, keep_x=keep, mask_drift=True)
    if len(ds) == 0:
        raise RuntimeError(f"No graphs loaded from {args.graph_root} (check paths/feature filters).")

    prot_counts = Counter(m["protein_id"] for m in ds.meta)
    for pid, n in prot_counts.items():
        log.info(f"{pid:22} → {n} graphs")

    proteins    = np.array([m["protein_id"] for m in ds.meta])
    all_indices = np.arange(len(ds))

    n_splits = min(args.cv, len(np.unique(proteins)))
    if n_splits < 2:
        raise ValueError(f"Need ≥2 proteins for CV, found {len(np.unique(proteins))}.")
    gkf = GroupKFold(n_splits=n_splits)
    splits = list(gkf.split(all_indices, groups=proteins))

    results = []
    MIN_ROC_THRESHOLD = 0.7
    start_time = time.time()

    for fold, (train_i, test_i) in enumerate(splits):
        mid = len(test_i) // 2
        val_i, test_i = test_i[:mid], test_i[mid:]
        for seed in args.seeds:
            MAX_RETRIES, retries = 2, 0
            while retries <= MAX_RETRIES:
                # clean old artifacts
                for f in args.attn_dir.glob(f"*f{fold}_s{seed}.csv"): f.unlink(missing_ok=True)
                ckpt = args.save_dir / f"f{fold}_s{seed}.pt"
                if ckpt.exists(): ckpt.unlink()

                mae, reg_roc, prot_contact_roc = run_fold(seed, fold, train_i, val_i, test_i, ds, args, lookup)
                test_pids = {proteins[i] for i in test_i}
                global_roc = summarise_attention(args.attn_dir, args.plip_csv, suffix=f"f{fold}_s{seed}", include_pids=test_pids)

                if not (np.isnan(global_roc) or reg_roc < MIN_ROC_THRESHOLD or prot_contact_roc < MIN_ROC_THRESHOLD):
                    results.append(dict(fold=fold, seed=seed, mae=mae, roc=reg_roc,
                                        prot_contact_roc=prot_contact_roc, global_roc=global_roc))
                    print("-" * 80)
                    print(f"fold{fold} seed{seed}  test_mae={mae:.3f}  test_reg_roc={reg_roc:.3f}  "
                          f"test_prot_contact_roc={prot_contact_roc:.3f}  global_roc={global_roc:.3f}")
                    print("-" * 80)
                    break
                log.warning(f"[RETRY] fold{fold} seed{seed}: reg_roc={reg_roc:.3f} prot_contact_roc={prot_contact_roc:.3f} — retry {retries+1}/{MAX_RETRIES}")
                retries += 1

    write_summary(args, results)
    log.info(f"Execution time: {(time.time() - start_time):.2f}s")

if __name__ == "__main__":
    main()

    # python train.py \
    # --graph-root path/to/graph_data \
    # --model egnn_rbf_att \
    # --cv 4 --seeds 0 1 --epochs 30 \
    # --batch 4 --huber --scheduler cosine --swa-start 0.6 \
    #  --nodecap --contacts path/to/contacts.csv \
    #  --plip-csv path/to/all_proteins_plip_summary.csv
