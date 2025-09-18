import torch
from torch import nn
from torch_scatter import scatter_sum, scatter_softmax
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import add_self_loops
from e3nn import o3
import torch._dynamo

from e3nn.o3 import Irreps, FullyConnectedTensorProduct, spherical_harmonics
from e3nn.nn import Activation, FullyConnectedNet


# ===== Full TFN building blocks =====
import torch
from torch import nn
import torch.nn.functional as F
from torch_scatter import scatter_sum
from e3nn.o3 import Irreps, spherical_harmonics, FullyConnectedTensorProduct
from e3nn.nn import Activation, FullyConnectedNet


class GaussianRBF(nn.Module):
    """Gaussian radial basis with optional learnable centers/widths."""
    def __init__(self, num_basis: int = 24, cutoff: float = 10.0,
                 learnable: bool = False):
        super().__init__()
        centers = torch.linspace(0.0, cutoff, num_basis)
        widths = torch.full((num_basis,), cutoff / num_basis)
        if learnable:
            self.centers = nn.Parameter(centers)
            self.widths = nn.Parameter(widths)
        else:
            self.register_buffer("centers", centers)
            self.register_buffer("widths", widths)

    def forward(self, r: torch.Tensor) -> torch.Tensor:
        # r: [E, 1] or [E]
        x = r.view(-1, 1)  # [E, 1]
        return torch.exp(-((x - self.centers) ** 2) / (2.0 * self.widths ** 2))  # [E, B]


class TFNFullLayer(nn.Module):
    """
    Full TFN-style message passing layer:
      - Inputs: x ∈ ℝ^{N × in_irreps.dim}, pos ∈ ℝ^{N×3}, edge_index (2×E)
      - For each edge: compute Y_lm(edge_vec), RBF(r), radial MLP -> TP weights
      - Apply tensor product: (x_src ⊗ Y) · W  → messages
      - Aggregate to targets (kept consistent with your earlier code's scatter over 'row')
      - Equivariant activation: silu on scalar (l=0) channels, identity otherwise
    """
    def __init__(
            self,
            in_irreps: Irreps,
            out_irreps: Irreps,
            sh_irreps: Irreps = Irreps("1x0e + 1x1o + 1x2e"),
            num_rbf: int = 24,
            r_cut: float = 10.0,
            rbf_hidden: int = 128,
            learnable_rbf: bool = False,
            aggregate_to: str = "row",  # "row" to match your earlier convention
            act_on_scalars: bool = True,
    ):
        super().__init__()
        self.in_irreps = in_irreps
        self.out_irreps = out_irreps
        self.sh_irreps = sh_irreps
        self.aggregate_to = aggregate_to
        self.edge_chunk = 16384

        # Tensor product (Clebsch–Gordan handled internally)
        self.tp = FullyConnectedTensorProduct(
            self.in_irreps, self.sh_irreps, self.out_irreps, shared_weights=False
        )

        # Radial path: RBF(r) -> MLP -> weights for TP
        self.rbf = GaussianRBF(num_basis=num_rbf, cutoff=r_cut, learnable=learnable_rbf)
        self.weight_mlp = FullyConnectedNet(
            [num_rbf, rbf_hidden, rbf_hidden, self.tp.weight_numel],
            F.silu
        )

        # Equivariant activation: apply silu only to l=0 channels
        if act_on_scalars:
            acts = [F.silu if mul_ir.ir.l == 0 else None for mul_ir in self.out_irreps]
        else:
            acts = [None for _ in self.out_irreps]
        self.act = Activation(self.out_irreps, acts)

    def forward(self, x, pos, edge_index):
        row, col = edge_index
        edge_vec = pos[row] - pos[col]                      # [E,3]
        r = edge_vec.norm(dim=-1, keepdim=True)             # [E,1]
        Y = spherical_harmonics(self.sh_irreps, edge_vec, normalize=True, normalization="component")
        E = edge_vec.size(0)

        accum_dtype = torch.float32
        out = torch.zeros(x.size(0), self.out_irreps.dim, device=x.device, dtype=accum_dtype)

        use_amp = x.is_cuda
        amp_ctx = torch.cuda.amp.autocast(enabled=use_amp)

        with amp_ctx:
            for s in range(0, E, self.edge_chunk):
                e = slice(s, min(s + self.edge_chunk, E))

                phi = self.rbf(r[e])                        # [e, num_rbf]
                w   = self.weight_mlp(phi)                  # [e, weight_numel]
                msg = self.tp(x[row[e]], Y[e], w)           # [e, out_dim]

                # >>> ensure dtype matches accumulator <<<
                if msg.dtype != out.dtype:
                    msg = msg.to(out.dtype)

                scatter_sum(msg, row[e], dim=0, out=out)    # or scatter_add_ if you prefer

        # Cast back to the model’s working dtype before activation / downstream
        y = out.to(x.dtype)
        y = self.act(y)
        return y


class TFNFullStack(nn.Module):
    """
    Stack L TFNFullLayer layers; feeds irreps_out of each layer as next layer's input.
    """
    def __init__(
            self,
            hidden_dim: int,
            n_layers: int = 2,
            out_irreps_str: str = "24x0e + 12x1o + 6x2e",  # moderate default; adjust if OOM
            sh_irreps_str: str = "1x0e + 1x1o + 1x2e",
            num_rbf: int = 24,
            r_cut: float = 10.0,
            rbf_hidden: int = 128,
            learnable_rbf: bool = False,
    ):
        super().__init__()
        in_irreps = Irreps(f"{hidden_dim}x0e")           # start from scalar node features
        out_irreps = Irreps(out_irreps_str)
        sh_irreps = Irreps(sh_irreps_str)

        layers = []
        curr_in = in_irreps
        for _ in range(n_layers):
            layer = TFNFullLayer(
                curr_in, out_irreps, sh_irreps=sh_irreps,
                num_rbf=num_rbf, r_cut=r_cut, rbf_hidden=rbf_hidden,
                learnable_rbf=learnable_rbf, aggregate_to="row",
                act_on_scalars=True
            )
            layers.append(layer)
            curr_in = out_irreps  # chain irreps across the stack

        self.layers = nn.ModuleList(layers)
        self.final_irreps = out_irreps

    def forward(self, x, pos, edge_index):
        for layer in self.layers:
            x = layer(x, pos, edge_index)
        return x  # irreps = self.final_irreps


# TFN Layer
def edge_length(edge_vec):
    return torch.norm(edge_vec, dim=-1, keepdim=True)

class RadialBasis(nn.Module):
    def __init__(self, num_basis: int = 16, cutoff: float = 10.0):
        super().__init__()
        self.register_buffer("centers", torch.linspace(0, cutoff, num_basis))  # [num_basis]
        self.register_buffer("widths", torch.full((num_basis,), cutoff / num_basis))  # [num_basis]

    def forward(self, edge_lengths: torch.Tensor) -> torch.Tensor:
        # edge_lengths: [num_edges, 1] or [num_edges]
        x = edge_lengths.view(-1, 1)  # [num_edges, 1]
        return torch.exp(-((x - self.centers) ** 2) / (2 * self.widths ** 2))  # [num_edges, num_basis]

class TFNLiteLayer(nn.Module):
    def __init__(self, in_irreps: Irreps, out_irreps: Irreps,
                 sh_irreps="2x0e + 2x1o + 1x2e"):
        super().__init__()
        self.in_irreps = in_irreps
        self.sh_irreps = Irreps(sh_irreps)

        # tensor product
        self.tp = FullyConnectedTensorProduct(
            in_irreps,                     # node features
            self.sh_irreps,                # spherical harmonics on edges
            out_irreps,
            shared_weights=False
        )

        # ❶  Radial basis → 16-d  →  fc → weight_numel
        self.rbf = RadialBasis(num_basis=16)
        self.fc  = FullyConnectedNet([16, 64, self.tp.weight_numel], torch.relu)

        # ❷  Activation on scalar channels only
        self.out_irreps = self.tp.irreps_out
        self.norm = Activation(
            self.out_irreps,
            [F.silu if mul_ir.ir.l == 0 else None for mul_ir in self.out_irreps]
        )

    def forward(self, x, pos, edge_index):
        row, col = edge_index                        # E edges
        # **keep x shape (N, hidden_dim)** ← matches in_irreps.dim already
        edge_vec = pos[row] - pos[col]               # [E,3]
        edge_len = torch.norm(edge_vec, dim=-1, keepdim=True)  # [E,1]

        # ❸  Encode distance with RBF, then fc → weights
        edge_rbf = self.rbf(edge_len)                # [E,16]
        w        = self.fc(edge_rbf)                 # [E, weight_numel]

        edge_attr = spherical_harmonics(
            self.sh_irreps, edge_vec,
            normalize=True, normalization="component"
        )                                            # [E, sh_dim]

        out = self.tp(x[row], edge_attr, w)          # [E, out_irreps.dim]
        y   = scatter_sum(out, row, dim=0, dim_size=x.size(0))
        return self.norm(y)


# --- EGNN with RBF distance features and optional edge attention ---
class EGNN_RBF_Att(MessagePassing):
    def __init__(self, hidden_dim: int, num_basis: int = 16, use_edge_attn: bool = True):
        super().__init__(aggr="add")
        self.use_edge_attn = use_edge_attn
        self.rbf = RadialBasis(num_basis=num_basis)

        edge_in = 2 * hidden_dim + num_basis  # x_i, x_j, RBF(dist)
        self.edge_mlp = nn.Sequential(
            nn.Linear(edge_in, hidden_dim), nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim), nn.SiLU(),
        )

        if use_edge_attn:
            self.attn = nn.Sequential(
                nn.Linear(edge_in, hidden_dim), nn.Tanh(),
                nn.Linear(hidden_dim, 1)
            )

        self.node_mlp = nn.Sequential(
            nn.Linear(2 * hidden_dim, hidden_dim), nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.coord_mlp = nn.Sequential(nn.Linear(hidden_dim, 1), nn.Tanh())

    def forward(self, x, pos, edge_index):
        edge_index, _ = add_self_loops(edge_index, num_nodes=x.size(0))
        row, col = edge_index

        diff = pos[row] - pos[col]
        dist = diff.norm(dim=-1, keepdim=True)  # [E,1]
        dfeat = self.rbf(dist)                  # [E,num_basis]

        edge_in = torch.cat([x[row], x[col], dfeat], dim=-1)  # [E, 2H+num_basis]
        m = self.edge_mlp(edge_in)  # [E,H]

        if self.use_edge_attn:
            alpha = scatter_softmax(self.attn(edge_in).squeeze(-1), row)  # softmax over in-edges of each node
            m = m * alpha.unsqueeze(-1)

        # coordinate update (same EGNN style)
        pos = pos.index_add(0, row, self.coord_mlp(m) * diff / dist.clamp_min(1e-6))

        # node update
        agg = self.propagate(edge_index, x=x, m=m, size=(x.size(0), x.size(0)))
        x = self.node_mlp(torch.cat([x, agg], dim=-1))
        return x, pos

    def message(self, m):  # MessagePassing API
        return m


# EGNN Layer
class EGNNLayer(MessagePassing):
    def __init__(self, hidden_dim: int):
        super().__init__(aggr="add")
        self.edge_mlp = nn.Sequential(
            nn.Linear(2 * hidden_dim + 1, hidden_dim), nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim), nn.SiLU(),
        )
        self.node_mlp = nn.Sequential(
            nn.Linear(2 * hidden_dim, hidden_dim), nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.coord_mlp = nn.Sequential(nn.Linear(hidden_dim, 1), nn.Tanh())

    def forward(self, x, pos, edge_index):
        edge_index, _ = add_self_loops(edge_index, num_nodes=x.size(0))
        row, col = edge_index
        diff = pos[row] - pos[col]
        dist2 = (diff ** 2).sum(dim=-1, keepdim=True)
        m = self.edge_mlp(torch.cat([x[row], x[col], dist2], dim=-1))
        coef = self.coord_mlp(m)
        norm = diff.norm(dim=-1, keepdim=True).clamp_min(1e-6)
        pos = pos.index_add(0, row, coef * diff / norm)
        agg = self.propagate(edge_index, x=x, m=m, size=(x.size(0), x.size(0)))
        x = self.node_mlp(torch.cat([x, agg], dim=-1))
        return x, pos

    def message(self, m): return m

# Attention Pooling
class AttnPool(nn.Module):
    def __init__(self, dim: int):
        super().__init__()
        self.q = nn.Linear(dim, 1)
    def forward(self, h: torch.Tensor, batch, *, return_weights: bool = False):
        w = scatter_softmax(self.q(h).squeeze(-1), batch)
        pooled = scatter_sum(w.unsqueeze(-1) * h, batch, dim=0)
        return (pooled, w) if return_weights else pooled

# PocketStabilityModel
MAX_PROTEINS = 20
EMBED_DIM = 8
ENC_PHYS_DIM = 80

class PocketStabilityModel(nn.Module):
    def __init__(
            self,
            in_dim: int,
            hidden_dim: int,
            global_dim: int,
            n_layers: int = 2,
            dropout: float = 0.3,
            model_type: str = "tfn",   # or
            use_dyn_attn: bool = True,
            use_u: bool = False,
    ):
        super().__init__()
        assert model_type in {"egnn", "tfn",  "tfn_full", "egnn_rbf", "egnn_rbf_att"}, \
            "model_type must be 'egnn', 'tfn', 'egnn_rbf', or 'egnn_rbf_att'"

        self.model_type = model_type
        self.use_u = use_u
        self.use_dyn_attn = use_dyn_attn
        self.hidden_dim = hidden_dim
        self.global_dim = global_dim

        self.x_norm = nn.LayerNorm(in_dim)
        self.lin_in = nn.Linear(in_dim, hidden_dim)
        self.gnn_layers = None          # placeholder for EGNN path
        self.gnn_tfn = None             # placeholder for TFN-full path
        self.to_scalars = None

        irreps_in = Irreps(f"{hidden_dim}x0e")
        irreps_out = Irreps("8x0e + 2x1o")
        irreps_out_str = "24x0e + 12x1o + 6x2e"  # adjust if you want a heavier/lighter TFN


        if model_type == "egnn":
            self.gnn_layers = nn.ModuleList([EGNNLayer(hidden_dim) for _ in range(n_layers)])
            self.readout = nn.Identity()
        elif model_type == "egnn_rbf":
            self.gnn_layers = nn.ModuleList([EGNN_RBF_Att(hidden_dim, num_basis=16, use_edge_attn=False) for _ in range(n_layers)])
            self.readout = nn.Identity()

        elif model_type == "egnn_rbf_att":
            self.gnn_layers = nn.ModuleList([EGNN_RBF_Att(hidden_dim, num_basis=16, use_edge_attn=True) for _ in range(n_layers)])
            self.readout = nn.Identity()
        elif model_type == "tfn_full":
            self.gnn_tfn = TFNFullStack(
                hidden_dim=hidden_dim,
                n_layers=n_layers,
                out_irreps_str=irreps_out_str,
                sh_irreps_str="1x0e + 1x1o + 1x2e",
                num_rbf=24, r_cut=10.0, rbf_hidden=128, learnable_rbf=False
            )
            # Map equivariant output back to hidden_dim scalars for the rest of the network
            self.to_scalars = o3.Linear(self.gnn_tfn.final_irreps, o3.Irreps(f"{hidden_dim}x0e"))
            self.readout = nn.Identity()
        else:
            self.gnn_layers = nn.ModuleList()
            if n_layers > 1:
                # First layer: from raw scalar to desired irreps
                self.gnn_layers.append(TFNLiteLayer(irreps_in, irreps_out))

                # Remaining layers: chain same irreps
                for _ in range(1, n_layers):
                    self.gnn_layers.append(TFNLiteLayer(irreps_out, irreps_out))
            else:
                self.gnn_layers = nn.ModuleList([
                 TFNLiteLayer(irreps_in, irreps_out) for _ in range(n_layers)
            ])
            self.to_scalars = o3.Linear(irreps_out, o3.Irreps(f"{hidden_dim}x0e"))
            self.readout = nn.Identity()

        self.pool = AttnPool(hidden_dim)
        self.pid_emb = nn.Embedding(MAX_PROTEINS, EMBED_DIM)
        self.u_encoder = nn.Sequential(nn.Linear(global_dim, ENC_PHYS_DIM), nn.ReLU(), nn.LayerNorm(ENC_PHYS_DIM))
        joint_dim = hidden_dim + ENC_PHYS_DIM + EMBED_DIM

        self.reg_head = nn.Sequential(nn.Linear(joint_dim, hidden_dim), nn.SiLU(), nn.Dropout(dropout), nn.Linear(hidden_dim, 1))
        self.contact_head = nn.Sequential(nn.Linear(joint_dim, 128), nn.ReLU(), nn.Linear(128, 1))
        self.node_attn_head = (
            nn.Sequential(nn.Linear(joint_dim, hidden_dim), nn.Tanh(), nn.Linear(hidden_dim, 1))
            if use_dyn_attn else nn.Linear(hidden_dim, 1)
        )
        self.node_contact_head = nn.Sequential(nn.Linear(hidden_dim + ENC_PHYS_DIM + EMBED_DIM, 128), nn.ReLU(), nn.Linear(128, 1))
        if use_u:
            self.u2x = nn.Linear(ENC_PHYS_DIM + EMBED_DIM, hidden_dim)

    @torch._dynamo.disable()
    def forward(self, x, edge_index, *, pos, u=None, batch=None, pid=None, return_node_attn=False):
        if batch is None:
            batch = torch.zeros(x.size(0), dtype=torch.long, device=x.device)
        x = self.lin_in(self.x_norm(x))

        if self.use_u and u is not None:
            if u.dim() == 1:
                u = u.unsqueeze(0)
            u_enc = self.u_encoder(u[:, : self.global_dim])
            if pid is not None:
                pid_vec = self.pid_emb(pid)
                if u_enc.size(0) == 1 and pid_vec.size(0) > 1:
                    u_enc = u_enc.expand(pid_vec.size(0), -1)
                u_enc = torch.cat([u_enc, pid_vec], dim=-1)
            x = x + self.u2x(u_enc[batch])


        if self.model_type in {"egnn", "egnn_rbf", "egnn_rbf_att"}:
            for layer in self.gnn_layers:
                x, pos = layer(x, pos, edge_index)
        elif self.model_type == "tfn":   # lite
            for layer in self.gnn_layers:   # or whatever you named it
                x = layer(x, pos, edge_index)
            x = self.to_scalars(x)
        else:
            if self.model_type == "tfn_full":
                x = self.gnn_tfn(x, pos, edge_index)
                x = self.to_scalars(x)

        x = self.readout(x)

        h_pool, w = (self.pool(x, batch, return_weights=True) if return_node_attn else (self.pool(x, batch), None))

        if u is None:
            u = torch.zeros(h_pool.size(0), self.global_dim, device=x.device)
            # ensure u is at least 2-D before slicing/encoding
        if u.dim() == 1:               # e.g. torch.Size([26])
            u = u.unsqueeze(0)         # → torch.Size([1, 26])

        u_enc = self.u_encoder(u[:, : self.global_dim])  # now safe

        # if pid is not None:
        #     u_enc = torch.cat([u_enc, self.pid_emb(pid)], dim=-1)
        if pid is not None:
            pid_vec = self.pid_emb(pid)
            if u_enc.size(0) == 1 and pid_vec.size(0) > 1:
                u_enc = u_enc.expand(pid_vec.size(0), -1)
            u_enc = torch.cat([u_enc, pid_vec], dim=-1)

        joint = torch.cat([h_pool, u_enc], dim=-1)

        y_reg = self.reg_head(joint).squeeze(-1)
        cls_out = self.contact_head(joint).squeeze(-1)
        node_logits = self.node_contact_head(torch.cat([x, u_enc[batch]], dim=-1)).squeeze(-1)

        return ((y_reg, cls_out, node_logits, {"h": x, "pos": pos}), w) if return_node_attn else (y_reg, cls_out, node_logits, {"h": x, "pos": pos})

    def set_phys_stats(self, rmsd_mu, rmsd_std, en_mu, en_std):
        self.rmsd_mu = rmsd_mu
        self.rmsd_std = rmsd_std
        self.en_mu = en_mu
        self.en_std = en_std

    def set_phys_stats_by_protein(self, stats_dict):
        self.phys_stats_by_protein = stats_dict
