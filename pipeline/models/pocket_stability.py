import torch
from torch import nn
from e3nn.o3 import Irreps
from .common import AttnPool, MAX_PROTEINS, EMBED_DIM, ENC_PHYS_DIM
from .egnn import EGNNLayer
from .egnn_rbf import EGNN_RBF_Att
from .tfn_lite import TFNLiteLayer
from .tfn_full import TFNFullStack

class PocketStabilityModel(nn.Module):
    def __init__(self, in_dim: int, hidden_dim: int, global_dim: int, n_layers: int = 2,
                 dropout: float = 0.3, model_type: str = "tfn",
                 use_dyn_attn: bool = True, use_u: bool = False):
        super().__init__()
        assert model_type in {"egnn", "tfn", "tfn_full", "egnn_rbf", "egnn_rbf_att"}
        self.model_type, self.use_u = model_type, use_u
        self.hidden_dim, self.global_dim = hidden_dim, global_dim

        self.x_norm = nn.LayerNorm(in_dim)
        self.lin_in = nn.Linear(in_dim, hidden_dim)

        irreps_in  = Irreps(f"{hidden_dim}x0e")
        irreps_out = Irreps("8x0e + 2x1o")

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
            self.gnn_tfn = TFNFullStack(hidden_dim=hidden_dim, n_layers=n_layers,
                                        out_irreps_str="24x0e + 12x1o + 6x2e",
                                        sh_irreps_str="1x0e + 1x1o + 1x2e")
            self.readout = nn.Linear(self.gnn_tfn.final_irreps.dim, hidden_dim)
        else:  # tfn (lite)
            self.gnn_layers = nn.ModuleList(
                [TFNLiteLayer(irreps_in if i == 0 else irreps_out, irreps_out) for i in range(n_layers)]
            )
            self.readout = nn.Linear(irreps_out.dim, hidden_dim)

        self.pool = AttnPool(hidden_dim)
        self.pid_emb = nn.Embedding(MAX_PROTEINS, EMBED_DIM)
        self.u_encoder = nn.Sequential(
            nn.Linear(global_dim, ENC_PHYS_DIM), nn.ReLU(), nn.LayerNorm(ENC_PHYS_DIM)
        )
        if use_u:
            self.u2x = nn.Linear(ENC_PHYS_DIM + EMBED_DIM, hidden_dim)

        joint_dim = hidden_dim + ENC_PHYS_DIM + EMBED_DIM
        self.reg_head = nn.Sequential(nn.Linear(joint_dim, hidden_dim), nn.SiLU(),
                                      nn.Dropout(dropout), nn.Linear(hidden_dim, 1))
        self.contact_head = nn.Sequential(nn.Linear(joint_dim, 128), nn.ReLU(), nn.Linear(128, 1))
        self.node_attn_head = (
            nn.Sequential(nn.Linear(joint_dim, hidden_dim), nn.Tanh(), nn.Linear(hidden_dim, 1))
            if use_dyn_attn else nn.Linear(hidden_dim, 1)
        )
        self.node_contact_head = nn.Sequential(nn.Linear(hidden_dim + ENC_PHYS_DIM + EMBED_DIM, 128),
                                               nn.ReLU(), nn.Linear(128, 1))

    def forward(self, x, edge_index, *, pos, u=None, batch=None, pid=None, return_node_attn=False):
        if batch is None:
            batch = torch.zeros(x.size(0), dtype=torch.long, device=x.device)
        x = self.lin_in(self.x_norm(x))

        if self.use_u and u is not None:
            if u.dim() == 1: u = u.unsqueeze(0)
            u_enc = self.u_encoder(u[:, : self.global_dim])
            if pid is not None:
                pid_vec = self.pid_emb(pid)
                if u_enc.size(0) == 1 and pid_vec.size(0) > 1:
                    u_enc = u_enc.expand(pid_vec.size(0), -1)
                u_enc = torch.cat([u_enc, pid_vec], dim=-1)
            x = x + self.u2x(u_enc[batch])

        if self.model_type in {"egnn", "egnn_rbf", "egnn_rbf_att"}:
            for layer in getattr(self, "gnn_layers"):
                x, pos = layer(x, pos, edge_index)
        elif self.model_type == "tfn":
            for layer in getattr(self, "gnn_layers"):
                x = layer(x, pos, edge_index)
        else:  # tfn_full
            x = self.gnn_tfn(x, pos, edge_index)

        x = self.readout(x)

        h_pool, w = (self.pool(x, batch, return_weights=True)
                     if return_node_attn else (self.pool(x, batch), None))

        if u is None:
            u = torch.zeros(h_pool.size(0), self.global_dim, device=x.device)
        if u.dim() == 1:
            u = u.unsqueeze(0)
        u_enc = self.u_encoder(u[:, : self.global_dim])
        if pid is not None:
            pid_vec = self.pid_emb(pid)
            if u_enc.size(0) == 1 and pid_vec.size(0) > 1:
                u_enc = u_enc.expand(pid_vec.size(0), -1)
            u_enc = torch.cat([u_enc, pid_vec], dim=-1)

        joint = torch.cat([h_pool, u_enc], dim=-1)
        y_reg = self.reg_head(joint).squeeze(-1)
        cls_out = self.contact_head(joint).squeeze(-1)
        node_logits = self.node_contact_head(torch.cat([x, u_enc[batch]], dim=-1)).squeeze(-1)

        return ((y_reg, cls_out, node_logits, {"h": x, "pos": pos}), w) if return_node_attn \
            else (y_reg, cls_out, node_logits, {"h": x, "pos": pos})
