import torch
from torch import nn
from torch_scatter import scatter_softmax
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import add_self_loops
from .common import RadialBasis

class EGNN_RBF_Att(MessagePassing):
    def __init__(self, hidden_dim: int, num_basis: int = 16, use_edge_attn: bool = True):
        super().__init__(aggr="add")
        self.use_edge_attn = use_edge_attn
        self.rbf = RadialBasis(num_basis=num_basis)
        edge_in = 2 * hidden_dim + num_basis
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
        dist = diff.norm(dim=-1, keepdim=True)
        dfeat = self.rbf(dist)
        edge_in = torch.cat([x[row], x[col], dfeat], dim=-1)
        m = self.edge_mlp(edge_in)
        if self.use_edge_attn:
            alpha = scatter_softmax(self.attn(edge_in).squeeze(-1), row)
            m = m * alpha.unsqueeze(-1)
        pos = pos.index_add(0, row, self.coord_mlp(m) * diff / dist.clamp_min(1e-6))
        agg = self.propagate(edge_index, x=x, m=m, size=(x.size(0), x.size(0)))
        x = self.node_mlp(torch.cat([x, agg], dim=-1))
        return x, pos

    def message(self, m): return m
