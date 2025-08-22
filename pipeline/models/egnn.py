import torch
from torch import nn
from torch_geometric.nn import MessagePassing
from torch_geometric.utils import add_self_loops

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

    def message(self, m):
        return m
