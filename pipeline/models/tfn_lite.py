import torch
from torch import nn
import torch.nn.functional as F
from torch_scatter import scatter_sum
from e3nn.o3 import Irreps, spherical_harmonics, FullyConnectedTensorProduct
from e3nn.nn import Activation, FullyConnectedNet
from .common import RadialBasis

class TFNLiteLayer(nn.Module):
    def __init__(self, in_irreps: Irreps, out_irreps: Irreps, sh_irreps="2x0e + 2x1o + 1x2e"):
        super().__init__()
        self.in_irreps = in_irreps
        self.sh_irreps = Irreps(sh_irreps)
        self.tp = FullyConnectedTensorProduct(in_irreps, self.sh_irreps, out_irreps, shared_weights=False)
        self.rbf = RadialBasis(num_basis=16)
        self.fc  = FullyConnectedNet([16, 64, self.tp.weight_numel], torch.relu)
        self.out_irreps = self.tp.irreps_out
        self.norm = Activation(
            self.out_irreps,
            [F.silu if mul_ir.ir.l == 0 else None for mul_ir in self.out_irreps]
        )

    def forward(self, x, pos, edge_index):
        row, col = edge_index
        edge_vec = pos[row] - pos[col]
        edge_len = edge_vec.norm(dim=-1, keepdim=True)
        edge_rbf = self.rbf(edge_len)
        w        = self.fc(edge_rbf)
        Y        = spherical_harmonics(self.sh_irreps, edge_vec, normalize=True, normalization="component")
        out = self.tp(x[row], Y, w)
        y   = scatter_sum(out, row, dim=0, dim_size=x.size(0))
        return self.norm(y)
