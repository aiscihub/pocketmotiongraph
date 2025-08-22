import torch
from torch import nn
import torch.nn.functional as F
from torch_scatter import scatter_sum
from e3nn.o3 import Irreps, spherical_harmonics, FullyConnectedTensorProduct
from e3nn.nn import Activation, FullyConnectedNet
from .common import GaussianRBF

class TFNFullLayer(nn.Module):
    def __init__(self, in_irreps: Irreps, out_irreps: Irreps,
                 sh_irreps: Irreps = Irreps("1x0e + 1x1o + 1x2e"),
                 num_rbf: int = 24, r_cut: float = 10.0, rbf_hidden: int = 128,
                 learnable_rbf: bool = False, aggregate_to: str = "row", act_on_scalars: bool = True):
        super().__init__()
        self.sh_irreps = sh_irreps
        self.tp = FullyConnectedTensorProduct(in_irreps, self.sh_irreps, out_irreps, shared_weights=False)
        self.rbf = GaussianRBF(num_basis=num_rbf, cutoff=r_cut, learnable=learnable_rbf)
        self.weight_mlp = FullyConnectedNet([num_rbf, rbf_hidden, rbf_hidden, self.tp.weight_numel], F.silu)
        acts = [F.silu if mul_ir.ir.l == 0 else None for mul_ir in out_irreps] if act_on_scalars else [None for _ in out_irreps]
        self.act = Activation(out_irreps, acts)
        self.out_irreps = out_irreps
        self.edge_chunk = 16384

    def forward(self, x, pos, edge_index):
        row, col = edge_index
        edge_vec = pos[row] - pos[col]
        r = edge_vec.norm(dim=-1, keepdim=True)
        Y = spherical_harmonics(self.sh_irreps, edge_vec, normalize=True, normalization="component")
        E = edge_vec.size(0)
        out = torch.zeros(x.size(0), self.out_irreps.dim, device=x.device, dtype=torch.float32)
        use_amp = x.is_cuda
        amp_ctx = torch.cuda.amp.autocast(enabled=use_amp)
        with amp_ctx:
            for s in range(0, E, self.edge_chunk):
                e = slice(s, min(s + self.edge_chunk, E))
                phi = self.rbf(r[e])
                w   = self.weight_mlp(phi)
                msg = self.tp(x[row[e]], Y[e], w)
                if msg.dtype != out.dtype:
                    msg = msg.to(out.dtype)
                scatter_sum(msg, row[e], dim=0, out=out)
        y = self.act(out.to(x.dtype))
        return y

class TFNFullStack(nn.Module):
    def __init__(self, hidden_dim: int, n_layers: int = 2,
                 out_irreps_str: str = "24x0e + 12x1o + 6x2e",
                 sh_irreps_str: str = "1x0e + 1x1o + 1x2e",
                 num_rbf: int = 24, r_cut: float = 10.0, rbf_hidden: int = 128, learnable_rbf: bool = False):
        super().__init__()
        in_irreps = Irreps(f"{hidden_dim}x0e")
        out_irreps = Irreps(out_irreps_str)
        sh_irreps  = Irreps(sh_irreps_str)
        layers = []
        curr = in_irreps
        for _ in range(n_layers):
            layers.append(TFNFullLayer(curr, out_irreps, sh_irreps=sh_irreps,
                                       num_rbf=num_rbf, r_cut=r_cut, rbf_hidden=rbf_hidden,
                                       learnable_rbf=learnable_rbf))
            curr = out_irreps
        self.layers = nn.ModuleList(layers)
        self.final_irreps = out_irreps

    def forward(self, x, pos, edge_index):
        for layer in self.layers:
            x = layer(x, pos, edge_index)
        return x
