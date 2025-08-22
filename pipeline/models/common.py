import torch
from torch import nn
from torch_scatter import scatter_sum, scatter_softmax

class RadialBasis(nn.Module):
    def __init__(self, num_basis: int = 16, cutoff: float = 10.0):
        super().__init__()
        self.register_buffer("centers", torch.linspace(0, cutoff, num_basis))
        self.register_buffer("widths", torch.full((num_basis,), cutoff / num_basis))
    def forward(self, r: torch.Tensor) -> torch.Tensor:
        x = r.view(-1, 1)
        return torch.exp(-((x - self.centers) ** 2) / (2 * self.widths ** 2))

class GaussianRBF(nn.Module):
    def __init__(self, num_basis: int = 24, cutoff: float = 10.0, learnable: bool = False):
        super().__init__()
        centers = torch.linspace(0.0, cutoff, num_basis)
        widths  = torch.full((num_basis,), cutoff / num_basis)
        if learnable:
            self.centers = nn.Parameter(centers); self.widths = nn.Parameter(widths)
        else:
            self.register_buffer("centers", centers); self.register_buffer("widths", widths)
    def forward(self, r: torch.Tensor) -> torch.Tensor:
        x = r.view(-1, 1)
        return torch.exp(-((x - self.centers) ** 2) / (2.0 * self.widths ** 2))

class AttnPool(nn.Module):
    def __init__(self, dim: int):
        super().__init__()
        self.q = nn.Linear(dim, 1)
    def forward(self, h: torch.Tensor, batch, *, return_weights: bool = False):
        w = scatter_softmax(self.q(h).squeeze(-1), batch)
        pooled = scatter_sum(w.unsqueeze(-1) * h, batch, dim=0)
        return (pooled, w) if return_weights else pooled

MAX_PROTEINS = 20
EMBED_DIM = 8
ENC_PHYS_DIM = 80
