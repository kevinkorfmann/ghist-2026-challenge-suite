from __future__ import annotations

import torch
from torch import nn
from torch.nn import functional as F


class ResidualConv(nn.Module):
    def __init__(self, width: int, dilation: int):
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv1d(width, width, 3, padding=dilation, dilation=dilation),
            nn.GELU(),
            nn.GroupNorm(8, width),
            nn.Conv1d(width, width, 1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.net(x)


class BottleneckNet(nn.Module):
    """Fast-context-style local convolutions plus global self-attention."""

    def __init__(
        self,
        window_features: int = 8,
        sfs_bins: int = 116,
        width: int = 96,
        mixture_components: int = 8,
    ):
        super().__init__()
        self.mixture_components = mixture_components
        self.window_in = nn.Linear(window_features, width)
        self.local = nn.Sequential(*[ResidualConv(width, d) for d in (1, 2, 4, 8)])
        layer = nn.TransformerEncoderLayer(
            d_model=width, nhead=8, dim_feedforward=width * 3,
            dropout=0.12, activation="gelu", batch_first=True, norm_first=True,
        )
        self.context = nn.TransformerEncoder(layer, num_layers=2)
        self.sfs = nn.Sequential(nn.Linear(sfs_bins, width), nn.GELU(), nn.LayerNorm(width))
        self.head = nn.Sequential(
            nn.Linear(width * 3, width * 2), nn.GELU(), nn.Dropout(0.12),
            # mixture logits + two means + three Cholesky parameters
            nn.Linear(width * 2, mixture_components * 6),
        )

    def mixture_parameters(
        self, windows: torch.Tensor, sfs: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        x = self.window_in(windows)
        x = self.local(x.transpose(1, 2)).transpose(1, 2)
        x = self.context(x)
        pooled = torch.cat((x.mean(1), x.amax(1), self.sfs(sfs)), dim=1)
        out = self.head(pooled)
        batch = out.shape[0]
        out = out.reshape(batch, self.mixture_components, 6)
        logits = out[..., 0]
        means = out[..., 1:3]
        raw = out[..., 3:6]
        scale = torch.zeros(
            batch, self.mixture_components, 2, 2,
            device=out.device, dtype=out.dtype,
        )
        scale[..., 0, 0] = F.softplus(raw[..., 0]) + 0.03
        scale[..., 1, 0] = raw[..., 1]
        scale[..., 1, 1] = F.softplus(raw[..., 2]) + 0.03
        return logits, means, scale

    def forward(self, windows: torch.Tensor, sfs: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        logits, means, scale = self.mixture_parameters(windows, sfs)
        weights = torch.softmax(logits, dim=-1)
        mean = torch.sum(weights[..., None] * means, dim=1)
        component_var = torch.sum(scale.square(), dim=-1)
        second = torch.sum(weights[..., None] * (component_var + means.square()), dim=1)
        var = (second - mean.square()).clamp_min(1e-7)
        return mean, torch.log(var)


def gaussian_nll(mean: torch.Tensor, log_var: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    return 0.5 * torch.mean(torch.exp(-log_var) * (target - mean).square() + log_var)


def mixture_nll(
    logits: torch.Tensor,
    means: torch.Tensor,
    scale: torch.Tensor,
    target: torch.Tensor,
) -> torch.Tensor:
    distribution = torch.distributions.MultivariateNormal(means, scale_tril=scale)
    component_log_prob = distribution.log_prob(target[:, None, :])
    log_prob = torch.logsumexp(torch.log_softmax(logits, dim=-1) + component_log_prob, dim=-1)
    return -log_prob.mean()
