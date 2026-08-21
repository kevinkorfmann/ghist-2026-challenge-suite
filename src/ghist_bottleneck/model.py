from __future__ import annotations

import torch
from torch import nn


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

    def __init__(self, window_features: int = 8, sfs_bins: int = 109, width: int = 96):
        super().__init__()
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
            nn.Linear(width * 2, 4),
        )

    def forward(self, windows: torch.Tensor, sfs: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        x = self.window_in(windows)
        x = self.local(x.transpose(1, 2)).transpose(1, 2)
        x = self.context(x)
        pooled = torch.cat((x.mean(1), x.amax(1), self.sfs(sfs)), dim=1)
        out = self.head(pooled)
        return out[:, :2], out[:, 2:].clamp(-7, 5)


def gaussian_nll(mean: torch.Tensor, log_var: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    return 0.5 * torch.mean(torch.exp(-log_var) * (target - mean).square() + log_var)
