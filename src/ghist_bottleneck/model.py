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


class BidirectionalMambaBlock(nn.Module):
    """FastCXT-style forward/backward state-space mixing with an FFN."""

    def __init__(
        self,
        width: int,
        *,
        state_size: int,
        convolution_width: int,
        expand: int,
        dropout: float,
    ):
        super().__init__()
        try:
            from mamba_ssm import Mamba2
        except ImportError as error:
            raise RuntimeError(
                "context_backend='mamba' requires the pinned mamba-ssm runtime"
            ) from error
        self.norm1 = nn.LayerNorm(width)
        self.forward_mixer = Mamba2(
            d_model=width,
            d_state=state_size,
            d_conv=convolution_width,
            expand=expand,
        )
        self.backward_mixer = Mamba2(
            d_model=width,
            d_state=state_size,
            d_conv=convolution_width,
            expand=expand,
        )
        self.merge = nn.Linear(width * 2, width)
        self.dropout = nn.Dropout(dropout)
        self.norm2 = nn.LayerNorm(width)
        self.ffn = nn.Sequential(
            nn.Linear(width, width * 4),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(width * 4, width),
            nn.Dropout(dropout),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        normalized = self.norm1(x)
        forward = self.forward_mixer(normalized)
        backward = self.backward_mixer(normalized.flip(1)).flip(1)
        x = x + self.dropout(self.merge(torch.cat((forward, backward), dim=-1)))
        return x + self.ffn(self.norm2(x))


class BottleneckNet(nn.Module):
    """Fast-context-style local convolutions plus global self-attention."""

    def __init__(
        self,
        window_features: int = 8,
        sfs_bins: int = 116,
        width: int = 96,
        mixture_components: int = 8,
        dilations: tuple[int, ...] = (1, 2, 4, 8),
        context_layers: int = 2,
        attention_heads: int = 8,
        dropout: float = 0.12,
        context_backend: str = "transformer",
        mamba_state: int = 64,
        mamba_convolution: int = 4,
        mamba_expand: int = 2,
        max_windows: int = 64,
    ):
        super().__init__()
        if width % 8 or width % attention_heads:
            raise ValueError("width must be divisible by 8 and attention_heads")
        self.mixture_components = mixture_components
        self.context_backend = context_backend
        self.window_in = nn.Linear(window_features, width)
        self.local = nn.Sequential(*[ResidualConv(width, d) for d in dilations])
        if context_backend == "transformer":
            layer = nn.TransformerEncoderLayer(
                d_model=width, nhead=attention_heads, dim_feedforward=width * 3,
                dropout=dropout, activation="gelu", batch_first=True, norm_first=True,
            )
            self.position = None
            self.context = nn.TransformerEncoder(layer, num_layers=context_layers)
        elif context_backend == "mamba":
            self.position = nn.Parameter(torch.zeros(1, max_windows, width))
            nn.init.trunc_normal_(self.position, std=0.02)
            self.context = nn.Sequential(
                *[
                    BidirectionalMambaBlock(
                        width,
                        state_size=mamba_state,
                        convolution_width=mamba_convolution,
                        expand=mamba_expand,
                        dropout=dropout,
                    )
                    for _ in range(context_layers)
                ]
            )
        else:
            raise ValueError(f"unknown context_backend: {context_backend}")
        self.sfs = nn.Sequential(nn.Linear(sfs_bins, width), nn.GELU(), nn.LayerNorm(width))
        self.head = nn.Sequential(
            nn.Linear(width * 3, width * 2), nn.GELU(), nn.Dropout(dropout),
            # mixture logits + two means + three Cholesky parameters
            nn.Linear(width * 2, mixture_components * 6),
        )

    def mixture_parameters(
        self, windows: torch.Tensor, sfs: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        x = self.window_in(windows)
        x = self.local(x.transpose(1, 2)).transpose(1, 2)
        if self.position is not None:
            if x.shape[1] > self.position.shape[1]:
                raise ValueError(
                    f"received {x.shape[1]} windows, maximum is {self.position.shape[1]}"
                )
            x = x + self.position[:, : x.shape[1]]
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
