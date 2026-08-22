from __future__ import annotations

import torch
from torch import nn
from torch.nn import functional as F

from ghist_bottleneck.model import BidirectionalMambaBlock


class ResidualConv1d(nn.Module):
    def __init__(self, width: int, dilation: int):
        super().__init__()
        self.block = nn.Sequential(
            nn.Conv1d(width, width, 3, padding=dilation, dilation=dilation),
            nn.GELU(),
            nn.GroupNorm(8, width),
            nn.Conv1d(width, width, 1),
        )

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        return values + self.block(values)


class ResidualConv2d(nn.Module):
    def __init__(self, width: int):
        super().__init__()
        self.block = nn.Sequential(
            nn.Conv2d(width, width, 3, padding=1),
            nn.GELU(),
            nn.GroupNorm(8, width),
            nn.Conv2d(width, width, 3, padding=1),
            nn.GroupNorm(8, width),
        )

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        return F.gelu(values + self.block(values))


class SecondaryContactNet(nn.Module):
    """Multiresolution joint-SFS and chromosome-context posterior network."""

    target_dimensions = 5

    def __init__(
        self,
        *,
        window_features: int,
        global_features: int,
        width: int = 128,
        attention_heads: int = 8,
        context_layers: int = 3,
        mixture_components: int = 8,
        dropout: float = 0.1,
        context_backend: str = "transformer",
        mamba_state: int = 64,
        mamba_convolution: int = 4,
        mamba_expand: int = 2,
        max_windows: int = 50,
    ):
        super().__init__()
        if width % attention_heads or width % 8:
            raise ValueError("width must be divisible by attention heads and eight")
        self.mixture_components = mixture_components
        self.window_input = nn.Linear(window_features, width)
        self.window_local = nn.Sequential(
            ResidualConv1d(width, 1),
            ResidualConv1d(width, 2),
            ResidualConv1d(width, 4),
            ResidualConv1d(width, 8),
        )
        if context_backend == "transformer":
            encoder_layer = nn.TransformerEncoderLayer(
                d_model=width,
                nhead=attention_heads,
                dim_feedforward=width * 4,
                dropout=dropout,
                activation="gelu",
                batch_first=True,
                norm_first=True,
            )
            self.window_position = None
            self.window_context = nn.TransformerEncoder(
                encoder_layer,
                num_layers=context_layers,
            )
        elif context_backend == "mamba":
            self.window_position = nn.Parameter(torch.zeros(1, max_windows, width))
            nn.init.trunc_normal_(self.window_position, std=0.02)
            self.window_context = nn.Sequential(
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
        self.sfs_input = nn.Sequential(
            nn.Conv2d(1, width // 2, 5, padding=2),
            nn.GELU(),
            nn.Conv2d(width // 2, width, 3, stride=2, padding=1),
            nn.GELU(),
            ResidualConv2d(width),
            nn.Conv2d(width, width, 3, stride=2, padding=1),
            nn.GELU(),
            ResidualConv2d(width),
            nn.AdaptiveAvgPool2d((6, 6)),
        )
        self.global_input = nn.Sequential(
            nn.Linear(global_features, width * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(width * 2, width),
            nn.LayerNorm(width),
        )
        self.modality = nn.Parameter(torch.zeros(1, 3, width))
        nn.init.trunc_normal_(self.modality, std=0.02)
        self.parameter_queries = nn.Parameter(
            torch.empty(1, self.target_dimensions, width)
        )
        nn.init.trunc_normal_(self.parameter_queries, std=0.02)
        self.cross_attention = nn.MultiheadAttention(
            width,
            attention_heads,
            dropout=dropout,
            batch_first=True,
        )
        self.query_norm = nn.LayerNorm(width)
        triangular = self.target_dimensions * (self.target_dimensions + 1) // 2
        component_width = 1 + self.target_dimensions + triangular
        self.posterior = nn.Sequential(
            nn.Linear(width * self.target_dimensions, width * 3),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(width * 3, mixture_components * component_width),
        )
        self.register_buffer(
            "tril_indices",
            torch.tril_indices(self.target_dimensions, self.target_dimensions),
            persistent=False,
        )

    def mixture_parameters(
        self,
        windows: torch.Tensor,
        joint_sfs: torch.Tensor,
        global_summary: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        window_tokens = self.window_input(windows)
        window_tokens = self.window_local(window_tokens.transpose(1, 2)).transpose(1, 2)
        if self.window_position is not None:
            if window_tokens.shape[1] > self.window_position.shape[1]:
                raise ValueError("window sequence exceeds configured maximum")
            window_tokens = window_tokens + self.window_position[:, : window_tokens.shape[1]]
        window_tokens = self.window_context(window_tokens)
        sfs_map = self.sfs_input(joint_sfs)
        sfs_tokens = sfs_map.flatten(2).transpose(1, 2)
        global_token = self.global_input(global_summary)[:, None, :]
        window_tokens = window_tokens + self.modality[:, 0:1]
        sfs_tokens = sfs_tokens + self.modality[:, 1:2]
        global_token = global_token + self.modality[:, 2:3]
        context = torch.cat((window_tokens, sfs_tokens, global_token), dim=1)
        queries = self.parameter_queries.expand(windows.shape[0], -1, -1)
        attended, _ = self.cross_attention(queries, context, context, need_weights=False)
        queries = self.query_norm(queries + attended)
        raw = self.posterior(queries.flatten(1))
        triangular = self.target_dimensions * (self.target_dimensions + 1) // 2
        raw = raw.reshape(
            windows.shape[0],
            self.mixture_components,
            1 + self.target_dimensions + triangular,
        )
        logits = raw[..., 0]
        means = raw[..., 1 : 1 + self.target_dimensions]
        packed_scale = raw[..., 1 + self.target_dimensions :]
        scale = torch.zeros(
            *packed_scale.shape[:-1],
            self.target_dimensions,
            self.target_dimensions,
            dtype=raw.dtype,
            device=raw.device,
        )
        scale[..., self.tril_indices[0], self.tril_indices[1]] = packed_scale
        diagonal = torch.arange(self.target_dimensions, device=raw.device)
        scale[..., diagonal, diagonal] = F.softplus(
            scale[..., diagonal, diagonal]
        ) + 0.025
        return logits, means, scale

    def forward(
        self,
        windows: torch.Tensor,
        joint_sfs: torch.Tensor,
        global_summary: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        logits, means, scale = self.mixture_parameters(
            windows,
            joint_sfs,
            global_summary,
        )
        weights = torch.softmax(logits, dim=-1)
        mean = torch.sum(weights[..., None] * means, dim=1)
        component_variance = torch.sum(scale.square(), dim=-1)
        second = torch.sum(
            weights[..., None] * (component_variance + means.square()),
            dim=1,
        )
        variance = (second - mean.square()).clamp_min(1e-7)
        return mean, torch.log(variance)


def mixture_nll(
    logits: torch.Tensor,
    means: torch.Tensor,
    scale: torch.Tensor,
    target: torch.Tensor,
) -> torch.Tensor:
    distribution = torch.distributions.MultivariateNormal(means, scale_tril=scale)
    component_log_probability = distribution.log_prob(target[:, None, :])
    log_probability = torch.logsumexp(
        torch.log_softmax(logits, dim=-1) + component_log_probability,
        dim=-1,
    )
    return -log_probability.mean()
