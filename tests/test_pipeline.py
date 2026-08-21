import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path

import numpy as np

from ghist_bottleneck.features import FeatureBundle, summarize_variants
from ghist_bottleneck.model import BottleneckNet, mixture_nll
from ghist_bottleneck.simulate import SimulationConfig, simulate_one
from ghist_bottleneck.sfs_baseline import projected_sfs_from_vcf
from ghist_bottleneck.train import _device


class PipelineTests(unittest.TestCase):
    @patch("ghist_bottleneck.train.torch.cuda.is_available", return_value=True)
    def test_cuda_is_preferred_when_available(self, _cuda_available):
        self.assertEqual(_device().type, "cuda")

    def test_missingness_projected_sfs(self):
        spectrum, metadata = projected_sfs_from_vcf("tests/fixtures/tiny.vcf")
        self.assertEqual(metadata["projection"], 4)
        self.assertAlmostEqual(float(spectrum.sum()), 2.0, places=8)

    def test_summary_is_fixed_shape_and_finite(self):
        b = summarize_variants(
            np.array([10, 20, 999]),
            np.array([1, 3, 5]),
            np.array([10, 10, 10]),
            np.array([0.2, 0.4, 0.6]),
            sequence_length=1000,
            n_haplotypes=10,
        )
        self.assertEqual(b.windows.shape, (64, 8))
        self.assertEqual(b.sfs.shape, (116,))
        self.assertTrue(np.isfinite(b.windows).all())
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "bundle.npz"
            b.save(path)
            loaded = FeatureBundle.load(path)
            np.testing.assert_allclose(loaded.windows, b.windows, rtol=1e-6)

    def test_coalescent_and_network_contract(self):
        b = simulate_one(
            12, 3000, 0.5, 10_000,
            SimulationConfig(sequence_length=50_000), seed=123,
        )
        self.assertEqual(b.sfs.shape, (116,))
        self.assertEqual(b.sfs[-1], 1.0)
        self.assertTrue(np.isfinite(b.sfs).all())
        model = BottleneckNet()
        import torch
        mean, log_var = model(
            torch.from_numpy(b.windows[None]).float(),
            torch.from_numpy(b.sfs[None]).float(),
        )
        self.assertEqual(tuple(mean.shape), (1, 2))
        self.assertEqual(tuple(log_var.shape), (1, 2))
        logits, means, scale = model.mixture_parameters(
            torch.from_numpy(b.windows[None]).float(),
            torch.from_numpy(b.sfs[None]).float(),
        )
        loss = mixture_nll(logits, means, scale, torch.zeros((1, 2)))
        self.assertTrue(torch.isfinite(loss))
        loss.backward()


if __name__ == "__main__":
    unittest.main()
