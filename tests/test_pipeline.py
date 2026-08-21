import tempfile
import unittest
from unittest.mock import patch
from pathlib import Path
import subprocess
import sys

import numpy as np

from ghist_bottleneck.features import LDAccumulator, FeatureBundle, summarize_variants
from ghist_bottleneck.model import BottleneckNet, mixture_nll
from ghist_bottleneck.simulate import SimulationConfig, simulate_one
from ghist_bottleneck.sfs_baseline import projected_sfs_from_vcf
from ghist_bottleneck.train import _device


class PipelineTests(unittest.TestCase):
    def test_simulation_merge_preserves_nuisance_provenance(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            names = np.asarray(["scenario", "missingness"])
            for index in range(2):
                np.savez(
                    root / f"part_{index}.npz",
                    windows=np.zeros((1, 64, 8)),
                    sfs=np.zeros((1, 116)),
                    targets=np.zeros((1, 2)),
                    ancestral_ne=np.ones(1),
                    nuisance=np.asarray([[index, 0.01 * index]]),
                    nuisance_names=names,
                    n_haplotypes=np.int64(12),
                    seed=np.int64(100 + index),
                )
            output = root / "merged.npz"
            subprocess.run(
                [
                    sys.executable,
                    "scripts/merge_simulations.py",
                    "--input",
                    str(root / "part_*.npz"),
                    "--output",
                    str(output),
                ],
                check=True,
            )
            merged = np.load(output)
            np.testing.assert_array_equal(merged["chunk_seeds"], [100, 101])
            np.testing.assert_array_equal(merged["nuisance_names"], names)
            self.assertEqual(merged["nuisance"].shape, (2, 2))

    def test_ld_accumulator_preserves_local_perfect_linkage(self):
        accumulator = LDAccumulator(pairs_per_bin=64)
        dosage = np.array([0, 0, 1, 1, 2, 2], dtype=np.float32)
        for position in np.arange(0, 2_000_000, 250):
            accumulator.add(float(position), dosage)
        ld = accumulator.result()
        self.assertGreater(np.count_nonzero(ld), 15)
        np.testing.assert_allclose(ld[ld > 0], 1.0, atol=1e-7)

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

    def test_confounder_domain_is_finite(self):
        b = simulate_one(
            12,
            3000,
            0.5,
            10_000,
            SimulationConfig(sequence_length=50_000),
            seed=321,
            mutation_multiplier=1.08,
            recombination_multiplier=0.9,
            recombination_cv=0.6,
            missingness=0.08,
            spatial_suppression=0.3,
        )
        self.assertEqual(b.windows.shape, (64, 8))
        self.assertTrue(np.isfinite(b.windows).all())
        self.assertTrue(np.isfinite(b.sfs).all())


if __name__ == "__main__":
    unittest.main()
