import tempfile
import unittest
import json
import os
import gzip
from unittest.mock import patch
from pathlib import Path
import subprocess
import sys

import numpy as np
import msprime

from ghist_bottleneck.features import (
    LDAccumulator,
    FeatureBundle,
    FeatureCollection,
    extract_vcf_blocks,
    summarize_variants,
)
from ghist_bottleneck.model import BottleneckNet, mixture_nll
from ghist_bottleneck.infer import apply_chromosome_calibration
from ghist_bottleneck.arg_features import (
    gaussian_tmrca_hazard_features,
    tmrca_hazard_features,
)
from ghist_bottleneck.simulate import (
    SimulationConfig,
    generate_chromosome_dataset,
    simulate_one,
)
from ghist_bottleneck.sfs_baseline import (
    complete_block_spectra_from_vcf,
    projected_sfs_from_vcf,
)
from ghist_bottleneck.train import _device
from ghist_bottleneck.genotype_matrix import export_genotype_matrix, strict_positions
from ghist_bottleneck.vcf_normalize import normalize_vcf
from ghist_bottleneck.fastcxt_dataset import grouped_splits, windowed_tmrca
from scripts.run_fastcxt_observed import select_checkpoint
from scripts.select_bottleneck_candidate import EXPECTED_CANDIDATES, select_candidate
from scripts.evaluate_bottleneck_chromosomes import (
    apply_affine,
    fit_affine,
    stratified_cell_split,
)


class PipelineTests(unittest.TestCase):
    def test_chromosome_calibration_applies_only_the_chosen_frozen_branch(self):
        calibration = {
            "chosen_branch": "summary",
            "branches": {
                "summary": {
                    "matrix": [[1.2, 0.1], [-0.3, 0.8]],
                    "intercept": [0.4, -0.2],
                },
                "genealogy_proxy": {
                    "matrix": [[99, 0], [0, 99]],
                    "intercept": [99, 99],
                },
            },
        }
        latent = np.array([2.0, -1.0])
        calibrated, provenance = apply_chromosome_calibration(latent, calibration)
        expected = latent @ calibration["branches"]["summary"]["matrix"] + [0.4, -0.2]
        np.testing.assert_allclose(calibrated, expected)
        self.assertEqual(provenance["branch"], "summary")

    def test_chromosome_calibration_is_fit_before_final_evaluation(self):
        cells = np.arange(40)
        scenarios = np.repeat(np.arange(4), 10)
        calibration, selection, evaluation = stratified_cell_split(cells, scenarios, seed=9)
        self.assertFalse(set(calibration) & set(selection))
        self.assertFalse(set(calibration) & set(evaluation))
        self.assertFalse(set(selection) & set(evaluation))
        self.assertEqual(set(np.concatenate((calibration, selection, evaluation))), set(cells))
        prediction = np.column_stack((np.arange(20), np.arange(20) ** 0.5))
        target = prediction @ np.array([[1.2, 0.1], [-0.3, 0.8]]) + [0.4, -0.2]
        matrix, intercept = fit_affine(prediction, target, ridge=1e-9)
        np.testing.assert_allclose(apply_affine(prediction, matrix, intercept), target, atol=1e-7)

    def test_chromosome_validation_repeats_parameters_within_cells(self):
        observed = FeatureBundle(
            windows=np.zeros((64, 8)),
            sfs=np.zeros(116),
            n_haplotypes=4,
            sequence_length=2_000_000,
            n_variants=1,
        )

        def fake_simulation(*args, **kwargs):
            return FeatureBundle(
                windows=np.zeros((64, 8)),
                sfs=np.zeros(116),
                n_haplotypes=4,
                sequence_length=2_000_000,
                n_variants=7,
            )

        with tempfile.TemporaryDirectory() as directory, patch(
            "ghist_bottleneck.simulate.simulate_one", side_effect=fake_simulation
        ):
            output = Path(directory) / "chromosomes.npz"
            generate_chromosome_dataset(
                observed, output, n_cells=2, blocks_per_cell=3, cell_start=10, seed=51
            )
            data = np.load(output)
            np.testing.assert_array_equal(data["cell_id"], [10, 10, 10, 11, 11, 11])
            np.testing.assert_array_equal(data["n_variants"], np.full(6, 7))
            np.testing.assert_allclose(
                data["targets"][:3], np.repeat(data["targets"][0:1], 3, axis=0)
            )
            np.testing.assert_allclose(
                data["targets"][3:], np.repeat(data["targets"][3:4], 3, axis=0)
            )
            self.assertTrue(np.isfinite(data["ancestral_ne"]).all())

    def test_complete_block_spectra_preserve_physical_blocks(self):
        spectra, metadata = complete_block_spectra_from_vcf(
            "tests/fixtures/tiny.vcf", sequence_length=1000, block_length=500
        )
        self.assertEqual(spectra.shape, (2, 5))
        self.assertEqual(spectra[0].sum(), 2)
        self.assertEqual(spectra[1].sum(), 0)
        self.assertEqual(metadata["records_per_block"], [2, 0])

    def test_screen_selection_requires_one_complete_fixed_protocol(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for index, candidate in enumerate(sorted(EXPECTED_CANDIDATES)):
                target = root / candidate
                target.mkdir()
                record = {
                    "candidate": candidate,
                    "selection_score": 1.0 + index / 10,
                    "elapsed_seconds": 10 + index,
                    "dataset_sha256": "a" * 64,
                    "dataset_bytes": 1234,
                    "git_revision": "b" * 40,
                    "protocol": {"split_seed": 991},
                }
                (target / "experiment.json").write_text(json.dumps(record))
            selected = select_candidate(root)
            self.assertEqual(selected["winner"], sorted(EXPECTED_CANDIDATES)[0])
            (root / "baseline" / "experiment.json").unlink()
            with self.assertRaisesRegex(ValueError, "incomplete screen"):
                select_candidate(root)

    def test_fastcxt_checkpoint_selection_uses_validation_only(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            log = root / "csv_metrics" / "version_0"
            checkpoints = log / "checkpoints"
            checkpoints.mkdir(parents=True)
            (log / "metrics.csv").write_text(
                "epoch,step,train_loss,val_loss\n"
                "0,10,3.0,2.0\n"
                "1,20,2.0,1.5\n"
                "2,30,1.0,1.8\n"
            )
            for epoch, step in enumerate((10, 20, 30)):
                (checkpoints / f"epoch={epoch}-step={step}.ckpt").write_bytes(b"x")
            selected, row = select_checkpoint(root)
            self.assertEqual(selected.name, "epoch=1-step=20.ckpt")
            self.assertEqual(float(row["val_loss"]), 1.5)

    def test_fastcxt_cells_are_grouped_across_replicas(self):
        first = grouped_splits(25, 0.2, 71)
        second = grouped_splits(25, 0.2, 71)
        np.testing.assert_array_equal(first, second)
        self.assertEqual(np.sum(first == "test"), 5)

    def test_windowed_tmrca_is_span_weighted(self):
        ts = msprime.sim_ancestry(
            samples=4,
            ploidy=1,
            population_size=10_000,
            sequence_length=10_000,
            recombination_rate=0,
            random_seed=19,
        )
        values = windowed_tmrca(ts, 0, 1, 2_000)
        self.assertEqual(values.shape, (5,))
        np.testing.assert_allclose(values, values[0])

    def test_simulations_use_continuous_mutation_coordinates(self):
        self.assertFalse(SimulationConfig().discrete_genome)

    def test_repeated_vcf_positions_receive_stable_sub_base_order(self):
        adjusted, groups, records = strict_positions(
            np.array([1, 2, 2, 2, 3, 8, 8, 10], dtype=float)
        )
        self.assertTrue(np.all(np.diff(adjusted) > 0))
        self.assertEqual(groups, 2)
        self.assertEqual(records, 5)
        np.testing.assert_allclose(adjusted[[0, 4, 7]], [1, 3, 10])

    def test_vcf_normalization_repairs_only_the_header(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "malformed.vcf.gz"
            destination = root / "normalized.vcf.gz"
            report_path = root / "normalization.json"
            body = (
                "1\t2\t.\tA\tG\t.\tPASS\tAA=A\tGT\t0|1\n"
                "1\t8\t.\tC\tT\t.\tPASS\tAA=C\tGT\t1|0\n"
            )
            malformed = (
                "##fileformat=VCFv4.2\n"
                "##contig=<ID=1,length=10>\n"
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral allele\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample_1\n"
                + body
            )
            with gzip.open(source, "wt") as handle:
                handle.write(malformed)
            report = normalize_vcf(source, destination, report_path)
            self.assertEqual(report["records"], 2)
            self.assertEqual(report["samples"], 1)
            self.assertTrue(Path(str(destination) + ".tbi").exists())
            with gzip.open(destination, "rt") as handle:
                normalized = handle.read()
            self.assertIn('Description="Ancestral Allele.', normalized)
            self.assertTrue(normalized.endswith(body))
            matrix_metadata = export_genotype_matrix(
                destination, root / "matrix", expected_records=2
            )
            self.assertEqual(matrix_metadata["haplotypes"], 2)
            matrix = np.load(root / "matrix" / "genotypes.npy")
            np.testing.assert_array_equal(matrix, [[0, 1], [1, 0]])

    def test_tmrca_posterior_becomes_conditional_hazards(self):
        hazards = tmrca_hazard_features(np.log([30.0, 30.0, 100.0]))
        self.assertEqual(hazards.shape, (32,))
        self.assertAlmostEqual(hazards[0], 0.0)
        self.assertAlmostEqual(hazards.sum(), 5 / 3, places=7)

    def test_gaussian_tmrca_posterior_is_integrated_exactly(self):
        values = np.log([30.0, 30.0, 100.0])
        sampled = tmrca_hazard_features(values)
        integrated = gaussian_tmrca_hazard_features(values, np.zeros_like(values))
        np.testing.assert_allclose(integrated, sampled, atol=1e-12)

        uncertain = gaussian_tmrca_hazard_features(
            np.full(20, np.log(3000.0)), np.full(20, 0.5)
        )
        self.assertEqual(uncertain.shape, (32,))
        self.assertTrue(np.isfinite(uncertain).all())
        self.assertTrue(np.all((uncertain >= 0) & (uncertain <= 1)))

    def test_fastcxt_bridge_maps_finer_translation_blocks(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            collection = FeatureCollection(
                windows=np.zeros((2, 64, 8), dtype=np.float32),
                sfs=np.zeros((2, 116), dtype=np.float32),
                n_haplotypes=12,
                sequence_length=4000.0,
                block_length=2000.0,
                block_starts=np.array([0.0, 2000.0]),
                block_ends=np.array([2000.0, 4000.0]),
                n_variants=np.array([10, 10]),
            )
            collection_path = root / "collection.npz"
            translation_path = root / "translation.npz"
            output_path = root / "enriched.npz"
            collection.save(collection_path)
            np.savez(
                translation_path,
                means=np.array(
                    [
                        [np.log(100.0)],
                        [np.log(120.0)],
                        [np.log(10_000.0)],
                        [np.log(12_000.0)],
                    ]
                ),
                variances=np.zeros((4, 1)),
                index_map=np.column_stack((np.arange(4), np.zeros(4, dtype=int))),
                block_starts=np.array([0.0, 1000.0, 2000.0, 3000.0]),
                block_ends=np.array([1000.0, 2000.0, 3000.0, 4000.0]),
            )
            subprocess.run(
                [
                    sys.executable,
                    "scripts/inject_fastcxt_hazards.py",
                    "--collection",
                    str(collection_path),
                    "--translation",
                    str(translation_path),
                    "--out",
                    str(output_path),
                ],
                check=True,
                env={**os.environ, "PYTHONPATH": "src"},
            )
            enriched = FeatureCollection.load(output_path)
            np.testing.assert_array_equal(enriched.sfs[:, -1], [1.0, 1.0])
            self.assertFalse(np.allclose(enriched.sfs[0, -33:-1], enriched.sfs[1, -33:-1]))

    def test_block_extraction_matches_training_scale_and_keeps_empty_blocks(self):
        collection = extract_vcf_blocks(
            "tests/fixtures/tiny.vcf", sequence_length=1000, block_length=500
        )
        self.assertEqual(collection.windows.shape, (2, 64, 8))
        self.assertEqual(collection.sfs.shape, (2, 116))
        np.testing.assert_array_equal(collection.n_variants, [2, 0])
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "blocks.npz"
            collection.save(path)
            loaded = FeatureCollection.load(path)
            np.testing.assert_allclose(loaded.windows, collection.windows)
            np.testing.assert_array_equal(loaded.block_starts, [0, 500])

    def test_validation_figure_exports_regular_weight_vector_and_raster(self):
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            rng = np.random.default_rng(7)
            target = np.column_stack(
                (rng.uniform(np.log(200), np.log(30_000), 40), rng.normal(0, 1, 40))
            )
            deep = target + rng.normal(0, 0.12, target.shape)
            reference = target + rng.normal(0, 0.18, target.shape)
            input_path = root / "validation.npz"
            np.savez(
                input_path,
                target_latent=target,
                deep_mean_latent=deep,
                reference_latent=reference,
                nuisance_scenario=np.repeat(np.arange(4), 10),
            )
            stem = root / "figure"
            subprocess.run(
                [
                    sys.executable,
                    "scripts/plot_bottleneck_validation.py",
                    "--input",
                    str(input_path),
                    "--output-stem",
                    str(stem),
                ],
                check=True,
            )
            for suffix in (".pdf", ".svg", ".png"):
                self.assertGreater(stem.with_suffix(suffix).stat().st_size, 1000)
            svg = stem.with_suffix(".svg").read_text()
            self.assertNotIn("font-weight: bold", svg.lower())

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
            manifest = json.loads(output.with_suffix(".npz.manifest.json").read_text())
            self.assertEqual(manifest["simulations"], 2)
            self.assertEqual(len(manifest["sha256"]), 64)

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
        polarized, polarized_metadata = projected_sfs_from_vcf(
            "tests/fixtures/tiny.vcf", polarized=True
        )
        self.assertAlmostEqual(float(polarized.sum()), 2.0, places=8)
        self.assertEqual(polarized_metadata["ancestral_reference_records"], 2)

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

    def test_unknown_context_backend_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "unknown context_backend"):
            BottleneckNet(context_backend="unknown")

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
