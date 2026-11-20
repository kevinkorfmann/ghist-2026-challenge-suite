from __future__ import annotations

from pathlib import Path
import gzip
import tempfile
import unittest

import numpy as np
import msprime
import torch

from ghist_secondary_contact.demography import (
    SecondaryContactParameters,
    SizeChange,
    build_demography,
)
from ghist_secondary_contact.coalescent import (
    secondary_contact_interval_probabilities,
    secondary_contact_pair_hazards,
)
from ghist_secondary_contact.calibration import (
    apply_log_correction,
    competition_scores,
    fit_log_correction,
)
from ghist_secondary_contact.spectrum import (
    generator_matched_secondary_contact_spectrum,
    joint_spectrum_from_vcf,
    secondary_contact_spectrum,
    transform_spectrum_counts,
)
from ghist_secondary_contact.model import SecondaryContactNet, mixture_nll
from ghist_secondary_contact.simulate import (
    SimulationNuisance,
    focused_parameter_design,
    nuisance_size_changes,
    parameter_design,
    simulate_block,
)
from ghist_secondary_contact.vcf_fast import scan_secondary_contact_vcf
from ghist_secondary_contact.features import (
    inverse_transform_targets,
    population_pair_hazards,
    summarize_tree_sequence,
    transform_targets,
)


class SecondaryContactTests(unittest.TestCase):
    def test_log_corrections_reduce_known_multiplicative_bias(self):
        rng = np.random.default_rng(90210)
        truth = np.exp(rng.normal(5.0, 0.5, (24, 5)))
        ratios = np.asarray((0.9, 1.1, 0.95, 1.2, 0.8))
        estimates = truth * ratios
        raw = competition_scores(estimates, truth)
        median = fit_log_correction("median_ratio", estimates, truth)
        corrected = apply_log_correction(median, estimates)
        self.assertLess(float(competition_scores(corrected, truth).max()), 1e-12)
        self.assertGreater(float(raw.min()), 0.0)

    def test_ridge_log_correction_is_positive_and_score_aligned(self):
        rng = np.random.default_rng(19073)
        truth = np.exp(rng.normal(5.0, 0.7, (40, 5)))
        logged = np.log(truth)
        estimates = truth * np.exp(0.12 + 0.04 * (logged - logged.mean(axis=0)))
        correction = fit_log_correction("ridge_0.1", estimates, truth)
        corrected = apply_log_correction(correction, estimates)
        self.assertTrue(np.all(corrected > 0))
        self.assertLess(
            float(competition_scores(corrected, truth).mean()),
            float(competition_scores(estimates, truth).mean()),
        )

    def test_analytic_pair_coalescent_is_normalized_and_ordered(self):
        parameters = SecondaryContactParameters(
            ancestral_size=16_000,
            population1_size=10_000,
            population2_size=50_000,
            split_time=5_000,
            contact_time=1_000,
            migration_rate=2e-4,
        )
        probabilities = secondary_contact_interval_probabilities(parameters)
        hazards = secondary_contact_pair_hazards(parameters)
        self.assertEqual(probabilities.shape, (3, 33))
        self.assertEqual(hazards.shape, (3, 32))
        np.testing.assert_allclose(probabilities.sum(axis=1), 1.0, rtol=1e-12)
        self.assertTrue(np.all(probabilities >= 0))
        self.assertTrue(np.all((hazards >= 0) & (hazards <= 1)))
        self.assertGreater(hazards[0, 0], hazards[1, 0])

    def test_multimodal_posterior_contract(self):
        model = SecondaryContactNet(
            window_features=18,
            global_features=40,
            width=64,
            attention_heads=8,
            context_layers=1,
            mixture_components=3,
        )
        windows = torch.randn(2, 50, 18)
        joint_sfs = torch.randn(2, 1, 87, 49)
        global_summary = torch.randn(2, 40)
        logits, means, scale = model.mixture_parameters(
            windows,
            joint_sfs,
            global_summary,
        )
        self.assertEqual(logits.shape, (2, 3))
        self.assertEqual(means.shape, (2, 3, 5))
        self.assertEqual(scale.shape, (2, 3, 5, 5))
        loss = mixture_nll(logits, means, scale, torch.randn(2, 5))
        self.assertTrue(torch.isfinite(loss))
        loss.backward()

    def test_demography_has_isolation_then_symmetric_contact(self):
        parameters = SecondaryContactParameters(
            ancestral_size=10_000,
            population1_size=12_000,
            population2_size=7_000,
            split_time=5_000,
            contact_time=1_000,
            migration_rate=1e-4,
        )
        demography = build_demography(
            parameters,
            size_changes=[SizeChange("ancestral", 5_500, 15_000)],
        )
        demography.validate()
        self.assertEqual(demography.num_populations, 3)
        self.assertAlmostEqual(demography.migration_matrix[1][2], 1e-4)
        self.assertAlmostEqual(demography.migration_matrix[2][1], 1e-4)
        event_times = [event.time for event in demography.events]
        self.assertEqual(event_times, sorted(event_times))
        self.assertIn(1_000, event_times)
        self.assertIn(5_000, event_times)

    def test_generator_matched_hidden_size_epochs(self):
        parameters = SecondaryContactParameters(
            ancestral_size=12_000,
            population1_size=24_000,
            population2_size=3_600,
            split_time=5_000,
            contact_time=700,
            migration_rate=1e-4,
        )
        nuisance = SimulationNuisance(
            ancestral_epoch_factors=(0.5, 1.5, 2.0),
            population1_epoch_factors=(0.8, 0.5),
            population2_epoch_factor=1.5,
            population2_contact_factor=1 / 3,
        )
        changes = nuisance_size_changes(parameters, nuisance)
        island = {
            change.time: change.size
            for change in changes
            if change.population == "population2"
        }
        self.assertAlmostEqual(island[315], 5_400)
        self.assertAlmostEqual(island[700], 1_200)

    def test_joint_spectrum_model_has_expected_dimensions(self):
        spectrum = secondary_contact_spectrum(
            [1.2, 0.7, 0.2, 0.05, 1.0],
            [8, 6],
        )
        self.assertEqual(spectrum.shape, (9, 7))
        self.assertTrue(np.isfinite(spectrum).all())
        self.assertGreater(spectrum.S(), 0)

    def test_generator_matched_spectrum_has_expected_dimensions(self):
        spectrum = generator_matched_secondary_contact_spectrum(
            [1.2, 0.7, 0.2, 0.25, 1.0],
            [8, 6],
        )
        self.assertEqual(spectrum.shape, (9, 7))
        self.assertTrue(np.isfinite(spectrum).all())
        self.assertGreater(spectrum.S(), 0)

    def test_spectrum_count_transforms_are_unit_norm(self):
        counts = np.asarray([[0.0, 3.0], [5.0, 11.0]])
        for method in ("log_l2", "sqrt_proportion"):
            transformed = transform_spectrum_counts(counts, method)
            self.assertAlmostEqual(float(np.linalg.norm(transformed)), 1.0)
            self.assertTrue(np.all(transformed >= 0))

    def test_joint_spectrum_respects_ancestral_allele(self):
        text = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=100>
##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral allele\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tpop1_A\tpop2_B
chr1\t10\t.\tA\tG\t.\tPASS\tAA=A\tGT\t0|1\t1|1
chr1\t20\t.\tC\tT\t.\tPASS\tAA=T\tGT\t0|0\t0|1
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tiny.vcf"
            path.write_text(text)
            spectrum, metadata = joint_spectrum_from_vcf(
                path,
                {"population1": ["pop1_A"], "population2": ["pop2_B"]},
            )
        self.assertEqual(spectrum.shape, (3, 3))
        self.assertEqual(spectrum[1, 2], 1)
        self.assertEqual(spectrum[2, 1], 1)
        self.assertEqual(metadata["ancestral_reference_records"], 1)
        self.assertEqual(metadata["ancestral_alternate_records"], 1)

    def test_population_pair_hazards_have_three_bounded_rows(self):
        parameters = SecondaryContactParameters(
            ancestral_size=20_000,
            population1_size=30_000,
            population2_size=10_000,
            split_time=3_000,
            contact_time=600,
            migration_rate=2e-4,
        )
        ts = msprime.sim_ancestry(
            samples={"population1": 4, "population2": 3},
            demography=build_demography(parameters),
            sequence_length=20_000,
            recombination_rate=5.39e-8,
            random_seed=18,
        )
        hazards = population_pair_hazards(
            ts,
            ts.samples(population=1),
            ts.samples(population=2),
        )
        self.assertEqual(hazards.shape, (3, 32))
        self.assertTrue(np.all((hazards >= 0) & (hazards <= 1)))

    def test_target_transform_round_trip(self):
        physical = np.array(
            [50_000.0, 12_000.0, 8_000.0, 1_500.0, 3e-4],
            dtype=np.float64,
        )
        recovered = inverse_transform_targets(transform_targets(physical))
        np.testing.assert_allclose(recovered, physical, rtol=1e-12)

    def test_tree_sequence_summary_contract(self):
        parameters = SecondaryContactParameters(
            ancestral_size=10_000,
            population1_size=12_000,
            population2_size=7_000,
            split_time=2_000,
            contact_time=500,
            migration_rate=5e-4,
        )
        ancestry = msprime.sim_ancestry(
            samples={"population1": 4, "population2": 3},
            demography=build_demography(parameters),
            sequence_length=100_000,
            recombination_rate=5.39e-8,
            random_seed=27,
        )
        ts = msprime.sim_mutations(
            ancestry,
            rate=8.4e-8,
            model=msprime.BinaryMutationModel(),
            random_seed=28,
        )
        bundle = summarize_tree_sequence(
            ts,
            ts.samples(population=1),
            ts.samples(population=2),
            n_windows=5,
        )
        self.assertEqual(bundle.windows.shape, (5, 12))
        self.assertEqual(bundle.joint_sfs.shape, (9, 7))
        self.assertEqual(bundle.global_summary.shape, (146,))
        self.assertTrue(np.isfinite(bundle.global_summary).all())

    def test_parameter_design_is_deterministic_and_grouped(self):
        first, split1 = parameter_design(40, seed=991)
        second, split2 = parameter_design(40, seed=991)
        np.testing.assert_array_equal(first, second)
        np.testing.assert_array_equal(split1, split2)
        self.assertTrue(np.all(first[:, 4] < first[:, 3]))
        self.assertEqual(set(split1), {0, 1, 2})

    def test_focused_design_has_core_guard_and_broad_domains(self):
        center = np.asarray((16_585, 18_880, 51_090, 4_622, 1_188, 2.43e-4))
        first, split1, domain1 = focused_parameter_design(100, center, seed=813)
        second, split2, domain2 = focused_parameter_design(100, center, seed=813)
        np.testing.assert_array_equal(first, second)
        np.testing.assert_array_equal(split1, split2)
        np.testing.assert_array_equal(domain1, domain2)
        self.assertEqual(dict(zip(*np.unique(domain1, return_counts=True))), {0: 70, 1: 20, 2: 10})
        self.assertTrue(np.all(first[:, 4] < first[:, 3]))
        self.assertEqual(set(split1), {0, 1, 2})

    def test_tiny_simulation_block(self):
        parameters = SecondaryContactParameters(
            ancestral_size=5_000,
            population1_size=6_000,
            population2_size=3_000,
            split_time=1_000,
            contact_time=200,
            migration_rate=5e-4,
        )
        bundle = simulate_block(
            parameters,
            SimulationNuisance(),
            seed=41,
            sequence_length=20_000,
            mutation_rate=2e-7,
            population1_diploids=4,
            population2_diploids=3,
        )
        self.assertEqual(bundle.windows.shape, (1, 12))
        self.assertEqual(bundle.joint_sfs.shape, (9, 7))

    def test_fast_registered_layout_scanner_orients_ancestral_alleles(self):
        names = [f"mainland_{index}" for index in range(1, 44)] + [
            f"island_{index}" for index in range(1, 25)
        ]
        header = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
            + "\t".join(names)
            + "\n"
        )
        genotypes1 = ["0|0"] * 67
        genotypes1[0] = "0|1"
        genotypes1[43] = "1|1"
        genotypes2 = ["0|0"] * 67
        genotypes2[0] = "0|0"
        genotypes2[43] = "0|1"
        records = (
            "1\t10\t.\tA\tG\t.\tPASS\tAA=A\tGT\t"
            + "\t".join(genotypes1)
            + "\n"
            + "1\t20\t.\tC\tT\t.\tPASS\tAA=T\tGT\t"
            + "\t".join(genotypes2)
            + "\n"
        )
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tiny.vcf.gz"
            with gzip.open(path, "wt") as handle:
                handle.write(header + records)
            spectrum, metadata = scan_secondary_contact_vcf(path)
        self.assertEqual(spectrum[1, 2], 1)
        self.assertEqual(spectrum[86, 47], 1)
        self.assertEqual(metadata["records"], 2)
        self.assertEqual(metadata["ancestral_alternate_records"], 1)


if __name__ == "__main__":
    unittest.main()


def test_fast_scanner_rejects_legacy_zero_position_and_malformed_gt(tmp_path):
    import pytest
    for position, gt, message in [(0, '0|1', '1-based'), (100000001, '0|1', '1-based'),
                                   (1, '01|', 'layout validation'),
                                   (1, '|||', 'layout validation')]:
        path = tmp_path / 'invalid.vcf.gz'
        record = f'1\t{position}\t.\tA\tG\t.\tPASS\tAA=A\tGT\t' + '\t'.join([gt] * 67) + '\n'
        with gzip.open(path, 'wt') as handle:
            handle.write(record)
        with pytest.raises(ValueError, match=message):
            scan_secondary_contact_vcf(path)
