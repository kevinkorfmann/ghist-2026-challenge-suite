from __future__ import annotations

from pathlib import Path
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
from ghist_secondary_contact.spectrum import (
    joint_spectrum_from_vcf,
    secondary_contact_spectrum,
)
from ghist_secondary_contact.model import SecondaryContactNet, mixture_nll
from ghist_secondary_contact.simulate import (
    SimulationNuisance,
    parameter_design,
    simulate_block,
)
from ghist_secondary_contact.features import (
    inverse_transform_targets,
    population_pair_hazards,
    summarize_tree_sequence,
    transform_targets,
)


class SecondaryContactTests(unittest.TestCase):
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

    def test_joint_spectrum_model_has_expected_dimensions(self):
        spectrum = secondary_contact_spectrum(
            [1.2, 0.7, 0.2, 0.05, 1.0],
            [8, 6],
        )
        self.assertEqual(spectrum.shape, (9, 7))
        self.assertTrue(np.isfinite(spectrum).all())
        self.assertGreater(spectrum.S(), 0)

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


if __name__ == "__main__":
    unittest.main()
