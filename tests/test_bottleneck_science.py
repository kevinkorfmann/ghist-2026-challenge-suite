from unittest.mock import patch
from types import SimpleNamespace

import msprime
import numpy as np
import pytest

from ghist_bottleneck.provenance import require_diploid_scale
from ghist_bottleneck.simulate import SimulationConfig, simulate_one
from ghist_bottleneck.train import train_ensemble
from ghist_bottleneck.dadi_reference import fit_two_epoch_dadi_spectrum
from ghist_bottleneck.sfs_baseline import spectrum_exposure


def test_gene_copy_sampling_preserves_diploid_coalescent_scale():
    original = msprime.sim_ancestry
    captured = {}

    def capture(**kwargs):
        captured['arguments'] = kwargs
        captured['tree_sequence'] = original(**kwargs)
        return captured['tree_sequence']

    with patch('ghist_bottleneck.simulate.msprime.sim_ancestry', side_effect=capture):
        simulate_one(
            8, 1e12, 1.0, 100.0,
            SimulationConfig(sequence_length=1000, recombination_rate=0,
                             mutation_rate=1e-4, n_windows=8),
            seed=20261120,
        )
    correct = captured['tree_sequence']
    legacy = original(**{**captured['arguments'], 'ploidy': 1})
    assert correct.num_samples == legacy.num_samples == 8
    # Same seed and haplotypes: a diploid population has twice the coalescence
    # times of the mistakenly haploid population, not twice the sample count.
    np.testing.assert_allclose(correct.tables.nodes.time, 2*legacy.tables.nodes.time)


@pytest.mark.parametrize('artifact', [{}, {'coalescent_ploidy': 1}, {'coalescent_ploidy': [2]}, {'coalescent_ploidy': np.nan}])
def test_legacy_population_scaling_is_rejected(artifact):
    with pytest.raises(ValueError, match='diploid coalescent scaling'):
        require_diploid_scale(artifact, context='test')


def test_verified_population_scale_survives_npz(tmp_path):
    path = tmp_path/'verified.npz'
    np.savez(path, coalescent_ploidy=np.int64(2))
    with np.load(path) as artifact:
        require_diploid_scale(artifact, context='test')


def test_variant_only_missing_calls_cannot_define_exposure():
    with pytest.raises(ValueError, match='variant-only'):
        spectrum_exposure({'incomplete_genotype_records': 1}, 1000, None)
    assert spectrum_exposure({'incomplete_genotype_records': 1}, 1000, 850) == 850
    assert spectrum_exposure({'incomplete_genotype_records': 0}, 1000, None) == 1000
    with pytest.raises(ValueError, match='within'):
        spectrum_exposure({}, 1000, 1001)


def test_training_does_not_consume_unverified_simulations(tmp_path):
    path = tmp_path/'legacy.npz'
    np.savez(path, coalescent_ploidy=1)
    with pytest.raises(ValueError, match='Regenerate'):
        train_ensemble(path, tmp_path/'checkpoint')
    assert not (tmp_path/'checkpoint').exists()


def test_retraining_cannot_overwrite_provenance_and_reuse_old_members(tmp_path):
    directory = tmp_path/'checkpoint'
    directory.mkdir()
    old = directory/'model_17.pt'
    old.write_bytes(b'legacy member')
    with pytest.raises(ValueError, match='empty output directory'):
        train_ensemble(tmp_path/'unused.npz', directory, quick=True)
    assert old.read_bytes() == b'legacy member'
    assert not (directory/'normalization.npz').exists()


@pytest.mark.parametrize('blocks', [False, True])
def test_individual_legacy_member_rejected_despite_verified_normalization(tmp_path, blocks):
    from ghist_bottleneck import infer
    np.savez(tmp_path/'normalization.npz', mean=0, std=1, coalescent_ploidy=2)
    (tmp_path/'model_17.pt').touch()
    observed = SimpleNamespace(windows=np.ones((2, 4, 8) if blocks else (4, 8)),
                               sfs=np.ones((2, 10) if blocks else (10,)))
    loader = 'FeatureCollection' if blocks else 'FeatureBundle'
    predictor = infer.predict_blocks if blocks else infer.predict
    with patch(f'ghist_bottleneck.infer.{loader}.load', return_value=observed), \
         patch('ghist_bottleneck.infer.torch.load', return_value={'coalescent_ploidy': 1}):
        with pytest.raises(ValueError, match='Model checkpoint model_17.pt'):
            predictor('unused', tmp_path, tmp_path/'prediction.yaml')
    assert not (tmp_path/'prediction.yaml').exists()


def test_legacy_chromosome_calibration_is_rejected():
    from ghist_bottleneck.infer import apply_chromosome_calibration
    with pytest.raises(ValueError, match='Chromosome calibration'):
        apply_chromosome_calibration(np.zeros(2), {})


def test_legacy_experiments_cannot_select_a_new_candidate(tmp_path):
    from scripts.select_bottleneck_candidate import select_candidate
    directory = tmp_path/'legacy'
    directory.mkdir()
    (directory/'experiment.json').write_text('{"candidate": "legacy", "selection_score": 0}')
    with pytest.raises(ValueError, match='Candidate experiment'):
        select_candidate(tmp_path, require_complete=False)


def test_large_recovery_panel_has_independent_simulation_seeds():
    from scripts.validate_bottleneck_science import panel_work, validate_reused_cases
    work = panel_work(101, 10_000, 202611201, 2)
    cases = [dict(point=p, replicate=r, scenario=s, seed=seed) for p, r, s, _, seed, _ in work]
    validate_reused_cases(cases, 101)
    cases[1]['seed'] = cases[0]['seed'] + 1
    with pytest.raises(ValueError, match='colliding'):
        validate_reused_cases(cases, 101)
    with pytest.raises(ValueError, match='32-bit'):
        panel_work(2, 10_000, 2**32-1, 2)


@pytest.mark.parametrize('kwargs', [{'mutation_rate': 0}, {'sequence_length': -1}, {'restarts': 0}])
def test_diffusion_reference_rejects_invalid_physical_scale(kwargs):
    with pytest.raises(ValueError):
        fit_two_epoch_dadi_spectrum(np.array([0, 10, 5, 3, 0]), **kwargs)


def test_diffusion_reference_recovers_known_expected_spectrum():
    import dadi
    model = dadi.Numerics.make_extrap_log_func(dadi.Demographics1D.two_epoch)
    spectrum = np.asarray(model([.6, .3], [30], [50, 60, 70])) * 10_000
    result = fit_two_epoch_dadi_spectrum(
        spectrum, mutation_rate=1e-7, sequence_length=1e6,
        restarts=4, seed=7319, folded=True,
    )
    assert result['generations'] == pytest.approx(15_000, rel=.04)
    assert result['post_decline_fraction'] == pytest.approx(.6, rel=.02)
    assert result['ancestral_effective_size'] == pytest.approx(25_000, rel=.02)
