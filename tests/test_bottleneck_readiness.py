import numpy as np
import pytest
from scipy.stats import hypergeom

from ghist_bottleneck.sfs_baseline import fit_two_epoch_spectrum, projected_sfs_from_vcf


def test_unpolarized_vcf_cannot_be_fit_as_ancestral_unfolded_spectrum():
    from ghist_bottleneck.sfs_baseline import fit_two_epoch
    with pytest.raises(ValueError, match='Unpolarized.*folded'):
        fit_two_epoch('unused.vcf', polarized=False, folded=False)


def test_grouped_projection_matches_per_site_hypergeometric(tmp_path):
    path = tmp_path / 'calls.vcf'
    calls = ['0|1\t0|0\t1|1'] * 11 + ['0|1\t0|0\t.|.'] * 3 + ['0|0\t0|1\t0|1'] * 2
    path.write_text('##fileformat=VCFv4.2\n##contig=<ID=1,length=1000>\n'
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\n'
                    + ''.join(f'1\t{i+1}\t.\tA\tC\t.\tPASS\t.\tGT\t{gt}\n' for i, gt in enumerate(calls)))
    spectrum, metadata = projected_sfs_from_vcf(path)
    expected = 11*hypergeom.pmf(np.arange(5), 6, 3, 4) + 3*hypergeom.pmf(np.arange(5), 4, 1, 4) + 2*hypergeom.pmf(np.arange(5), 6, 2, 4)
    np.testing.assert_allclose(spectrum, expected)
    assert metadata['projection'] == 4
    assert metadata['segregating_records'] == metadata['projected_records'] == 16
    assert metadata['incomplete_genotype_records'] == 3


def test_finite_nonconverged_fit_is_not_reported_as_converged(monkeypatch):
    import moments
    spectrum = np.asarray(moments.Demographics1D.two_epoch([.4, .2], [10])) * 200
    monkeypatch.setattr(moments.Inference, 'optimize_log_lbfgsb', lambda *a, **k: (
        np.array([.4, .2]), 0., dict(warnflag=1, task='ITERATION LIMIT', funcalls=99, nit=50),
    ))
    fit = fit_two_epoch_spectrum(spectrum, restarts=2, optimizer_strategy="legacy_coarse_v1")
    assert fit['successful_restarts'] == 2
    assert fit['converged_restarts'] == 0
    assert fit['optimizer']['converged'] is False
    from scripts.bootstrap_bottleneck_sfs import record_replica
    record = record_replica(fit)
    assert record['fit'] == fit
    assert 'failure' in record


def test_boundary_reference_can_initialize_bootstrap():
    import moments
    spectrum = np.asarray(moments.Demographics1D.two_epoch([.999, .2], [10])) * 200
    fit = fit_two_epoch_spectrum(spectrum, initial=(.999, .2), restarts=1)
    assert np.isfinite(fit['generations'])
    with pytest.raises(ValueError, match='within'):
        fit_two_epoch_spectrum(spectrum, initial=(1.01, .2), restarts=1)


@pytest.mark.parametrize('changed', [dict(mutation_rate=2e-8), dict(sequence_length=2000), dict(folded=False)])
def test_bootstrap_rejects_changed_physical_scale(changed, monkeypatch):
    # This test isolates scale checks; complete optimizer contracts are tested separately.
    monkeypatch.setattr("scripts.bootstrap_bottleneck_sfs.validate_optimizer_contract", lambda *a, **k: {})
    from scripts.bootstrap_bottleneck_sfs import validate_reference
    from ghist_bottleneck.optimization import optimizer_contract
    reference = dict(folded=True, sequence_length=1000, theta=.4, ancestral_effective_size=10000,
                     engine="moments", optimizer_contract=optimizer_contract("moments"),
                     optimizer={"converged": True, "warnflag": 0}, best_restart_converged=True,
                     optimizer_diagnostics={"restarts_within_0_1_log_likelihood": 2},
                     optimizer_neighborhood={"passes_finite_neighborhood_check": True})
    settings = dict(mutation_rate=1e-8, sequence_length=1000, folded=True)
    validate_reference(reference, **settings)
    with pytest.raises(ValueError, match='disagree'):
        validate_reference(reference, **{**settings, **changed})


def test_bootstrap_rejects_unrelated_input_and_records_legacy_status(tmp_path):
    import hashlib
    from ghist_inference.checkpoint import spectrum_digest
    from scripts.bootstrap_bottleneck_sfs import validate_input_reference
    path = tmp_path / 'input.vcf'
    path.write_text('first input')
    blocks = np.array([[0, 5, 2, 1, 0], [0, 8, 1, 2, 0]])
    reference = {'sample_haplotypes': 4, 'sfs': {'projection': 4, 'spectrum_sha256': spectrum_digest(blocks.sum(axis=0))},
                 'input_vcf_sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    check = validate_input_reference(reference, path, blocks, {'projection': 4})
    assert check['reference_provenance'] == 'verified_input_and_spectrum'
    with pytest.raises(ValueError, match='spectrum_sha256'):
        validate_input_reference(reference, path, blocks + 1, {'projection': 4})
    path.write_text('different input')
    with pytest.raises(ValueError, match='input_vcf_sha256'):
        validate_input_reference(reference, path, blocks, {'projection': 4}, allow_unverified=True)
    legacy = {'sample_haplotypes': 4, 'sfs': {'projection': 4}}
    with pytest.raises(ValueError, match='lacks input provenance'):
        validate_input_reference(legacy, path, blocks, {'projection': 4})
    check = validate_input_reference(legacy, path, blocks, {'projection': 4}, allow_unverified=True)
    assert check['reference_provenance'] == 'unverified_legacy_reference'


def test_vcf_contract_requires_actual_spectrum_checkpoint(tmp_path):
    import json
    from scripts.summarize_bottleneck_readiness import summarize
    (tmp_path / 'protocol.json').write_text(json.dumps({'sequence_length': 1000, 'seed': 1729, 'restarts': 2}))
    case = tmp_path / 'case_00'
    case.mkdir()
    (case / 'vcf_contract.json').write_text('{"passed": true}')
    with pytest.raises(FileNotFoundError):
        summarize(tmp_path)


@pytest.mark.parametrize('corruption', ['error', 'spectrum', 'truth'])
def test_summary_rejects_forged_fit_evidence(tmp_path, corruption):
    import hashlib
    import json
    from ghist_inference.checkpoint import save_spectrum_checkpoint
    from scripts.summarize_bottleneck_readiness import summarize
    from scripts.validate_bottleneck_science import POINTS
    config = dict(schema='bottleneck_target_length_v2', sequence_length=1000,
                  seed=1729, restarts=2, parameter_points=POINTS)
    protocol = tmp_path / 'protocol.json'
    protocol.write_text(json.dumps(config))
    identity = dict(case=1, point=0, replicate=0, scenario='ancestral_epoch',
                    ancestry_seed=1732, mutation_seed=1733, optimizer_seed=1734,
                    protocol_sha256=hashlib.sha256(protocol.read_bytes()).hexdigest())
    truth = dict(zip(('generations', 'post_decline_fraction', 'ancestral_effective_size'), POINTS[0]))
    case = tmp_path / 'case_01'
    case.mkdir()
    metadata = save_spectrum_checkpoint(case / 'spectrum.npz', np.ones((50, 5)),
                                        {**identity, 'truth': truth})
    record = dict(identity=identity, truth=truth, spectrum_sha256=metadata['spectrum_sha256'],
                  fit={**truth, 'optimizer': {'converged': True}},
                  relative_error={key: 0. for key in truth})
    if corruption == 'error':
        record['relative_error']['generations'] = -.5
    elif corruption == 'spectrum':
        record['spectrum_sha256'] = 'wrong'
    else:
        record['truth'] = {**truth, 'generations': 1}
    (case / 'moments.json').write_text(json.dumps(record))
    with pytest.raises(ValueError, match='error|spectrum|truth'):
        summarize(tmp_path)


def test_real_vcf_export_has_exact_length_and_binds_spectrum(tmp_path):
    import msprime
    from ghist_inference.checkpoint import spectrum_digest
    from scripts.validate_bottleneck_readiness import vcf_contract
    ancestry = msprime.sim_ancestry(samples=3, population_size=1000, sequence_length=100,
                                    random_seed=12345)
    ts = msprime.sim_mutations(ancestry, rate=1e-4, model=msprime.BinaryMutationModel(),
                              discrete_genome=False, random_seed=12346)
    blocks = ts.allele_frequency_spectrum(windows=np.linspace(0, 100, 51), mode='site',
                                        polarised=True, span_normalise=False)
    identity = {'case': 'software_test'}
    digest = spectrum_digest(blocks)
    result = vcf_contract(ts, 100, blocks, tmp_path, identity, digest)
    assert result['passed'] is True
    assert result['contig_length'] == 100
    assert result['identity'] == identity
    assert result['spectrum_sha256'] == digest
    assert result['diploid_individuals'] == 3
    assert result['haplotypes'] == 6
