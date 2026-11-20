"""Scientific folding/scaling and BAM input-contract regression tests."""
import numpy as np
import pysam
import pytest

from ghist_lowcoverage.pipeline import audit_inputs, folded_spectrum, physical_parameters


def make_inputs(tmp_path, task="bottleneck", length=1000, contig="1"):
    ref = tmp_path / "reference.fasta"
    ref.write_text(">" + contig + "\n" + "A" * length + "\n")
    pysam.faidx(str(ref))
    samples = {f"S{i}": "pop1" if task == "bottleneck" or i < 10 else "pop2" for i in range(20)}
    (tmp_path / "popfile.txt").write_text("".join(f"{s} {p}\n" for s, p in samples.items()))
    for name in samples:
        with pysam.AlignmentFile(str(tmp_path / f"{name}.bam"), "wb", header={
            "HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": contig, "LN": length}],
            "RG": [{"ID": name, "SM": name}],
        }):
            pass
        pysam.index(str(tmp_path / f"{name}.bam"))
    return ref


def test_fold_keeps_fixed_differences_and_symmetrizes_boundary():
    values = np.zeros((5, 5))
    values[0, 0], values[0, 1], values[4, 0], values[0, 4] = 900, 50, 30, 10
    spectrum = folded_spectrum(values, [4, 4])
    assert not spectrum.mask[4, 0] and not spectrum.mask[0, 4]
    assert spectrum[4, 0] == spectrum[0, 4] == 20
    assert spectrum[0, 1] == 50 and spectrum.sum() == 90 and spectrum.folded


def test_fold_matches_moments_without_double_folding():
    import moments
    expected = moments.Spectrum(np.arange(25, dtype=float).reshape(5, 5)).fold()
    got = folded_spectrum(expected.data, [4, 4])
    np.testing.assert_allclose(got.compressed(), expected.compressed())
    np.testing.assert_equal(got.mask, expected.mask)


@pytest.mark.parametrize("bad", [[1, -1, 0, 0, 0], [1, float("nan"), 0, 0, 0], [1, 0, 0, 0, 2], [1, 2]])
def test_fold_rejects_invalid(bad):
    with pytest.raises(ValueError):
        folded_spectrum(bad, [4])


def test_time_uses_ancestral_size():
    result = physical_parameters("bottleneck", [.2, .1], 4000, 1_000_000)
    assert result["ancestral_population_size"] == pytest.approx(10000)
    assert result["targets"] == pytest.approx({"generations": 2000, "post_decline_fraction": .2})
    split = physical_parameters("split_isolation", [.5, 1.5, .25], 4000, 1_000_000)
    assert split["targets"] == pytest.approx({"pop1_population_size": 5000, "pop2_population_size": 15000,
                                             "generations_since_split": 5000})


@pytest.mark.parametrize("task", ["bottleneck", "split_isolation"])
def test_full_roster(tmp_path, task):
    ref = make_inputs(tmp_path, task)
    result = audit_inputs(tmp_path, ref, task, 1000)
    assert len(result["samples"]) == 20 and result["length"] == 1000


def test_missing_bam(tmp_path):
    ref = make_inputs(tmp_path)
    (tmp_path / "S0.bam").unlink()
    with pytest.raises(ValueError, match="incomplete"):
        audit_inputs(tmp_path, ref, "bottleneck", 1000)


def test_coherent_non_numeric_contig(tmp_path):
    ref = make_inputs(tmp_path, contig="chrSynthetic")
    assert audit_inputs(tmp_path, ref, "bottleneck", 1000)["contig"] == "chrSynthetic"


def test_duplicate_roster(tmp_path):
    ref = make_inputs(tmp_path)
    with (tmp_path / "popfile.txt").open("a") as out:
        out.write("S0 pop1\n")
    with pytest.raises(ValueError, match="unique samples"):
        audit_inputs(tmp_path, ref, "bottleneck", 1000)


def test_mismatched_length(tmp_path):
    ref = make_inputs(tmp_path)
    with pytest.raises(ValueError, match="expected length"):
        audit_inputs(tmp_path, ref, "bottleneck", 10_000_000)


def test_unknown_population(tmp_path):
    ref = make_inputs(tmp_path, "split_isolation")
    roster = tmp_path / "popfile.txt"
    roster.write_text(roster.read_text().replace("pop2", "unknown"))
    with pytest.raises(ValueError, match="ten pop1"):
        audit_inputs(tmp_path, ref, "split_isolation", 1000)


def test_unindexed_bam(tmp_path):
    ref = make_inputs(tmp_path)
    (tmp_path / "S0.bam.bai").unlink()
    with pytest.raises(ValueError, match="indexed"):
        audit_inputs(tmp_path, ref, "bottleneck", 1000)


def test_independent_simulated_reads_pass_contract(tmp_path):
    from scripts.lowcoverage_validate import simulate_bams
    output = tmp_path / "simulated"
    record = simulate_bams("bottleneck", output, 1000)
    audit = audit_inputs(output, output / "reference.fasta", "bottleneck", 1000)
    assert record["truth"] == {"generations": 2000, "post_decline_fraction": .2}
    assert len(audit["samples"]) == 20
    starts = []
    for sample in ["synthetic_00", "synthetic_01"]:
        with pysam.AlignmentFile(str(output / f"{sample}.bam"), "rb") as bam:
            reads = list(bam)
            starts.append([r.reference_start for r in reads])
            assert all(r.query_length == 126 and r.is_proper_pair for r in reads)
    assert starts[0] != starts[1]


def test_realsfs_convergence_uses_absolute_change_and_both_stops():
    from ghist_lowcoverage.pipeline import realsfs_diagnostics
    assert not realsfs_diagnostics('lik[50]=-100 diff=-100')['convergence_observed']
    assert realsfs_diagnostics('lik[50]=-100 diff=-1e-10')['convergence_observed']
    assert realsfs_diagnostics('-> Breaking EM(sr2) at iter:501, sqrt(sr2):9e-9')['convergence_observed']
    assert realsfs_diagnostics('-> Breaking EM(sq2) at iter:502, sqrt(sq2):9e-9')['convergence_observed']
    assert not realsfs_diagnostics('Please check that it has converged!')['convergence_observed']


def test_explicit_confirmation_noise_and_design(tmp_path):
    from scripts.lowcoverage_validate import simulate_bams
    design = {'ancestral_size': 8000, 'fraction': .05, 'generations': 500, 'seed': 842002}
    output = tmp_path / 'confirmation'
    record = simulate_bams('bottleneck', output, 1000, design_override=design, depth_sigma=.5, quality=20)
    assert record['design'] == design
    assert record['error_probability'] == .01 and record['individual_depth_lognormal_sigma'] == .5
    with pysam.AlignmentFile(str(output/'synthetic_00.bam'), 'rb') as bam:
        assert all(q == 20 for read in bam for q in read.query_qualities)


def test_better_boundary_fit_cannot_be_hidden(monkeypatch):
    from types import SimpleNamespace
    import scipy.optimize
    from ghist_lowcoverage.pipeline import fit_spectrum
    calls = []
    def optimizer(*args, **kwargs):
        calls.append(1)
        boundary = len(calls) == 1
        return SimpleNamespace(x=np.log([.001 if boundary else .2, .1]), fun=1 if boundary else 10,
                               success=True, message="mock", nfev=1)
    monkeypatch.setattr(scipy.optimize, 'minimize', optimizer)
    values = np.zeros(41)
    values[0], values[1:21] = 100000, 10
    result = fit_spectrum(values, [40], 'bottleneck', restarts=2)
    assert result['best_fit']['at_boundary']
    assert result['best_fit']['negative_log_likelihood'] == 1
    assert not result['optimizer_agreement']['passed']
