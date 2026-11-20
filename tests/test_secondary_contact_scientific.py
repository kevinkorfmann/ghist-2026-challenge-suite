import numpy as np
import pytest
import tskit

from ghist_secondary_contact.contract import (
    SPECTRUM_CONTRACT,
    require_spectrum_contract,
)
from ghist_secondary_contact.demography import SecondaryContactParameters
from ghist_secondary_contact.features import summarize_tree_sequence
from ghist_secondary_contact.simulate import SimulationNuisance, nuisance_size_changes


def test_fixed_differences_survive_training_feature_extraction():
    tables = tskit.TableCollection(sequence_length=100)
    for _ in range(4):
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    left = tables.nodes.add_row(time=2)
    right = tables.nodes.add_row(time=3)
    root = tables.nodes.add_row(time=10)
    for parent, child in (
        (left, 0),
        (left, 1),
        (right, 2),
        (right, 3),
        (root, left),
        (root, right),
    ):
        tables.edges.add_row(0, 100, parent, child)
    for position, node in ((10, left), (20, right), (30, root)):
        site = tables.sites.add_row(position, ancestral_state="0")
        tables.mutations.add_row(site, node=node, derived_state="1")
    tables.sort()
    bundle = summarize_tree_sequence(
        tables.tree_sequence(), np.array([0, 1]), np.array([2, 3]), n_windows=1
    )
    assert bundle.n_variants == 2
    assert bundle.joint_sfs[0, 2] == bundle.joint_sfs[2, 0] == pytest.approx(2**-0.5)
    assert bundle.joint_sfs[0, 0] == bundle.joint_sfs[2, 2] == 0
    assert bundle.windows[0, 8] == 1  # fixed-difference fraction
    assert bundle.windows[0, 3] == 0.02  # dxy per base


def test_hidden_mainland_epochs_apply_even_after_contact():
    parameters = SecondaryContactParameters(10000, 15000, 6000, 12000, 10000, 1e-4)
    changes = nuisance_size_changes(
        parameters, SimulationNuisance(population1_epoch_factors=(0.8, 0.5))
    )
    assert [(change.time, change.size) for change in changes] == [
        (9600, 8000),
        (6000, 5000),
    ]


def test_old_secondary_contact_evidence_requires_regeneration(tmp_path):
    from scripts.train_secondary_contact_tabular import _load

    path = tmp_path / "old.npz"
    np.savez(path, truth=np.ones((2, 5)))
    with pytest.raises(ValueError, match="fixed-difference"):
        _load(str(path))
    with pytest.raises(ValueError, match="fixed-difference"):
        require_spectrum_contract({})
    assert (
        require_spectrum_contract({"spectrum_contract": SPECTRUM_CONTRACT})
        == SPECTRUM_CONTRACT
    )


def test_surrogate_sampling_excludes_masked_monomorphic_bins(monkeypatch):
    import moments
    from scripts import generate_secondary_contact_sfs_dataset as generate

    model = moments.Spectrum(np.ones((87, 49)))
    monkeypatch.setattr(
        generate, "secondary_contact_spectrum", lambda *args, **kwargs: model
    )
    row = generate._one_point(
        (
            0,
            np.array([10000, 15000, 6000, 12000, 4000, 1e-4]),
            1,
            1e-7,
            100000,
            42,
            "simple",
            "log_l2",
        )
    )
    assert row["spectra"][0, 0, 0] == row["spectra"][0, -1, -1] == 0
    assert row["spectra"][0, 0, -1] > 0
    assert row["spectra"][0, -1, 0] > 0
    model[1, 1] = -1
    with pytest.raises(ValueError, match="likelihood bins"):
        generate._one_point(
            (
                0,
                np.array([10000, 15000, 6000, 12000, 4000, 1e-4]),
                1,
                1e-7,
                100000,
                42,
                "simple",
                "log_l2",
            )
        )
