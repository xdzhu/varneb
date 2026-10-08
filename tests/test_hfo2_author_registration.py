from pathlib import Path
import numpy as np
from ase.io import read
import pytest

from scripts.register_hfo2_author_path import register


ROOT = Path(__file__).resolve().parents[1]


def test_one_species_assignment_is_retained_for_entire_author_path(tmp_path):
    source = ROOT / "benchmarks/hfo2_channels/20261008/author_path_registration/published_seed.traj"
    target = ROOT / "benchmarks/hfo2_channels/20261008/reference_variants/PO.vasp"
    report = register(source, target, tmp_path / "registered")
    assert report["assignment_max_distance_A"] < .013
    assert report["minimum_assignment_margin_A"] > .1
    original, shifted = read(source, index=":"), read(tmp_path / "registered/registered_author_path.traj", index=":")
    r = np.asarray(report["rotation_fractional"])
    shift = np.asarray(report["translation_fractional"])
    mapping = report["source_to_target_indices"]
    for before, after in zip(original, shifted):
        expected = before.get_scaled_positions() @ r.T + shift
        actual = after.get_scaled_positions()[mapping]
        delta = actual - expected
        delta -= np.rint(delta)
        np.testing.assert_allclose(delta, 0, atol=1e-12)
        np.testing.assert_allclose(after.cell.array, r @ before.cell.array @ r.T, atol=1e-12)
        np.testing.assert_allclose(after.get_volume(), before.get_volume())
        np.testing.assert_allclose(after.get_all_distances(mic=True)[np.ix_(mapping, mapping)],
                                   before.get_all_distances(mic=True), atol=1e-10)
    assert all(v["number"] == 14 for v in report["endpoint_audits"]["M_seed_registered"]["symmetry"])
    with pytest.raises(FileExistsError):
        register(source, target, tmp_path / "registered")
