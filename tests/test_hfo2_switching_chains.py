from pathlib import Path
import json

from ase.io import read
import numpy as np
import pytest

from scripts.prepare_hfo2_switching_chains import ordered_seed, prepare
from vcneb import endpoint_structure_record


def test_real_ordered_variants_are_distinct_safe_seeds_and_exact_endpoints():
    root = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008/reference_variants"
    initial = read(root / "PO.vasp")
    chains = []
    for name in ("PO_minus_T_preserving_inversion", "PO_minus_T_reversing_inversion"):
        final = read(root / f"{name}.vasp")
        chain = ordered_seed(initial, final)
        assert len(chain) == 9
        for expected, actual in ((initial, chain[0]), (final, chain[-1])):
            assert endpoint_structure_record(expected)["sha256"] == endpoint_structure_record(actual)["sha256"]
        for image in chain:
            assert image.get_chemical_symbols() == initial.get_chemical_symbols()
            assert image.calc is None
            distances = image.get_all_distances(mic=True)
            np.fill_diagonal(distances, np.inf)
            assert distances.min() > 2.0
        chains.append(chain)
    assert np.linalg.norm(chains[0][4].positions - chains[1][4].positions) > .5
    changed = initial.copy()
    changed.cell[0, 0] += .01
    with pytest.raises(ValueError, match="same-cell"):
        ordered_seed(initial, changed)


def test_no_switching_seed_before_real_electronic_gate(tmp_path):
    polar = tmp_path / "polar"
    polar.mkdir()
    (polar / "manifest.json").write_text("{}")
    (polar / "summary.json").write_text(json.dumps({"status": "review_required"}))
    with pytest.raises(ValueError, match="not passed"):
        prepare(tmp_path, polar, tmp_path / "not_created")
    assert not (tmp_path / "not_created").exists()
