from pathlib import Path

from ase.io import read
import numpy as np

from scripts.prepare_hfo2_PO_M_chain import uniform_seed
from vcneb import endpoint_structure_record, validate_path_geometry


def test_registered_seed_keeps_order_and_exact_cached_endpoints():
    root = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008"
    controls = read(root / "author_path_registration/registered_author_path.traj", index=":")
    po = read(root / "reference_variants/PO.vasp")
    m = read(root / "author_path_registration/M_seed_registered.vasp")
    m.positions[0, 0] += .001  # synthetic endpoint change; no physical relaxation claim
    images = uniform_seed(controls, po, m)
    assert len(images) == 9
    for expected, actual in ((po, images[0]), (m, images[-1])):
        assert endpoint_structure_record(expected)["sha256"] == endpoint_structure_record(actual)["sha256"]
    assert all(a.get_chemical_symbols() == po.get_chemical_symbols() for a in images)
    assert all(a.calc is None for a in images)
    validate_path_geometry(images, minimum_distance=1.6, maximum_deformation=.25)
    np.testing.assert_array_equal(images[0].positions, po.positions)
