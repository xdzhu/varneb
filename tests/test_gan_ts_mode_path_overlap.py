import numpy as np
import pytest
from ase import Atoms

from scripts.analyze_gan_ts_mode_path_overlap import joint_delta
from vcneb.joint_curvature import JointCurvatureCoordinates


def test_joint_delta_recovers_local_symmetric_displacement():
    reference = Atoms(
        "Ga2N2",
        scaled_positions=[[0.1, 0.2, 0.3], [0.3, 0.4, 0.5],
                          [0.6, 0.7, 0.1], [0.8, 0.1, 0.6]],
        cell=[[3.0, 0.0, 0.0], [0.1, 3.2, 0.0], [0.2, 0.1, 4.0]], pbc=True,
    )
    chart = JointCurvatureCoordinates(reference, 4.0)
    delta = np.linspace(-0.015, 0.015, chart.size)
    target = chart.displaced(delta)
    recovered, rotation = joint_delta(reference, target, chart.cell_scale_A)
    assert recovered == pytest.approx(delta, abs=1e-12)
    assert rotation == pytest.approx(0.0, abs=1e-12)


def test_joint_delta_rejects_atom_order_change():
    reference = Atoms("GaN", positions=[[0, 0, 0], [1, 1, 1]], cell=[3, 3, 3], pbc=True)
    target = Atoms("NGa", positions=[[0, 0, 0], [1, 1, 1]], cell=[3, 3, 3], pbc=True)
    with pytest.raises(ValueError, match="atoms/order"):
        joint_delta(reference, target, 3.0)
