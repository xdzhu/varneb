"""The GaN 17x17 increment must nest the existing 9x9 DFT coordinates."""

import numpy as np
import pytest

from scripts.prepare_gan_600eV_ts_2d_dense17 import missing_17x17
from scripts.prepare_gan_600eV_ts_2d_refinement import coordinate_key
from scripts.plot_gan_600eV_local_joint_cut17 import surface


def test_17x17_adds_only_208_unique_midpoints():
    prior = {coordinate_key(float(u), float(v))
             for u in np.linspace(-0.02, 0.02, 9)
             for v in np.linspace(-0.0125, 0.0125, 9)}
    u, v, missing = missing_17x17(prior)
    assert len(u) == len(v) == 17
    assert len(missing) == 208
    assert {coordinate_key(q_u, q_v) for _, _, q_u, q_v in missing}.isdisjoint(prior)
    assert len(prior | {coordinate_key(q_u, q_v) for _, _, q_u, q_v in missing}) == 289


def test_17x17_rejects_missing_or_off_grid_prior():
    prior = {coordinate_key(float(u), float(v))
             for u in np.linspace(-0.02, 0.02, 9)
             for v in np.linspace(-0.0125, 0.0125, 9)}
    with pytest.raises(ValueError, match="audited 9x9"):
        missing_17x17(prior - {(0.0, 0.0)})
    with pytest.raises(ValueError, match="audited 9x9"):
        missing_17x17((prior - {(0.0, 0.0)}) | {(0.001, 0.0)})


def test_17x17_plot_interpolation_is_finite_and_bounded_to_samples():
    rows = [{"q_u_A": float(u), "q_v_A": float(v),
             "measured_delta_H_meV_per_GaN": float(-u * u + v * v)}
            for v in np.linspace(-0.0125, 0.0125, 17)
            for u in np.linspace(-0.02, 0.02, 17)]
    uu, vv, energy = surface(rows)
    assert uu.shape == vv.shape == energy.shape == (251, 401)
    assert np.isfinite(energy).all()
    assert energy[125, 200] == pytest.approx(0.0, abs=1e-12)
