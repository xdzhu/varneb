from __future__ import annotations

import json
import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator

from scripts.prepare_gan_45p7_two_segment_vcneb import (
    _json_default,
    _static_summary,
    split_frames,
)


def _frame(x: float) -> Atoms:
    atoms = Atoms("Ga2N2", positions=[(x, 0, 0), (1.8, 0, 0),
                                     (0, 1.8, 0), (0, 0, 1.8)],
                  cell=[4, 4, 4], pbc=True)
    return atoms


def test_split_keeps_original_order_and_one_shared_center():
    frames = [_frame(index * 0.001) for index in range(29)]
    center = frames[15].copy()
    center.positions[0, 1] += 0.004
    left, right = split_frames(frames, center)
    assert len(left) == 16 and len(right) == 14
    assert np.array_equal(left[-1].positions, right[0].positions)
    assert np.array_equal(left[0].positions, frames[0].positions)
    assert np.array_equal(right[-1].positions, frames[-1].positions)
    assert not np.array_equal(left[-1].positions, frames[15].positions)
    assert np.array_equal(frames[15].positions, _frame(0.015).positions)


def test_split_rejects_nonlocal_center():
    frames = [_frame(index * 0.001) for index in range(29)]
    center = frames[15].copy()
    center.positions[0, 1] += 0.1
    with pytest.raises(ValueError, match="not a local"):
        split_frames(frames, center)


def test_static_summary_preserves_raw_results_for_fixed_endpoint(tmp_path):
    target = _frame(0.015)
    raw = target.copy()
    raw.calc = SinglePointCalculator(raw, energy=-23.1,
                                     forces=np.full((4, 3), 0.001),
                                     stress=np.full(6, -0.29))
    outcar = tmp_path / "OUTCAR"
    outcar.write_text("test raw output", encoding="utf-8")
    report = _static_summary(raw, target, role="final", n_images=16,
                             parameters={"encut": 600, "kpts": np.array([8, 8, 6])}, fingerprints={},
                             outcar=outcar)
    loaded = json.loads(json.dumps(report, default=_json_default))
    assert report["evaluated_image_index"] == 15
    assert report["n_images"] == 16
    assert report["potential_energy_eV"] == -23.1
    assert report["forces_eV_per_A"] == [[0.001] * 3] * 4
    assert report["stress_eV_per_A3_voigt"] == [-0.29] * 6
    assert loaded["calculator_parameters"]["kpts"] == [8, 8, 6]
