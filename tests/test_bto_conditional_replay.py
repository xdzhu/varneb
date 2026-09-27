"""Read-only branch replay must reproduce a path or stop at a cache miss."""

import hashlib
from argparse import Namespace
from pathlib import Path

import numpy as np
import pytest

from scripts.replay_bto_conditional_branches import (
    _load_replay_optimizer, cache_only_evaluator, replay,
)
from vcneb.mode_surface import ModePlane, relax_orthogonal_at_q
import vcneb.mode_surface as mode_surface


def test_cache_only_replay_preserves_three_branches_and_rejects_missing_point() -> None:
    plane = ModePlane(
        reference=np.zeros(3), basis=np.eye(3), metric_weights=np.ones(3),
        axis_weights=np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]),
        axis_labels=("Q1", "Q2"), amplitude_unit="toy", reference_id="toy/replay",
    )
    q = np.array([0.6, 0.0])
    base = plane.frozen_coordinates(q)
    starts = [base + np.array([0.0, 0.0, 0.8]),
              base + np.array([0.0, 0.0, -0.8])]
    cache: dict[str, tuple[dict, str]] = {}

    def record(coordinates: np.ndarray) -> tuple[float, np.ndarray]:
        z = float(coordinates[2])
        energy = (z * z - 1.0)**2
        gradient = np.array([0.0, 0.0, 4.0 * z * (z * z - 1.0)])
        key = hashlib.sha256(np.asarray(coordinates, dtype=float).tobytes()).hexdigest()
        cache[key] = ({"enthalpy_eV": energy, "gradient": gradient.tolist()}, key)
        return energy, gradient

    options = dict(
        starts=starts, optimizer="safeguarded_bfgs", trial_validator=lambda _: True,
        gradient_tolerance=1e-8, max_iterations=80, initial_trust_radius=0.2,
    )
    original = relax_orthogonal_at_q(plane, q, record, **options)
    replay_optimizer = _load_replay_optimizer(Path(mode_surface.__file__))
    cached_only, hits = cache_only_evaluator(cache, reference_energy=0.0)
    repeated = replay_optimizer(plane, q, cached_only, **options)
    assert hits
    assert repeated.selected_start == original.selected_start
    assert np.array_equal(repeated.coordinates, original.coordinates)
    assert [row.converged for row in repeated.branch_outcomes] == [True, True, True]
    assert np.allclose([row.energy for row in repeated.branch_outcomes], [1.0, 0.0, 0.0])

    incomplete = dict(cache)
    incomplete.pop(hashlib.sha256(base.tobytes()).hexdigest())
    cached_only, _ = cache_only_evaluator(incomplete, reference_energy=0.0)
    with pytest.raises(RuntimeError, match="absent from audited DFT cache"):
        replay_optimizer(plane, q, cached_only, **options)


def test_replay_summary_name_rejects_path_traversal(tmp_path) -> None:
    args = Namespace(workdir=tmp_path, output=tmp_path / "result.json",
                     summary_name="../another-point.json")
    with pytest.raises(ValueError, match="simple visible"):
        replay(args)
