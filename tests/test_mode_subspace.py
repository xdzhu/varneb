"""Exact endpoint-gated two-mode lift into VCNEB's 3N+9 coordinates."""

from __future__ import annotations

import json
import numpy as np
import pytest
from ase import Atoms

from vcneb import (
    ModePlane, ReferenceCellCoordinates, VCNEB,
    endpoint_completed_mode_subspace_for_vcneb,
    interpolate_vcneb, load_vcneb_subspace_artifact,
    make_vcneb_subspace_artifact, strict_mode_subspace_for_vcneb,
)
from vcneb.mode_subspace import _vcneb_x


def _setup():
    reference = Atoms(
        "He", scaled_positions=[[0.21, 0.34, 0.43]],
        cell=[[4.0, 0.0, 0.0], [0.4, 4.8, 0.0], [0.1, 0.2, 5.1]], pbc=True,
    )
    chart = ReferenceCellCoordinates(reference)
    axes = np.zeros((9, 2))
    axes[0, 0] = 1.0 / np.sqrt(2.0)
    axes[3, 0] = 1.0 / np.sqrt(2.0)  # eta_xx coupled to atomic x
    axes[1, 1] = 1.0
    plane = ModePlane(
        reference=np.zeros(9), basis=np.eye(9), metric_weights=np.ones(9),
        axis_weights=axes, axis_labels=("atomic+strain", "atomic y"),
        amplitude_unit="toy metric", reference_id="toy/nonorthogonal-cell",
    )
    initial = chart.to_atoms(plane.frozen_coordinates([-0.12, 0.05]))
    final = chart.to_atoms(plane.frozen_coordinates([0.12, -0.05]))
    return chart, plane, initial, final


def test_strict_subspace_lift_handles_nonorthogonal_cell_and_vcneb() -> None:
    chart, plane, initial, final = _setup()
    scale = 4.7
    subspace = strict_mode_subspace_for_vcneb(
        plane, chart, initial, final, cell_scale_A=scale,
    )
    assert subspace.mode_basis.shape == (12, 2)
    assert np.allclose(subspace.initial_q, [-0.12, 0.05])
    assert np.allclose(subspace.final_q, [0.12, -0.05])
    assert subspace.endpoint_residual_vcneb_A < 1e-10
    x_initial = _vcneb_x(initial, initial.cell.array, scale)
    coordinates = chart.from_atoms(initial)
    for index in range(2):
        trial = chart.to_atoms(coordinates + 1e-5 * plane.axis_vectors[:, index])
        numerical = (_vcneb_x(trial, initial.cell.array, scale) - x_initial) / 1e-5
        assert np.allclose(numerical, subspace.mode_basis[:, index], atol=1e-9)
    images = interpolate_vcneb(initial, final, 5, align_cells=False)
    chain = VCNEB(
        images, cell_scale=scale, mode_basis=subspace.mode_basis,
        constraint_mode="subspace", climb=False,
    )
    assert chain.n_images == 5


def test_strict_subspace_rejects_off_plane_endpoint_and_periodic_gauge() -> None:
    chart, plane, initial, final = _setup()
    off_plane = final.copy()
    off_plane.positions[0, 2] += 0.04
    with pytest.raises(ValueError, match="does not contain both endpoints"):
        strict_mode_subspace_for_vcneb(
            plane, chart, initial, off_plane, cell_scale_A=4.7,
        )
    different_gauge = final.copy()
    different_gauge.set_scaled_positions(
        different_gauge.get_scaled_positions(wrap=False) + np.array([[1.0, 0.0, 0.0]])
    )
    with pytest.raises(ValueError, match="periodic atom gauge"):
        strict_mode_subspace_for_vcneb(
            plane, chart, initial, different_gauge, cell_scale_A=4.7,
        )


def test_endpoint_completion_adds_explicit_strain_direction_and_connects_vcneb() -> None:
    chart, plane, initial, final = _setup()
    final_coordinates = chart.from_atoms(final)
    final_coordinates[-4] += 0.015  # eta_zz: not included by Q1/Q2
    completed_final = chart.to_atoms(final_coordinates)
    with pytest.raises(ValueError, match="does not contain both endpoints"):
        strict_mode_subspace_for_vcneb(
            plane, chart, initial, completed_final, cell_scale_A=4.7,
        )
    completed = endpoint_completed_mode_subspace_for_vcneb(
        plane, chart, initial, completed_final, cell_scale_A=4.7,
    )
    assert completed.mode_basis.shape == (12, 3)
    assert completed.completion_amplitude == pytest.approx(0.015)
    assert completed.initial_offplane_norm < 1e-12
    assert completed.final_offplane_norm == pytest.approx(0.015)
    assert completed.endpoint_residual_vcneb_A < 1e-10
    assert np.allclose(completed.completion_direction_chart[-6:], [0, 0, 1, 0, 0, 0])
    assert np.allclose(completed.initial_q, [-0.12, 0.05])
    assert np.allclose(completed.final_q, [0.12, -0.05])
    chain = VCNEB(
        interpolate_vcneb(initial, completed_final, 5, align_cells=False),
        cell_scale=4.7, mode_basis=completed.mode_basis,
        constraint_mode="subspace", climb=False,
    )
    assert chain.n_images == 5
    original = endpoint_completed_mode_subspace_for_vcneb(
        plane, chart, initial, final, cell_scale_A=4.7,
    )
    assert original.mode_basis.shape == (12, 2)
    assert original.completion_direction_chart is None


def test_endpoint_completion_refuses_rigid_translation_as_mode() -> None:
    chart, plane, initial, final = _setup()
    translated = final.copy()
    translated.positions[:, 0] += 0.02
    with pytest.raises(ValueError, match="rigid translation"):
        endpoint_completed_mode_subspace_for_vcneb(
            plane, chart, initial, translated, cell_scale_A=4.7,
        )


def test_subspace_artifact_checks_effective_endpoints_and_unwrapped_gauge(tmp_path) -> None:
    chart, plane, initial, final = _setup()
    scale = abs(np.linalg.det(initial.cell.array)) ** (1.0 / 3.0)
    completed = endpoint_completed_mode_subspace_for_vcneb(
        plane, chart, initial, final, cell_scale_A=scale,
    )
    artifact = make_vcneb_subspace_artifact(
        completed, initial, final, source_sha256={"toy_reference": "a" * 64},
    )
    path = tmp_path / "mode_subspace.json"
    path.write_text(json.dumps(artifact), encoding="utf-8")
    basis, recovered = load_vcneb_subspace_artifact(path, initial, final)
    assert recovered["n_directions"] == 2
    assert np.allclose(basis, completed.mode_basis)
    assert not basis.flags.writeable

    periodic_image = final.copy()
    periodic_image.set_scaled_positions(
        periodic_image.get_scaled_positions(wrap=False) + [[1.0, 0.0, 0.0]]
    )
    with pytest.raises(ValueError, match="periodic image gauge"):
        load_vcneb_subspace_artifact(path, initial, periodic_image)
    altered = dict(artifact)
    altered["mode_basis"] = np.asarray(artifact["mode_basis"], dtype=float).copy().tolist()
    altered["mode_basis"][0][0] += 0.01
    path.write_text(json.dumps(altered), encoding="utf-8")
    with pytest.raises(ValueError, match="malformed or changed"):
        load_vcneb_subspace_artifact(path, initial, final)
    for field, replacement in (
        ("initial_q", [float("nan"), 0.0]),
        ("source_sha256", {"toy_reference": "invalid"}),
        ("subspace_kind", "unspecified"),
        ("reference_id", ""),
    ):
        altered = dict(artifact)
        altered[field] = replacement
        path.write_text(json.dumps(altered), encoding="utf-8")
        with pytest.raises(ValueError, match="metadata is malformed"):
            load_vcneb_subspace_artifact(path, initial, final)
    with pytest.raises(ValueError, match="source SHA256"):
        make_vcneb_subspace_artifact(
            completed, initial, final, source_sha256={"toy_reference": "invalid"},
        )
