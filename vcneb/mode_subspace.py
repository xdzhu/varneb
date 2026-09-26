"""Lift an audited two-mode plane into VCNEB's optimizer coordinates.

Only an endpoint-representable *global* two-axis plane is a strict-subspace
path. An off-plane endpoint must not be silently projected: that would change
the requested physical transition and its barrier.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import numpy as np
from ase import Atoms

from .mode_surface import ModePlane
from .provenance import endpoint_structure_record
from .reference_cell import ReferenceCellCoordinates


Array = np.ndarray


def _symmetric_strain_matrix(voigt: Array) -> Array:
    xx, yy, zz, yz, xz, xy = voigt
    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]], dtype=float)


def _vcneb_x(atoms: Atoms, reference_cell: Array, cell_scale: float) -> Array:
    fractional = atoms.get_scaled_positions(wrap=False)
    deformation = np.linalg.solve(reference_cell, atoms.cell.array).T
    return np.r_[
        (fractional @ reference_cell).reshape(-1),
        (cell_scale * (deformation - np.eye(3))).reshape(-1),
    ]


def _lift_directions_to_vcneb(
    directions: Array, chart: ReferenceCellCoordinates, initial_cell: Array, cell_scale: float,
) -> Array:
    """Exact linear chart-to-VCNEB map for atomic/strain direction columns."""

    atomic_count = 3 * chart.n_atoms
    chart_cell = chart.reference_cell
    basis = np.empty((atomic_count + 9, directions.shape[1]), dtype=float)
    for axis in range(directions.shape[1]):
        displacement = directions[:atomic_count, axis].reshape(chart.n_atoms, 3)
        basis[:atomic_count, axis] = (
            displacement @ np.linalg.inv(chart_cell) @ initial_cell
        ).reshape(-1)
        strain = _symmetric_strain_matrix(directions[atomic_count:, axis])
        deformation_direction = np.linalg.solve(initial_cell, chart_cell @ strain).T
        basis[atomic_count:, axis] = (cell_scale * deformation_direction).reshape(-1)
    if not np.all(np.isfinite(basis)) or np.linalg.matrix_rank(basis) != directions.shape[1]:
        raise ValueError("lifted VCNEB mode basis is non-finite or rank deficient")
    return basis


def _require_matching_periodic_gauge(
    chart: ReferenceCellCoordinates, atoms: Atoms, coordinates: Array, name: str,
    tolerance: float,
) -> None:
    aligned = chart.to_atoms(coordinates)
    if not np.allclose(
        atoms.get_scaled_positions(wrap=False),
        aligned.get_scaled_positions(wrap=False), rtol=0.0, atol=tolerance,
    ):
        raise ValueError(f"{name} endpoint uses a different periodic atom gauge; unwrap it first")


@dataclass(frozen=True)
class StrictModeSubspace:
    """Two VCNEB direction columns plus exact endpoint order parameters."""

    mode_basis: Array
    initial_q: Array
    final_q: Array
    reference_id: str
    cell_scale_A: float
    endpoint_residual_vcneb_A: float


@dataclass(frozen=True)
class EndpointCompletedModeSubspace:
    """Q1/Q2 plus explicit endpoint-residual difference, never a phonon label.

    The VCNEB affine-subspace origin is the supplied initial endpoint. The
    third column, when needed, is the metric-normalized off-plane component
    of the final-minus-initial displacement. It guarantees representability
    but does not imply that the completed path is the unrestricted MEP.
    """

    mode_basis: Array
    initial_q: Array
    final_q: Array
    completion_direction_chart: Array | None
    completion_amplitude: float
    amplitude_unit: str
    initial_offplane_norm: float
    final_offplane_norm: float
    reference_id: str
    cell_scale_A: float
    endpoint_residual_vcneb_A: float


def strict_mode_subspace_for_vcneb(
    plane: ModePlane,
    chart: ReferenceCellCoordinates,
    initial: Atoms,
    final: Atoms,
    *,
    cell_scale_A: float,
    endpoint_tolerance: float = 1e-8,
) -> StrictModeSubspace:
    """Convert a globally affine Q1/Q2 plane to ``VCNEB(mode_basis=...)``.

    The reference-cell chart uses ``u`` in Angstrom and symmetric ``eta``;
    VCNEB uses fractional atoms multiplied by the *initial* cell and nine
    length-scaled deformation entries. This exact linear conversion handles
    a nonorthogonal reference cell and a different initial cell. Passing the
    returned basis to ``VCNEB(..., constraint_mode='subspace')`` is valid only
    after this endpoint gate succeeds and all images use the same atom gauge.
    """

    if plane.reference.size != chart.coordinate_count:
        raise ValueError("mode plane and reference-cell chart have different coordinate counts")
    scale = float(cell_scale_A)
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError("cell_scale_A must be finite and positive")
    initial_coordinates = chart.from_atoms(initial)
    final_coordinates = chart.from_atoms(final)
    initial_q, final_q = plane.require_representable_endpoints(
        initial_coordinates, final_coordinates, tolerance=endpoint_tolerance,
    )
    for name, atoms, coordinates in (
        ("initial", initial, initial_coordinates), ("final", final, final_coordinates)
    ):
        _require_matching_periodic_gauge(chart, atoms, coordinates, name, endpoint_tolerance)
    first_cell = np.asarray(initial.cell.array, dtype=float)
    axes = plane.axis_vectors
    if axes.shape != (chart.coordinate_count, 2):
        raise ValueError("mode plane has an invalid pair of direction vectors")
    basis = _lift_directions_to_vcneb(axes, chart, first_cell, scale)
    delta = _vcneb_x(final, first_cell, scale) - _vcneb_x(initial, first_cell, scale)
    residual = float(np.linalg.norm(delta - basis @ np.linalg.lstsq(basis, delta, rcond=None)[0]))
    threshold = endpoint_tolerance * max(1.0, float(np.linalg.norm(delta)))
    if residual > threshold:
        raise ValueError(
            f"lifted strict Q1/Q2 subspace cannot connect endpoints: {residual:.3e} > {threshold:.3e} A"
        )
    for value in (basis, initial_q, final_q):
        value.flags.writeable = False
    return StrictModeSubspace(
        mode_basis=basis, initial_q=initial_q, final_q=final_q,
        reference_id=plane.reference_id, cell_scale_A=scale,
        endpoint_residual_vcneb_A=residual,
    )


def endpoint_completed_mode_subspace_for_vcneb(
    plane: ModePlane,
    chart: ReferenceCellCoordinates,
    initial: Atoms,
    final: Atoms,
    *,
    cell_scale_A: float,
    endpoint_tolerance: float = 1e-8,
) -> EndpointCompletedModeSubspace:
    """Lift Q1/Q2 and the missing endpoint motion into one global subspace.

    Unlike the strict *two*-axis gate, the VCNEB affine origin is ``initial``
    and a third, explicitly labeled residual direction is admitted if the
    endpoint difference has a nonzero metric-orthogonal part. Whole-crystal
    translation is rejected as a mapping error, not silently made a mode.
    This is still a constrained path and its barrier must be compared to an
    independent unrestricted VCNEB path under identical DFT settings.
    """

    if plane.reference.size != chart.coordinate_count:
        raise ValueError("mode plane and reference-cell chart have different coordinate counts")
    scale = float(cell_scale_A)
    tolerance = float(endpoint_tolerance)
    if not np.isfinite(scale) or scale <= 0.0 or not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("cell_scale_A must be positive and endpoint_tolerance nonnegative")
    coordinates = [chart.from_atoms(atoms) for atoms in (initial, final)]
    for name, atoms, value in zip(("initial", "final"), (initial, final), coordinates):
        _require_matching_periodic_gauge(chart, atoms, value, name, tolerance)
    q = [plane.project(value) for value in coordinates]
    residuals = [value - plane.frozen_coordinates(projected)
                 for value, projected in zip(coordinates, q)]
    residual_norms = [float(np.sqrt(np.dot(plane.metric_weights * value, value)))
                      for value in residuals]
    completion = residuals[1] - residuals[0]
    metric = plane.metric_weights
    amplitude = float(np.sqrt(np.dot(metric * completion, completion)))
    directions = plane.axis_vectors
    if np.max(np.abs(directions.T @ (metric * completion))) > tolerance:
        raise ValueError("endpoint completion is not metric-orthogonal to Q1/Q2")
    translation = chart.rigid_translation_directions()
    translation_gram = translation.T @ (metric[:, None] * translation)
    translation_coefficients = np.linalg.solve(
        translation_gram, translation.T @ (metric * completion),
    )
    translation_norm = float(np.sqrt(translation_coefficients @ translation_gram @ translation_coefficients))
    if translation_norm > tolerance:
        raise ValueError(
            f"endpoint completion contains rigid translation {translation_norm:.3e} {plane.amplitude_unit}; "
            "align endpoint origins before defining the mode subspace"
        )
    completion_direction = None
    if amplitude > tolerance:
        completion_direction = completion / amplitude
        directions = np.column_stack([directions, completion_direction])
    first_cell = np.asarray(initial.cell.array, dtype=float)
    basis = _lift_directions_to_vcneb(directions, chart, first_cell, scale)
    delta = _vcneb_x(final, first_cell, scale) - _vcneb_x(initial, first_cell, scale)
    endpoint_error = float(np.linalg.norm(delta - basis @ np.linalg.lstsq(basis, delta, rcond=None)[0]))
    threshold = tolerance * max(1.0, float(np.linalg.norm(delta)))
    if endpoint_error > threshold:
        raise ValueError(
            f"endpoint-completed subspace still cannot connect endpoints: {endpoint_error:.3e} > {threshold:.3e} A"
        )
    for value in (basis, *q):
        value.flags.writeable = False
    if completion_direction is not None:
        completion_direction.flags.writeable = False
    return EndpointCompletedModeSubspace(
        mode_basis=basis, initial_q=q[0], final_q=q[1],
        completion_direction_chart=completion_direction,
        completion_amplitude=amplitude, amplitude_unit=plane.amplitude_unit,
        initial_offplane_norm=residual_norms[0], final_offplane_norm=residual_norms[1],
        reference_id=plane.reference_id, cell_scale_A=scale,
        endpoint_residual_vcneb_A=endpoint_error,
    )


def _array_sha256(values: Array) -> str:
    return hashlib.sha256(np.ascontiguousarray(values, dtype="<f8").tobytes()).hexdigest()


def _valid_source_sha256(value: object) -> bool:
    return isinstance(value, dict) and bool(value) and all(
        isinstance(key, str) and bool(key)
        and isinstance(digest, str) and len(digest) == 64
        and all(char in "0123456789abcdef" for char in digest.lower())
        for key, digest in value.items()
    )


def _valid_q_pair(value: object) -> bool:
    try:
        coordinates = np.asarray(value, dtype=float)
    except (TypeError, ValueError):
        return False
    return coordinates.shape == (2,) and bool(np.all(np.isfinite(coordinates)))


def make_vcneb_subspace_artifact(
    subspace: StrictModeSubspace | EndpointCompletedModeSubspace,
    initial: Atoms,
    final: Atoms,
    *,
    source_sha256: dict[str, str],
) -> dict:
    """Serialize a global VCNEB subspace with exact endpoint/gauge identity.

    ``endpoint_structure_record`` identifies the physical ordered endpoints;
    a second hash of unwrapped VCNEB coordinates identifies the periodic image
    gauge, which the wrapped structural hash deliberately ignores.
    """

    if not _valid_source_sha256(source_sha256):
        raise ValueError("subspace artifact requires named source SHA256 digests")
    basis = np.asarray(subspace.mode_basis, dtype=float)
    scale = float(subspace.cell_scale_A)
    if (basis.ndim != 2 or basis.shape[0] != 3 * len(initial) + 9
            or len(initial) != len(final) or basis.shape[1] < 1
            or not np.all(np.isfinite(basis)) or np.linalg.matrix_rank(basis) != basis.shape[1]
            or not np.isfinite(scale) or scale <= 0.0):
        raise ValueError("subspace basis/cell scale is invalid for the supplied endpoints")
    first_x = _vcneb_x(initial, initial.cell.array, scale)
    final_x = _vcneb_x(final, initial.cell.array, scale)
    delta = final_x - first_x
    residual = float(np.linalg.norm(delta - basis @ np.linalg.lstsq(basis, delta, rcond=None)[0]))
    if residual > 1e-8 * max(1.0, float(np.linalg.norm(delta))):
        raise ValueError("subspace artifact cannot represent both supplied endpoints")
    subspace_kind = ("endpoint_completed" if isinstance(subspace, EndpointCompletedModeSubspace)
                     else "strict_mode_plane")
    if (not isinstance(subspace.reference_id, str) or not subspace.reference_id
            or not _valid_q_pair(subspace.initial_q) or not _valid_q_pair(subspace.final_q)
            or (subspace_kind == "strict_mode_plane" and basis.shape[1] != 2)
            or (subspace_kind == "endpoint_completed" and basis.shape[1] not in (2, 3))):
        raise ValueError("subspace order-parameter metadata is malformed")
    return {
        "schema_version": 1,
        "kind": "varneb_global_vcneb_mode_subspace",
        "constraint_mode": "subspace",
        "subspace_kind": subspace_kind,
        "reference_id": subspace.reference_id,
        "cell_scale_A": scale,
        "coordinate_convention": "vcneb_initial_cell_fractional_atoms_and_3x3_deformation_scaled_A",
        "mode_basis": basis.tolist(),
        "mode_basis_sha256": _array_sha256(basis),
        "n_directions": basis.shape[1],
        "initial_q": subspace.initial_q.tolist(),
        "final_q": subspace.final_q.tolist(),
        "endpoint_structure_sha256": {
            "initial": endpoint_structure_record(initial)["sha256"],
            "final": endpoint_structure_record(final)["sha256"],
        },
        "endpoint_unwrapped_coordinate_sha256": {
            "initial": _array_sha256(np.round(first_x, 12)),
            "final": _array_sha256(np.round(final_x, 12)),
        },
        "endpoint_residual_vcneb_A": residual,
        "source_sha256": dict(source_sha256),
        "interpretation": "global constrained VCNEB subspace; not an unrestricted MEP or barrier",
    }


def load_vcneb_subspace_artifact(
    path: str | Path, initial: Atoms, final: Atoms,
) -> tuple[Array, dict]:
    """Reject stale or different-gauge artifacts before any calculator call."""

    artifact = json.loads(Path(path).read_text(encoding="utf-8"))
    if (artifact.get("schema_version") != 1
            or artifact.get("kind") != "varneb_global_vcneb_mode_subspace"
            or artifact.get("constraint_mode") != "subspace"
            or artifact.get("coordinate_convention") !=
            "vcneb_initial_cell_fractional_atoms_and_3x3_deformation_scaled_A"):
        raise ValueError("unrecognized global VCNEB mode-subspace artifact")
    kind = artifact.get("subspace_kind")
    direction_count = artifact.get("n_directions")
    if (kind not in ("strict_mode_plane", "endpoint_completed")
            or type(direction_count) is not int
            or (kind == "strict_mode_plane" and direction_count != 2)
            or (kind == "endpoint_completed" and direction_count not in (2, 3))
            or not isinstance(artifact.get("reference_id"), str)
            or not artifact["reference_id"]
            or not _valid_q_pair(artifact.get("initial_q"))
            or not _valid_q_pair(artifact.get("final_q"))
            or not _valid_source_sha256(artifact.get("source_sha256"))):
        raise ValueError("mode-subspace order-parameter/provenance metadata is malformed")
    basis = np.asarray(artifact["mode_basis"], dtype=float)
    scale = float(artifact["cell_scale_A"])
    if (basis.ndim != 2 or basis.shape != (3 * len(initial) + 9, artifact["n_directions"])
            or not np.all(np.isfinite(basis)) or np.linalg.matrix_rank(basis) != basis.shape[1]
            or _array_sha256(basis) != artifact.get("mode_basis_sha256")
            or not np.isfinite(scale) or scale <= 0.0):
        raise ValueError("mode-subspace basis is malformed or changed")
    records = {
        "initial": endpoint_structure_record(initial)["sha256"],
        "final": endpoint_structure_record(final)["sha256"],
    }
    if records != artifact.get("endpoint_structure_sha256"):
        raise ValueError("mode-subspace physical endpoints differ from the artifact")
    first_x = _vcneb_x(initial, initial.cell.array, scale)
    final_x = _vcneb_x(final, initial.cell.array, scale)
    gauge = {
        "initial": _array_sha256(np.round(first_x, 12)),
        "final": _array_sha256(np.round(final_x, 12)),
    }
    if gauge != artifact.get("endpoint_unwrapped_coordinate_sha256"):
        raise ValueError("mode-subspace endpoints use a different periodic image gauge")
    delta = final_x - first_x
    residual = float(np.linalg.norm(delta - basis @ np.linalg.lstsq(basis, delta, rcond=None)[0]))
    if residual > 1e-8 * max(1.0, float(np.linalg.norm(delta))):
        raise ValueError("mode-subspace cannot represent the effective endpoints")
    stored_residual = float(artifact["endpoint_residual_vcneb_A"])
    if not np.isfinite(stored_residual) or not np.isclose(
        residual, stored_residual, atol=1e-8, rtol=0.0,
    ):
        raise ValueError("mode-subspace endpoint residual changed")
    basis.flags.writeable = False
    return basis, artifact


__all__ = [
    "StrictModeSubspace", "strict_mode_subspace_for_vcneb",
    "EndpointCompletedModeSubspace", "endpoint_completed_mode_subspace_for_vcneb",
    "make_vcneb_subspace_artifact", "load_vcneb_subspace_artifact",
]
