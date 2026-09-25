"""Variable-cell NEB in an ASE-compatible extended coordinate space.

The implementation optimizes fractional atomic coordinates together with a
cell deformation gradient.  It deliberately stays calculator-agnostic: any ASE
calculator that can provide energy, forces, and stress can be used.
"""

from __future__ import annotations

from dataclasses import dataclass
import inspect
import json
import os
from pathlib import Path
import tempfile
from typing import Callable, Iterable, Iterator, Mapping, Optional, Sequence

import numpy as np

from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write
from ase.io.trajectory import Trajectory
from ase.optimize import BFGS, FIRE, LBFGS
from ase.optimize.bfgslinesearch import BFGSLineSearch
from ase.parallel import world
from ase.units import GPa

from .calculator import calculator_context, classify_calculator_failure, validate_image_calculators
from .executor import ImageEvaluation
from .step_control import CandidateStepRejected, CheckedFIRE
from .block_fire import BlockFIRE
from .scaled_fire import ImageScaledFIRE, StagedFIRE


Array = np.ndarray


def _calculator_diagnostic_paths(images: Sequence[Atoms]) -> list[Path]:
    """Return existing tails worth reading when an external calculator fails."""

    relative_names = (
        Path("OUT.ABACUS") / "running_scf.log",
        Path("running_scf.log"),
        Path("OUTCAR"),
        Path("vasp.out"),
        Path("vasp.err"),
        Path("stdout"),
        Path("stderr"),
    )
    paths: list[Path] = []
    seen: set[Path] = set()
    for image in images:
        directory = getattr(getattr(image, "calc", None), "directory", None)
        if directory is None:
            continue
        root = Path(directory)
        for relative in relative_names:
            path = root / relative
            if path.is_file() and path not in seen:
                seen.add(path)
                paths.append(path)
    return paths


@dataclass
class VCNEBState:
    """Extended coordinates for one image.

    ``q`` stores fractional coordinates without wrapping.  ``deform`` is the
    ASE/UnitCellFilter-style deformation gradient, with the current cell given
    by ``cell0 @ deform.T``.
    """

    q: Array
    deform: Array


def cell_matrix(atoms: Atoms) -> Array:
    return atoms.cell.array.copy()


def _validate_cell_matrix(cell: Array, *, context: str = "cell") -> None:
    matrix = np.asarray(cell, dtype=float)
    determinant = float(np.linalg.det(matrix)) if matrix.shape == (3, 3) else 0.0
    if matrix.shape != (3, 3) or not np.all(np.isfinite(matrix)) or determinant <= 1e-12:
        raise ValueError(f"{context} must be a finite 3x3 cell with positive determinant")


def deformation_from_cell(cell: Array, reference_cell: Array) -> Array:
    """Return deformation gradient F for ASE row-vector cells.

    ASE stores the three cell vectors as rows and Cartesian positions are
    ``scaled @ cell``.  The UnitCellFilter convention is
    ``F = solve(cell0, cell).T`` and ``cell = cell0 @ F.T``.
    """

    return np.linalg.solve(reference_cell, cell).T


def cell_from_deformation(deform: Array, reference_cell: Array) -> Array:
    return reference_cell @ deform.T


def state_from_atoms(atoms: Atoms, reference_cell: Array) -> VCNEBState:
    return VCNEBState(
        q=atoms.get_scaled_positions(wrap=False).copy(),
        deform=deformation_from_cell(cell_matrix(atoms), reference_cell),
    )


def apply_state(
    atoms: Atoms,
    state: VCNEBState,
    reference_cell: Array,
    *,
    wrap_positions: bool = False,
) -> None:
    atoms.set_cell(cell_from_deformation(state.deform, reference_cell), scale_atoms=False)
    atoms.set_scaled_positions(state.q)
    if wrap_positions:
        atoms.wrap(eps=1e-12)


def polar_rotation(matrix: Array) -> Array:
    u, _, vt = np.linalg.svd(matrix)
    rot = u @ vt
    if np.linalg.det(rot) < 0.0:
        vt[-1] *= -1.0
        rot = u @ vt
    return rot


def remove_global_rotation(reference: Atoms, atoms: Atoms) -> None:
    """Rotate ``atoms`` into the gauge closest to ``reference``."""

    ref_cell = cell_matrix(reference)
    deform = deformation_from_cell(cell_matrix(atoms), ref_cell)
    rot = polar_rotation(deform)
    atoms.set_cell(cell_matrix(atoms) @ rot, scale_atoms=False)
    atoms.set_positions(atoms.positions @ rot)


def _fractional_delta(a: Array, b: Array, pbc: Sequence[bool], mic: bool) -> Array:
    delta = np.asarray(a) - np.asarray(b)
    if mic:
        pbc_arr = np.asarray(pbc, dtype=bool)
        delta[:, pbc_arr] -= np.rint(delta[:, pbc_arr])
    return delta


def _minimum_cost_assignment(cost: Array) -> list[int]:
    """Solve a finite square assignment problem without a SciPy dependency."""

    matrix = np.asarray(cost, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or matrix.shape[0] == 0:
        raise ValueError("assignment cost must be a non-empty square matrix")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("assignment cost contains non-finite entries")

    size = matrix.shape[0]
    row_potentials = np.zeros(size + 1)
    column_potentials = np.zeros(size + 1)
    matched_row = np.zeros(size + 1, dtype=int)
    predecessor = np.zeros(size + 1, dtype=int)
    for row in range(1, size + 1):
        matched_row[0] = row
        column = 0
        minimum = np.full(size + 1, np.inf)
        used = np.zeros(size + 1, dtype=bool)
        while True:
            used[column] = True
            current_row = matched_row[column]
            delta = np.inf
            next_column = 0
            for candidate in range(1, size + 1):
                if used[candidate]:
                    continue
                reduced = (
                    matrix[current_row - 1, candidate - 1]
                    - row_potentials[current_row]
                    - column_potentials[candidate]
                )
                if reduced < minimum[candidate]:
                    minimum[candidate] = reduced
                    predecessor[candidate] = column
                if minimum[candidate] < delta:
                    delta = minimum[candidate]
                    next_column = candidate
            if not np.isfinite(delta):
                raise ValueError("atom mapping assignment has no finite solution")
            for candidate in range(size + 1):
                if used[candidate]:
                    row_potentials[matched_row[candidate]] += delta
                    column_potentials[candidate] -= delta
                else:
                    minimum[candidate] -= delta
            column = next_column
            if matched_row[column] == 0:
                break
        while True:
            previous = predecessor[column]
            matched_row[column] = matched_row[previous]
            column = previous
            if column == 0:
                break

    assignment = [0] * size
    for column in range(1, size + 1):
        assignment[matched_row[column] - 1] = column - 1
    return assignment


def _mapping_cost_matrix(
    initial: Atoms,
    final: Atoms,
    initial_indices: Sequence[int],
    final_indices: Sequence[int],
    *,
    mic: bool,
    translation: Array | None = None,
) -> Array:
    q0 = initial.get_scaled_positions(wrap=False)
    q1 = final.get_scaled_positions(wrap=False)
    if translation is not None:
        q1 = q1 + np.asarray(translation, dtype=float)
    distance_cell = 0.5 * (cell_matrix(initial) + cell_matrix(final))
    costs = np.zeros((len(initial_indices), len(final_indices)), dtype=float)
    for row, initial_index in enumerate(initial_indices):
        for column, final_index in enumerate(final_indices):
            delta = _fractional_delta(
                q1[final_index : final_index + 1],
                q0[initial_index : initial_index + 1],
                initial.pbc,
                mic,
            )[0]
            # A tiny stable tie-breaker makes repeated coordinates deterministic.
            costs[row, column] = float(np.linalg.norm(delta @ distance_cell)) + 1e-12 * column
    return costs


def _translation_candidates(
    initial: Atoms,
    final: Atoms,
    pairs: Sequence[tuple[int, int]],
) -> list[Array]:
    """Generate deterministic fractional translations from matched atom pairs."""

    q0 = initial.get_scaled_positions(wrap=False)
    q1 = final.get_scaled_positions(wrap=False)
    pbc = np.asarray(initial.pbc, dtype=bool)
    candidates = [np.zeros(3, dtype=float)]
    for initial_index, final_index in pairs:
        shift = q0[initial_index] - q1[final_index]
        shift = np.asarray(shift, dtype=float)
        shift[pbc] -= np.floor(shift[pbc])
        shift[~pbc] = 0.0
        if not any(np.allclose(shift, previous, rtol=0.0, atol=1e-12) for previous in candidates):
            candidates.append(shift)
    return candidates


def _translation_cost(
    initial: Atoms,
    final: Atoms,
    mapping: Sequence[int],
    translation: Array,
    *,
    mic: bool,
) -> float:
    q0 = initial.get_scaled_positions(wrap=False)
    q1 = final.get_scaled_positions(wrap=False) + np.asarray(translation, dtype=float)
    distance_cell = 0.5 * (cell_matrix(initial) + cell_matrix(final))
    total = 0.0
    for initial_index, final_index in enumerate(mapping):
        delta = _fractional_delta(
            q1[final_index : final_index + 1],
            q0[initial_index : initial_index + 1],
            initial.pbc,
            mic,
        )[0]
        total += float(np.linalg.norm(delta @ distance_cell))
    return total


def _infer_mapping_and_translation(
    initial: Atoms,
    final: Atoms,
    *,
    mic: bool,
) -> tuple[list[int], Array]:
    """Infer atom assignment and a common periodic translation together."""

    symbols0 = initial.get_chemical_symbols()
    symbols1 = final.get_chemical_symbols()
    if sorted(symbols0) != sorted(symbols1):
        raise ValueError("Initial and final structures do not contain the same elements")
    pairs = [
        (initial_index, final_index)
        for initial_index, symbol0 in enumerate(symbols0)
        for final_index, symbol1 in enumerate(symbols1)
        if symbol0 == symbol1
    ]
    best_cost = np.inf
    best_mapping: list[int] | None = None
    best_translation: Array | None = None
    for translation in _translation_candidates(initial, final, pairs):
        mapping = [-1] * len(initial)
        total = 0.0
        for symbol in sorted(set(symbols0)):
            initial_indices = [index for index, value in enumerate(symbols0) if value == symbol]
            final_indices = [index for index, value in enumerate(symbols1) if value == symbol]
            costs = _mapping_cost_matrix(
                initial,
                final,
                initial_indices,
                final_indices,
                mic=mic,
                translation=translation,
            )
            assignment = _minimum_cost_assignment(costs)
            for row, final_column in enumerate(assignment):
                mapping[initial_indices[row]] = final_indices[final_column]
                total += float(costs[row, final_column])
        if total < best_cost - 1e-12:
            best_cost = total
            best_mapping = mapping
            best_translation = translation.copy()
    if best_mapping is None or best_translation is None:
        raise ValueError("could not infer a finite atom mapping and translation")
    return best_mapping, best_translation


def _infer_translation_for_mapping(
    initial: Atoms,
    final: Atoms,
    mapping: Sequence[int],
    *,
    mic: bool,
) -> Array:
    pairs = list(enumerate(mapping))
    candidates = _translation_candidates(initial, final, pairs)
    return min(
        candidates,
        key=lambda translation: _translation_cost(initial, final, mapping, translation, mic=mic),
    ).copy()


def infer_atom_mapping(initial: Atoms, final: Atoms, *, mic: bool = True) -> list[int]:
    """Infer a final-atom permutation by element and periodic Cartesian distance.

    The result is indexed by the initial atom order: ``mapping[i]`` is the
    corresponding atom index in ``final``.  This is a geometry-based heuristic,
    not a proof of chemical identity; inspect the returned report or provide an
    explicit permutation when a phase transition contains substantial atom
    rearrangement.
    """

    if len(initial) != len(final):
        raise ValueError("Initial and final structures have different atom counts")
    initial_symbols = initial.get_chemical_symbols()
    final_symbols = final.get_chemical_symbols()
    if sorted(initial_symbols) != sorted(final_symbols):
        raise ValueError("Initial and final structures do not contain the same elements")
    mapping = [-1] * len(initial)
    for symbol in sorted(set(initial_symbols)):
        initial_indices = [index for index, value in enumerate(initial_symbols) if value == symbol]
        final_indices = [index for index, value in enumerate(final_symbols) if value == symbol]
        assignment = _minimum_cost_assignment(
            _mapping_cost_matrix(initial, final, initial_indices, final_indices, mic=mic)
        )
        for row, final_column in enumerate(assignment):
            mapping[initial_indices[row]] = final_indices[final_column]
    return mapping


def validate_atom_mapping(
    initial: Atoms,
    final: Atoms,
    mapping: Sequence[int] | str | None = None,
    *,
    mic: bool = True,
    maximum_displacement: float | None = None,
) -> dict:
    """Validate and report an explicit or geometry-inferred atom mapping."""

    if mapping is None:
        resolved = list(range(len(initial)))
        source = "identity"
    elif isinstance(mapping, str):
        if mapping.lower() != "auto":
            raise ValueError("mapping string must be 'auto'")
        resolved = infer_atom_mapping(initial, final, mic=mic)
        source = "auto"
    else:
        resolved = [int(value) for value in mapping]
        source = "explicit"
    if len(resolved) != len(initial) or sorted(resolved) != list(range(len(final))):
        raise ValueError("mapping must be a permutation of all final atom indices")
    if initial.get_chemical_symbols() != [final.get_chemical_symbols()[index] for index in resolved]:
        raise ValueError("mapping pairs atoms with different chemical elements")
    if maximum_displacement is not None and maximum_displacement <= 0.0:
        raise ValueError("maximum_displacement must be positive when provided")

    q0 = initial.get_scaled_positions(wrap=False)
    q1 = final.get_scaled_positions(wrap=False)
    distance_cell = 0.5 * (cell_matrix(initial) + cell_matrix(final))
    records = []
    for initial_index, final_index in enumerate(resolved):
        delta = _fractional_delta(
            q1[final_index : final_index + 1],
            q0[initial_index : initial_index + 1],
            initial.pbc,
            mic,
        )[0]
        distance = float(np.linalg.norm(delta @ distance_cell))
        records.append(
            {
                "initial_index": initial_index,
                "final_index": final_index,
                "element": initial.get_chemical_symbols()[initial_index],
                "fractional_delta": delta.tolist(),
                "distance_A": distance,
            }
        )
    maximum = max((record["distance_A"] for record in records), default=0.0)
    issues = []
    if maximum_displacement is not None and maximum > maximum_displacement:
        issues.append(
            f"maximum mapped displacement {maximum:.6g} A exceeds "
            f"{maximum_displacement:.6g} A"
        )
    return {
        "source": source,
        "mapping": resolved,
        "maximum_displacement_A": maximum,
        "maximum_displacement_threshold_A": maximum_displacement,
        "records": records,
        "issues": issues,
        "valid": not issues,
    }


def _symmetric_matrix_log(matrix: Array, *, context: str) -> Array:
    """Return the real matrix logarithm for a symmetric positive matrix."""

    value = np.asarray(matrix, dtype=float)
    if value.shape != (3, 3) or not np.all(np.isfinite(value)):
        raise ValueError(f"{context} must be a finite 3x3 matrix")
    if not np.allclose(value, value.T, rtol=1e-10, atol=1e-10):
        raise ValueError(
            f"{context} must be symmetric for log_strain interpolation; "
            "use align_cells=True or provide a custom interpolator"
        )
    eigenvalues, eigenvectors = np.linalg.eigh(value)
    if np.any(eigenvalues <= 1e-12):
        raise ValueError(
            f"{context} must be positive definite for log_strain interpolation"
        )
    return (eigenvectors * np.log(eigenvalues)) @ eigenvectors.T


def _symmetric_matrix_exp(matrix: Array) -> Array:
    eigenvalues, eigenvectors = np.linalg.eigh(np.asarray(matrix, dtype=float))
    return (eigenvectors * np.exp(eigenvalues)) @ eigenvectors.T


def _interpolate_deformation(
    lam: float,
    deform0: Array,
    deform1: Array,
    *,
    strategy: str | Callable[[float, Array, Array], Array],
    log_deform0: Array | None = None,
    log_deform1: Array | None = None,
) -> Array:
    """Interpolate two deformation gradients with an explicit strategy."""

    if lam <= 0.0:
        return np.asarray(deform0, dtype=float).copy()
    if lam >= 1.0:
        return np.asarray(deform1, dtype=float).copy()
    if callable(strategy):
        value = strategy(
            float(lam),
            np.asarray(deform0, dtype=float).copy(),
            np.asarray(deform1, dtype=float).copy(),
        )
        result = np.asarray(value, dtype=float)
        if result.shape != (3, 3) or not np.all(np.isfinite(result)):
            raise ValueError("custom cell interpolation must return a finite 3x3 matrix")
        return result

    name = str(strategy).lower().replace("-", "_")
    if name in {"linear", "affine"}:
        return (1.0 - lam) * deform0 + lam * deform1
    if name in {"log", "log_strain", "logarithmic"}:
        if log_deform0 is None or log_deform1 is None:
            raise RuntimeError("log_strain interpolation was not initialized")
        return _symmetric_matrix_exp((1.0 - lam) * log_deform0 + lam * log_deform1)
    raise ValueError(
        "cell_interpolation must be 'linear', 'log_strain', or a callable"
    )


def interpolate_vcneb(
    initial: Atoms,
    final: Atoms,
    n_images: int,
    *,
    align_cells: bool = True,
    align_translation: bool = False,
    mic: bool = False,
    wrap_positions: bool = False,
    cell_interpolation: str | Callable[[float, Array, Array], Array] = "linear",
    mapping: Sequence[int] | str | None = None,
    minimum_distance: float | None = None,
    maximum_deformation: float | None = None,
) -> list[Atoms]:
    """Create an initial variable-cell band, including both endpoints.

    ``cell_interpolation`` selects the path in deformation-gradient space.
    ``linear`` is the historical default.  ``log_strain`` interpolates the
    matrix logarithm of symmetric-positive deformation gradients; this is
    well-defined after the default global-rotation alignment.  A callable is
    given ``(lambda, deform0, deform1)`` and must return a finite 3x3 matrix.
    ``mapping`` is ``None`` for the backward-compatible identity mapping,
    ``"auto"`` for an element-grouped geometry assignment, or a permutation
    indexed by the initial atom order.  With ``align_translation=True``, a
    common periodic endpoint translation is optimized together with an
    automatic mapping (or separately for an explicit mapping).
    ``minimum_distance`` and ``maximum_deformation`` are optional preflight
    thresholds.  They are useful for rejecting an unphysical initial path
    before any calculator is invoked; the default keeps the historical
    interpolation behavior unchanged.
    """

    if n_images < 2:
        raise ValueError("n_images must include endpoints and be at least 2")
    if len(initial) != len(final):
        raise ValueError("Initial and final structures have different atom counts")
    _validate_cell_matrix(cell_matrix(initial), context="initial cell")
    _validate_cell_matrix(cell_matrix(final), context="final cell")

    first = initial.copy()
    last_source = final.copy()
    if align_cells:
        remove_global_rotation(first, last_source)

    if align_translation and (mapping is None or str(mapping).lower() == "auto"):
        resolved_mapping, translation = _infer_mapping_and_translation(
            first,
            last_source,
            mic=mic,
        )
    else:
        mapping_report = validate_atom_mapping(first, last_source, mapping, mic=mic)
        resolved_mapping = mapping_report["mapping"]
        translation = (
            _infer_translation_for_mapping(
                first,
                last_source,
                resolved_mapping,
                mic=mic,
            )
            if align_translation
            else np.zeros(3, dtype=float)
        )
    last = last_source[resolved_mapping].copy()
    if align_translation:
        last.set_scaled_positions(last.get_scaled_positions(wrap=False) + translation)

    reference_cell = cell_matrix(first)
    q0 = first.get_scaled_positions(wrap=False)
    q1 = last.get_scaled_positions(wrap=False)
    dq = _fractional_delta(q1, q0, first.pbc, mic)
    deform0 = deformation_from_cell(cell_matrix(first), reference_cell)
    deform1 = deformation_from_cell(cell_matrix(last), reference_cell)
    log_deform0 = None
    log_deform1 = None
    if not callable(cell_interpolation) and str(cell_interpolation).lower().replace("-", "_") in {
        "log",
        "log_strain",
        "logarithmic",
    }:
        log_deform0 = _symmetric_matrix_log(deform0, context="initial deformation")
        log_deform1 = _symmetric_matrix_log(deform1, context="final deformation")

    path_metadata = {
        "mapping": [int(value) for value in resolved_mapping],
        "translation_fractional_final_cell": [float(value) for value in translation],
        "align_translation": bool(align_translation),
        "cell_interpolation": "callable" if callable(cell_interpolation) else str(cell_interpolation),
        "mic": bool(mic),
    }
    images = []
    for index in range(n_images):
        lam = index / (n_images - 1)
        image = first.copy()
        state = VCNEBState(
            q=q0 + lam * dq,
            deform=_interpolate_deformation(
                lam,
                deform0,
                deform1,
                strategy=cell_interpolation,
                log_deform0=log_deform0,
                log_deform1=log_deform1,
            ),
        )
        _validate_cell_matrix(
            cell_from_deformation(state.deform, reference_cell),
            context=f"interpolated cell at image {index}",
        )
        apply_state(image, state, reference_cell, wrap_positions=wrap_positions)
        image.info["vcneb_path_metadata"] = path_metadata.copy()
        images.append(image)
    if minimum_distance is not None or maximum_deformation is not None:
        validate_path_geometry(
            images,
            minimum_distance=minimum_distance,
            maximum_deformation=maximum_deformation,
        )
    return images


def path_geometry_diagnostics(
    images: Sequence[Atoms],
    *,
    reference_cell: Array | None = None,
    minimum_distance: float | None = None,
    maximum_deformation: float | None = None,
    cell_scale: float | None = None,
    fold_cosine_threshold: float | None = None,
    minimum_endpoint_separation: float | None = None,
) -> dict:
    """Report image volumes, cell deformation, and MIC atom separations.

    The report is calculator-independent.  ``minimum_distance`` is compared
    with the shortest distinct-atom distance under periodic MIC.  For a path
    with fewer than two atoms, the minimum distance is reported as infinity.
    ``maximum_deformation`` is a Frobenius-norm threshold on ``F-I`` relative
    to ``reference_cell``.  If ``fold_cosine_threshold`` is supplied, adjacent
    extended-coordinate segments with a cosine below that threshold are
    reported as path folds and make validation fail.  The optional endpoint
    separation is measured after the caller's atom mapping and translation
    alignment, in the same atom-plus-cell coordinates used by VC-NEB.  It
    catches collapsed endpoints but does not by itself establish phase identity.
    """

    if len(images) < 2:
        raise ValueError("At least two images are required for path geometry diagnostics")
    if minimum_distance is not None and minimum_distance <= 0.0:
        raise ValueError("minimum_distance must be positive when provided")
    if maximum_deformation is not None and maximum_deformation <= 0.0:
        raise ValueError("maximum_deformation must be positive when provided")
    if cell_scale is not None and cell_scale <= 0.0:
        raise ValueError("cell_scale must be positive when provided")
    if fold_cosine_threshold is not None and not -1.0 <= fold_cosine_threshold <= 1.0:
        raise ValueError("fold_cosine_threshold must be between -1 and 1")
    if minimum_endpoint_separation is not None and minimum_endpoint_separation <= 0.0:
        raise ValueError("minimum_endpoint_separation must be positive when provided")
    reference = (
        cell_matrix(images[0]) if reference_cell is None else np.asarray(reference_cell, dtype=float)
    )
    _validate_cell_matrix(reference, context="reference cell")
    coordinate_scale = (
        float(abs(np.linalg.det(reference)) ** (1.0 / 3.0))
        if cell_scale is None
        else float(cell_scale)
    )

    records = []
    issues = []
    extended_coordinates = []
    for image_index, image in enumerate(images):
        cell = cell_matrix(image)
        _validate_cell_matrix(cell, context=f"image {image_index} cell")
        deform = deformation_from_cell(cell, reference)
        q = image.get_scaled_positions(wrap=False)
        extended_coordinates.append(
            np.concatenate(
                [
                    (q @ reference).reshape(-1),
                    (coordinate_scale * (deform - np.eye(3))).reshape(-1),
                ]
            )
        )
        minimum = float("inf")
        for atom_index in range(len(image)):
            for other_index in range(atom_index + 1, len(image)):
                distance = float(image.get_distance(atom_index, other_index, mic=True))
                minimum = min(minimum, distance)
        deformation_norm = float(np.linalg.norm(deform - np.eye(3)))
        record = {
            "image_index": image_index,
            "volume_A3": float(image.get_volume()),
            "minimum_interatomic_distance_A": minimum,
            "deformation_from_reference_frobenius": deformation_norm,
        }
        if minimum_distance is not None and minimum < minimum_distance:
            issues.append(
                f"image {image_index} minimum interatomic distance {minimum:.6g} A "
                f"is below {minimum_distance:.6g} A"
            )
        if maximum_deformation is not None and deformation_norm > maximum_deformation:
            issues.append(
                f"image {image_index} deformation norm {deformation_norm:.6g} "
                f"exceeds {maximum_deformation:.6g}"
            )
        records.append(record)
    segment_lengths = []
    for left, right in zip(extended_coordinates[:-1], extended_coordinates[1:]):
        segment_lengths.append(float(np.linalg.norm(right - left)))
    endpoint_separation = float(np.linalg.norm(extended_coordinates[-1] - extended_coordinates[0]))
    if (
        minimum_endpoint_separation is not None
        and endpoint_separation < minimum_endpoint_separation
    ):
        issues.append(
            f"endpoint extended-coordinate separation {endpoint_separation:.6g} A "
            f"is below {minimum_endpoint_separation:.6g} A; "
            "check that the requested phases did not collapse to the same structure"
        )
    adjacent_cosines = []
    zero_length_segments = []
    folded_junctions = []
    for junction in range(1, len(extended_coordinates) - 1):
        left = extended_coordinates[junction] - extended_coordinates[junction - 1]
        right = extended_coordinates[junction + 1] - extended_coordinates[junction]
        left_norm = float(np.linalg.norm(left))
        right_norm = float(np.linalg.norm(right))
        if left_norm <= 1e-14 or right_norm <= 1e-14:
            adjacent_cosines.append(None)
            if left_norm <= 1e-14:
                zero_length_segments.append(junction - 1)
            if right_norm <= 1e-14:
                zero_length_segments.append(junction)
            continue
        cosine = float(np.dot(left, right) / (left_norm * right_norm))
        adjacent_cosines.append(cosine)
        if fold_cosine_threshold is not None and cosine < fold_cosine_threshold:
            folded_junctions.append(junction)
    zero_length_segments = sorted(set(zero_length_segments))
    if fold_cosine_threshold is not None:
        for junction in folded_junctions:
            issues.append(
                f"path fold at image {junction}: adjacent extended-coordinate "
                f"segment cosine {adjacent_cosines[junction - 1]:.6g} is below "
                f"{fold_cosine_threshold:.6g}"
            )
        for segment in zero_length_segments:
            issues.append(f"zero-length extended-coordinate segment between images {segment} and {segment + 1}")
    return {
        "n_images": len(images),
        "minimum_distance_threshold_A": minimum_distance,
        "maximum_deformation_threshold": maximum_deformation,
        "cell_scale_A": coordinate_scale,
        "images": records,
        "segment_lengths_A": segment_lengths,
        "endpoint_separation_A": endpoint_separation,
        "minimum_endpoint_separation_A": minimum_endpoint_separation,
        "adjacent_segment_cosines": adjacent_cosines,
        "zero_length_segments": zero_length_segments,
        "fold_cosine_threshold": fold_cosine_threshold,
        "folded_junctions": folded_junctions,
        "issues": issues,
        "valid": not issues,
    }


def validate_path_geometry(
    images: Sequence[Atoms],
    *,
    reference_cell: Array | None = None,
    minimum_distance: float | None = None,
    maximum_deformation: float | None = None,
    cell_scale: float | None = None,
    fold_cosine_threshold: float | None = None,
    minimum_endpoint_separation: float | None = None,
) -> dict:
    """Validate a calculator-independent initial path and return its report."""

    report = path_geometry_diagnostics(
        images,
        reference_cell=reference_cell,
        minimum_distance=minimum_distance,
        maximum_deformation=maximum_deformation,
        cell_scale=cell_scale,
        fold_cosine_threshold=fold_cosine_threshold,
        minimum_endpoint_separation=minimum_endpoint_separation,
    )
    if report["issues"]:
        detail = "\n".join(f"- {issue}" for issue in report["issues"])
        raise ValueError(f"VC-NEB path geometry preflight failed:\n{detail}")
    return report


def validate_candidate_cell_step(
    previous_images: Sequence[Atoms],
    candidate_images: Sequence[Atoms],
    *,
    maximum_cell_step: float,
) -> None:
    """Reject a candidate with an excessive one-step cell deformation.

    ``maximum_deformation`` validates the absolute initial path.  This check
    is complementary: it compares each candidate cell with the cell that is
    currently committed, so one image cannot jump to a different volume while
    the rest of the band remains near its previous geometry.  The relative
    deformation is ``solve(old_cell, new_cell).T - I`` in ASE's row-vector
    convention.
    """

    try:
        limit = float(maximum_cell_step)
    except (TypeError, ValueError) as exc:
        raise ValueError("maximum_cell_step must be a finite positive number") from exc
    if not np.isfinite(limit) or limit <= 0.0:
        raise ValueError("maximum_cell_step must be a finite positive number")
    if len(previous_images) != len(candidate_images):
        raise ValueError("previous and candidate paths must contain the same number of images")
    for image_index in range(1, len(previous_images) - 1):
        old_cell = cell_matrix(previous_images[image_index])
        new_cell = cell_matrix(candidate_images[image_index])
        try:
            _validate_cell_matrix(old_cell, context=f"previous image {image_index} cell")
            _validate_cell_matrix(new_cell, context=f"candidate image {image_index} cell")
            relative = np.linalg.solve(old_cell, new_cell).T
        except (np.linalg.LinAlgError, ValueError) as exc:
            raise CandidateStepRejected(
                f"candidate image {image_index} has an invalid cell step",
                details=[{"image_index": image_index, "category": "candidate_cell_step"}],
            ) from exc
        norm = float(np.linalg.norm(relative - np.eye(3)))
        volume_ratio = float(np.linalg.det(new_cell) / np.linalg.det(old_cell))
        if not np.isfinite(norm) or not np.isfinite(volume_ratio) or norm > limit:
            raise CandidateStepRejected(
                f"candidate image {image_index} cell step {norm:.6g} exceeds {limit:.6g}",
                details=[
                    {
                        "image_index": image_index,
                        "category": "candidate_cell_step",
                        "relative_cell_step_frobenius": norm,
                        "volume_ratio": volume_ratio,
                        "maximum_cell_step": limit,
                    }
                ],
            )


def fractional_force(atoms: Atoms) -> Array:
    """Convert Cartesian forces to forces conjugate to fractional coordinates."""

    return fractional_force_from_arrays(atoms.get_forces(), cell_matrix(atoms))


def fractional_force_from_arrays(forces: Array, cell: Array) -> Array:
    """Convert a Cartesian force array to fractional-coordinate forces."""

    return np.asarray(forces, dtype=float) @ np.asarray(cell, dtype=float).T


def cell_force(
    atoms: Atoms,
    reference_cell: Array,
    *,
    pressure: float = 0.0,
    mask: Optional[Array] = None,
) -> Array:
    """Force conjugate to the deformation gradient.

    This follows ASE's UnitCellFilter row-vector convention.  Stress and
    pressure are in eV/A^3; pressure is positive for compression.
    """

    return cell_force_from_arrays(
        atoms.get_stress(voigt=False),
        atoms,
        reference_cell,
        pressure=pressure,
        mask=mask,
    )


def cell_force_from_arrays(
    stress: Array,
    atoms: Atoms,
    reference_cell: Array,
    *,
    pressure: float = 0.0,
    mask: Optional[Array] = None,
) -> Array:
    """Convert a full stress array to deformation-gradient forces."""

    stress = np.asarray(stress, dtype=float)
    current_cell = cell_matrix(atoms)
    _validate_cell_matrix(current_cell, context="current cell")
    if stress.shape != (3, 3):
        raise ValueError(f"stress must have shape (3, 3), got {stress.shape}")
    volume = atoms.get_volume()
    virial = -volume * (stress + np.eye(3) * pressure)
    deform = deformation_from_cell(current_cell, reference_cell)
    force = np.linalg.solve(deform, virial.T).T
    if mask is not None:
        force = force * np.asarray(mask, dtype=float)
    return force


def _sum_array(value: Array) -> Array:
    summed = world.sum(value)
    return value if summed is None else summed


def _sum_scalar(value: float) -> float:
    try:
        return float(world.sum_scalar(value))
    except AttributeError:
        summed = world.sum(value)
        return float(value if summed is None else summed)


class VCNEB:
    """ASE optimizer target for variable-cell NEB.

    Coordinates used by the optimizer are length-like:

    - fractional coordinates are mapped through the reference cell,
    - deformation-gradient components are multiplied by ``cell_scale``.

    This keeps atomic and cell steps on comparable numerical footing and makes
    ``fmax`` roughly an eV/A criterion for both blocks.
    """

    def __init__(
        self,
        images: Sequence[Atoms],
        *,
        pressure: float = 0.0,
        k: float | Iterable[float] = 0.2,
        climb: bool = True,
        cell_scale: Optional[float] = None,
        atom_mask: Optional[Array] = None,
        cell_mask: Optional[Array] = None,
        mode_basis: Optional[Array] = None,
        constraint_mode: Optional[str] = None,
        mic: bool = False,
        wrap_positions: bool = False,
        parallel: bool = False,
        image_executor: object | None = None,
        candidate_validator: Optional[Callable[[Sequence[Atoms]], None]] = None,
        dynamic_relaxation: float = 1.0,
        dynamic_energy_scale: float = 0.5,
        log: Optional[Callable[[str], None]] = None,
    ) -> None:
        if len(images) < 2:
            raise ValueError("VCNEB needs at least two images")
        if any(len(image) != len(images[0]) for image in images):
            raise ValueError("All images must have the same atom count")
        symbols = images[0].get_chemical_symbols()
        if any(image.get_chemical_symbols() != symbols for image in images):
            raise ValueError("All images must have the same atom order")

        self.images = list(images)
        self.n_images = len(self.images)
        for image_index, image in enumerate(self.images):
            _validate_cell_matrix(cell_matrix(image), context=f"image {image_index} cell")
        self.reference_cell = cell_matrix(self.images[0])
        self.pressure = float(pressure)
        self.climb = bool(climb)
        self.cell_scale = float(
            cell_scale if cell_scale is not None else abs(np.linalg.det(self.reference_cell)) ** (1.0 / 3.0)
        )
        if self.cell_scale <= 0.0:
            raise ValueError("cell_scale must be positive")
        self.mic = bool(mic)
        self.wrap_positions = bool(wrap_positions)
        self.parallel = bool(parallel)
        if image_executor is not None and self.parallel and world.size > 1:
            raise ValueError("image_executor cannot be combined with MPI parallel=True")
        if image_executor is not None and not callable(getattr(image_executor, "evaluate", None)):
            raise TypeError("image_executor must provide evaluate(images)")
        self.image_executor = image_executor
        if candidate_validator is not None and not callable(candidate_validator):
            raise TypeError("candidate_validator must be callable or None")
        self.candidate_validator = candidate_validator
        self.atom_mask = None
        if atom_mask is not None:
            atom_mask_array = np.asarray(atom_mask, dtype=float)
            if atom_mask_array.size != 3 * self.n_atoms:
                raise ValueError("atom_mask must have shape (n_atoms, 3)")
            if not np.all(np.isfinite(atom_mask_array)) or not np.all(
                np.isclose(atom_mask_array, 0.0) | np.isclose(atom_mask_array, 1.0)
            ):
                raise ValueError("atom_mask must contain only finite 0/1 values")
            self.atom_mask = atom_mask_array.reshape(self.n_atoms, 3).copy()
        self.cell_mask = None if cell_mask is None else np.asarray(cell_mask, dtype=float).reshape(3, 3)
        if self.cell_mask is not None and not np.all(np.isfinite(self.cell_mask)):
            raise ValueError("cell_mask contains non-finite values")
        if self.cell_mask is not None and not np.all(
            np.isclose(self.cell_mask, 0.0) | np.isclose(self.cell_mask, 1.0)
        ):
            raise ValueError("cell_mask must contain only finite 0/1 values")
        if constraint_mode is None:
            constraint_mode = "subspace" if mode_basis is not None else "none"
        constraint_mode = str(constraint_mode).lower()
        if constraint_mode not in {"none", "subspace", "projected"}:
            raise ValueError("constraint_mode must be 'none', 'subspace', or 'projected'")
        if constraint_mode != "none" and mode_basis is None:
            raise ValueError("mode_basis is required when constraint_mode is enabled")
        if constraint_mode == "none" and mode_basis is not None:
            raise ValueError("constraint_mode must be enabled when mode_basis is supplied")
        self.constraint_mode = constraint_mode
        self.mode_basis = self._prepare_mode_basis(mode_basis)
        self._constraint_origin = self._image_x(0) if self.mode_basis is not None else None
        self._constraint_references = (
            [self._image_x(index) for index in range(self.n_images)]
            if self.mode_basis is not None
            else None
        )
        if self.mode_basis is not None and self.constraint_mode == "subspace":
            self._validate_endpoint_subspace()
            self._project_initial_interior_images()
        self.dynamic_relaxation = float(dynamic_relaxation)
        self.dynamic_energy_scale = float(dynamic_energy_scale)
        if not 0.0 <= self.dynamic_relaxation <= 1.0:
            raise ValueError("dynamic_relaxation must be between 0 and 1")
        if self.dynamic_energy_scale <= 0.0:
            raise ValueError("dynamic_energy_scale must be positive")
        self.log = log

        k_arr = np.atleast_1d(np.asarray(k, dtype=float))
        if k_arr.size == 1:
            self.k = np.full(self.n_images - 1, float(k_arr[0]))
        elif k_arr.size == self.n_images - 1:
            self.k = k_arr.copy()
        else:
            raise ValueError("k must be a scalar or have length n_images - 1")

        self._last_enthalpies: Optional[Array] = None
        self._last_forces_x: Optional[Array] = None
        self._last_true_forces_x: Optional[Array] = None
        self._last_evaluations: Optional[list[ImageEvaluation]] = None
        self._last_evaluation_geometry = None
        # Endpoints are fixed during VC-NEB; cache their evaluations and only
        # dispatch interior images to the worker executor on later iterations.
        self._endpoint_evaluations: dict[int, ImageEvaluation] = {}
        self._owners = self._assign_image_owners()
        self._validate_candidate_images([image.copy() for image in self.images])

    def __ase_optimizable__(self) -> "VCNEB":
        return self

    def __len__(self) -> int:
        return self.ndofs() // 3

    @property
    def n_atoms(self) -> int:
        return len(self.images[0])

    @property
    def image_ndofs(self) -> int:
        return 3 * self.n_atoms + 9

    def ndofs(self) -> int:
        return max(0, self.n_images - 2) * self.image_ndofs

    def _assign_image_owners(self) -> list[int]:
        owners = [0] * self.n_images
        if not self.parallel or world.size == 1:
            for image_index in range(1, self.n_images - 1):
                owners[image_index] = 0
            return owners
        for rank_index, image_index in enumerate(range(1, self.n_images - 1)):
            owners[image_index] = rank_index % world.size
        return owners

    def _own_image(self, image_index: int) -> bool:
        if image_index in (0, self.n_images - 1):
            return (not self.parallel) or world.size == 1 or world.rank == 0
        if not self.parallel or world.size == 1:
            return True
        return self._owners[image_index] == world.rank

    def _state(self, image_index: int) -> VCNEBState:
        return state_from_atoms(self.images[image_index], self.reference_cell)

    def _state_to_x(self, state: VCNEBState) -> Array:
        x_atoms = state.q @ self.reference_cell
        x_cell = self.cell_scale * (state.deform - np.eye(3))
        return np.concatenate([x_atoms.reshape(-1), x_cell.reshape(-1)])

    def _active_x_mask(self) -> Array:
        atom_mask = (
            np.ones((self.n_atoms, 3), dtype=float)
            if self.atom_mask is None
            else self.atom_mask
        ).reshape(-1)
        cell_mask = (
            np.ones((3, 3), dtype=float)
            if self.cell_mask is None
            else self.cell_mask
        ).reshape(-1)
        return np.concatenate([atom_mask, cell_mask])

    def _prepare_mode_basis(self, mode_basis: Optional[Array]) -> Optional[Array]:
        if mode_basis is None:
            return None
        basis = np.asarray(mode_basis, dtype=float)
        if basis.ndim == 1:
            basis = basis.reshape((-1, 1))
        if basis.ndim != 2 or basis.shape[0] != self.image_ndofs:
            raise ValueError(f"mode_basis must have shape ({self.image_ndofs}, n_modes)")
        if not np.all(np.isfinite(basis)):
            raise ValueError("mode_basis contains non-finite values")
        basis = basis * self._active_x_mask()[:, None]
        if not np.any(np.abs(basis) > 0.0):
            raise ValueError("mode_basis has no active components after masks")
        left_vectors, singular_values, _ = np.linalg.svd(basis, full_matrices=False)
        if singular_values.size == 0:
            raise ValueError("mode_basis is empty")
        tolerance = max(basis.shape) * np.finfo(float).eps * singular_values[0] * 100.0
        rank = singular_values > tolerance
        if not np.any(rank):
            raise ValueError("mode_basis is numerically rank deficient")
        return left_vectors[:, rank].copy()

    def _project_constraint(self, x: Array, image_index: int) -> Array:
        if self.mode_basis is None:
            return np.asarray(x, dtype=float)
        if self.constraint_mode == "subspace":
            assert self._constraint_origin is not None
            reference = self._constraint_origin
        else:
            assert self._constraint_references is not None
            reference = self._constraint_references[image_index]
        delta = np.asarray(x, dtype=float) - reference
        return reference + self.mode_basis @ (self.mode_basis.T @ delta)

    def _project_constraint_force(self, force: Array) -> Array:
        if self.mode_basis is None:
            return np.asarray(force, dtype=float)
        return self.mode_basis @ (self.mode_basis.T @ np.asarray(force, dtype=float))

    def _validate_endpoint_subspace(self) -> None:
        assert self.mode_basis is not None
        endpoint_delta = self._image_x(self.n_images - 1) - self._image_x(0)
        residual = endpoint_delta - self.mode_basis @ (self.mode_basis.T @ endpoint_delta)
        tolerance = 1e-8 * max(1.0, float(np.linalg.norm(endpoint_delta)))
        if np.linalg.norm(residual) > tolerance:
            raise ValueError(
                "VCNEB subspace constraint cannot connect the endpoints: "
                f"residual={np.linalg.norm(residual):.3e} > {tolerance:.3e}"
            )

    def _project_initial_interior_images(self) -> None:
        for image_index in range(1, self.n_images - 1):
            projected = self._project_constraint(self._image_x(image_index), image_index)
            if np.allclose(projected, self._image_x(image_index), rtol=0.0, atol=1e-14):
                continue
            apply_state(
                self.images[image_index],
                self._x_to_state(projected, image_index),
                self.reference_cell,
                wrap_positions=self.wrap_positions,
            )

    def _x_to_state(self, x: Array, image_index: int) -> VCNEBState:
        atom_size = 3 * self.n_atoms
        x_atoms = x[:atom_size].reshape(self.n_atoms, 3)
        x_cell = x[atom_size:].reshape(3, 3)
        q = x_atoms @ np.linalg.inv(self.reference_cell)
        if self.atom_mask is not None:
            current = state_from_atoms(self.images[image_index], self.reference_cell)
            current_x_atoms = current.q @ self.reference_cell
            q = (
                (x_atoms * self.atom_mask + current_x_atoms * (1.0 - self.atom_mask))
                @ np.linalg.inv(self.reference_cell)
            )
        deform = np.eye(3) + x_cell / self.cell_scale
        if self.cell_mask is not None:
            current = deformation_from_cell(cell_matrix(self.images[image_index]), self.reference_cell)
            deform = np.eye(3) + (deform - np.eye(3)) * self.cell_mask + (current - np.eye(3)) * (1.0 - self.cell_mask)
        _validate_cell_matrix(
            cell_from_deformation(deform, self.reference_cell),
            context=f"proposed cell for image {image_index}",
        )
        return VCNEBState(q=q, deform=deform)

    def _force_to_x(self, force: VCNEBState) -> Array:
        f_atoms = force.q @ np.linalg.inv(self.reference_cell.T)
        f_cell = force.deform / self.cell_scale
        return np.concatenate([f_atoms.reshape(-1), f_cell.reshape(-1)])

    def _x_to_force(self, x_force: Array) -> VCNEBState:
        atom_size = 3 * self.n_atoms
        f_atoms_x = x_force[:atom_size].reshape(self.n_atoms, 3)
        f_cell_x = x_force[atom_size:].reshape(3, 3)
        return VCNEBState(
            q=f_atoms_x @ self.reference_cell.T,
            deform=f_cell_x * self.cell_scale,
        )

    def _image_x(self, image_index: int) -> Array:
        return self._state_to_x(self._state(image_index))

    def get_x(self) -> Array:
        if self.ndofs() == 0:
            return np.zeros(0)
        return np.concatenate([self._image_x(i) for i in range(1, self.n_images - 1)])

    def candidate_images_from_x(self, x: Array) -> list[Atoms]:
        """Build, but do not validate or commit, a candidate image chain."""

        x = np.asarray(x, dtype=float).reshape(-1)
        if x.size != self.ndofs():
            raise ValueError(f"Expected {self.ndofs()} coordinates, got {x.size}")
        if not np.isfinite(x).all():
            raise ValueError("Candidate coordinates must be finite")
        candidates = [image.copy() for image in self.images]
        for offset, image_index in enumerate(range(1, self.n_images - 1)):
            lo = offset * self.image_ndofs
            hi = lo + self.image_ndofs
            x_image = self._project_constraint(x[lo:hi], image_index)
            try:
                apply_state(
                    candidates[image_index],
                    self._x_to_state(x_image, image_index),
                    self.reference_cell,
                    wrap_positions=self.wrap_positions,
                )
            except ValueError as error:
                raise ValueError(f"image {image_index}: {error}") from error
        return candidates

    def set_x(self, x: Array) -> None:
        x = np.asarray(x, dtype=float).reshape(-1)
        try:
            candidates = self.candidate_images_from_x(x)
        except ValueError as error:
            if self.candidate_validator is None:
                raise
            rejected = CandidateStepRejected(str(error))
            rejected.candidate_coordinates = x.copy()
            raise rejected from error
        try:
            self._validate_candidate_images(candidates)
        except CandidateStepRejected as error:
            error.candidate_coordinates = x.copy()
            raise
        # Only commit once every image passes; preserve calculator identities,
        # endpoints and constraints. No calculator was attached to preview copies.
        for image_index in range(1, self.n_images - 1):
            self.images[image_index].set_cell(candidates[image_index].cell, scale_atoms=False, apply_constraint=False)
            self.images[image_index].set_positions(candidates[image_index].positions, apply_constraint=False)
        self._last_enthalpies = None
        self._last_forces_x = None
        self._last_true_forces_x = None
        self._last_evaluations = None
        self._last_evaluation_geometry = None

    def _validate_candidate_images(self, candidates: Sequence[Atoms]) -> None:
        if self.candidate_validator is None:
            return
        snapshots = [(image.cell.array.copy(), image.positions.copy(), image.numbers.copy(), image.pbc.copy())
                     for image in candidates]
        self.candidate_validator(candidates)
        if len(candidates) != len(snapshots):
            raise RuntimeError("candidate_validator must not alter the number of images")
        for image, snapshot in zip(candidates, snapshots):
            if not all(np.array_equal(value, expected) for value, expected in
                       zip((image.cell.array, image.positions, image.numbers, image.pbc), snapshot)):
                raise RuntimeError("candidate_validator must not modify the candidate structures")

    def get_positions(self) -> Array:
        return self.get_x().reshape((-1, 3))

    def set_positions(self, positions: Array) -> None:
        self.set_x(np.asarray(positions, dtype=float).reshape(-1))

    def get_masses(self) -> Array:
        return np.ones(len(self))

    def iterimages(self) -> Iterator[Atoms]:
        yield from self.images

    def _enthalpy_and_force(self, image_index: int) -> tuple[float, VCNEBState]:
        if self.image_executor is not None:
            if self._last_evaluations is None:
                self._ensure_executor_evaluations()
            evaluation = self._last_evaluations[image_index]
            atoms = self.images[image_index]
            enthalpy = float(evaluation.energy + self.pressure * atoms.get_volume())
            f_q = fractional_force_from_arrays(evaluation.forces, cell_matrix(atoms))
            f_cell = cell_force_from_arrays(
                evaluation.stress,
                atoms,
                self.reference_cell,
                pressure=self.pressure,
                mask=self.cell_mask,
            )
            return enthalpy, VCNEBState(q=f_q, deform=f_cell)

        enthalpy = 0.0
        f_q: Optional[Array] = None
        f_cell: Optional[Array] = None
        if self._own_image(image_index):
            atoms = self.images[image_index]
            try:
                energy = atoms.get_potential_energy()
                enthalpy = float(energy + self.pressure * atoms.get_volume())
                f_q = fractional_force(atoms)
                f_cell = cell_force(
                    atoms,
                    self.reference_cell,
                    pressure=self.pressure,
                    mask=self.cell_mask,
                )
            except Exception as exc:
                raise RuntimeError(
                    f"VC-NEB evaluation failed for image {image_index} "
                    f"({calculator_context(atoms.calc)}): {exc}"
                ) from exc
        enthalpy = _sum_scalar(enthalpy)
        if f_q is None:
            f_q = np.zeros((self.n_atoms, 3), dtype=float)
        if f_cell is None:
            f_cell = np.zeros((3, 3), dtype=float)
        return enthalpy, VCNEBState(q=_sum_array(f_q), deform=_sum_array(f_cell))

    def _evaluate_images(self) -> list[ImageEvaluation]:
        """Evaluate all images through the optional image-level executor."""

        return self._evaluate_image_indices(range(self.n_images))

    def _evaluate_image_indices(self, indices: Sequence[int]) -> list[ImageEvaluation]:
        """Evaluate a selected subset through the image-level executor."""

        if self.image_executor is None:
            raise RuntimeError("_evaluate_image_indices called without an image executor")
        indices = [int(index) for index in indices]
        if any(index < 0 or index >= self.n_images for index in indices):
            raise IndexError("image index is outside the VC-NEB chain")
        selected = [self.images[index] for index in indices]
        evaluate = self.image_executor.evaluate
        try:
            parameters = inspect.signature(evaluate).parameters
            supports_indices = "indices" in parameters or any(
                parameter.kind is inspect.Parameter.VAR_KEYWORD
                for parameter in parameters.values()
            )
        except (TypeError, ValueError):
            supports_indices = True
        raw = evaluate(selected, indices=indices) if supports_indices else evaluate(selected)
        if len(raw) != len(indices):
            raise ValueError(
                f"image executor returned {len(raw)} evaluations for {len(indices)} images"
            )
        evaluations: list[ImageEvaluation] = []
        for image_index, value in zip(indices, raw):
            if isinstance(value, ImageEvaluation):
                evaluation = value
            else:
                try:
                    evaluation = ImageEvaluation(
                        float(value.energy),
                        np.asarray(value.forces, dtype=float),
                        np.asarray(value.stress, dtype=float),
                    )
                except AttributeError as exc:
                    raise TypeError(
                        f"image executor result {image_index} is not an ImageEvaluation"
                    ) from exc
            evaluations.append(evaluation.validate(image_index=image_index, n_atoms=self.n_atoms))
        return evaluations

    def _evaluate_fixed_endpoint(self, image_index: int) -> ImageEvaluation:
        """Evaluate one fixed endpoint once and retain its validated result."""

        cached = self._endpoint_evaluations.get(image_index)
        if cached is not None:
            return cached
        atoms = self.images[image_index]
        try:
            evaluation = ImageEvaluation(
                float(atoms.get_potential_energy()),
                np.asarray(atoms.get_forces(), dtype=float),
                np.asarray(atoms.get_stress(voigt=False), dtype=float),
            ).validate(image_index=image_index, n_atoms=self.n_atoms)
        except Exception as exc:
            raise RuntimeError(
                f"VC-NEB endpoint evaluation failed for image {image_index} "
                f"({calculator_context(atoms.calc)}): {exc}"
            ) from exc
        self._endpoint_evaluations[image_index] = evaluation
        return evaluation

    def _ensure_executor_evaluations(self) -> None:
        """Refresh interiors while reusing the fixed endpoint evaluations."""

        first = self._evaluate_fixed_endpoint(0)
        last = self._evaluate_fixed_endpoint(self.n_images - 1)
        interior = self._evaluate_image_indices(range(1, self.n_images - 1))
        self._last_evaluations = [first, *interior, last]
        self._last_evaluation_geometry = self._geometry_signature()

    def _geometry_signature(self):
        return tuple((image.cell.array.tobytes(), image.positions.tobytes(),
                      image.numbers.tobytes(), image.pbc.tobytes()) for image in self.images)

    def evaluated_snapshot_images(self) -> list[Atoms]:
        """Serialize executor results, not potentially empty calculator caches.

        No evaluation is initiated here and live calculators are never replaced.
        Serial calculators retain their existing ASE writer behavior.
        """
        if self.image_executor is None:
            return self.images
        if self._last_evaluations is None or len(self._last_evaluations) != self.n_images:
            raise RuntimeError("A complete executor evaluation is required before saving result snapshots")
        if self._last_evaluation_geometry != self._geometry_signature():
            raise RuntimeError("Snapshot geometry changed since the complete executor evaluation")
        snapshots = []
        for index, (image, value) in enumerate(zip(self.images, self._last_evaluations)):
            evaluation = value.validate(image_index=index, n_atoms=len(image))
            snapshot = image.copy()
            snapshot.calc = SinglePointCalculator(snapshot, energy=evaluation.energy,
                                                   forces=evaluation.forces, stress=evaluation.stress)
            snapshots.append(snapshot)
        return snapshots

    def _tangent(self, image_index: int, enthalpies: Array, image_x: list[Array]) -> Array:
        d_minus = image_x[image_index] - image_x[image_index - 1]
        d_plus = image_x[image_index + 1] - image_x[image_index]
        if enthalpies[image_index + 1] > enthalpies[image_index] > enthalpies[image_index - 1]:
            tangent = d_plus
        elif enthalpies[image_index + 1] < enthalpies[image_index] < enthalpies[image_index - 1]:
            tangent = d_minus
        else:
            tangent = (
                abs(enthalpies[image_index + 1] - enthalpies[image_index]) * d_plus
                + abs(enthalpies[image_index - 1] - enthalpies[image_index]) * d_minus
            )
        norm = np.linalg.norm(tangent)
        if norm < 1e-14:
            tangent = d_plus
            norm = np.linalg.norm(tangent)
        if norm < 1e-14:
            return np.zeros_like(tangent)
        return tangent / norm

    def _dynamic_weights(self, enthalpies: Array) -> Array:
        weights = np.ones(self.n_images)
        alpha = self.dynamic_relaxation
        if alpha >= 1.0 or self.n_images <= 2:
            return weights
        scale = max(abs(self.dynamic_energy_scale), 1e-12)
        emax = float(enthalpies.max())
        for image_index in range(1, self.n_images - 1):
            weights[image_index] = alpha + (1.0 - alpha) * np.exp(-(emax - enthalpies[image_index]) / scale)
        return weights

    def get_forces(self) -> Array:
        return self._compute_forces().reshape((-1, 3))

    def get_gradient(self) -> Array:
        return -self._compute_forces()

    def _compute_forces(self) -> Array:
        if self.image_executor is not None:
            self._ensure_executor_evaluations()
        else:
            self._last_evaluations = None
        enthalpies = np.zeros(self.n_images)
        true_forces = []
        image_x = []
        image_x_active = []
        active_mask = self._active_x_mask()
        for image_index in range(self.n_images):
            h_i, f_i = self._enthalpy_and_force(image_index)
            enthalpies[image_index] = h_i
            true_forces.append(self._force_to_x(f_i))
            x_i = self._image_x(image_index)
            image_x.append(x_i)
            image_x_active.append(x_i * active_mask)

        self._last_enthalpies = enthalpies.copy()
        self._last_true_forces_x = np.asarray(true_forces, dtype=float)
        force_x = np.zeros(self.ndofs())
        if self.n_images <= 2:
            self._last_forces_x = force_x
            return force_x

        climbing_image = int(1 + np.argmax(enthalpies[1:-1]))
        weights = self._dynamic_weights(enthalpies)

        for offset, image_index in enumerate(range(1, self.n_images - 1)):
            tangent = self._tangent(image_index, enthalpies, image_x_active)
            true_force = true_forces[image_index] * active_mask
            true_force = self._project_constraint_force(true_force) * active_mask
            true_parallel = np.dot(true_force, tangent) * tangent
            force_perp = true_force - true_parallel

            d_plus = np.linalg.norm(image_x_active[image_index + 1] - image_x_active[image_index])
            d_minus = np.linalg.norm(image_x_active[image_index] - image_x_active[image_index - 1])
            spring = self.k[image_index - 1] * (d_plus - d_minus) * tangent
            neb_force = force_perp + spring

            if self.climb and image_index == climbing_image:
                neb_force = true_force - 2.0 * true_parallel
            neb_force = self._project_constraint_force(neb_force)
            neb_force *= weights[image_index] * active_mask

            lo = offset * self.image_ndofs
            hi = lo + self.image_ndofs
            force_x[lo:hi] = neb_force

        self._last_forces_x = force_x.copy()
        if self.log is not None:
            max_force = 0.0 if force_x.size == 0 else float(np.linalg.norm(force_x.reshape(-1, 3), axis=1).max())
            self.log(f"max_force={max_force:.6g} enthalpies=" + " ".join(f"{e:.8f}" for e in enthalpies))
        return force_x

    def get_potential_energy(self, force_consistent: bool = False) -> float:
        return self.get_value()

    def get_value(self) -> float:
        if self._last_enthalpies is None:
            if self.image_executor is not None:
                self._compute_forces()
                assert self._last_enthalpies is not None
                return float(self._last_enthalpies.max())
            enthalpies = [self._enthalpy_and_force(i)[0] for i in range(self.n_images)]
            self._last_enthalpies = np.asarray(enthalpies, dtype=float)
        return float(self._last_enthalpies.max())

    def converged(self, gradient: Array, fmax: float) -> bool:
        max_force = self.gradient_norm(gradient)
        return bool(np.isfinite(max_force) and max_force < fmax)

    def gradient_norm(self, gradient: Array) -> float:
        if gradient.size == 0:
            return 0.0
        forces = (-gradient).reshape(-1, 3)
        return float(np.linalg.norm(forces, axis=1).max())

    @property
    def enthalpies(self) -> Array:
        if self._last_enthalpies is None:
            self.get_forces()
        assert self._last_enthalpies is not None
        return self._last_enthalpies.copy()

    def highest_image_index(self) -> int | None:
        """Return the highest-energy interior image, or ``None`` for no interior."""

        if self.n_images <= 2:
            return None
        return int(1 + np.argmax(self.enthalpies[1:-1]))

    def _interior_peak_indices(self, enthalpies: Array, tolerance: float = 1e-8) -> list[int]:
        """Return strict-enough local interior peaks in the current band."""

        peaks: list[int] = []
        for image_index in range(1, self.n_images - 1):
            energy = float(enthalpies[image_index])
            left = float(enthalpies[image_index - 1])
            right = float(enthalpies[image_index + 1])
            is_not_lower = energy >= left - tolerance and energy >= right - tolerance
            is_higher_on_one_side = energy > left + tolerance or energy > right + tolerance
            if is_not_lower and is_higher_on_one_side:
                peaks.append(image_index)
        return peaks

    def saddle_diagnostics(self) -> dict:
        """Summarize the current highest-image/saddle quality.

        The curvature is a local finite-difference estimate along the current
        extended-coordinate path.  It is a diagnostic, not a replacement for
        a Hessian calculation. Legacy ``tangential_force_eV_per_A`` refers to
        the *NEB residual*, whose ordinary-NEB tangent is only a spring term.
        Use ``true_tangential_force_eV_per_A`` to assess physical stationarity.
        """

        if self.n_images <= 2:
            return {
                "image_index": None,
                "is_climbing": False,
                "is_local_peak": False,
                "has_interior_barrier": False,
                "interior_peak_indices": [],
                "interior_barrier_indices": [],
                "ci_warning": "At least three images are required for an interior saddle candidate",
                "relative_enthalpy_eV": None,
                "residual_generalized_force_eV_per_A": 0.0,
                "tangential_force_eV_per_A": 0.0,
                "perpendicular_force_eV_per_A": 0.0,
                "true_tangential_force_eV_per_A": 0.0,
                "true_perpendicular_force_eV_per_A": 0.0,
                "true_generalized_force_max_vector_eV_per_A": 0.0,
                "tangent_curvature_eV_per_A2": None,
            }

        self._compute_forces()
        assert self._last_enthalpies is not None
        assert self._last_forces_x is not None
        assert self._last_true_forces_x is not None
        image_index = self.highest_image_index()
        assert image_index is not None
        interior_peak_indices = self._interior_peak_indices(self._last_enthalpies)
        interior_barrier_indices = [
            index
            for index in interior_peak_indices
            if self._last_enthalpies[index]
            > max(self._last_enthalpies[0], self._last_enthalpies[-1]) + 1e-8
        ]
        active_mask = self._active_x_mask()
        image_x = [self._image_x(index) for index in range(self.n_images)]
        image_x_active = [x * active_mask for x in image_x]
        tangent = self._tangent(image_index, self._last_enthalpies, image_x_active)
        offset = (image_index - 1) * self.image_ndofs
        residual = self._last_forces_x[offset : offset + self.image_ndofs]
        tangential = float(np.dot(residual, tangent))
        perpendicular = residual - tangential * tangent
        residual_norm = float(np.linalg.norm(residual.reshape(-1, 3), axis=1).max())
        perpendicular_norm = float(np.linalg.norm(perpendicular))
        physical_force = self._last_true_forces_x[image_index] * active_mask
        physical_tangential = float(np.dot(physical_force, tangent))
        physical_perpendicular = physical_force - physical_tangential * tangent

        d_minus = float(np.linalg.norm(image_x_active[image_index] - image_x_active[image_index - 1]))
        d_plus = float(np.linalg.norm(image_x_active[image_index + 1] - image_x_active[image_index]))
        curvature = None
        if d_minus > 1e-14 and d_plus > 1e-14:
            slope_plus = (self._last_enthalpies[image_index + 1] - self._last_enthalpies[image_index]) / d_plus
            slope_minus = (self._last_enthalpies[image_index] - self._last_enthalpies[image_index - 1]) / d_minus
            curvature = float(2.0 * (slope_plus - slope_minus) / (d_plus + d_minus))

        return {
            "image_index": image_index,
            "is_climbing": bool(self.climb),
            "is_local_peak": bool(image_index in interior_peak_indices),
            "has_interior_barrier": bool(interior_barrier_indices),
            "interior_peak_indices": interior_peak_indices,
            "interior_barrier_indices": interior_barrier_indices,
            "ci_warning": (
                None
                if interior_barrier_indices
                else "No interior energy peak rises above both endpoints; CI is not a validated transition-state search"
            ),
            "relative_enthalpy_eV": float(self._last_enthalpies[image_index] - self._last_enthalpies[0]),
            "residual_generalized_force_eV_per_A": residual_norm,
            "tangential_force_eV_per_A": tangential,
            "perpendicular_force_eV_per_A": perpendicular_norm,
            "true_tangential_force_eV_per_A": physical_tangential,
            "true_perpendicular_force_eV_per_A": float(np.linalg.norm(physical_perpendicular)),
            "true_generalized_force_max_vector_eV_per_A": float(
                np.linalg.norm(physical_force.reshape(-1, 3), axis=1).max()
            ),
            "tangent_curvature_eV_per_A2": curvature,
        }

    def path_diagnostics(self) -> dict:
        """Return per-image physical and NEB force diagnostics.

        The physical force fields are reported separately from the projected
        NEB residual.  This distinction is important for variable-cell runs:
        a small NEB residual does not imply small raw atomic forces or stress,
        and a large residual can be caused by the cell block rather than atoms.
        In parallel mode, image-local arrays are reduced so every rank receives
        the same JSON-serializable result.
        """

        self._compute_forces()
        assert self._last_enthalpies is not None
        assert self._last_forces_x is not None
        active_mask = self._active_x_mask()
        image_x = [self._image_x(index) for index in range(self.n_images)]
        image_x_active = [x * active_mask for x in image_x]
        physical: list[tuple[VCNEBState, Array, Array]] = []
        for image_index in range(self.n_images):
            _, force_state = self._enthalpy_and_force(image_index)
            atom_forces = np.zeros((self.n_atoms, 3), dtype=float)
            stress = np.zeros((3, 3), dtype=float)
            if self.image_executor is not None:
                assert self._last_evaluations is not None
                evaluation = self._last_evaluations[image_index]
                atom_forces = np.asarray(evaluation.forces, dtype=float)
                stress = np.asarray(evaluation.stress, dtype=float)
            elif self._own_image(image_index):
                atoms = self.images[image_index]
                atom_forces = np.asarray(atoms.get_forces(), dtype=float)
                stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
            physical.append((force_state, _sum_array(atom_forces), _sum_array(stress)))

        climbing_image = self.highest_image_index()
        interior_peak_indices = self._interior_peak_indices(self._last_enthalpies)
        interior_barrier_indices = [
            index
            for index in interior_peak_indices
            if self._last_enthalpies[index]
            > max(self._last_enthalpies[0], self._last_enthalpies[-1]) + 1e-8
        ]
        image_records = []
        for image_index, (force_state, atom_forces, stress) in enumerate(physical):
            true_force_x = self._force_to_x(force_state)
            true_force_norm = float(np.linalg.norm(true_force_x.reshape(-1, 3), axis=1).max())
            cell_force_x = true_force_x[3 * self.n_atoms :]
            active_true_force_x = true_force_x * active_mask
            constrained_true_force_x = self._project_constraint_force(active_true_force_x) * active_mask
            released_force_x = active_true_force_x - constrained_true_force_x
            record = {
                "image_index": image_index,
                "interior": 0 < image_index < self.n_images - 1,
                "enthalpy_eV": float(self._last_enthalpies[image_index]),
                "relative_enthalpy_eV": float(self._last_enthalpies[image_index] - self._last_enthalpies[0]),
                "volume_A3": float(self.images[image_index].get_volume()),
                "cell_lengths_A": [float(value) for value in self.images[image_index].cell.lengths()],
                "cell_angles_deg": [float(value) for value in self.images[image_index].cell.angles()],
                "max_atom_force_eV_per_A": float(np.linalg.norm(atom_forces, axis=1).max()),
                "max_stress_eV_per_A3": float(np.max(np.abs(stress))),
                "stress_eV_per_A3": stress.tolist(),
                "cell_force_eV": force_state.deform.tolist(),
                "max_cell_force_eV": float(np.max(np.abs(force_state.deform))),
                "max_true_generalized_force_eV_per_A": true_force_norm,
                "cell_generalized_force_norm_eV_per_A": float(np.linalg.norm(cell_force_x)),
            }
            if self.mode_basis is not None:
                record.update({
                    "raw_force_orthogonal_to_mode_subspace_norm_eV_per_A": float(
                        np.linalg.norm(released_force_x)
                    ),
                    "raw_force_orthogonal_to_mode_subspace_max_vector_eV_per_A": float(
                        np.linalg.norm(released_force_x.reshape(-1, 3), axis=1).max()
                    ),
                })
            if 0 < image_index < self.n_images - 1:
                tangent = self._tangent(image_index, self._last_enthalpies, image_x_active)
                constrained_true = constrained_true_force_x
                true_parallel = float(np.dot(constrained_true, tangent))
                true_perpendicular = constrained_true - true_parallel * tangent
                d_plus = float(np.linalg.norm(image_x_active[image_index + 1] - image_x_active[image_index]))
                d_minus = float(np.linalg.norm(image_x_active[image_index] - image_x_active[image_index - 1]))
                spring = self.k[image_index - 1] * (d_plus - d_minus) * tangent
                offset = (image_index - 1) * self.image_ndofs
                residual = self._last_forces_x[offset : offset + self.image_ndofs]
                record.update(
                    {
                        "is_climbing_image": bool(self.climb and image_index == climbing_image),
                        "spacing_minus_A": d_minus,
                        "spacing_plus_A": d_plus,
                        "true_tangential_force_eV_per_A": true_parallel,
                        "true_perpendicular_force_eV_per_A": float(np.linalg.norm(true_perpendicular)),
                        "true_perpendicular_force_max_vector_eV_per_A": float(
                            np.linalg.norm(true_perpendicular.reshape(-1, 3), axis=1).max()
                        ),
                        "spring_force_eV_per_A": float(np.linalg.norm(spring)),
                        "spring_force_max_vector_eV_per_A": float(
                            np.linalg.norm(spring.reshape(-1, 3), axis=1).max()
                        ),
                        "neb_residual_generalized_force_eV_per_A": float(
                            np.linalg.norm(residual.reshape(-1, 3), axis=1).max()
                        ),
                        "neb_residual_force_euclidean_eV_per_A": float(np.linalg.norm(residual)),
                    }
                )
            else:
                record.update(
                    {
                        "is_climbing_image": False,
                        "spacing_minus_A": None,
                        "spacing_plus_A": None,
                        "true_tangential_force_eV_per_A": None,
                        "true_perpendicular_force_eV_per_A": None,
                        "true_perpendicular_force_max_vector_eV_per_A": None,
                        "spring_force_eV_per_A": None,
                        "spring_force_max_vector_eV_per_A": None,
                        "neb_residual_generalized_force_eV_per_A": None,
                        "neb_residual_force_euclidean_eV_per_A": None,
                    }
                )
            image_records.append(record)

        geometry = path_geometry_diagnostics(self.images, cell_scale=self.cell_scale)
        return {
            "n_images": self.n_images,
            "n_atoms": self.n_atoms,
            "constraint_mode": self.constraint_mode,
            "n_mode_directions": None if self.mode_basis is None else int(self.mode_basis.shape[1]),
            "pressure_eV_per_A3": self.pressure,
            "cell_scale_A": self.cell_scale,
            "geometry": geometry,
            "highest_image_index": climbing_image,
            "interior_peak_indices": interior_peak_indices,
            "interior_barrier_indices": interior_barrier_indices,
            "has_interior_barrier": bool(interior_barrier_indices),
            "ci_warning": (
                None
                if interior_barrier_indices
                else "No interior energy peak rises above both endpoints; inspect the path before interpreting CI output"
            ),
            "images": image_records,
        }

    def reaction_coordinate(self) -> Array:
        coords = [0.0]
        active_mask = self._active_x_mask()
        for image_index in range(1, self.n_images):
            dx = (self._image_x(image_index) - self._image_x(image_index - 1)) * active_mask
            coords.append(coords[-1] + float(np.linalg.norm(dx)))
        return np.asarray(coords)

    def barrier(self) -> tuple[float, float]:
        enthalpies = self.enthalpies
        return float(enthalpies.max() - enthalpies[0]), float(enthalpies[-1] - enthalpies[0])

    def write_chain(self, path: str | Path) -> None:
        write(str(path), self.images)

    def write_step_directory(self, directory: str | Path, step: int) -> None:
        snapshots = (self.evaluated_snapshot_images()
                     if self.image_executor is not None and self._last_evaluations is not None else self.images)
        root = Path(directory)
        root.mkdir(parents=True, exist_ok=True)
        chain_path = root / f"chain_step_{step:04d}.traj"
        chain_fd, chain_tmp = tempfile.mkstemp(prefix=f".{chain_path.name}.", suffix=".tmp", dir=root)
        os.close(chain_fd)
        try:
            write(chain_tmp, snapshots, format="traj")
            os.replace(chain_tmp, chain_path)
        finally:
            try:
                os.unlink(chain_tmp)
            except FileNotFoundError:
                pass
        step_dir = root / f"step_{step:04d}"
        step_dir.mkdir(exist_ok=True)
        for image_index, image in enumerate(snapshots):
            poscar_path = step_dir / f"POSCAR_{image_index:02d}"
            fd, temporary = tempfile.mkstemp(prefix=f".{poscar_path.name}.", suffix=".tmp", dir=step_dir)
            os.close(fd)
            try:
                write(temporary, image, format="vasp", direct=True, vasp5=True)
                os.replace(temporary, poscar_path)
            finally:
                try:
                    os.unlink(temporary)
                except FileNotFoundError:
                    pass

    def plot_band(self, filename: str | Path) -> None:
        import matplotlib.pyplot as plt

        s = self.reaction_coordinate()
        h = self.enthalpies
        fig, ax = plt.subplots()
        ax.plot(s, h - h[0], marker="o")
        ax.set_xlabel("Reaction coordinate (A, extended)")
        ax.set_ylabel("Relative enthalpy (eV)")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(filename, dpi=200)
        plt.close(fig)


def read_chain_trajectory(
    path: str | Path,
    *,
    n_images: int,
    step: int = -1,
) -> list[Atoms]:
    """Read one complete VC-NEB chain from a flat ASE trajectory.

    ``run_vcneb`` stores each optimizer snapshot as ``n_images`` consecutive
    trajectory frames.  Interrupted writes can leave a partial tail; negative
    steps count over complete snapshots only, so ``step=-1`` returns the latest
    complete chain and ignores any partial trailing frames.
    """

    if n_images < 2:
        raise ValueError("n_images must include endpoints and be at least 2")
    frames = list(Trajectory(str(path), "r"))
    complete_steps = len(frames) // n_images
    if complete_steps == 0:
        raise ValueError(f"No complete {n_images}-image chain found in {path}")
    step_index = complete_steps + step if step < 0 else step
    if step_index < 0 or step_index >= complete_steps:
        raise IndexError(f"Trajectory {path} has {complete_steps} complete chain snapshots; got step {step}")
    chain = [frame.copy() for frame in frames[step_index * n_images : (step_index + 1) * n_images]]
    _validate_chain_compatibility(chain)
    return chain


def apply_chain_state(target: Sequence[Atoms], source: Sequence[Atoms]) -> None:
    """Copy cells and positions from ``source`` images while preserving calculators."""

    if len(target) != len(source):
        raise ValueError(f"Expected {len(target)} source images, got {len(source)}")
    _validate_chain_compatibility(source)
    target_symbols = target[0].get_chemical_symbols()
    for image_index, (dst, src) in enumerate(zip(target, source)):
        if dst.get_chemical_symbols() != target_symbols or src.get_chemical_symbols() != target_symbols:
            raise ValueError(f"Image {image_index} atom symbols/order are incompatible with the target chain")
        dst.set_cell(src.cell.array, scale_atoms=False)
        dst.set_positions(src.positions)
        dst.pbc = src.pbc


def _validate_chain_compatibility(images: Sequence[Atoms]) -> None:
    if not images:
        raise ValueError("Chain is empty")
    symbols = images[0].get_chemical_symbols()
    for image_index, image in enumerate(images):
        if image.get_chemical_symbols() != symbols:
            raise ValueError(f"Image {image_index} atom symbols/order differ from image 0")
        if image.cell.rank != 3 or image.get_volume() <= 0.0:
            raise ValueError(f"Image {image_index} has an invalid 3D cell")


def _make_optimizer(
    name: str,
    chain: VCNEB,
    logfile: str | Path | None,
    optimizer_kwargs: Optional[Mapping[str, object]] = None,
    *,
    candidate_step_retries: int = 0,
    candidate_step_retry_factor: float = 0.5,
    candidate_step_manifest: str | Path | None = None,
):
    kwargs = {} if optimizer_kwargs is None else dict(optimizer_kwargs)
    key = name.upper()
    if key == "FIRE":
        if chain.candidate_validator is not None:
            return CheckedFIRE(chain, logfile=logfile, max_candidate_retries=candidate_step_retries,
                               candidate_retry_factor=candidate_step_retry_factor,
                               candidate_manifest=candidate_step_manifest, **kwargs)
        return FIRE(chain, logfile=logfile, **kwargs)
    if key in {"BLOCKFIRE", "BLOCK_FIRE", "BLOCK-FIRE"}:
        return BlockFIRE(
            chain,
            logfile=logfile,
            block_mode="image",
            max_candidate_retries=candidate_step_retries,
            candidate_retry_factor=candidate_step_retry_factor,
            candidate_manifest=candidate_step_manifest,
            **kwargs,
        )
    if key in {"SPLITFIRE", "SPLIT_FIRE", "SPLIT-FIRE"}:
        return BlockFIRE(
            chain,
            logfile=logfile,
            block_mode="atomic-cell",
            max_candidate_retries=candidate_step_retries,
            candidate_retry_factor=candidate_step_retry_factor,
            candidate_manifest=candidate_step_manifest,
            **kwargs,
        )
    if key in {"IMAGESCALEDFIRE", "IMAGE_SCALED_FIRE", "IMAGE-SCALED-FIRE"}:
        return ImageScaledFIRE(
            chain,
            logfile=logfile,
            max_candidate_retries=candidate_step_retries,
            candidate_retry_factor=candidate_step_retry_factor,
            candidate_manifest=candidate_step_manifest,
            **kwargs,
        )
    if key in {"STAGEDFIRE", "STAGED_FIRE", "STAGED-FIRE"}:
        return StagedFIRE(
            chain,
            logfile=logfile,
            max_candidate_retries=candidate_step_retries,
            candidate_retry_factor=candidate_step_retry_factor,
            candidate_manifest=candidate_step_manifest,
            **kwargs,
        )
    if key == "LBFGS":
        return LBFGS(chain, logfile=logfile, **kwargs)
    if key == "BFGS":
        return BFGS(chain, logfile=logfile, **kwargs)
    if key in {"BFGSLINESEARCH", "BFGS_LINESEARCH", "BFGS-LINESEARCH"}:
        return BFGSLineSearch(chain, logfile=logfile, **kwargs)
    raise ValueError(
        f"Unknown optimizer {name!r}; choose FIRE, ImageScaledFIRE, StagedFIRE, "
        "BlockFIRE, SplitFIRE, BFGS, LBFGS, or BFGSLineSearch"
    )


def _is_line_search_failure(error: BaseException) -> bool:
    """Return whether an optimizer error is the recoverable ASE line-search failure."""

    message = str(error).lower().replace("-", " ")
    return "line search failed" in message or "linesearch failed" in message


def _optimizer_uses_line_search(name: str, optimizer_kwargs: Mapping[str, object]) -> bool:
    """Return whether the selected ASE optimizer owns an explicit line search."""

    key = str(name).upper().replace("-", "_")
    if key in {"BFGSLINESEARCH", "BFGS_LINESEARCH"}:
        return True
    return key == "LBFGS" and bool(optimizer_kwargs.get("use_line_search", False))


def _line_search_retry_kwargs(
    optimizer_kwargs: Mapping[str, object],
    *,
    retry_index: int,
    factor: float,
    optimizer_key: str,
) -> dict[str, object]:
    """Reduce the ASE line-search step cap for one fresh optimizer attempt."""

    kwargs = dict(optimizer_kwargs)
    scale = factor ** int(retry_index)
    maxstep = kwargs.get("maxstep")
    if maxstep is None:
        # ASE's BFGSLineSearch default is 0.2 Angstrom.  Materialize the
        # reduced value so a retry is deterministic rather than relying on a
        # version-dependent default.
        kwargs["maxstep"] = 0.2 * scale
    else:
        kwargs["maxstep"] = float(maxstep) * scale
    if "stpmax" in kwargs:
        kwargs["stpmax"] = float(kwargs["stpmax"]) * scale
    elif optimizer_key in {"BFGSLINESEARCH", "BFGS_LINESEARCH", "BFGS-LINESEARCH"}:
        # BFGSLineSearch defaults to 50.0; make the retry cap explicit.
        kwargs["stpmax"] = 50.0 * scale
    return kwargs


def _next_snapshot_step(directory: str | Path) -> int:
    root = Path(directory)
    if not root.exists():
        return 0
    next_step = 0
    for path in root.iterdir():
        name = path.name
        if name.startswith("chain_step_") and name.endswith(".traj"):
            token = name.removeprefix("chain_step_").removesuffix(".traj")
        elif name.startswith("step_"):
            token = name.removeprefix("step_")
        else:
            continue
        try:
            next_step = max(next_step, int(token) + 1)
        except ValueError:
            continue
    return next_step


def _write_json_atomic(path: str | Path, payload: Mapping[str, object]) -> None:
    """Persist a small machine-readable report without exposing partial JSON."""

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, destination)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def run_vcneb(
    images: Sequence[Atoms],
    *,
    pressure_gpa: float = 0.0,
    k: float | Iterable[float] = 0.2,
    climb: bool = True,
    climb_after: int | None = None,
    cell_scale: Optional[float] = None,
    atom_mask: Optional[Array] = None,
    cell_mask: Optional[Array] = None,
    mode_basis: Optional[Array] = None,
    constraint_mode: Optional[str] = None,
    mic: bool = False,
    wrap_positions: bool = False,
    parallel: bool = False,
    image_executor: object | None = None,
    candidate_validator: Optional[Callable[[Sequence[Atoms]], None]] = None,
    candidate_step_retries: int = 0,
    candidate_step_retry_factor: float = 0.5,
    candidate_step_manifest: str | Path | None = None,
    maximum_cell_step: float | None = None,
    optimizer: str = "FIRE",
    optimizer_kwargs: Optional[Mapping[str, object]] = None,
    line_search_retries: int = 0,
    line_search_retry_factor: float = 0.5,
    fmax: float = 0.10,
    steps: int = 300,
    logfile: str | Path | None = "vcneb-opt.log",
    trajectory: str | Path | None = "vcneb.traj",
    trajectory_mode: str = "w",
    snapshot_dir: str | Path | None = None,
    snapshot_start: Optional[int] = None,
    failure_report: str | Path | None = None,
    validate_calculators: bool = True,
    log: Optional[Callable[[str], None]] = None,
) -> tuple[VCNEB, object]:
    """Run a VC-NEB optimization with an ASE optimizer.

    When ``climb_after`` is provided with ``climb=True``, ordinary NEB is
    used for that many completed optimizer steps before the climbing-image
    force is enabled.  This prevents a noisy initial path from selecting a
    wrong climbing image too early.

    If ``failure_report`` is supplied, optimizer/calculator exceptions are
    recorded atomically with the completed-step count and recovery locations
    before the original exception is propagated.

    ``line_search_retries`` enables a bounded recovery for ASE's explicit
    line-search optimizers.  On a ``LineSearch failed!`` exception, the
    optimizer is recreated at the last complete image state with a reduced
    ``maxstep``/``stpmax`` and only the remaining steps are attempted.  Other
    optimizer or calculator errors are never retried automatically.

    ``candidate_validator`` previews every image before an atomic geometry
    commit. With FIRE, ``candidate_step_retries`` permits fixed-direction
    geometry backtracking before DFT; accepted shortened moves reset momentum.
    ``maximum_cell_step`` adds a calculator-independent per-iteration cell
    trust radius and can be combined with those retries.
    This is distinct from retrying a failed electronic calculation.
    """

    if climb_after is not None:
        if isinstance(climb_after, bool) or int(climb_after) != climb_after or climb_after < 0:
            raise ValueError("climb_after must be a non-negative integer or None")
        climb_after = int(climb_after)
    if not climb:
        climb_after = None
    if (isinstance(candidate_step_retries, bool) or int(candidate_step_retries) != candidate_step_retries
            or candidate_step_retries < 0):
        raise ValueError("candidate_step_retries must be a nonnegative integer")
    if not np.isfinite(candidate_step_retry_factor) or not 0 < candidate_step_retry_factor < 1:
        raise ValueError("candidate_step_retry_factor must be in (0, 1)")
    if maximum_cell_step is not None:
        try:
            maximum_cell_step = float(maximum_cell_step)
        except (TypeError, ValueError) as exc:
            raise ValueError("maximum_cell_step must be a finite positive number") from exc
        if not np.isfinite(maximum_cell_step) or maximum_cell_step <= 0.0:
            raise ValueError("maximum_cell_step must be a finite positive number")
    effective_candidate_validator = candidate_validator
    if maximum_cell_step is not None:
        def validate_candidate(candidates: Sequence[Atoms]) -> None:
            if candidate_validator is not None:
                candidate_validator(candidates)
            validate_candidate_cell_step(images, candidates, maximum_cell_step=maximum_cell_step)

        effective_candidate_validator = validate_candidate
    candidate_optimizers = {
        "FIRE", "BLOCKFIRE", "BLOCK_FIRE", "BLOCK-FIRE",
        "SPLITFIRE", "SPLIT_FIRE", "SPLIT-FIRE",
        "IMAGESCALEDFIRE", "IMAGE_SCALED_FIRE", "IMAGE-SCALED-FIRE",
        "STAGEDFIRE", "STAGED_FIRE", "STAGED-FIRE",
    }
    if candidate_step_retries and (
        effective_candidate_validator is None or optimizer.upper() not in candidate_optimizers
    ):
        raise ValueError("candidate step backtracking requires a candidate validator and a FIRE optimizer")
    if isinstance(line_search_retries, bool) or int(line_search_retries) != line_search_retries or line_search_retries < 0:
        raise ValueError("line_search_retries must be a non-negative integer")
    line_search_retries = int(line_search_retries)
    try:
        line_search_retry_factor = float(line_search_retry_factor)
    except (TypeError, ValueError) as exc:
        raise ValueError("line_search_retry_factor must be a finite number in (0, 1)") from exc
    if not np.isfinite(line_search_retry_factor) or not 0.0 < line_search_retry_factor < 1.0:
        raise ValueError("line_search_retry_factor must be a finite number in (0, 1)")

    if validate_calculators:
        validate_image_calculators(
            images,
            require_stress=True,
            require_variable_cell=True,
        )

    chain = VCNEB(
        images,
        pressure=pressure_gpa * GPa,
        k=k,
        climb=bool(climb and (climb_after is None or climb_after == 0)),
        cell_scale=cell_scale,
        atom_mask=atom_mask,
        cell_mask=cell_mask,
        mode_basis=mode_basis,
        constraint_mode=constraint_mode,
        mic=mic,
        wrap_positions=wrap_positions,
        parallel=parallel,
        image_executor=image_executor,
        candidate_validator=effective_candidate_validator,
        log=log,
    )
    base_optimizer_kwargs = {} if optimizer_kwargs is None else dict(optimizer_kwargs)
    line_search_enabled = _optimizer_uses_line_search(optimizer, base_optimizer_kwargs)
    opt = _make_optimizer(optimizer, chain, logfile, base_optimizer_kwargs,
                          candidate_step_retries=candidate_step_retries,
                          candidate_step_retry_factor=candidate_step_retry_factor,
                          candidate_step_manifest=candidate_step_manifest)

    if trajectory_mode not in {"w", "a"}:
        raise ValueError("trajectory_mode must be 'w' or 'a'")
    traj = Trajectory(str(trajectory), trajectory_mode) if trajectory else None
    if snapshot_start is None:
        snapshot_start = _next_snapshot_step(snapshot_dir) if trajectory_mode == "a" and snapshot_dir is not None else 0
    if snapshot_start < 0:
        raise ValueError("snapshot_start must be non-negative")
    step_counter = {"value": int(snapshot_start)}

    def save_snapshot() -> None:
        snapshots = chain.evaluated_snapshot_images()
        if traj is not None:
            for image in snapshots:
                traj.write(image)
        if snapshot_dir is not None:
            chain.write_step_directory(snapshot_dir, step_counter["value"])
        step_counter["value"] += 1

    def attach_observers(current_optimizer: object) -> None:
        current_optimizer.attach(save_snapshot, interval=1)
        if climb_after is not None and climb_after > 0:
            def enable_climbing_image() -> None:
                if getattr(opt, "nsteps", 0) >= climb_after:
                    chain.climb = True

            current_optimizer.attach(enable_climbing_image, interval=1)

    attach_observers(opt)
    retry_history: list[dict[str, object]] = []
    completed_steps = 0
    try:
        while True:
            try:
                opt.run(fmax=fmax, steps=max(0, int(steps) - completed_steps))
                break
            except Exception as exc:
                completed_steps = int(getattr(opt, "nsteps", completed_steps))
                if (
                    len(retry_history) >= line_search_retries
                    or not line_search_enabled
                    or not _is_line_search_failure(exc)
                ):
                    raise
                retry_index = len(retry_history) + 1
                retry_kwargs = _line_search_retry_kwargs(
                    base_optimizer_kwargs,
                    retry_index=retry_index,
                    factor=line_search_retry_factor,
                    optimizer_key=str(optimizer).upper(),
                )
                retry_history.append(
                    {
                        "retry_index": retry_index,
                        "failed_step": completed_steps,
                        "error": str(exc),
                        "optimizer_kwargs": retry_kwargs,
                    }
                )
                if log is not None:
                    log(
                        f"line_search_retry={retry_index} failed_step={completed_steps} "
                        f"maxstep={retry_kwargs.get('maxstep')} stpmax={retry_kwargs.get('stpmax')}"
                    )
                close = getattr(opt, "close", None)
                if callable(close):
                    close()
                opt = _make_optimizer(optimizer, chain, logfile, retry_kwargs)
                # ASE counts successful optimizer steps in ``nsteps``.  Carry
                # that count into the fresh instance so observers, staged CI,
                # and the remaining-step budget retain global semantics.
                opt.nsteps = completed_steps
                attach_observers(opt)
    except Exception as exc:
        if traj is not None:
            traj.close()
            traj = None
        if failure_report is not None:
            diagnostic_paths = _calculator_diagnostic_paths(chain.images)
            payload = {
                "status": "failed",
                "error_type": type(exc).__name__,
                "failure_category": "candidate_step_rejected" if isinstance(exc, CandidateStepRejected) else classify_calculator_failure(
                    exc,
                    diagnostic_paths=diagnostic_paths,
                ),
                "error": str(exc),
                "optimizer": str(optimizer),
                "optimizer_steps_completed": int(getattr(opt, "nsteps", 0)),
                "candidate_step_retries_allowed": candidate_step_retries,
                "candidate_step_history": getattr(opt, "candidate_step_history", []),
                "candidate_rejection_details": exc.details if isinstance(exc, CandidateStepRejected) else [],
                "line_search_retries_allowed": line_search_retries,
                "line_search_retries_used": len(retry_history),
                "line_search_retry_history": retry_history,
                "snapshot_steps_written": max(0, int(step_counter["value"])),
                "n_images": len(chain.images),
                "trajectory": str(trajectory) if trajectory is not None else None,
                "snapshot_dir": str(snapshot_dir) if snapshot_dir is not None else None,
                "diagnostic_paths": [str(path) for path in diagnostic_paths],
                "recovery_hint": (
                    "Inspect the preserved trajectory/snapshots, correct the calculator or "
                    "optimizer settings, then resume from the latest complete chain snapshot."
                ),
            }
            try:
                _write_json_atomic(failure_report, payload)
            except Exception:
                # Preserve the original calculator/optimizer exception if the
                # optional diagnostic path itself is unavailable.
                pass
        raise
    finally:
        if traj is not None:
            traj.close()
    # Expose recovery provenance on the returned ASE optimizer without
    # changing ASE's public state model.  A zero-length list is the common
    # no-retry case and is safe for callers to serialize directly.
    setattr(opt, "line_search_retries_used", len(retry_history))
    setattr(opt, "line_search_retry_history", list(retry_history))
    return chain, opt
