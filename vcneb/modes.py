"""Mode-guided path construction and modal diagnostics for VC-NEB.

The mode feature is deliberately split into two operations:

* ``mode_guided_path`` bends an otherwise ordinary endpoint interpolation and
  is used as an initial-path generator.  The subsequent NEB optimization is
  still unconstrained, so the converged path remains a physical MEP.
* ``project_path_onto_mode`` measures how much of a path lies along one or
  more supplied modes.  This is a diagnostic, not a replacement for a
  transition-state calculation.

Atomic modes are Cartesian displacements in Angstrom after normalization.  A
cell mode is a dimensionless deformation-gradient direction and is scaled by
``cell_scale`` when combined with atomic coordinates.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from ase import Atoms

from .core import (
    VCNEBState,
    apply_state,
    cell_matrix,
    interpolate_vcneb,
    state_from_atoms,
)


Array = np.ndarray


@dataclass
class Mode:
    """One real vibrational or collective mode.

    Parameters
    ----------
    atomic
        Array with shape ``(n_atoms, 3)`` or a flattened ``(3*n_atoms,)``
        vector.  It is interpreted as a Cartesian displacement pattern.
    cell
        Optional ``(3, 3)`` deformation-gradient pattern.  This is useful for
        a strain-coupled collective mode, but ordinary fixed-cell phonons do
        not need it.
    """

    atomic: Array
    cell: Array | None = None
    label: str = ""
    frequency: float | None = None

    def __post_init__(self) -> None:
        atomic = np.asarray(self.atomic, dtype=float)
        if atomic.ndim == 1:
            if atomic.size % 3:
                raise ValueError("Flattened mode length must be divisible by 3")
            atomic = atomic.reshape((-1, 3))
        if atomic.shape[1:] != (3,):
            raise ValueError("Mode.atomic must have shape (n_atoms, 3)")
        if not np.all(np.isfinite(atomic)):
            raise ValueError("Mode.atomic contains non-finite values")
        self.atomic = atomic.copy()

        if self.cell is not None:
            cell = np.asarray(self.cell, dtype=float)
            if cell.shape != (3, 3):
                raise ValueError("Mode.cell must have shape (3, 3)")
            if not np.all(np.isfinite(cell)):
                raise ValueError("Mode.cell contains non-finite values")
            self.cell = cell.copy()

    @classmethod
    def from_file(cls, path: str | Path, n_atoms: int | None = None) -> "Mode":
        """Read a simple text, JSON, or NPZ mode file.

        Text files contain either ``3*n_atoms`` or ``n_atoms`` rows of three
        numbers.  JSON accepts ``{"atomic": ..., "cell": ...}``; NPZ uses
        arrays named ``atomic`` and optionally ``cell``.
        """

        source = Path(path)
        suffix = source.suffix.lower()
        if suffix == ".json":
            data = json.loads(source.read_text(encoding="utf-8"))
            return cls(
                data["atomic"],
                cell=data.get("cell"),
                label=str(data.get("label", "")),
                frequency=data.get("frequency"),
            )
        if suffix == ".npz":
            with np.load(source) as data:
                return cls(
                    data["atomic"],
                    cell=data["cell"] if "cell" in data else None,
                    label=str(data["label"]) if "label" in data else "",
                    frequency=float(data["frequency"]) if "frequency" in data else None,
                )

        values = np.loadtxt(source, dtype=float)
        if values.ndim == 1:
            values = values.reshape(-1)
        if n_atoms is not None and values.size == 3 * n_atoms:
            values = values.reshape((n_atoms, 3))
        return cls(values)

    def normalized(
        self,
        *,
        masses: Array | None = None,
        mass_weighted_input: bool = False,
        remove_translation: bool = False,
        cell_scale: float = 1.0,
    ) -> "Mode":
        """Return a unit-norm mode in the Cartesian path metric.

        If ``mass_weighted_input`` is true, the supplied atomic eigenvector is
        converted to a Cartesian displacement by dividing each atom by
        ``sqrt(mass)``.  This makes the convention explicit for phonon files
        whose eigenvectors are mass weighted.
        """

        if cell_scale <= 0.0:
            raise ValueError("cell_scale must be positive")
        atomic = self.atomic.copy()
        if mass_weighted_input:
            if masses is None:
                raise ValueError("masses are required for mass_weighted_input")
            masses = np.asarray(masses, dtype=float).reshape(-1)
            if masses.size != len(atomic) or np.any(masses <= 0.0):
                raise ValueError("masses must contain one positive value per atom")
            atomic = atomic / np.sqrt(masses)[:, None]
        if remove_translation:
            if masses is None:
                atomic = atomic - atomic.mean(axis=0, keepdims=True)
            else:
                masses = np.asarray(masses, dtype=float).reshape(-1)
                if masses.size != len(atomic):
                    raise ValueError("masses must contain one value per atom")
                atomic = atomic - np.average(atomic, axis=0, weights=masses)[None, :]

        cell = None if self.cell is None else self.cell.copy()
        norm = np.linalg.norm(atomic)
        if cell is not None:
            norm = float(np.sqrt(norm**2 + (cell_scale * np.linalg.norm(cell)) ** 2))
        if norm < 1e-14:
            raise ValueError("Cannot normalize a zero mode")
        return Mode(
            atomic=atomic / norm,
            cell=None if cell is None else cell / norm,
            label=self.label,
            frequency=self.frequency,
        )


def _envelope(lam: float, name: str) -> float:
    key = name.lower()
    if key == "sin":
        return float(np.sin(np.pi * lam))
    if key == "bell":
        return float(4.0 * lam * (1.0 - lam))
    if key == "linear":
        return float(lam)
    raise ValueError("envelope must be 'sin', 'bell', or 'linear'")


def mode_guided_path(
    initial: Atoms,
    final: Atoms,
    n_images: int,
    mode: Mode | Array,
    *,
    amplitude: float = 0.25,
    envelope: str = "sin",
    align_cells: bool = True,
    mic: bool = False,
    wrap_positions: bool = False,
    normalize: bool = True,
    masses: Array | None = None,
    mass_weighted_input: bool = False,
    remove_translation: bool = False,
) -> list[Atoms]:
    """Create a mode-guided VC-NEB initial path with exact endpoints.

    The mode perturbation is multiplied by an endpoint-zero envelope.  Thus
    the endpoints are unchanged and the mode only selects an initial route;
    a later unconstrained NEB calculation can relax away from it.
    """

    mode_obj = mode if isinstance(mode, Mode) else Mode(mode)
    if len(initial) != len(mode_obj.atomic):
        raise ValueError("Mode atom count does not match the initial structure")

    images = interpolate_vcneb(
        initial,
        final,
        n_images,
        align_cells=align_cells,
        mic=mic,
        wrap_positions=False,
    )
    reference_cell = cell_matrix(images[0])
    cell_scale = abs(np.linalg.det(reference_cell)) ** (1.0 / 3.0)
    if normalize:
        mode_obj = mode_obj.normalized(
            masses=masses,
            mass_weighted_input=mass_weighted_input,
            remove_translation=remove_translation,
            cell_scale=cell_scale,
        )
    for image_index in range(1, n_images - 1):
        lam = image_index / (n_images - 1)
        factor = float(amplitude) * _envelope(lam, envelope)
        state = state_from_atoms(images[image_index], reference_cell)
        x_atoms = state.q @ reference_cell + factor * mode_obj.atomic
        deform = state.deform.copy()
        if mode_obj.cell is not None:
            deform = deform + (factor / cell_scale) * mode_obj.cell
        apply_state(
            images[image_index],
            VCNEBState(q=x_atoms @ np.linalg.inv(reference_cell), deform=deform),
            reference_cell,
            wrap_positions=wrap_positions,
        )
    return images


def project_path_onto_modes(
    images: Sequence[Atoms],
    reference: Atoms,
    modes: Mode | Iterable[Mode],
    *,
    cell_scale: float | None = None,
    normalize: bool = True,
) -> Array:
    """Project every path image onto one or more mode vectors.

    Returns an array with shape ``(n_images, n_modes)``.  Atomic components use
    Angstrom and cell components use ``cell_scale * (F-I)`` so the result is
    compatible with the extended VC-NEB coordinate metric.
    """

    if isinstance(modes, Mode):
        mode_list = [modes]
    else:
        mode_list = list(modes)
    if not mode_list:
        raise ValueError("At least one mode is required")
    if any(len(mode.atomic) != len(reference) for mode in mode_list):
        raise ValueError("All modes must match the reference atom count")

    ref_cell = cell_matrix(reference)
    if cell_scale is None:
        cell_scale = abs(np.linalg.det(ref_cell)) ** (1.0 / 3.0)
    if cell_scale <= 0.0:
        raise ValueError("cell_scale must be positive")

    prepared = [
        mode.normalized(cell_scale=cell_scale) if normalize else mode
        for mode in mode_list
    ]
    vectors = []
    for mode in prepared:
        cell = np.zeros((3, 3)) if mode.cell is None else mode.cell
        vectors.append(np.concatenate([mode.atomic.reshape(-1), (cell * cell_scale).reshape(-1)]))
    basis = np.asarray(vectors, dtype=float)
    result = np.zeros((len(images), len(prepared)), dtype=float)
    ref_state = state_from_atoms(reference, ref_cell)
    ref_x = np.concatenate([
        (ref_state.q @ ref_cell).reshape(-1),
        (ref_state.deform - np.eye(3)).reshape(-1) * cell_scale,
    ])
    for image_index, image in enumerate(images):
        state = state_from_atoms(image, ref_cell)
        x = np.concatenate([
            (state.q @ ref_cell).reshape(-1),
            (state.deform - np.eye(3)).reshape(-1) * cell_scale,
        ])
        result[image_index] = basis @ (x - ref_x)
    return result


def build_mode_basis(
    modes: Mode | Iterable[Mode | Array],
    reference: Atoms,
    *,
    cell_scale: float | None = None,
    normalize: bool = True,
    masses: Array | None = None,
    mass_weighted_input: bool = False,
    remove_translation: bool = False,
) -> Array:
    """Build an extended-coordinate basis for strict mode constraints.

    The returned array has shape ``(3*n_atoms + 9, n_modes)``.  Each column
    uses the same coordinate convention as :class:`vcneb.core.VCNEB`: atomic
    Cartesian components followed by ``cell_scale * deformation``
    components.  ``VCNEB`` orthonormalizes and combines this basis with any
    atom/cell masks before projecting updates.
    """

    if isinstance(modes, Mode):
        mode_list = [modes]
    else:
        mode_list = [mode if isinstance(mode, Mode) else Mode(mode) for mode in modes]
    if not mode_list:
        raise ValueError("At least one mode is required")
    if any(len(mode.atomic) != len(reference) for mode in mode_list):
        raise ValueError("All modes must match the reference atom count")

    ref_cell = cell_matrix(reference)
    if cell_scale is None:
        cell_scale = abs(np.linalg.det(ref_cell)) ** (1.0 / 3.0)
    if cell_scale <= 0.0:
        raise ValueError("cell_scale must be positive")

    prepared = [
        mode.normalized(
            masses=masses,
            mass_weighted_input=mass_weighted_input,
            remove_translation=remove_translation,
            cell_scale=cell_scale,
        )
        if normalize
        else mode
        for mode in mode_list
    ]
    vectors = []
    for mode in prepared:
        cell = np.zeros((3, 3)) if mode.cell is None else mode.cell
        vectors.append(
            np.concatenate([mode.atomic.reshape(-1), (cell * cell_scale).reshape(-1)])
        )
    return np.asarray(vectors, dtype=float).T


def build_direction_basis(directions: Array, *, cell_basis: Array | None = None) -> Array:
    """Build a generalized-coordinate basis from atomic directions.

    ``directions`` may have shape ``(n_atoms, 3)`` for one allowed direction
    per atom, or ``(n_atoms, 3, n_modes)`` for arbitrary linear combinations.
    Optional ``cell_basis`` supplies cell deformation directions with shape
    ``(3, 3)``, ``(9,)`` or ``(3, 3, n_modes)``.  The result can be passed to
    ``VCNEB(mode_basis=..., constraint_mode=...)``.
    """

    atomic = np.asarray(directions, dtype=float)
    if atomic.ndim == 2:
        if atomic.shape[1] != 3:
            raise ValueError("directions must have shape (n_atoms, 3)")
        atomic_basis = np.zeros((3 * len(atomic), len(atomic)), dtype=float)
        for index, direction in enumerate(atomic):
            atomic_basis[3 * index : 3 * index + 3, index] = direction
        atomic = atomic_basis
    elif atomic.ndim == 3 and atomic.shape[1] == 3:
        atomic = atomic.reshape((3 * atomic.shape[0], atomic.shape[2]))
    else:
        raise ValueError("directions must have shape (n_atoms, 3) or (n_atoms, 3, n_modes)")

    if not np.all(np.isfinite(atomic)):
        raise ValueError("directions contains non-finite values")

    if cell_basis is None:
        cell = np.zeros((9, atomic.shape[1]), dtype=float)
    else:
        raw_cell = np.asarray(cell_basis, dtype=float)
        if raw_cell.shape == (3, 3):
            raw_cell = raw_cell.reshape(9, 1)
        elif raw_cell.shape == (9,):
            raw_cell = raw_cell.reshape(9, 1)
        elif raw_cell.ndim == 3 and raw_cell.shape[:2] == (3, 3):
            raw_cell = raw_cell.reshape(9, raw_cell.shape[2])
        elif raw_cell.ndim == 2 and raw_cell.shape[0] == 9:
            pass
        else:
            raise ValueError("cell_basis must have shape (3, 3), (9,), (9, n_modes), or (3, 3, n_modes)")
        if not np.all(np.isfinite(raw_cell)):
            raise ValueError("cell_basis contains non-finite values")
        if raw_cell.shape[1] not in (1, atomic.shape[1]):
            raise ValueError("cell_basis mode count must be one or match directions")
        cell = np.repeat(raw_cell, atomic.shape[1], axis=1) if raw_cell.shape[1] == 1 else raw_cell

    return np.vstack([atomic, cell])


__all__ = [
    "Mode",
    "mode_guided_path",
    "project_path_onto_modes",
    "build_mode_basis",
    "build_direction_basis",
]
