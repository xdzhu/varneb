"""Immutable VASP input policy and non-mutating, per-evaluation geometry gates.

These checks do not emulate VASP's Bravais classifier.  A calculator-free
preflight must never be advertised as a successful electronic calculation.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from ase import Atoms


# Project default: a looser symmetry-recognition tolerance, not an SCF or
# force-convergence target. It does not certify VASP Bravais consistency.
DEFAULT_VASP_SYMMETRY_PARAMETERS = {"isym": -1, "symprec": 1e-4}
POSCAR_LATTICE_SIGNIFICANT_DIGITS = 17


class VaspInputContractError(ValueError):
    """Input changed or became invalid before launching VASP."""

    def __init__(self, message):
        super().__init__(message if message.startswith("VASP input contract:") else "VASP input contract: " + message)


def poscar_lattice_roundtrip(cell) -> np.ndarray:
    """Return the exact array emitted by VARNEB's binary64 POSCAR policy."""
    array = np.asarray(cell, dtype=float)
    if array.shape != (3, 3) or not np.isfinite(array).all():
        raise VaspInputContractError("finite 3x3 cell required for POSCAR serialization")
    serialized = np.asarray(
        [[float(format(float(value), f".{POSCAR_LATTICE_SIGNIFICANT_DIGITS}g")) for value in row]
         for row in array], dtype=float,
    )
    if not np.array_equal(serialized, array):
        raise VaspInputContractError("17-significant-digit POSCAR lattice did not round-trip exactly")
    return serialized


def rewrite_poscar_lattice_exact(path: str | Path, cell) -> None:
    """Replace only ASE's POSCAR lattice lines with exact-round-trip decimals."""
    target = Path(path)
    lines = target.read_text(encoding="utf-8").splitlines()
    if len(lines) < 5:
        raise VaspInputContractError("serialized POSCAR is truncated")
    try:
        scale = float(lines[1].split()[0])
        existing = np.asarray([[float(value) for value in lines[index].split()]
                               for index in range(2, 5)], dtype=float) * scale
    except (ValueError, IndexError) as exc:
        raise VaspInputContractError("cannot parse serialized POSCAR lattice") from exc
    expected = np.asarray(cell, dtype=float)
    if scale != 1.0:
        raise VaspInputContractError("ASE POSCAR writer emitted an unexpected non-unit scale")
    if existing.shape != (3, 3) or not np.allclose(existing, expected, rtol=0, atol=1e-12):
        raise VaspInputContractError("ASE POSCAR lattice differs materially before precision rewrite")
    roundtripped = poscar_lattice_roundtrip(expected)
    lines[2:5] = ["  " + "  ".join(format(float(value), f".{POSCAR_LATTICE_SIGNIFICANT_DIGITS}g")
                              for value in row) for row in expected]
    temporary = target.with_name(target.name + ".precision.tmp")
    temporary.write_text("\n".join(lines) + "\n", encoding="utf-8")
    parsed = np.asarray([[float(value) for value in lines[index].split()]
                         for index in range(2, 5)], dtype=float)
    if not np.array_equal(parsed, roundtripped):
        temporary.unlink(missing_ok=True)
        raise VaspInputContractError("rewritten POSCAR lattice failed exact round-trip verification")
    temporary.replace(target)


def canonical_parameters(parameters: dict) -> str:
    def encode(value):
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        raise TypeError(f"unsupported VASP parameter: {type(value).__name__}")

    return json.dumps(parameters, sort_keys=True, default=encode, allow_nan=False)


def parameter_digest(parameters: dict) -> str:
    return hashlib.sha256(canonical_parameters(parameters).encode()).hexdigest()


def validate_vasp_image_geometry(
    atoms: Atoms, *, expected_symbols=None, minimum_distance: float | None = None,
    maximum_cell_condition: float = 1e8,
) -> dict:
    """Check an image without wrapping, symmetrizing, rotating or modifying it.

    ``minimum_distance=None`` deliberately permits the coincident VCA
    representation.  Its physical sites must be checked by the VCA adapter.
    This is a safety gate, not a chemistry-specific bonding criterion.
    """
    cell = np.asarray(atoms.cell.array, dtype=float)
    positions = np.asarray(atoms.positions, dtype=float)
    if not len(atoms) or not np.all(atoms.pbc):
        raise VaspInputContractError("VASP image requires atoms and three periodic directions")
    if not np.isfinite(cell).all() or not np.isfinite(positions).all():
        raise VaspInputContractError("VASP image contains nonfinite cell or coordinates")
    if expected_symbols is not None and list(expected_symbols) != atoms.get_chemical_symbols():
        raise VaspInputContractError("VASP image species/order changed from its input contract")
    if not np.isfinite(maximum_cell_condition) or maximum_cell_condition < 1:
        raise VaspInputContractError("maximum_cell_condition must be finite and >= 1")
    volume = float(np.linalg.det(cell))
    condition = float(np.linalg.cond(cell))
    if volume <= 1e-12 or not np.isfinite(condition) or condition > maximum_cell_condition:
        raise VaspInputContractError("VASP image has invalid cell: non-positive volume or ill-conditioned lattice")
    closest = None
    if minimum_distance is not None:
        if not np.isfinite(minimum_distance) or minimum_distance <= 0:
            raise VaspInputContractError("minimum_distance must be finite and positive")
        if len(atoms) > 1:
            distances = atoms.get_all_distances(mic=True)
            closest = float(distances[np.triu_indices(len(atoms), 1)].min())
        # Include collisions with each atom's own periodic replicas.  A Minkowski
        # reduced cell supplies the shortest lattice vector, including skew cells.
        from ase.geometry import minkowski_reduce
        reduced, _ = minkowski_reduce(cell, pbc=True)
        replica = float(np.linalg.norm(reduced, axis=1).min())
        closest = replica if closest is None else min(closest, replica)
        if closest < minimum_distance:
            raise VaspInputContractError(
                f"VASP image minimum distance {closest:.8g} A < {minimum_distance:.8g} A"
            )
    return {"volume_A3": volume, "cell_condition": condition, "minimum_distance_A": closest}
