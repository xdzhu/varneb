"""Small calculator-independent provenance records for VCNEB inputs."""

from __future__ import annotations

import hashlib
import json

from ase import Atoms
import numpy as np


def _canonical_coordinates(values: np.ndarray) -> np.ndarray:
    """Remove insignificant file-I/O roundoff before a structural hash."""

    canonical = np.round(np.asarray(values, dtype=float), decimals=12)
    canonical[np.abs(canonical) < 0.5e-12] = 0.0
    return canonical


def _digest_payload(symbols: list[str], pbc: np.ndarray, cell: np.ndarray, scaled_positions: np.ndarray) -> str:
    """Return a stable binary digest without locale- or text-format dependence."""

    digest = hashlib.sha256()
    digest.update(json.dumps(symbols, ensure_ascii=True, separators=(",", ":")).encode("ascii"))
    digest.update(np.ascontiguousarray(pbc, dtype=np.uint8).tobytes())
    digest.update(np.ascontiguousarray(cell, dtype="<f8").tobytes())
    digest.update(np.ascontiguousarray(scaled_positions, dtype="<f8").tobytes())
    return digest.hexdigest()


def endpoint_structure_record(atoms: Atoms) -> dict:
    """Record an exact ordered periodic endpoint identity for backend comparison.

    Fractional coordinates are wrapped into ``[0, 1)`` before hashing, so an
    atom written one lattice translation away has the same identity.  Atom
    order is deliberately retained: NEB mapping is order-sensitive and a
    cross-calculator validation must expose an unreviewed reorder rather than
    silently accepting it.
    """

    if len(atoms) == 0:
        raise ValueError("endpoint structure must contain at least one atom")
    if atoms.cell.rank != 3 or atoms.get_volume() <= 0.0:
        raise ValueError("endpoint structure requires a non-singular 3D cell")
    symbols = atoms.get_chemical_symbols()
    cell = _canonical_coordinates(np.asarray(atoms.cell.array, dtype=float))
    pbc = np.asarray(atoms.pbc, dtype=bool)
    scaled_positions = _canonical_coordinates(np.asarray(atoms.get_scaled_positions(wrap=True), dtype=float))
    if not np.all(np.isfinite(cell)) or not np.all(np.isfinite(scaled_positions)):
        raise ValueError("endpoint structure contains non-finite cell or positions")
    composition: dict[str, int] = {}
    for symbol in symbols:
        composition[symbol] = composition.get(symbol, 0) + 1
    return {
        "format_version": 1,
        "sha256": _digest_payload(symbols, pbc, cell, scaled_positions),
        "n_atoms": len(atoms),
        "species_order": symbols,
        "composition": composition,
        "pbc": pbc.astype(bool).tolist(),
        "volume_A3": float(atoms.get_volume()),
        "cell_A": cell.tolist(),
        "fractional_positions_wrapped": scaled_positions.tolist(),
        "identity_convention": "ordered species, PBC, cell_A and wrapped fractional positions rounded to 1e-12",
    }


def compare_endpoint_records(reference: dict, candidate: dict) -> dict:
    """Compare two preflight endpoint-record bundles without reading structures.

    The comparison is intentionally strict: an atom reorder is observable and
    blocks cross-calculator production until the user records an audited common
    mapping. This prevents two physically similar but differently ordered
    inputs from being reported as the same NEB validation case.
    """

    result = {"matches": True, "endpoints": {}}
    for label in ("initial", "final"):
        left = reference.get(label)
        right = candidate.get(label)
        if not isinstance(left, dict) or not isinstance(right, dict):
            raise ValueError("both endpoint bundles require initial and final records")
        left_digest, right_digest = left.get("sha256"), right.get("sha256")
        matches = isinstance(left_digest, str) and left_digest == right_digest
        result["endpoints"][label] = {
            "matches": matches,
            "reference_sha256": left_digest,
            "candidate_sha256": right_digest,
        }
        result["matches"] = bool(result["matches"] and matches)
    return result


__all__ = ["compare_endpoint_records", "endpoint_structure_record"]
