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


def validate_static_endpoint_identity(summary: dict, images: list[Atoms]) -> dict:
    """Require a static summary to match the effective band endpoints exactly.

    The effective endpoints are checked after atom mapping, cell alignment and
    optional translation alignment.  This matters for grid-based calculators:
    a numerically non-invariant translation must not inherit a static gate from
    the untranslated source structure.
    """

    identities = summary.get("endpoint_structures")
    if not isinstance(identities, dict):
        raise ValueError("static endpoint summary has no endpoint_structures")
    if len(images) < 2:
        raise ValueError("endpoint identity validation requires at least two images")
    actual = {
        "initial": endpoint_structure_record(images[0]),
        "final": endpoint_structure_record(images[-1]),
    }
    comparison = compare_endpoint_records(identities, actual)
    issues = [
        f"{label}: static endpoint does not match the mapped/aligned VCNEB endpoint"
        for label, record in comparison["endpoints"].items()
        if not record["matches"]
    ]
    report = {
        "status": "passed" if not issues else "failed",
        "endpoints": comparison["endpoints"],
        "issues": issues,
    }
    if issues:
        raise ValueError("; ".join(issues))
    return report


__all__ = [
    "compare_endpoint_records",
    "endpoint_structure_record",
    "validate_static_endpoint_identity",
]
