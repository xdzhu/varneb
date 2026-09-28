"""Audit the BTO cubic endpoint obstruction to a two-soft-mode lower envelope.

This is a read-only analysis of archived Γ eigenpairs and measured branches.
It neither invokes ABACUS nor constructs a new conditional surface.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def anchor_projection(eigenvectors: np.ndarray, masses: np.ndarray) -> dict:
    """Project Ti-Ba z/x/y polar anchors into the unstable triplet."""
    if eigenvectors.shape != (15, 15) or masses.shape != (5,):
        raise ValueError("expected five-atom BTO Γ eigenvectors and masses")
    if not np.allclose(eigenvectors.T @ eigenvectors, np.eye(15), atol=1e-8):
        raise ValueError("mass-weighted eigenvectors are not orthonormal")
    anchors = np.zeros((15, 3))
    for axis in range(3):
        anchors[3 + axis, axis] = 1.0  # Ti
        anchors[axis, axis] = -1.0     # Ba
    coefficients = eigenvectors[:, :3].T @ (
        np.repeat(np.sqrt(masses), 3)[:, None] * anchors
    )
    fixed, _ = np.linalg.qr(coefficients[:, :2])
    y = coefficients[:, 2]
    residual = y - fixed @ (fixed.T @ y)
    if np.linalg.norm(residual) < 1e-8:
        raise ValueError("omitted Ti_y-Ba_y direction is not independent")
    return {
        "anchor_order": ["Ti_z-Ba_z", "Ti_x-Ba_x", "Ti_y-Ba_y"],
        "projected_anchor_gram": (coefficients.T @ coefficients).tolist(),
        "projected_anchor_singular_values": np.linalg.svd(coefficients, compute_uv=False).tolist(),
        "omitted_y_residual_norm": float(np.linalg.norm(residual)),
        "omitted_y_fraction_outside_fixed_zx": float(np.linalg.norm(residual) / np.linalg.norm(y)),
    }


def audit(eigenpairs: Path, force_constants: Path, provenance: Path,
          force_sets: Path, branches: Path, path_projection: Path) -> dict:
    source = json.loads(provenance.read_text(encoding="utf-8"))
    if (source.get("frequencies_unit")
            != "THz; negative values denote imaginary harmonic modes"
            or source.get("force_sets_sha256") != sha256(force_sets)
            or source.get("n_modes") != 15):
        raise ValueError("phonon eigenpair provenance is inconsistent")
    with np.load(eigenpairs, allow_pickle=False) as data:
        frequencies = np.asarray(data["frequencies_thz"], dtype=float)
        eigenvectors = np.asarray(data["eigenvectors_mass_weighted"], dtype=float)
    with np.load(force_constants, allow_pickle=False) as data:
        masses = np.asarray(data["masses_amu"], dtype=float)
    if (frequencies.shape != (15,) or not np.isfinite(frequencies).all()
            or np.any(frequencies[:3] >= 0)
            or np.ptp(frequencies[:3]) > 1e-5):
        raise ValueError("first three cubic Γ modes are not a degenerate imaginary triplet")
    projection = anchor_projection(eigenvectors, masses)
    with branches.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 4:
        raise ValueError("expected the four independently measured branch points")
    lowerings = []
    branch_points = []
    for row in rows:
        low = float(row["selected_E_minus_C_eV_per_BTO"])
        qy_zero = float(row["relaxed_Qy_zero_E_minus_C_eV_per_BTO"])
        delta = 1000 * (qy_zero - low)
        if (not np.isfinite([low, qy_zero, delta]).all() or delta <= 0
                or abs(delta - float(row["selected_branch_lowering_meV_per_BTO"])) > 1e-6):
            raise ValueError("measured lower branch is not consistently below Qy=0")
        lowerings.append(delta)
        branch_points.append({
            "fixed_Qz_Qx_sqrt_amu_A": [float(row["Q_z_sqrt_amu_A"]),
                                           float(row["Q_x_sqrt_amu_A"])],
            "Qy_zero_minus_low_branch_meV_per_BTO": delta,
        })
    with path_projection.open(newline="", encoding="utf-8") as stream:
        images = [row for row in csv.DictReader(stream) if row["record"] == "VCNEB_image"]
    if len(images) != 7:
        raise ValueError("expected the archived seven-image BTO path projection")
    t = [float(images[0][key]) for key in ("q_parallel_sqrt_amu_A", "q_transverse_sqrt_amu_A")]
    c = [float(images[-1][key]) for key in ("q_parallel_sqrt_amu_A", "q_transverse_sqrt_amu_A")]
    if max(abs(value) for value in c) > 1e-10 or max(abs(float(row["q_transverse_sqrt_amu_A"])) for row in images) > 1e-10:
        raise ValueError("archived T-to-C path is not on the Qx=0 line ending at cubic C")
    return {
        "kind": "bto_cubic_gamma_triplet_obstructs_full_two_mode_conditional_lower_envelope_at_C",
        "calculator_contract": "ABACUS/PBE/100 Ry/10 au DZP; Γ phonon 1x1x1; electronic 4x4x4",
        "imaginary_triplet_frequencies_THz_signed": frequencies[:3].tolist(),
        "anchor_projection": projection,
        "archived_T_projection_Qz_Qx_sqrt_amu_A": t,
        "archived_C_projection_Qz_Qx_sqrt_amu_A": c,
        "four_point_branch_lowering_meV_per_BTO": branch_points,
        "branch_lowering_range_meV_per_BTO": [min(lowerings), max(lowerings)],
        "deduction": (
            "At fixed Qz=Qx=0 and the cubic cell, the independent omitted Qy direction "
            "has negative harmonic curvature. Thus C is not a local minimum when Qy is "
            "released, and a fully minimized two-coordinate lower envelope cannot pass "
            "through the C endpoint at E_C. A T-to-C depiction must retain a constrained "
            "branch, include all three soft coordinates, or use a path-adapted cut."
        ),
        "limitations": [
            "The negative-curvature implication is local to the cubic reference, not a global branch map.",
            "The four off-center branch splittings do not establish a continuous conditional surface.",
            "No new DFT point, finite-temperature free energy, or transition barrier is inferred.",
        ],
        "source_sha256": {str(path.name): sha256(path) for path in (
            eigenpairs, force_constants, provenance, force_sets, branches, path_projection
        )},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("eigenpairs", "force-constants", "provenance", "force-sets",
                 "branches", "path-projection", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.eigenpairs, args.force_constants, args.provenance,
                   args.force_sets, args.branches, args.path_projection)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"frequencies_THz": report["imaginary_triplet_frequencies_THz_signed"],
                      "omitted_y_fraction": report["anchor_projection"]["omitted_y_fraction_outside_fixed_zx"],
                      "branch_lowering_meV": report["branch_lowering_range_meV_per_BTO"]}))


if __name__ == "__main__":
    main()
