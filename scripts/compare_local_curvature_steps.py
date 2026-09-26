"""Compare two independently audited finite-difference Hessian step sizes.

This checks common physical provenance and numerical sensitivity. It never
certifies a BTO conditional minimum or GaN transition state on its own.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare_bto_energy_force_gradients(first: dict, second: dict) -> dict | None:
    """Check whether central-energy derivative error is consistent with O(h²)."""

    keys = ("energy_gradient_per_open_direction_eV_per_sqrt_amu_A",
            "force_gradient_per_open_direction_eV_per_sqrt_amu_A")
    if not all(key in report for report in (first, second) for key in keys):
        return None
    h1, h2 = (float(report["step_sqrt_amu_A"]) for report in (first, second))
    energy1, energy2 = (np.asarray(report[keys[0]], dtype=float)
                        for report in (first, second))
    force1, force2 = (np.asarray(report[keys[1]], dtype=float)
                      for report in (first, second))
    if (any(values.shape != (16,) or not np.isfinite(values).all()
            for values in (energy1, energy2, force1, force2))
            or not np.allclose(force1, force2, atol=1e-10, rtol=0)):
        raise ValueError("BTO energy/force derivative vectors or shared center invalid")
    extrapolated = (h2**2 * energy1 - h1**2 * energy2) / (h2**2 - h1**2)
    error_first = np.abs(energy1 - force1)
    error_second = np.abs(energy2 - force1)
    error_extrapolated = np.abs(extrapolated - force1)
    return {
        "method": "two_step_O_h_squared_central_energy_derivative_extrapolation",
        "max_abs_error_first_eV_per_sqrt_amu_A": float(np.max(error_first)),
        "max_abs_error_second_eV_per_sqrt_amu_A": float(np.max(error_second)),
        "max_abs_error_extrapolated_eV_per_sqrt_amu_A": float(np.max(error_extrapolated)),
        "worst_direction_after_extrapolation": int(np.argmax(error_extrapolated)),
        "interpretation": (
            "A decrease supports finite-step truncation as one contributor; "
            "it is not an independent force-consistency or minimum certificate."
        ),
    }


def compare(first: dict, second: dict, kind: str) -> dict:
    if kind == "bto":
        for key in ("preflight", "refined_result", "branch_replay"):
            if (not first["source_sha256"].get(key)
                    or first["source_sha256"][key] != second["source_sha256"].get(key)):
                raise ValueError(f"BTO physical source differs: {key}")
        if (first["q1_q2_sqrt_amu_A"] != second["q1_q2_sqrt_amu_A"]
                or first["axis_kind"] != second["axis_kind"]
                or first["n_matched_signed_probes"] != 32
                or second["n_matched_signed_probes"] != 32):
            raise ValueError("BTO fixed-Q point or signed probe count differs")
        step_key = "step_sqrt_amu_A"
        eig_key = "eigenvalues_eV_per_amu_A2"
        curvature_key = "energy_hessian_diagonal_max_abs_difference_eV_per_amu_A2"
        gradient_key = "energy_gradient_max_abs_difference_eV_per_sqrt_amu_A"
        unit = "eV/(amu Å²)"
        expected_count = 16
    elif kind == "gan":
        if (not first["source_sha256"].get("center_OUTCAR")
                or first["source_sha256"]["center_OUTCAR"]
                != second["source_sha256"].get("center_OUTCAR")
                or first["pressure_GPa"] != second["pressure_GPa"]
                or first["encut_eV"] != second["encut_eV"]
                or first["n_static_displacements"] != 36
                or second["n_static_displacements"] != 36
                or not np.isclose(
                    first["center_enthalpy_eV_per_cell"],
                    second["center_enthalpy_eV_per_cell"], atol=1e-9, rtol=0,
                )):
            raise ValueError("GaN Hessians do not share an audited static center and contract")
        step_key = "step_A"
        eig_key = "translation_free_eigenvalues_eV_per_A2"
        curvature_key = "energy_hessian_diagonal_max_abs_difference_eV_per_A2"
        gradient_key = "energy_gradient_max_abs_difference_eV_per_A"
        unit = "eV/Å²"
        expected_count = 15
    else:
        raise ValueError("kind must be bto or gan")
    steps = [float(report[step_key]) for report in (first, second)]
    spectra = [np.asarray(report[eig_key], dtype=float) for report in (first, second)]
    if (not all(np.isfinite(step) and step > 0 for step in steps)
            or np.isclose(steps[0], steps[1], atol=0, rtol=1e-12)
            or any(values.shape != (expected_count,) or not np.isfinite(values).all()
                   or not np.all(np.diff(values) >= 0) for values in spectra)):
        raise ValueError("step sizes or ordered eigenspectra are invalid")
    cutoffs = (0.0, 0.01, 0.05, 0.1) if kind == "gan" else (0.0, 0.002, 0.005, 0.01)
    negative_counts = {
        str(cutoff): [int(np.count_nonzero(values < -cutoff)) for values in spectra]
        for cutoff in cutoffs
    }
    smallest_magnitude = min(float(np.min(np.abs(values))) for values in spectra)
    mismatch = [float(report[curvature_key]) for report in (first, second)]
    if any(not np.isfinite(value) or value < 0 for value in mismatch):
        raise ValueError("energy-gradient curvature diagnostics are invalid")
    result = {
        "status": "two_step_local_curvature_comparison_not_minimum_or_TS_certificate",
        "kind": kind,
        "step_sizes": steps,
        "curvature_unit": unit,
        "lowest_eigenvalues": [float(values[0]) for values in spectra],
        "next_eigenvalues": [float(values[1]) for values in spectra],
        "maximum_sorted_eigenvalue_change": float(np.max(np.abs(spectra[0] - spectra[1]))),
        "lowest_eigenvalue_change": float(abs(spectra[0][0] - spectra[1][0])),
        "negative_counts_by_cutoff": negative_counts,
        "negative_index_stable_by_cutoff": {
            cutoff: counts[0] == counts[1] for cutoff, counts in negative_counts.items()
        },
        "energy_gradient_max_differences": [float(report[gradient_key]) for report in (first, second)],
        "energy_vs_force_diagonal_curvature_max_differences": mismatch,
        "smallest_absolute_curvature_below_observed_diagonal_mismatch": (
            smallest_magnitude <= max(mismatch)
        ),
        "limitations": [
            "The largest diagonal energy-force mismatch is a diagnostic scale, not a rigorous eigenvalue error bound.",
            "Sorted eigenvalue differences do not track eigenvector identity across crossings.",
            "A second step alone cannot prove global branch selection, basin connections or finite-q stability.",
        ],
    }
    if kind == "bto":
        diagnostic = compare_bto_energy_force_gradients(first, second)
        if diagnostic is not None:
            result["energy_force_gradient_step_extrapolation"] = diagnostic
    return result


def compare_archived_modes(first_npz: Path, second_npz: Path,
                           first: dict, second: dict) -> dict:
    """Compare GaN joint eigendirections when both Hessian archives exist."""

    with np.load(first_npz) as archive:
        h1 = np.asarray(archive["hessian"], dtype=float)
        e1 = np.asarray(archive["eigenvalues"], dtype=float)
        v1 = np.asarray(archive["eigenvectors"], dtype=float)
    with np.load(second_npz) as archive:
        h2 = np.asarray(archive["hessian"], dtype=float)
        e2 = np.asarray(archive["eigenvalues"], dtype=float)
        v2 = np.asarray(archive["eigenvectors"], dtype=float)
    if (h1.shape != (18, 18) or h2.shape != (18, 18)
            or v1.shape != (18, 15) or v2.shape != (18, 15)
            or not np.allclose(e1, first["translation_free_eigenvalues_eV_per_A2"], atol=1e-8, rtol=0)
            or not np.allclose(e2, second["translation_free_eigenvalues_eV_per_A2"], atol=1e-8, rtol=0)
            or not np.allclose(v1.T @ v1, np.eye(15), atol=1e-8, rtol=0)
            or not np.allclose(v2.T @ v2, np.eye(15), atol=1e-8, rtol=0)):
        raise ValueError("GaN Hessian archives disagree with audited spectra")
    return {
        "lowest_mode_absolute_overlap": float(abs(v1[:, 0] @ v2[:, 0])),
        "full_hessian_relative_frobenius_change": float(
            np.linalg.norm(h1 - h2) / max(np.linalg.norm(h1), 1e-30)
        ),
        "first_hessian_npz_sha256": sha256(first_npz),
        "second_hessian_npz_sha256": sha256(second_npz),
        "metric_note": "Mode overlap uses the same declared atomic/symmetric-strain coordinate scale.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kind", choices=("bto", "gan"), required=True)
    for name in ("first", "second", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--first-npz", type=Path)
    parser.add_argument("--second-npz", type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    first = json.loads(args.first.read_text(encoding="utf-8"))
    second = json.loads(args.second.read_text(encoding="utf-8"))
    result = compare(first, second, args.kind)
    if (args.first_npz is None) != (args.second_npz is None):
        raise ValueError("both Hessian archives are required together")
    if args.first_npz is not None:
        if args.kind != "gan":
            raise ValueError("archived eigenvectors are currently supported for GaN only")
        result["archived_mode_comparison"] = compare_archived_modes(
            args.first_npz, args.second_npz, first, second,
        )
    result["source_sha256"] = {"first_audit": sha256(args.first),
                               "second_audit": sha256(args.second)}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "lowest_eigenvalues", "negative_counts_by_cutoff",
        "smallest_absolute_curvature_below_observed_diagonal_mismatch",
    )}))


if __name__ == "__main__":
    main()
