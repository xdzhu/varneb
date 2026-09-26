"""Compare GaN endpoint-Gamma path projections across two FD step sizes.

Reports invariant degenerate-subspace amplitudes, not arbitrary individual
eigenvector labels or an energy decomposition. The 1x1x1 Gamma force
constants cover atomic coordinates only; cell strain is separate.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


THZ_TO_CM1 = 33.35640951981521


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_family(directory: Path, label: str) -> tuple[dict, np.ndarray, np.ndarray]:
    projection = json.loads((directory / f"{label}_path_projection.json").read_text(encoding="utf-8"))
    with np.load(directory / f"{label}.npz", allow_pickle=False) as data:
        frequency_thz = np.asarray(data["frequencies_thz"], dtype=float)
        basis = np.asarray(data["eigenvectors_mass_weighted"], dtype=float)
    coordinates = np.asarray(projection["normal_coordinates_sqrt_amu_A"], dtype=float)
    if (projection["n_images"] != 29 or projection["n_atoms"] != 4
            or coordinates.shape != (29, 12) or basis.shape != (12, 12)
            or frequency_thz.shape != (12,)):
        raise ValueError(f"incomplete GaN Gamma path projection: {label}")
    np.testing.assert_allclose(
        np.asarray(projection["frequencies_cm1"]),
        frequency_thz * THZ_TO_CM1, rtol=0, atol=1e-8,
    )
    np.testing.assert_allclose(basis.T @ basis, np.eye(12), rtol=0, atol=1e-8)
    return projection, coordinates, basis


def _dominant_groups(projection: dict, coordinates: np.ndarray) -> list[dict]:
    groups = []
    for record in projection["degenerate_mode_subspaces"]:
        indices = [int(index) for index in record["mode_indices"]]
        frequency = float(record["frequency_cm1"])
        if min(indices) < 3 or abs(frequency) < 10.0:
            continue
        groups.append({
            "mode_indices": indices,
            "frequency_cm1": frequency,
            "squared_path_norm": float(np.sum(coordinates[:, indices] ** 2)),
        })
    if len(groups) < 3:
        raise ValueError("fewer than three resolved optical Gamma subspaces")
    return sorted(groups, key=lambda item: -item["squared_path_norm"])[:3]


def audit(directory: Path, path_csv: Path) -> tuple[dict, list[dict]]:
    dft_audit = json.loads((directory / "audit.json").read_text(encoding="utf-8"))
    if dft_audit["status"] != "computed_GaN_endpoint_atomic_Gamma_1x1x1_requires_numerical_interpretation":
        raise ValueError("GaN DFT force-constant audit status is not recognized")
    with path_csv.open(newline="", encoding="utf-8") as handle:
        path_rows = list(csv.DictReader(handle))
    if len(path_rows) != 29 or [int(row["image_index"]) for row in path_rows] != list(range(29)):
        raise ValueError("audited GaN path source data must have 29 ordered images")
    all_rows = []
    phase_reports = {}
    for phase in ("B4", "B1"):
        small, q_small, u_small = _read_family(directory, f"{phase}_d0.01")
        large, q_large, u_large = _read_family(directory, f"{phase}_d0.02")
        groups = _dominant_groups(small, q_small)
        selected = sorted(index for group in groups for index in group["mode_indices"])
        optical = list(range(3, 12))
        total = float(np.sum(q_small[:, optical] ** 2))
        captured = float(np.sum(q_small[:, selected] ** 2) / total)
        residual = np.linalg.norm(np.delete(q_small, selected, axis=1), axis=1)
        if not 0.0 <= captured <= 1.0 + 1e-8:
            raise ValueError("invalid captured atomic-path norm")
        group_reports = []
        for group in groups:
            indices = group["mode_indices"]
            if not any(record["mode_indices"] == indices
                       for record in large["degenerate_mode_subspaces"]):
                raise ValueError("finite-difference step changed a selected degeneracy group")
            norm_small = np.linalg.norm(q_small[:, indices], axis=1)
            norm_large = np.linalg.norm(q_large[:, indices], axis=1)
            overlaps = np.linalg.svd(u_small[:, indices].T @ u_large[:, indices],
                                     compute_uv=False)
            group_reports.append({
                "mode_indices": indices,
                "frequency_cm1_at_0p01": group["frequency_cm1"],
                "mean_frequency_cm1_at_0p02": next(
                    float(record["frequency_cm1"])
                    for record in large["degenerate_mode_subspaces"]
                    if record["mode_indices"] == indices
                ),
                "fraction_of_optical_path_squared_norm": group["squared_path_norm"] / total,
                "maximum_Q_norm_sqrt_amu_A": float(np.max(norm_small)),
                "max_Q_norm_step_size_difference_sqrt_amu_A": float(np.max(np.abs(norm_small - norm_large))),
                "minimum_principal_subspace_overlap": float(np.min(overlaps)),
            })
        for index, path_row in enumerate(path_rows):
            row = {
                "reference_phase": phase,
                "image_index": index,
                "reaction_coordinate_normalized": float(path_row["reaction_coordinate_normalized"]),
                "relative_enthalpy_eV_per_GaN": float(path_row["relative_enthalpy_eV_per_GaN"]),
                "symmetric_strain_xx": float(path_row["symmetric_strain_xx"]),
                "symmetric_strain_yy": float(path_row["symmetric_strain_yy"]),
                "symmetric_strain_zz": float(path_row["symmetric_strain_zz"]),
                "three_group_atomic_residual_sqrt_amu_A": float(residual[index]),
            }
            for rank, group in enumerate(groups, start=1):
                indices = group["mode_indices"]
                row[f"group{rank}_Q_norm_d0p01_sqrt_amu_A"] = float(np.linalg.norm(q_small[index, indices]))
                row[f"group{rank}_Q_norm_d0p02_sqrt_amu_A"] = float(np.linalg.norm(q_large[index, indices]))
            all_rows.append(row)
        phase_reports[phase] = {
            "three_group_captured_atomic_path_squared_norm_fraction": captured,
            "max_three_group_atomic_reconstruction_residual_sqrt_amu_A": float(np.max(residual)),
            "groups_ranked_by_path_amplitude": group_reports,
            "max_complete_basis_reconstruction_residual_sqrt_amu_A": float(
                max(small["reconstruction_residual_mass_weighted_sqrt_amu_A"])
            ),
        }
    source_paths = [directory / "audit.json", path_csv]
    for phase in ("B4", "B1"):
        for suffix in ("d0.01", "d0.02"):
            source_paths.extend((directory / f"{phase}_{suffix}.npz",
                                 directory / f"{phase}_{suffix}_path_projection.json"))
    report = {
        "status": "audited_GaN_endpoint_Gamma_atomic_subspaces_not_TS_modes",
        "route": "GaN B4-to-B1 tetragonal, 45.7 GPa",
        "n_images_total": 29,
        "reference_cell_rule": "B4 or B1 four-atom path endpoint; 1x1x1 q=0",
        "phases": phase_reports,
        "source_sha256": {path.name: sha256(path) for path in source_paths},
        "limitations": (
            "Atomic Gamma subspace amplitudes are not energy contributions; "
            "the cell strain is reported separately. Exact individual modes "
            "within degenerate groups are gauge-dependent. These endpoint "
            "fixed-cell modes do not certify the full-variable-cell saddle, "
            "finite-q instabilities, or polar non-analytic LO-TO splitting."
        ),
    }
    return report, all_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--path-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or args.csv.exists():
        raise FileExistsError("refusing to overwrite GaN Gamma subspace analysis")
    report, rows = audit(args.directory, args.path_csv)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    with args.csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({phase: data["three_group_captured_atomic_path_squared_norm_fraction"]
                      for phase, data in report["phases"].items()}))


if __name__ == "__main__":
    main()
