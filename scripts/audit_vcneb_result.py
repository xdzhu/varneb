"""Audit a completed VC-NEB run summary without rerunning a calculator."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import tempfile


EV_A3_TO_KBAR = 1602.176634


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workdir", type=Path, help="VC-NEB work directory containing vcneb_summary.json")
    parser.add_argument(
        "--max-min-distance",
        type=float,
        default=1.6,
        help="Minimum acceptable periodic interatomic distance in Angstrom",
    )
    parser.add_argument(
        "--max-deformation",
        type=float,
        default=None,
        help="Optional maximum acceptable path deformation norm",
    )
    parser.add_argument(
        "--max-stress-kbar",
        type=float,
        default=None,
        help="Optional maximum absolute per-image stress in kbar",
    )
    return parser.parse_args()


def audit(
    summary: dict,
    *,
    minimum_distance: float,
    maximum_deformation: float | None,
    maximum_stress_kbar: float | None = None,
) -> dict:
    diagnostics = summary.get("path_diagnostics") or {}
    image_records = diagnostics.get("images") or []
    geometry = diagnostics.get("geometry") or {}
    geometry_records = geometry.get("images") or []
    issues: list[str] = []
    status = summary.get("status")
    if status not in (None, "completed"):
        issues.append(f"summary status is {status!r}, not completed")
    expected_images = summary.get("n_images", diagnostics.get("n_images"))
    if expected_images is not None and len(image_records) != int(expected_images):
        issues.append(
            f"summary has {len(image_records)} image diagnostics, expected {int(expected_images)}"
        )
    final_force = summary.get("final_max_generalized_force_eV_per_A")
    target = summary.get("fmax_target_eV_per_A")
    if final_force is None or target is None:
        issues.append("summary lacks final generalized force or target")
    elif final_force > target:
        issues.append(f"final generalized force {final_force:.6g} exceeds target {target:.6g}")

    if not image_records:
        issues.append("summary contains no per-image path diagnostics")
    if geometry.get("valid") is False:
        issues.extend(f"geometry: {issue}" for issue in geometry.get("issues", []))
    min_distance = None
    max_deformation = None
    max_stress_kbar = None
    for record in geometry_records:
        distance = record.get("minimum_interatomic_distance_A")
        if distance is not None:
            min_distance = distance if min_distance is None else min(min_distance, distance)
        deformation = record.get("deformation_from_reference_frobenius")
        if deformation is not None:
            max_deformation = deformation if max_deformation is None else max(max_deformation, deformation)
    for record in image_records:
        if record.get("volume_A3", 0.0) <= 0.0:
            issues.append(f"image {record.get('image_index')} has non-positive volume")
        if record.get("max_atom_force_eV_per_A") is None:
            issues.append(f"image {record.get('image_index')} lacks atom-force diagnostic")
        stress_eV_per_A3 = record.get("max_stress_eV_per_A3")
        if stress_eV_per_A3 is None:
            if maximum_stress_kbar is not None:
                issues.append(f"image {record.get('image_index')} lacks stress diagnostic")
        else:
            stress_kbar = float(stress_eV_per_A3) * EV_A3_TO_KBAR
            max_stress_kbar = stress_kbar if max_stress_kbar is None else max(max_stress_kbar, stress_kbar)
    if min_distance is not None and min_distance < minimum_distance:
        issues.append(f"minimum path distance {min_distance:.6g} A is below {minimum_distance:.6g} A")
    if maximum_deformation is not None and max_deformation is not None and max_deformation > maximum_deformation:
        issues.append(f"maximum deformation {max_deformation:.6g} exceeds {maximum_deformation:.6g}")
    if maximum_stress_kbar is not None and max_stress_kbar is not None and max_stress_kbar > maximum_stress_kbar:
        issues.append(f"maximum path stress {max_stress_kbar:.6g} kbar exceeds {maximum_stress_kbar:.6g} kbar")

    return {
        "status": "ok" if not issues else "failed",
        "issues": issues,
        "n_images": diagnostics.get("n_images", summary.get("n_images")),
        "summary_status": status,
        "final_max_generalized_force_eV_per_A": final_force,
        "fmax_target_eV_per_A": target,
        "barrier_enthalpy_eV": summary.get("barrier_enthalpy_eV"),
        "reaction_enthalpy_eV": summary.get("reaction_enthalpy_eV"),
        "has_interior_barrier": diagnostics.get("has_interior_barrier"),
        "highest_image_index": diagnostics.get("highest_image_index"),
        "minimum_path_distance_A": min_distance,
        "maximum_deformation": max_deformation,
        "maximum_stress_kbar": max_stress_kbar,
        "maximum_stress_threshold_kbar": maximum_stress_kbar,
        "geometry_valid": geometry.get("valid"),
        "geometry_issues": geometry.get("issues", []),
    }


def write_report(path: Path, report: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(report, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass


def main() -> int:
    args = parse_args()
    summary_path = args.workdir / "vcneb_summary.json"
    try:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        print(f"[ERROR] cannot read {summary_path}: {exc}", file=sys.stderr)
        return 2
    report = audit(
        summary,
        minimum_distance=args.max_min_distance,
        maximum_deformation=args.max_deformation,
        maximum_stress_kbar=args.max_stress_kbar,
    )
    output = args.workdir / "vcneb_audit.json"
    write_report(output, report)
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0 if report["status"] == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
