"""Recompute the GaN highest-image physical force from an audited trajectory.

This is an offline evidence gate, not a TS optimization or Hessian calculation.
The ordinary-NEB tangent residual is deliberately reported separately from
the raw enthalpy-force tangent.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from ase import units

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.analyze_gan_qian_path import latest_evaluated_chain  # noqa: E402
from vcneb import VCNEB  # noqa: E402


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _close(actual: float, expected: float, label: str, tolerance: float = 1e-7) -> None:
    if not np.isclose(actual, expected, atol=tolerance, rtol=0.0):
        raise ValueError(f"{label} disagrees with the original audited result")


def audit(summary_path: Path, trajectory_path: Path, generic_audit_path: Path,
          route_analysis_path: Path) -> dict:
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    generic = json.loads(generic_audit_path.read_text(encoding="utf-8"))
    analysis = json.loads(route_analysis_path.read_text(encoding="utf-8"))
    hashes = {
        "summary": _sha256(summary_path),
        "trajectory": _sha256(trajectory_path),
        "audit": _sha256(generic_audit_path),
        "route_analysis": _sha256(route_analysis_path),
        "auditor": _sha256(Path(__file__)),
    }
    if (summary.get("status") != "completed" or not summary.get("converged")
            or generic.get("status") != "ok" or generic.get("issues")
            or analysis.get("route_requested") != "tetragonal"
            or analysis.get("status") != "numerically_audited_structure_diagnostics_require_interpretation"
            or analysis.get("input_sha256", {}).get("summary") != hashes["summary"]
            or analysis.get("input_sha256", {}).get("trajectory") != hashes["trajectory"]
            or analysis.get("input_sha256", {}).get("audit") != hashes["audit"]
            or summary.get("n_images") != 29
            or summary.get("endpoint_evaluation_policy") != "fixed_cached_once"):
        raise ValueError("GaN source, complete-chain or provenance gate failed")
    pressure = float(summary["path_diagnostics"]["pressure_eV_per_A3"])
    _close(pressure, 45.7 * units.GPa, "external pressure", tolerance=1e-9)
    images = latest_evaluated_chain(trajectory_path, 29)
    if any(image.calc is None for image in images):
        raise ValueError("the saved trajectory lacks static image results")
    chain = VCNEB(
        images, pressure=pressure,
        cell_scale=float(summary["path_diagnostics"]["cell_scale_A"]),
        k=0.2, climb=False,
    )
    candidate = chain.saddle_diagnostics()
    index = int(candidate["image_index"])
    original = summary["path_diagnostics"]["images"][index]
    if (index != 15 or index != analysis["highest_image_index"]
            or not candidate["has_interior_barrier"]
            or candidate["tangent_curvature_eV_per_A2"] >= 0
            or original["is_climbing_image"]):
        raise ValueError("unexpected GaN highest-image or CI state")
    for label, actual, expected in (
        ("true tangent force", candidate["true_tangential_force_eV_per_A"],
         original["true_tangential_force_eV_per_A"]),
        ("true perpendicular force", candidate["true_perpendicular_force_eV_per_A"],
         original["true_perpendicular_force_eV_per_A"]),
        ("true max-vector force", candidate["true_generalized_force_max_vector_eV_per_A"],
         original["max_true_generalized_force_eV_per_A"]),
        ("barrier", candidate["relative_enthalpy_eV"], summary["barrier_enthalpy_eV"]),
        ("path curvature", candidate["tangent_curvature_eV_per_A2"],
         summary["saddle_diagnostics"]["tangent_curvature_eV_per_A2"]),
    ):
        _close(float(actual), float(expected), label)
    true_parallel = float(candidate["true_tangential_force_eV_per_A"])
    true_perpendicular = float(candidate["true_perpendicular_force_eV_per_A"])
    return {
        "status": "audited_highest_image_TS_candidate_not_full_saddle_certificate",
        "route": "GaN B4-to-B1 tetragonal",
        "pressure_GPa": 45.7,
        "n_images_total": 29,
        "n_images_interior": 27,
        "n_formula_units": 2,
        "highest_image_index": index,
        "barrier_eV_per_GaN": float(analysis["barrier_eV_per_GaN"]),
        "ordinary_NEB_residual_tangent_eV_per_A": float(
            summary["saddle_diagnostics"]["tangential_force_eV_per_A"]
        ),
        "physical_enthalpy_force_tangent_eV_per_A": true_parallel,
        "physical_enthalpy_force_perpendicular_norm_eV_per_A": true_perpendicular,
        "physical_enthalpy_force_euclidean_norm_eV_per_A": float(
            np.hypot(true_parallel, true_perpendicular)
        ),
        "physical_enthalpy_force_max_vector_eV_per_A": float(
            candidate["true_generalized_force_max_vector_eV_per_A"]
        ),
        "discrete_path_tangent_curvature_eV_per_A2": float(
            candidate["tangent_curvature_eV_per_A2"]
        ),
        "full_variable_cell_Hessian_index_checked": False,
        "TS_stationarity_refinement_checked": False,
        "interpretation": (
            "A converged ordinary-NEB internal maximum with negative discrete "
            "path curvature is a TS candidate. The near-zero ordinary-NEB tangent "
            "residual is a spring term, not the physical tangent force. Full "
            "variable-cell saddle certification still requires local refinement "
            "and a one-negative-mode Hessian after removing gauge zero modes."
        ),
        "source_sha256": hashes,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--generic-audit", type=Path, required=True)
    parser.add_argument("--route-analysis", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    result = audit(args.summary, args.trajectory, args.generic_audit, args.route_analysis)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "highest_image_index", "ordinary_NEB_residual_tangent_eV_per_A",
        "physical_enthalpy_force_tangent_eV_per_A",
        "full_variable_cell_Hessian_index_checked",
    )}, indent=2))


if __name__ == "__main__":
    main()
