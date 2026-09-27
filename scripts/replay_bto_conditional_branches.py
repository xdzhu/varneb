"""Replay a completed BTO soft/soft conditional optimizer from audited DFT cache.

The original job's source selected only its lowest converged start. This
read-only replay uses the same optimizer equations but records every start's
terminal state. It refuses any geometry absent from the immutable DFT cache;
it never starts a calculator or silently extrapolates a missing branch.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.bto_q1q2_reference import bto_transverse_soft_plane, load_bto_q1q2_reference, sha256


def _load_replay_optimizer(path: Path):
    spec = importlib.util.spec_from_file_location("vcneb._mode_surface_replay", path)
    if spec is None or spec.loader is None:
        raise ValueError("cannot load explicitly supplied replay optimizer")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.relax_orthogonal_at_q


def cache_only_evaluator(cache: dict[str, tuple[dict, str]], reference_energy: float):
    """Return a gradient oracle that cannot launch DFT or interpolate a cache miss."""

    hits: list[str] = []

    def energy_and_gradient(coordinates: np.ndarray) -> tuple[float, np.ndarray]:
        key = hashlib.sha256(np.asarray(coordinates, dtype=float).tobytes()).hexdigest()
        if key not in cache:
            raise RuntimeError("replay requested a geometry absent from audited DFT cache")
        item, name = cache[key]
        hits.append(name)
        return float(item["enthalpy_eV"] - reference_energy), np.asarray(item["gradient"], dtype=float)

    return energy_and_gradient, hits


def replay(args: argparse.Namespace) -> dict:
    if args.output.exists():
        raise FileExistsError(f"refusing to replace replay report: {args.output}")
    summary_name = getattr(args, "summary_name", "conditional_q060_result.json")
    if (Path(summary_name).name != summary_name
            or not summary_name.endswith(".json") or summary_name.startswith(".")):
        raise ValueError("summary-name must be a simple visible .json basename")
    summary_path = args.workdir / summary_name
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    preflight = json.loads(args.preflight.read_text(encoding="utf-8"))
    audit = json.loads(args.raw_audit.read_text(encoding="utf-8"))
    grid = json.loads(args.grid_result_manifest.read_text(encoding="utf-8"))
    if (preflight.get("kind") != "bto_transverse_soft_conditional_preflight_no_dft"
            or summary.get("kind") != "bto_transverse_soft_variable_cell_conditional_local_candidate_not_PES_or_barrier"
            or summary.get("preflight_sha256") != sha256(args.preflight)
            or summary.get("runner_sha256") != sha256(args.frozen_runner)
            or audit.get("source_sha256", {}).get("summary") != sha256(summary_path)
            or audit.get("n_individually_audited_DFT_points") < 3
            or audit.get("status") != "verified_gradient_stationary_candidate_curvature_and_branches_unchecked"
            or summary.get("optimizer") != "safeguarded_bfgs"
            or summary.get("n_starts") != 3
            or summary.get("start_from_result_sha256") is not None
            or summary.get("branch_start_labels") != ["frozen", "+Q_y", "-Q_y"]):
        raise ValueError("completed pilot, frozen source, or independent raw audit is inconsistent")
    cubic = [point for point in grid["points"]
             if abs(point["q1"]) < 1e-12 and abs(point["q2"]) < 1e-12]
    if len(cubic) != 1:
        raise ValueError("no unique cubic energy reference")
    reference_energy = float(cubic[0]["energy_eV"])
    weights = np.asarray(preflight["strain_metric_weights_amu_A2"], dtype=float)
    loaded = load_bto_q1q2_reference(
        args.report, args.reference, args.force_constants, args.phonopy_eigenpairs,
        strain_metric_weights_amu_A2=weights,
    )
    chart = loaded.chart
    plane = bto_transverse_soft_plane(loaded, strain_metric_weights_amu_A2=weights)
    q = np.asarray(preflight["q_parallel_q_transverse_sqrt_amu_A"], dtype=float)
    starts = [np.asarray(row["coordinates_u_A_eta_voigt"], dtype=float)
              for row in preflight["branch_starts"]]
    if (not np.allclose(plane.project(starts), q, rtol=0, atol=1e-10)
            or not np.allclose(summary["q_parallel_q_transverse_sqrt_amu_A"], q, rtol=0, atol=1e-10)):
        raise ValueError("replay starts changed fixed-Q coordinates")
    audited_directories = {row["directory"] for row in audit["evaluations"]}
    cache = {}
    for path in sorted(args.workdir.glob("eval-*-*/result.json")):
        if path.parent.name not in audited_directories:
            raise ValueError(f"unaudited DFT point in cache: {path}")
        item = json.loads(path.read_text(encoding="utf-8"))
        coordinates = np.asarray(item["coordinates"], dtype=float)
        digest = hashlib.sha256(coordinates.tobytes()).hexdigest()
        if (item.get("coordinate_sha256") != digest
                or item.get("contract_sha256") != summary["evaluator_contract_sha256"]
                or digest in cache):
            raise ValueError(f"invalid or duplicate immutable DFT cache point: {path}")
        cache[digest] = (item, path.parent.name)
    if len(cache) != len(audited_directories):
        raise ValueError("the independent raw audit and cache point counts differ")
    energy_and_gradient, hits = cache_only_evaluator(cache, reference_energy)

    def trial_geometry_valid(coordinates: np.ndarray) -> bool:
        try:
            atoms = chart.to_atoms(coordinates)
            distances = atoms.get_all_distances(mic=True)
            np.fill_diagonal(distances, np.inf)
            return float(np.min(distances)) >= preflight["minimum_allowed_atomic_distance_A"]
        except (ValueError, np.linalg.LinAlgError):
            return False

    optimize = _load_replay_optimizer(args.mode_surface_code)
    result = optimize(
        plane, q, energy_and_gradient, starts=starts[1:], include_frozen_start=True,
        frozen_directions=chart.rigid_translation_directions(),
        orthogonal_amplitude_bound=summary["orthogonal_amplitude_bound_sqrt_amu_A"],
        gradient_tolerance=summary["gradient_tolerance_eV_per_sqrt_amu_A"],
        max_iterations=summary["max_iterations"], optimizer="safeguarded_bfgs",
        trial_validator=trial_geometry_valid,
        initial_trust_radius=args.initial_trust_radius,
    )
    if (len(result.branch_outcomes) != 3 or result.selected_start != summary["selected_start"]
            or not np.array_equal(result.coordinates,
                                  np.asarray(summary["coordinates_u_A_eta_voigt"], dtype=float))
            or abs(result.energy - summary["energy_minus_c_eV_per_BTO"]) > 1e-8):
        raise ValueError("cached replay did not reproduce the completed frozen-source optimizer")
    third = np.asarray(preflight["remaining_soft_y_metric_unit_direction"], dtype=float)
    base = plane.frozen_coordinates(q)
    outcomes = []
    for label, item in zip(summary["branch_start_labels"], result.branch_outcomes):
        key = hashlib.sha256(item.coordinates.tobytes()).hexdigest()
        cached, directory = cache[key]
        outcomes.append({
            "label": label,
            "start_index": item.start_index,
            "converged": item.converged,
            "energy_minus_c_eV_per_BTO": item.energy,
            "orthogonal_gradient_norm_eV_per_sqrt_amu_A": item.orthogonal_gradient_norm,
            "third_soft_y_amplitude_sqrt_amu_A": float(third @ (plane.metric_weights * (item.coordinates - base))),
            "n_result_evaluations": item.n_evaluations,
            "final_evaluation_directory": directory,
            "maximum_atomic_force_eV_per_A": cached["maximum_atomic_force_eV_per_A"],
            "maximum_absolute_stress_kbar": float(
                np.max(np.abs(cached["stress_eV_per_A3"])) * 1602.176634
            ),
            "coordinates_u_A_eta_voigt": item.coordinates.tolist(),
        })
    report = {
        "kind": "read_only_replay_of_BTO_three_soft_branch_outcomes_not_curvature_or_PES_certificate",
        "status": "all_three_branch_outcomes_reproduced_from_raw_audited_cache",
        "selected_start": result.selected_start,
        "n_cached_DFT_points": len(cache),
        "n_replay_cache_hits": len(hits),
        "branch_outcomes": outcomes,
        "source_sha256": {
            "summary": sha256(summary_path), "raw_audit": sha256(args.raw_audit),
            "preflight": sha256(args.preflight), "frozen_runner": sha256(args.frozen_runner),
            "frozen_mode_surface": sha256(args.frozen_mode_surface),
            "mode_surface_replay": sha256(args.mode_surface_code),
            "replay_script": sha256(Path(__file__)),
        },
        "limitations": "Deterministic cache replay only; individual branch curvature and global optimality remain unverified.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("workdir", "preflight", "report", "reference", "force_constants",
                 "phonopy_eigenpairs", "grid_result_manifest", "raw_audit",
                 "frozen_runner", "frozen_mode_surface", "mode_surface_code", "output"):
        parser.add_argument("--" + name.replace("_", "-"), type=Path, required=True)
    parser.add_argument("--summary-name", default="conditional_q060_result.json")
    parser.add_argument("--initial-trust-radius", type=float, default=0.2)
    args = parser.parse_args()
    report = replay(args)
    print(json.dumps({key: report[key] for key in (
        "status", "selected_start", "n_cached_DFT_points", "n_replay_cache_hits",
    )}, indent=2))


if __name__ == "__main__":
    main()
