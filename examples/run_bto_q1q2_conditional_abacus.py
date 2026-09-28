"""Run one fixed-Q1/Q2, variable-cell BTO conditional-relaxation pilot.

The frozen cubic reference is only a starting geometry. Q1/Q2 and the rigid
translation gauge remain fixed while all other atomic and six symmetric cell
coordinates may relax. This one point is NOT a T-to-C barrier or a complete
2D PES. Every ABACUS call receives a new recoverable directory and is audited
before the optimizer can consume its energy/gradient.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.bto_q1q2_reference import (
    bto_transverse_soft_plane, load_bto_q1q2_reference, sha256,
    validate_bto_gamma_source,
)
from examples.run_bto_q1q2_frozen_abacus import BASIS, PARAMETERS, PP, _validate_written_case
from vcneb import (CalculatorModeEvaluator, audit_orthogonal_curvature,
                   orthogonal_branch_seeds, relax_orthogonal_at_q)
from vcneb.abacus import make_ase_abacus_factory


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--force-constants", type=Path, required=True)
    parser.add_argument("--phonopy-eigenpairs", type=Path, required=True)
    parser.add_argument("--gamma-provenance", type=Path)
    parser.add_argument("--force-sets", type=Path)
    parser.add_argument("--eigenpairs-provenance", type=Path)
    parser.add_argument("--grid-result-manifest", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--pseudo-dir", type=Path, required=True)
    parser.add_argument("--basis-dir", type=Path, required=True)
    parser.add_argument("--abacus-bin", type=Path, required=True)
    parser.add_argument("--mpi-ranks", type=int, default=32)
    parser.add_argument("--gradient-tolerance", type=float, required=True,
                        help="eV per declared metric amplitude, not the NEB fmax")
    parser.add_argument("--orthogonal-amplitude-bound", type=float, required=True)
    parser.add_argument("--max-iterations", type=int, default=16)
    parser.add_argument("--max-new-evaluations", type=int, default=80)
    parser.add_argument("--result-name", default="conditional_result.json",
                        help="new, nonexisting JSON basename for this optimizer stage")
    parser.add_argument("--start-from-result", type=Path,
                        help="audited result JSON in the same workdir for warm continuation")
    parser.add_argument("--max-absolute-stress-kbar", type=float,
                        help="optional raw Cartesian stress gate, e.g. 1.0 kbar")
    parser.add_argument("--curvature-only", action="store_true",
                        help="finite-difference the full fixed-Q orthogonal Hessian at an audited prior result")
    parser.add_argument("--curvature-step", type=float,
                        help="metric amplitude for each +/- orthogonal-gradient probe")
    parser.add_argument("--curvature-result", type=Path,
                        help="completed curvature summary used to authenticate branch seeds")
    parser.add_argument("--curvature-audit", type=Path,
                        help="independently reconstructed audited Hessian with signed unstable vectors")
    parser.add_argument("--branch-direction-index", type=int,
                        help="index into the audited negative directions")
    parser.add_argument("--branch-amplitude", type=float,
                        help="signed seed magnitude in sqrt(amu)*angstrom; both signs are tested")
    parser.add_argument("--branch-trust-radius", type=float, default=0.2,
                        help="maximum initial optimization step in metric amplitude")
    parser.add_argument("--canary-starts-only", action="store_true",
                        help="evaluate the three reviewed soft/soft starting structures, without optimization")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--run-dft", action="store_true", help="otherwise perform read-only input checks")
    return parser.parse_args()


def _require_slurm_allocation(ranks: int) -> None:
    if (os.environ.get("SLURM_JOB_PARTITION") != "hfacnormal01"
            or int(os.environ.get("SLURM_JOB_NUM_NODES", "0")) != 1
            or int(os.environ.get("SLURM_NTASKS", "0")) != 1
            or int(os.environ.get("SLURM_CPUS_PER_TASK", "0")) < ranks):
        raise RuntimeError("DFT requires one Slurm task owning >=MPI-ranks CPUs on hfacnormal01")


def main() -> None:
    args = parse_args()
    if (not np.isfinite(args.gradient_tolerance) or args.gradient_tolerance <= 0
            or not np.isfinite(args.orthogonal_amplitude_bound)
            or args.orthogonal_amplitude_bound <= 0
            or args.max_iterations < 1 or args.max_new_evaluations < 1
            or args.mpi_ranks != 32):
        raise ValueError("pilot tolerances/bounds/iterations must be positive and MPI ranks must be 32")
    if (args.max_absolute_stress_kbar is not None
            and (not np.isfinite(args.max_absolute_stress_kbar) or args.max_absolute_stress_kbar <= 0)):
        raise ValueError("maximum absolute stress must be finite and positive")
    if args.curvature_only and (args.start_from_result is None or args.curvature_step is None):
        raise ValueError("curvature-only requires an audited --start-from-result and --curvature-step")
    if args.curvature_step is not None and (not np.isfinite(args.curvature_step) or args.curvature_step <= 0):
        raise ValueError("curvature-step must be finite and positive")
    branch_options = (args.curvature_result, args.curvature_audit,
                      args.branch_direction_index, args.branch_amplitude)
    branch_requested = any(value is not None for value in branch_options)
    if branch_requested and (any(value is None for value in branch_options)
                             or args.start_from_result is None or args.curvature_only):
        raise ValueError("branch search requires prior result, curvature result/audit, direction index and amplitude")
    if not np.isfinite(args.branch_trust_radius) or args.branch_trust_radius <= 0:
        raise ValueError("branch-trust-radius must be finite and positive")
    if (Path(args.result_name).name != args.result_name or not args.result_name.endswith(".json")
            or args.result_name.startswith(".")):
        raise ValueError("result-name must be a simple visible .json basename")
    if args.workdir.exists() and not args.resume:
        raise FileExistsError(f"refusing to overwrite existing workdir: {args.workdir}")
    if args.resume and not args.workdir.is_dir():
        raise FileNotFoundError("resume requires an existing evaluator workdir")
    result_path = args.workdir / args.result_name
    if result_path.exists():
        raise FileExistsError(f"this optimizer-stage result already exists: {result_path}")
    if args.start_from_result is not None and not args.resume:
        raise ValueError("start-from-result requires --resume")
    preflight = json.loads(args.preflight.read_text(encoding="utf-8"))
    lock_soft_y = preflight.get("kind") == "bto_symmetry_restricted_soft_qy_zero_preflight_no_dft"
    physical_soft_plane = (preflight.get("kind") == "bto_transverse_soft_conditional_preflight_no_dft"
                           or lock_soft_y)
    if not physical_soft_plane and preflight.get("kind") != "bto_q1q2_conditional_preflight_no_dft":
        raise ValueError("missing reviewed BTO conditional preflight")
    if args.canary_starts_only and (not physical_soft_plane or lock_soft_y or args.curvature_only
                                   or branch_requested or args.start_from_result is not None):
        raise ValueError("starts-only canary requires the transverse-soft preflight without continuation")
    paths = {
        "report": args.report,
        "reference": args.reference,
        "force_constants": args.force_constants,
        "phonopy_eigenpairs": args.phonopy_eigenpairs,
        "grid_result_manifest": args.grid_result_manifest,
    }
    if physical_soft_plane:
        if any(path is None for path in (args.gamma_provenance, args.force_sets,
                                         args.eigenpairs_provenance)):
            raise ValueError("transverse-soft preflight requires the 1x1x1 Gamma source files")
        paths.update({
            "gamma_provenance": args.gamma_provenance,
            "force_sets": args.force_sets,
            "eigenpairs_provenance": args.eigenpairs_provenance,
        })
        validate_bto_gamma_source(args.gamma_provenance, args.force_sets,
                                  args.eigenpairs_provenance)
    for key, path in paths.items():
        if sha256(path) != preflight["source_sha256"][key]:
            raise ValueError(f"conditional source hash changed: {key}")
    loaded = load_bto_q1q2_reference(
        args.report, args.reference, args.force_constants, args.phonopy_eigenpairs,
        strain_metric_weights_amu_A2=np.asarray(preflight["strain_metric_weights_amu_A2"], dtype=float),
    )
    chart = loaded.chart
    plane = (bto_transverse_soft_plane(
        loaded, strain_metric_weights_amu_A2=np.asarray(preflight["strain_metric_weights_amu_A2"], dtype=float),
    ) if physical_soft_plane else loaded.plane)
    q_key = "q_parallel_q_transverse_sqrt_amu_A" if physical_soft_plane else "q1_q2_sqrt_amu_A"
    q = np.asarray(preflight[q_key], dtype=float)
    if q.shape != (2,) or not np.all(np.isfinite(q)) or plane.reference_id != preflight["reference_id"]:
        raise ValueError("conditional Q coordinates or mode reference changed")
    expected_open = chart.coordinate_count - 5 - int(lock_soft_y)
    if (chart.coordinate_count != preflight["n_total_coordinates"]
            or preflight["n_relaxed_orthogonal_coordinates"] != expected_open):
        raise ValueError("coordinate or gauge dimensions changed")
    frozen_start = (preflight["branch_starts"][0]["coordinates_u_A_eta_voigt"] if physical_soft_plane
                    else preflight["starting_coordinates_u_A_eta_voigt"])
    if not np.allclose(plane.frozen_coordinates(q), frozen_start,
                       atol=1e-10, rtol=0.0):
        raise ValueError("starting coordinate chart differs from preflight")
    starts = ()
    frozen_directions = chart.rigid_translation_directions()
    if physical_soft_plane:
        branch_starts = preflight.get("branch_starts", [])
        expected_labels = (["frozen"] if lock_soft_y else ["frozen", "+Q_y", "-Q_y"])
        if [row.get("label") for row in branch_starts] != expected_labels:
            raise ValueError("transverse-soft preflight has the wrong third-soft-mode seed contract")
        starts = tuple(np.asarray(row["coordinates_u_A_eta_voigt"], dtype=float)
                       for row in branch_starts[1:])
        third = np.asarray(preflight["remaining_soft_y_metric_unit_direction"], dtype=float)
        amplitude = float(preflight["seed_amplitude_sqrt_amu_A"])
        if (third.shape != (chart.coordinate_count,) or not np.all(np.isfinite(third))
                or not np.isfinite(amplitude) or amplitude <= 0.0
                or abs(np.sqrt(np.dot(plane.metric_weights * third, third)) - 1.0) > 1e-8
                or np.max(np.abs(plane.axis_vectors.T @ (plane.metric_weights * third))) > 1e-8
                or np.max(np.abs(frozen_directions.T @ (plane.metric_weights * third))) > 1e-8):
            raise ValueError("transverse-soft branch seeds differ from the reviewed signed direction")
        if lock_soft_y:
            if (preflight.get("n_fixed_third_soft_axes") != 1
                    or preflight.get("third_soft_mode_constraint")
                    != "Q_y=0 at every evaluation; transverse stability is not implied"
                    or not np.allclose(branch_starts[0]["coordinates_u_A_eta_voigt"],
                                       plane.frozen_coordinates(q), rtol=0.0, atol=1e-10)):
                raise ValueError("the symmetry-restricted branch must start at Q_y=0")
            frozen_directions = np.column_stack([frozen_directions, third])
        elif (not np.allclose(starts[0], plane.frozen_coordinates(q) + amplitude * third,
                              rtol=0.0, atol=1e-10)
              or not np.allclose(starts[1], plane.frozen_coordinates(q) - amplitude * third,
                                  rtol=0.0, atol=1e-10)):
            raise ValueError("transverse-soft branch seeds differ from the reviewed signed direction")
        if any(start.shape != (chart.coordinate_count,)
               or not np.allclose(plane.project(start), q, rtol=0.0, atol=1e-10)
               for start in starts):
            raise ValueError("transverse-soft branch starts changed fixed Q")
    old_summary = None
    if args.start_from_result is not None:
        if args.start_from_result.resolve().parent != args.workdir.resolve():
            raise ValueError("warm-start result must belong to the same evaluator workdir")
        old_summary = json.loads(args.start_from_result.read_text(encoding="utf-8"))
        expected_kind = ("bto_symmetry_restricted_qy_zero_variable_cell_local_candidate_not_PES_or_barrier"
                         if lock_soft_y else
                         "bto_transverse_soft_variable_cell_conditional_local_candidate_not_PES_or_barrier"
                         if physical_soft_plane else
                         "bto_fixed_q1q2_variable_cell_conditional_local_candidate_not_T_to_C_barrier")
        if (old_summary.get("kind") != expected_kind
                or old_summary.get("preflight_sha256") != sha256(args.preflight)
                or not np.allclose(old_summary.get("q1_q2_sqrt_amu_A"), q, atol=1e-12, rtol=0.0)):
            raise ValueError("warm-start result has a different BTO/Q/preflight contract")
        start = np.asarray(old_summary["coordinates_u_A_eta_voigt"], dtype=float)
        if (start.shape != (chart.coordinate_count,) or not np.all(np.isfinite(start))
                or not np.allclose(plane.project(start), q, atol=1e-8, rtol=0.0)):
            raise ValueError("warm-start coordinates are invalid or change Q")
        if lock_soft_y and abs(float(third @ (plane.metric_weights * start))) > 1e-8:
            raise ValueError("warm-start coordinates leave the Q_y=0 restricted sheet")
        if old_summary.get("evaluator_contract_sha256") is None:
            raise ValueError("warm-start result lacks an evaluator contract hash")
        cached_matches = []
        for directory in args.workdir.glob("eval-*-*"):
            cached_path = directory / "result.json"
            if cached_path.is_file():
                cached = json.loads(cached_path.read_text(encoding="utf-8"))
                if np.array_equal(np.asarray(cached["coordinates"], dtype=float), start):
                    cached_matches.append((directory, cached))
        if (len(cached_matches) != 1 or cached_matches[0][1]["contract_sha256"] !=
                old_summary["evaluator_contract_sha256"]):
            raise ValueError("warm-start coordinates must match one completed cached DFT point")
        starts = (start,)
    grid = json.loads(args.grid_result_manifest.read_text(encoding="utf-8"))
    if grid.get("status") != "converged" or grid.get("mpi_ranks") != args.mpi_ranks:
        raise ValueError("the reviewed frozen grid is incomplete or used different MPI ranks")
    parameters = dict(grid["calculator_parameters"])
    if any(parameters.get(key) != value for key, value in PARAMETERS.items() if key != "kpts"):
        raise ValueError("calculator settings differ from the reviewed BTO frozen grid")
    if list(parameters.get("kpts", [])) != list(PARAMETERS["kpts"]):
        raise ValueError("calculator k mesh differs from reviewed BTO frozen grid")
    if parameters.get("pp") != PP or parameters.get("basis") != BASIS:
        raise ValueError("pseudopotential or orbital map changed")
    if parameters.get("pseudo_dir") != str(args.pseudo_dir) or parameters.get("basis_dir") != str(args.basis_dir):
        raise ValueError("asset directories differ from the reviewed grid")
    for label, filename in [*[(f"pp_{k}", v) for k, v in PP.items()],
                            *[(f"basis_{k}", v) for k, v in BASIS.items()]]:
        directory = args.pseudo_dir if label.startswith("pp_") else args.basis_dir
        if sha256(directory / filename) != grid["asset_sha256"][label]:
            raise ValueError(f"reviewed ABACUS asset hash changed: {label}")
    if sha256(args.abacus_bin) != grid["abacus_binary_sha256"]:
        raise ValueError("ABACUS executable changed from the reviewed frozen grid")
    frozen_input = {key: grid["points"][0]["input_sha256"][key] for key in ("INPUT", "KPT")}
    cubic = [item for item in grid["points"] if abs(item["q1"]) < 1e-12 and abs(item["q2"]) < 1e-12]
    if len(cubic) != 1:
        raise ValueError("reviewed grid lacks a unique cubic reference")
    energy_reference = float(cubic[0]["energy_eV"])
    if old_summary is not None and abs(
        cached_matches[0][1]["enthalpy_eV"] - energy_reference - old_summary["energy_minus_c_eV_per_BTO"]
    ) > 1e-8:
        raise ValueError("warm-start energy differs from its completed cached DFT point")
    if args.curvature_only and (
        old_summary.get("status") != "orthogonal_gradient_and_stress_converged_curvature_unchecked"
        or not old_summary.get("stress_target_passed")
    ):
        raise ValueError("curvature-only requires a previously audited gradient/stress-converged result")
    branch_provenance = None
    if branch_requested:
        if physical_soft_plane:
            raise ValueError("curvature-guided branch continuation is not yet supported for the transverse-soft pilot")
        curvature = json.loads(args.curvature_result.read_text(encoding="utf-8"))
        audit = json.loads(args.curvature_audit.read_text(encoding="utf-8"))
        if (audit.get("status") != "negative_orthogonal_curvature_rejects_conditional_local_minimum"
                or audit.get("source_sha256", {}).get("curvature_result") != sha256(args.curvature_result)
                or audit.get("source_sha256", {}).get("refined_result") != sha256(args.start_from_result)
                or curvature.get("start_from_result_sha256") != sha256(args.start_from_result)
                or curvature.get("evaluator_contract_sha256") != old_summary["evaluator_contract_sha256"]
                or audit.get("n_audited_DFT_points") != 48
                or audit.get("n_matched_signed_probes") != 32
                or not np.allclose(audit.get("q1_q2_sqrt_amu_A"), q, atol=1e-12, rtol=0.0)):
            raise ValueError("branch seed lacks the matching audited BTO fixed-Q curvature chain")
        index = args.branch_direction_index
        negative = audit["negative_directions"]
        if index < 0 or index >= len(negative) or negative[index]["eigenvalue_eV_per_amu_A2"] >= 0:
            raise ValueError("branch-direction-index must select a negative curvature direction")
        starts = orthogonal_branch_seeds(
            plane, starts[0], negative[index]["metric_normalized_chart_direction"],
            amplitude=args.branch_amplitude,
            frozen_directions=frozen_directions,
            orthogonal_amplitude_bound=args.orthogonal_amplitude_bound,
        )
        seed_distances = []
        for seed in starts:
            atoms = chart.to_atoms(seed)
            distances = atoms.get_all_distances(mic=True)
            np.fill_diagonal(distances, np.inf)
            minimum_distance = float(np.min(distances))
            seed_distances.append(minimum_distance)
            if minimum_distance < preflight["minimum_allowed_atomic_distance_A"]:
                raise ValueError("curvature branch seed violates minimum atomic distance")
        branch_provenance = {
            "curvature_result_sha256": sha256(args.curvature_result),
            "curvature_audit_sha256": sha256(args.curvature_audit),
            "negative_direction_index": index,
            "negative_eigenvalue_eV_per_amu_A2": negative[index]["eigenvalue_eV_per_amu_A2"],
            "seed_amplitude_sqrt_amu_A": args.branch_amplitude,
            "seed_signs": [1, -1],
            "seed_minimum_atomic_distances_A": seed_distances,
            "initial_trust_radius_sqrt_amu_A": args.branch_trust_radius,
        }
    validation_id = f"bto-abacus-converged-static-v1:{sha256(Path(__file__))}"
    if args.resume:
        recorded = json.loads((args.workdir / "contract.json").read_text(encoding="utf-8"))
        validation_id = recorded.get("validation_id", "")
        if not isinstance(validation_id, str) or not validation_id.startswith("bto-abacus-converged-static-v1:"):
            raise ValueError("existing evaluator has an unrecognized validation contract")
        if old_summary is not None:
            origin_summary = old_summary
            seen_hashes = set()
            while origin_summary.get("start_from_result_sha256") is not None:
                parent_hash = origin_summary["start_from_result_sha256"]
                if parent_hash in seen_hashes:
                    raise ValueError("warm-start provenance contains a cycle")
                seen_hashes.add(parent_hash)
                parent_paths = [path for path in args.workdir.glob("*.json") if sha256(path) == parent_hash]
                if len(parent_paths) != 1:
                    raise ValueError("warm-start parent result hash has no unique file")
                origin_summary = json.loads(parent_paths[0].read_text(encoding="utf-8"))
                if (origin_summary.get("kind") != old_summary["kind"]
                        or origin_summary.get("evaluator_contract_sha256") !=
                        old_summary["evaluator_contract_sha256"]):
                    raise ValueError("warm-start provenance changed physical contract")
            if (validation_id != f"bto-abacus-converged-static-v1:{origin_summary.get('runner_sha256')}"
                    or recorded.get("contract_sha256") != old_summary.get("evaluator_contract_sha256")):
                raise ValueError("warm-start source or evaluator contract differs from cached DFT results")
    if not args.run_dft:
        print(json.dumps({"status": "read_only_checks_passed_no_dft", "q": q.tolist(),
                          "calculator_id": preflight["calculator_id"]}, indent=2))
        return
    _require_slurm_allocation(args.mpi_ranks)
    command = f"mpirun -np {args.mpi_ranks} {args.abacus_bin}"
    factory = make_ase_abacus_factory(parameters=parameters, command=command)

    def validate(index, atoms, directory):
        hashes = _validate_written_case(directory)
        if hashes["INPUT"] != frozen_input["INPUT"] or hashes["KPT"] != frozen_input["KPT"]:
            raise ValueError("conditional point changed the reviewed BTO INPUT/KPT contract")
        log_path = directory / "OUT.ABACUS" / "running_scf.log"
        log = log_path.read_text(encoding="utf-8", errors="replace")
        if ("charge density convergence is achieved" not in log
                or "!FINAL_ETOT_IS" not in log or "PMI server not found" in log):
            raise ValueError(f"ABACUS SCF/energy/MPI markers failed: {log_path}")
        sizes = re.findall(r"\bDSIZE\s*=\s*(\d+)", log)
        if len(sizes) != 1 or int(sizes[0]) != args.mpi_ranks:
            raise ValueError(f"ABACUS MPI size differs from requested {args.mpi_ranks}: {log_path}")
        return {"log_sha256": sha256(log_path), "input_sha256": hashes,
                "mpi_dsize": int(sizes[0]), "slurm_job_id": os.environ["SLURM_JOB_ID"]}

    evaluator = CalculatorModeEvaluator(
        chart, factory, workdir=args.workdir,
        calculator_id=preflight["calculator_id"],
        minimum_distance_A=float(preflight["minimum_allowed_atomic_distance_A"]),
        result_validator=validate,
        validation_id=validation_id,
        resume=args.resume,
    )

    def relative_energy_and_gradient(coordinates):
        if lock_soft_y and abs(float(third @ (plane.metric_weights * coordinates))) > 1e-8:
            raise ValueError("candidate left the Q_y=0 restricted sheet before ABACUS")
        if evaluator.n_new_evaluations >= args.max_new_evaluations:
            raise RuntimeError("maximum new DFT evaluations reached; existing point directories are preserved")
        energy, gradient = evaluator(coordinates)
        return energy - energy_reference, gradient

    if args.canary_starts_only:
        if args.max_new_evaluations < 3:
            raise ValueError("three-start canary requires a budget of at least three new evaluations")
        rows = []
        for entry in preflight["branch_starts"]:
            coordinates = np.asarray(entry["coordinates_u_A_eta_voigt"], dtype=float)
            relative_energy, gradient = relative_energy_and_gradient(coordinates)
            matched = []
            for directory in args.workdir.glob("eval-*-*"):
                cached_path = directory / "result.json"
                if cached_path.is_file():
                    cached = json.loads(cached_path.read_text(encoding="utf-8"))
                    if np.array_equal(np.asarray(cached["coordinates"], dtype=float), coordinates):
                        matched.append((directory.name, cached_path, cached))
            if len(matched) != 1 or matched[0][2]["contract_sha256"] != evaluator.contract["contract_sha256"]:
                raise RuntimeError("canary structure must match one validated cached DFT point")
            rows.append({
                "label": entry["label"],
                "energy_minus_c_eV_per_BTO": relative_energy,
                "full_chart_gradient_norm": float(np.linalg.norm(gradient)),
                "evaluation_directory": matched[0][0],
                "cached_result_sha256": sha256(matched[0][1]),
                "coordinates_u_A_eta_voigt": coordinates.tolist(),
            })
        summary = {
            "kind": "bto_transverse_soft_three_start_canary_not_PES_or_barrier",
            "status": "three_start_statics_complete_optimizer_not_run",
            "q_parallel_q_transverse_sqrt_amu_A": q.tolist(),
            "preflight_sha256": sha256(args.preflight),
            "evaluator_contract_sha256": evaluator.contract["contract_sha256"],
            "runner_sha256": sha256(Path(__file__)),
            "slurm_job_id": os.environ["SLURM_JOB_ID"],
            "n_new_evaluations": evaluator.n_new_evaluations,
            "points": rows,
        }
        result_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(json.dumps({"status": summary["status"], "points": [
            {"label": row["label"], "energy_minus_c_eV_per_BTO": row["energy_minus_c_eV_per_BTO"]}
            for row in rows], "n_new_evaluations": evaluator.n_new_evaluations}, indent=2))
        return

    if args.curvature_only:
        if args.max_new_evaluations < 2 * expected_open:
            raise ValueError("curvature budget must cover two probes per open orthogonal direction")
        curvature = audit_orthogonal_curvature(
            plane, starts[0], relative_energy_and_gradient,
            step_amplitude=args.curvature_step,
            frozen_directions=frozen_directions,
        )
        summary = {
            "kind": "bto_fixed_q1q2_single_point_orthogonal_curvature_screen_not_PES_or_barrier",
            "status": "orthogonal_curvature_screen_complete_not_branch_certification",
            "q1_q2_sqrt_amu_A": q.tolist(),
            "curvature_step_sqrt_amu_A": args.curvature_step,
            "orthogonal_hessian_eigenvalues_eV_per_amu_A2": curvature.eigenvalues.tolist(),
            "minimum_eigenvalue_eV_per_amu_A2": curvature.minimum_eigenvalue,
            "n_gradient_evaluations": curvature.n_gradient_evaluations,
            "n_new_evaluations": evaluator.n_new_evaluations,
            "start_from_result_sha256": sha256(args.start_from_result),
            "evaluator_contract_sha256": evaluator.contract["contract_sha256"],
            "runner_sha256": sha256(Path(__file__)),
            "slurm_job_id": os.environ["SLURM_JOB_ID"],
            "limitations": "one finite-difference step and one branch; repeat a second step before claiming stable sign",
        }
        result_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(json.dumps({key: summary[key] for key in (
            "status", "minimum_eigenvalue_eV_per_amu_A2", "n_new_evaluations",
        )}, indent=2))
        return

    def trial_geometry_valid(coordinates):
        try:
            atoms = chart.to_atoms(coordinates)
            distances = atoms.get_all_distances(mic=True)
            np.fill_diagonal(distances, np.inf)
            return float(np.min(distances)) >= preflight["minimum_allowed_atomic_distance_A"]
        except (ValueError, np.linalg.LinAlgError):
            return False

    safeguarded = branch_requested or physical_soft_plane
    result = relax_orthogonal_at_q(
        plane, q, relative_energy_and_gradient,
        starts=starts, include_frozen_start=((physical_soft_plane and old_summary is None)
                                             or not bool(starts)),
        frozen_directions=frozen_directions,
        orthogonal_amplitude_bound=args.orthogonal_amplitude_bound,
        gradient_tolerance=args.gradient_tolerance,
        max_iterations=args.max_iterations,
        optimizer="safeguarded_bfgs" if safeguarded else "L-BFGS-B",
        trial_validator=trial_geometry_valid if safeguarded else None,
        initial_trust_radius=args.branch_trust_radius,
    )
    matched = []
    for directory in args.workdir.glob("eval-*-*"):
        cached = directory / "result.json"
        if cached.is_file():
            record = json.loads(cached.read_text(encoding="utf-8"))
            if np.array_equal(np.asarray(record["coordinates"], dtype=float), result.coordinates):
                matched.append((directory, record))
    if len(matched) != 1:
        raise RuntimeError("optimized geometry must match exactly one complete cached DFT evaluation")
    final_directory, final_record = matched[0]
    if lock_soft_y and old_summary is None:
        start_labels = ["Q_y=0_frozen"]
    elif physical_soft_plane and old_summary is None:
        start_labels = ["frozen", "+Q_y", "-Q_y"]
    elif branch_requested:
        start_labels = ["+curvature", "-curvature"]
    elif old_summary is not None:
        start_labels = ["warm_start"]
    else:
        start_labels = ["frozen"]
    if len(start_labels) != result.n_starts or len(result.branch_outcomes) != result.n_starts:
        raise RuntimeError("optimizer did not retain exactly one outcome per declared start")
    branch_outcomes = []
    for label, outcome in zip(start_labels, result.branch_outcomes):
        matches = []
        for directory in args.workdir.glob("eval-*-*"):
            cached = directory / "result.json"
            if cached.is_file():
                record = json.loads(cached.read_text(encoding="utf-8"))
                if np.array_equal(np.asarray(record["coordinates"], dtype=float), outcome.coordinates):
                    matches.append(directory)
        if len(matches) != 1:
            raise RuntimeError(f"{label} terminal geometry is not one complete cached DFT evaluation")
        branch_outcomes.append({
            "label": label,
            "start_index": outcome.start_index,
            "converged": outcome.converged,
            "energy_minus_c_eV_per_BTO": outcome.energy,
            "orthogonal_gradient_norm_eV_per_sqrt_amu_A": outcome.orthogonal_gradient_norm,
            "n_result_evaluations": outcome.n_evaluations,
            "coordinates_u_A_eta_voigt": outcome.coordinates.tolist(),
            "final_evaluation_directory": matches[0].name,
        })
    maximum_stress_kbar = float(np.max(np.abs(final_record["stress_eV_per_A3"])) * 1602.176634)
    stress_passed = (args.max_absolute_stress_kbar is None
                     or maximum_stress_kbar <= args.max_absolute_stress_kbar)
    status = ("orthogonal_gradient_and_stress_converged_curvature_unchecked" if stress_passed
              and args.max_absolute_stress_kbar is not None else
              "orthogonal_gradient_converged_curvature_unchecked" if stress_passed else
              "orthogonal_gradient_converged_stress_target_failed_curvature_unchecked")
    summary = {
        "kind": ("bto_symmetry_restricted_qy_zero_variable_cell_local_candidate_not_PES_or_barrier"
                 if lock_soft_y else
                 "bto_transverse_soft_variable_cell_conditional_local_candidate_not_PES_or_barrier"
                 if physical_soft_plane else
                 "bto_fixed_q1q2_variable_cell_conditional_local_candidate_not_T_to_C_barrier"),
        "axis_kind": preflight.get("axis_kind", "archived_soft_stable"),
        "status": status,
        "q1_q2_sqrt_amu_A": q.tolist(),
        "q_parallel_q_transverse_sqrt_amu_A": q.tolist() if physical_soft_plane else None,
        "phonon_supercell": preflight.get("phonon_supercell") if physical_soft_plane else None,
        "third_soft_mode_constraint": preflight.get("third_soft_mode_constraint") if lock_soft_y else None,
        "q_y_sqrt_amu_A": (float(third @ (plane.metric_weights * result.coordinates))
                             if lock_soft_y else None),
        "electronic_kpoints": parameters["kpts"],
        "branch_start_labels": start_labels,
        "branch_outcomes": branch_outcomes,
        "energy_minus_c_eV_per_BTO": result.energy,
        "orthogonal_gradient_norm_eV_per_sqrt_amu_A": result.orthogonal_gradient_norm,
        "n_new_evaluations": evaluator.n_new_evaluations,
        "n_result_evaluations": result.n_evaluations,
        "n_starts": result.n_starts,
        "selected_start": result.selected_start,
        "max_iterations": args.max_iterations,
        "gradient_tolerance_eV_per_sqrt_amu_A": args.gradient_tolerance,
        "orthogonal_amplitude_bound_sqrt_amu_A": args.orthogonal_amplitude_bound,
        "maximum_absolute_stress_kbar": maximum_stress_kbar,
        "stress_target_kbar": args.max_absolute_stress_kbar,
        "stress_target_passed": stress_passed,
        "final_evaluation_directory": final_directory.name,
        "start_from_result_sha256": (None if args.start_from_result is None
                                     else sha256(args.start_from_result)),
        "coordinates_u_A_eta_voigt": result.coordinates.tolist(),
        "preflight_sha256": sha256(args.preflight),
        "evaluator_contract_sha256": evaluator.contract["contract_sha256"],
        "evaluator_validation_id": validation_id,
        "runner_sha256": sha256(Path(__file__)),
        "slurm_job_id": os.environ["SLURM_JOB_ID"],
        "branch_seed": branch_provenance,
        "optimizer": "safeguarded_bfgs" if safeguarded else "L-BFGS-B",
        "limitations": "gradient/stress candidate only; local curvature and competing branches require independent audit",
    }
    result_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": summary["status"], "q": q.tolist(),
                      "energy_minus_c_eV_per_BTO": result.energy,
                      "orthogonal_gradient_norm": result.orthogonal_gradient_norm,
                      "maximum_absolute_stress_kbar": maximum_stress_kbar,
                      "n_new_evaluations": evaluator.n_new_evaluations}, indent=2))
    if not stress_passed:
        raise RuntimeError("conditional gradient converged but the declared raw-stress target did not")


if __name__ == "__main__":
    main()
