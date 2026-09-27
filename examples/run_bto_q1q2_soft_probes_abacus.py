"""Evaluate four audited BTO fixed-Q soft-direction probes with ABACUS.

This is a selected local-curvature cross-check, not a 2D conditional PES or a
T-to-C barrier. Existing point directories and physical calculator settings
are reused only under the evaluator's exact resume contract.
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
    bto_transverse_soft_plane,
    load_bto_q1q2_reference,
    sha256,
)
from examples.run_bto_q1q2_frozen_abacus import BASIS, PARAMETERS, PP, _validate_written_case
from vcneb import CalculatorModeEvaluator
from vcneb.abacus import make_ase_abacus_factory


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("preflight", "soft_preflight", "branch_result", "report", "reference",
                 "force_constants", "phonopy_eigenpairs", "grid_result_manifest",
                 "workdir", "pseudo_dir", "basis_dir", "abacus_bin", "output"):
        parser.add_argument("--" + name.replace("_", "-"), type=Path, required=True)
    parser.add_argument("--mpi-ranks", type=int, default=32)
    parser.add_argument("--run-dft", action="store_true")
    return parser.parse_args()


def _check_slurm(ranks: int) -> None:
    if (os.environ.get("SLURM_JOB_PARTITION") != "hfacnormal01"
            or int(os.environ.get("SLURM_JOB_NUM_NODES", "0")) != 1
            or int(os.environ.get("SLURM_NTASKS", "0")) != 1
            or int(os.environ.get("SLURM_CPUS_PER_TASK", "0")) < ranks):
        raise RuntimeError("four DFT probes require one 32-CPU task on hfacnormal01")


def main() -> None:
    args = parse_args()
    if args.mpi_ranks != 32 or args.output.exists() or not args.workdir.is_dir():
        raise ValueError("requires 32 MPI ranks, existing evaluator cache and a new output path")
    prior = json.loads(args.preflight.read_text(encoding="utf-8"))
    soft = json.loads(args.soft_preflight.read_text(encoding="utf-8"))
    branch = json.loads(args.branch_result.read_text(encoding="utf-8"))
    physical_soft_plane = prior.get("kind") == "bto_transverse_soft_conditional_preflight_no_dft"
    if (prior.get("kind") not in {
            "bto_q1q2_conditional_preflight_no_dft",
            "bto_transverse_soft_conditional_preflight_no_dft",
        }
            or soft.get("kind") != (
                "bto_transverse_soft_mixed_eigenvector_probe_preflight_no_dft"
                if physical_soft_plane else
                "bto_fixed_q1q2_soft_mixed_eigenvector_probe_preflight_no_dft")
            or branch.get("kind") != (
                "bto_transverse_soft_variable_cell_conditional_local_candidate_not_PES_or_barrier"
                if physical_soft_plane else
                "bto_fixed_q1q2_variable_cell_conditional_local_candidate_not_T_to_C_barrier")
            or soft.get("status") != "four_soft_direction_geometries_passed_no_dft"
            or soft.get("source_sha256", {}).get("preflight") != sha256(args.preflight)
            or soft.get("source_sha256", {}).get("branch_result") != sha256(args.branch_result)
            or branch.get("evaluator_contract_sha256") is None
            or branch.get("preflight_sha256") != sha256(args.preflight)
            or branch.get("status") != "orthogonal_gradient_and_stress_converged_curvature_unchecked"):
        raise ValueError("soft-direction preflight is not tied to the audited branch and reference")
    loaded = load_bto_q1q2_reference(
        args.report, args.reference, args.force_constants, args.phonopy_eigenpairs,
        strain_metric_weights_amu_A2=np.asarray(prior["strain_metric_weights_amu_A2"], dtype=float),
    )
    if any(soft.get("source_sha256", {}).get(key) != digest
           for key, digest in loaded.source_hashes.items()):
        raise ValueError("soft-direction mode/reference sources changed")
    chart = loaded.chart
    plane = (bto_transverse_soft_plane(
        loaded,
        strain_metric_weights_amu_A2=np.asarray(prior["strain_metric_weights_amu_A2"], dtype=float),
    ) if physical_soft_plane else loaded.plane)
    q = np.asarray(soft["q1_q2_sqrt_amu_A"], dtype=float)
    center = np.asarray(soft["center_coordinates_u_A_eta_voigt"], dtype=float)
    direction = np.asarray(soft["soft_direction_metric_unit_chart"], dtype=float)
    if (q.shape != (2,) or center.shape != (chart.coordinate_count,)
            or direction.shape != center.shape
            or not np.all(np.isfinite(np.r_[q, center, direction]))
            or not np.array_equal(center, np.asarray(branch["coordinates_u_A_eta_voigt"], dtype=float))
            or abs(float(direction @ (plane.metric_weights * direction)) - 1.0) > 1e-8
            or not np.allclose(plane.project(center + direction), q, rtol=0.0, atol=1e-8)):
        raise ValueError("soft direction has a changed center, metric or fixed Q")
    entries = soft["probes"]
    steps = sorted({float(item["step_sqrt_amu_A"]) for item in entries})
    if len(entries) != 4 or len(steps) != 2 or sorted((float(item["step_sqrt_amu_A"]), int(item["sign"]))
                                                 for item in entries) != [(steps[0], -1), (steps[0], 1),
                                                                          (steps[1], -1), (steps[1], 1)]:
        raise ValueError("exactly four signed probes at two distinct steps are required")
    minimum_allowed = float(prior["minimum_allowed_atomic_distance_A"])
    for item in entries:
        coordinates = center + int(item["sign"]) * float(item["step_sqrt_amu_A"]) * direction
        if not np.array_equal(coordinates, np.asarray(item["coordinates_u_A_eta_voigt"], dtype=float)):
            raise ValueError("probe coordinates changed since the no-DFT geometry preflight")
        atoms = chart.to_atoms(coordinates)
        distances = atoms.get_all_distances(mic=True)
        np.fill_diagonal(distances, np.inf)
        if (float(np.min(distances)) < minimum_allowed
                or abs(float(np.min(distances)) - item["minimum_atomic_distance_A"]) > 1e-8
                or abs(float(atoms.get_volume()) - item["volume_A3"]) > 1e-8):
            raise ValueError("probe geometry differs from the no-DFT preflight")
    grid = json.loads(args.grid_result_manifest.read_text(encoding="utf-8"))
    if grid.get("status") != "converged" or grid.get("mpi_ranks") != args.mpi_ranks:
        raise ValueError("frozen-grid calculator reference is not complete")
    parameters = dict(grid["calculator_parameters"])
    if (any(parameters.get(key) != value for key, value in PARAMETERS.items() if key != "kpts")
            or list(parameters.get("kpts", [])) != list(PARAMETERS["kpts"])
            or parameters.get("pp") != PP or parameters.get("basis") != BASIS
            or parameters.get("pseudo_dir") != str(args.pseudo_dir)
            or parameters.get("basis_dir") != str(args.basis_dir)
            or sha256(args.abacus_bin) != grid["abacus_binary_sha256"]):
        raise ValueError("ABACUS/PBE/100 Ry/10 au DZP/4x4x4 calculation contract changed")
    for label, filename in [*[(f"pp_{key}", value) for key, value in PP.items()],
                            *[(f"basis_{key}", value) for key, value in BASIS.items()]]:
        folder = args.pseudo_dir if label.startswith("pp_") else args.basis_dir
        if sha256(folder / filename) != grid["asset_sha256"][label]:
            raise ValueError(f"pseudopotential/orbital asset changed: {label}")
    contract = json.loads((args.workdir / "contract.json").read_text(encoding="utf-8"))
    validation_id = contract.get("validation_id")
    if (contract.get("contract_sha256") != branch["evaluator_contract_sha256"]
            or not isinstance(validation_id, str)
            or not validation_id.startswith("bto-abacus-converged-static-v1:")):
        raise ValueError("existing evaluator cache has a different physical/validation contract")
    center_path = args.workdir / branch["final_evaluation_directory"] / "result.json"
    center_record = json.loads(center_path.read_text(encoding="utf-8"))
    if (center_record.get("contract_sha256") != contract["contract_sha256"]
            or not np.array_equal(np.asarray(center_record["coordinates"], dtype=float), center)):
        raise ValueError("branch center does not match its completed DFT cache point")
    if not args.run_dft:
        print(json.dumps({"status": "read_only_four_probe_input_checks_passed_no_dft",
                          "steps_sqrt_amu_A": steps, "evaluator_contract_sha256": contract["contract_sha256"]}, indent=2))
        return
    _check_slurm(args.mpi_ranks)
    frozen_input = {key: grid["points"][0]["input_sha256"][key] for key in ("INPUT", "KPT")}
    factory = make_ase_abacus_factory(
        parameters=parameters, command=f"mpirun -np {args.mpi_ranks} {args.abacus_bin}",
    )

    def validate(index, atoms, directory):
        hashes = _validate_written_case(directory)
        if hashes["INPUT"] != frozen_input["INPUT"] or hashes["KPT"] != frozen_input["KPT"]:
            raise ValueError("soft probe changed reviewed BTO INPUT/KPT")
        log_path = directory / "OUT.ABACUS" / "running_scf.log"
        log = log_path.read_text(encoding="utf-8", errors="replace")
        if ("charge density convergence is achieved" not in log
                or "!FINAL_ETOT_IS" not in log or "PMI server not found" in log):
            raise ValueError("ABACUS soft probe lacks valid SCF/energy/MPI markers")
        sizes = re.findall(r"\bDSIZE\s*=\s*(\d+)", log)
        if len(sizes) != 1 or int(sizes[0]) != args.mpi_ranks:
            raise ValueError("ABACUS soft probe MPI rank count changed")
        return {"log_sha256": sha256(log_path), "input_sha256": hashes,
                "mpi_dsize": int(sizes[0]), "slurm_job_id": os.environ["SLURM_JOB_ID"]}

    evaluator = CalculatorModeEvaluator(
        chart, factory, workdir=args.workdir, calculator_id=prior["calculator_id"],
        minimum_distance_A=minimum_allowed, result_validator=validate,
        validation_id=validation_id, resume=True,
    )
    results = []
    for item in entries:
        coordinates = np.asarray(item["coordinates_u_A_eta_voigt"], dtype=float)
        energy, gradient = evaluator(coordinates)
        results.append({"step_sqrt_amu_A": float(item["step_sqrt_amu_A"]),
                        "sign": int(item["sign"]), "energy_eV": float(energy),
                        "directional_gradient_eV_per_sqrt_amu_A": float(direction @ gradient)})
    if evaluator.n_new_evaluations != 4:
        raise RuntimeError("expected four new, distinct static DFT evaluations")
    center_energy = float(center_record["enthalpy_eV"])
    curvatures = []
    for step in steps:
        plus = next(item for item in results if item["step_sqrt_amu_A"] == step and item["sign"] == 1)
        minus = next(item for item in results if item["step_sqrt_amu_A"] == step and item["sign"] == -1)
        curvatures.append({
            "step_sqrt_amu_A": step,
            "energy_curvature_eV_per_amu_A2": (plus["energy_eV"] + minus["energy_eV"] - 2 * center_energy) / step**2,
            "gradient_curvature_eV_per_amu_A2": (
                plus["directional_gradient_eV_per_sqrt_amu_A"]
                - minus["directional_gradient_eV_per_sqrt_amu_A"]
            ) / (2 * step),
        })
    output = {
        "kind": ("bto_transverse_soft_mixed_direction_four_DFT_probe_crosscheck_not_PES_or_barrier"
                 if physical_soft_plane else
                 "bto_fixed_q1q2_soft_mixed_direction_four_DFT_probe_crosscheck_not_PES_or_barrier"),
        "status": "four_static_probes_complete_pending_independent_SCF_audit",
        "q1_q2_sqrt_amu_A": q.tolist(),
        "center_energy_eV": center_energy,
        "center_result_sha256": sha256(center_path),
        "probes": results,
        "curvatures": curvatures,
        "n_new_evaluations": evaluator.n_new_evaluations,
        "evaluator_contract_sha256": evaluator.contract["contract_sha256"],
        "preflight_sha256": sha256(args.soft_preflight),
        "runner_sha256": sha256(Path(__file__)),
        "slurm_job_id": os.environ["SLURM_JOB_ID"],
        "limitations": "requires independent per-probe INPUT/KPT/STRU/SCF/force/stress audit; "
                       "one direction at one fixed Q does not establish a continuous conditional PES",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": output["status"], "curvatures": curvatures}, indent=2))


if __name__ == "__main__":
    main()
