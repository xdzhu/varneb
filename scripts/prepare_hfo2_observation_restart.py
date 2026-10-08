"""Promote a complete terminal-run observation to an exact-cache restart.

No atom mapping, periodic relifting, DFT, physical-input changes or submission.
The previous allocation must be terminal. A fresh FIRE state is explicit:
geometry/cache continuation is not an optimizer-state or acceleration claim.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np
from ase.io import read, write

from examples.hfo2_fixed_input_factory import CONTRACT, same_ordered_geometry
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.export_hfo2_chain_observation import replay


TERMINAL = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL", "OUT_OF_MEMORY", "PREEMPTED"}


def terminal_job_record(job_id):
    result = subprocess.run(["sacct", "-j", str(job_id), "-n", "-P",
                             "--format=JobIDRaw,State,ExitCode,End"],
                            check=True, capture_output=True, text=True)
    rows = [line.split("|") for line in result.stdout.splitlines()
            if line.split("|")[0] == str(job_id)]
    if len(rows) != 1 or len(rows[0]) < 4:
        raise ValueError("one authoritative main-allocation sacct record required")
    _, state, exit_code, end = rows[0][:4]
    state = state.split()[0].rstrip("+")
    if state not in TERMINAL or not end or end == "Unknown":
        raise ValueError("source allocation must be confirmed terminal; do not fork a live chain")
    return {"job_id": str(job_id), "state": state, "exit_code": exit_code, "end": end}


def prepare(observation, output, *, job_reader=terminal_job_record):
    if output.exists():
        raise FileExistsError("refusing existing restart namespace")
    record_path = observation / "observation.json"
    report = json.loads(record_path.read_text())
    if (report["status"] != "complete_observation_not_final_result" or report["climb"]
            or report["pressure_GPa"] != 0 or report["k_eV_A2"] != .2
            or report["fmax_target_eV_A"] != .10 or report["n_total_images"] not in (9, 10)
            or report["formula_units"] != 4 or not np.isfinite(report["replayed_fmax_eV_A"])):
        raise ValueError("ordinary fixed-contract Hf4O8 observation required")
    state = job_reader(report["source_job_id"])
    if state["job_id"] != report["source_job_id"] or state["state"] not in TERMINAL:
        raise ValueError("terminal source job identity mismatch")
    images = read(observation / "evaluated_chain.traj", index=":")
    if sha256(observation / "evaluated_chain.traj") != report["evaluated_chain_sha256"]:
        raise ValueError("observation trajectory changed")
    points = report["raw_image_evaluations"]
    if len(images) != report["n_total_images"] or len(points) != len(images):
        raise ValueError("incomplete restart observation")
    workdir = Path(report["source_workdir"])
    names = [f"POSCAR_{i:02d}" for i in range(len(images))]
    complete = [p for p in (workdir / "snapshots").glob("step_*")
                if p.name[5:].isdigit() and sorted(q.name for q in p.glob("POSCAR_*")) == names]
    if not complete or max(int(p.name[5:]) for p in complete) != report["snapshot_step"]:
        raise ValueError("restart must use the latest complete source snapshot")
    for i, (image, point) in enumerate(zip(images, points)):
        source = Path(point["raw_source"])
        snapshot = observation / names[i]
        original = workdir / "snapshots" / f'step_{report["snapshot_step"]:04d}' / names[i]
        if (point["image_index"] != i or image.get_chemical_symbols() != ["Hf"] * 4 + ["O"] * 8
                or sha256(snapshot) != point["snapshot_POSCAR_sha256"]
                or sha256(original) != point["snapshot_POSCAR_sha256"]
                or not same_ordered_geometry(image, read(snapshot, format="vasp"))
                or not same_ordered_geometry(image, read(source / "STRU", format="abacus"))):
            raise ValueError("snapshot/raw ordered geometry evidence differs")
        if (any(sha256(source / name) != expected for name, expected in CONTRACT.items())
                or any(sha256(source / name) != expected for name, expected in point["input_sha256"].items())
                or sha256(source / "OUT.ABACUS/running_scf.log") != point["raw_log_sha256"]):
            raise ValueError("raw SCF evidence or fixed physical contract changed")
        raw = audited_results(source)
        if (not np.isclose(raw["energy"], image.get_potential_energy(), rtol=0, atol=1e-10)
                or not np.allclose(raw["forces"], image.get_forces(), rtol=0, atol=1e-12)
                or not np.allclose(raw["stress"], image.get_stress(), rtol=0, atol=1e-12)):
            raise ValueError("fresh raw E/F/stress differs from the observation")
    _, forces, _, geometry = replay(images)
    fmax = float(np.linalg.norm(forces, axis=1).max())
    if abs(fmax - report["replayed_fmax_eV_A"]) > 1e-10:
        raise ValueError("ordinary force replay differs from the observation")
    if fmax <= .10:
        raise ValueError("ordinary residual already converged; audit the result instead of resubmitting")
    output.mkdir(parents=True, exist_ok=False)
    for image in images:
        image.calc = None
    write(output / "seed.traj", images)
    shutil.copyfile(observation / names[0], output / "initial.vasp")
    shutil.copyfile(observation / names[-1], output / "final.vasp")
    parameters = {"source_directory": points[0]["raw_source"],
                  "seed_static_directories": [p["raw_source"] for p in points]}
    (output / "factory_parameters.json").write_text(json.dumps(parameters, indent=2) + "\n")
    manifest = {
        "purpose": "ordinary_observation_restart", "source_job": state,
        "source_observation_sha256": sha256(record_path), "source_snapshot_step": report["snapshot_step"],
        "source_workdir": str(workdir), "prepare_script_sha256": sha256(Path(__file__)),
        "n_total_images": len(images), "n_active_images": len(images) - 2,
        "fmax_eV_A": .10, "climb": False, "pressure_GPa": 0, "spring_eV_A2": .2,
        "optimizer": "FIRE", "maxstep_A": .02, "restart_fmax_eV_A": fmax,
        "next_optimizer_step_cap": 80, "wall_cap_hours": 24,
        "physical_inputs_changed": False, "old_optimizer_state_reused": False,
        "all_seed_SCFs_reused": True, "new_DFT_calls": 0, "job_submitted": False,
        "initial_geometry_gate": geometry, "raw_image_evaluations": points,
        "seed_file_sha256": {n: sha256(output / n) for n in
                             ("seed.traj", "initial.vasp", "final.vasp", "factory_parameters.json")},
        "limitations": "fresh-state geometry continuation, not acceleration or full FIRE-state restoration",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observation", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = prepare(args.observation, args.output)
    print(json.dumps({k: report[k] for k in ("purpose", "source_job", "restart_fmax_eV_A", "new_DFT_calls")}))


if __name__ == "__main__":
    main()
