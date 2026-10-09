"""Bounded ordinary refinement of the audited twelve-image PO/M reconstruction."""
import argparse
import json
from pathlib import Path
import tempfile

import numpy as np
from ase.io import read, write

from scripts.analyze_hfo2_G1_peak_sampling import analyze
from scripts.analyze_hfo2_channel_network import read_evaluated_observation
from scripts.prepare_hfo2_G1_peak_sampling import BASE, SOURCES, load_manifest, save
from scripts.prepare_hfo2_observation_restart import terminal_job_record
from scripts.audit_hfo2_static_replica import sha256
from vcneb import validate_periodic_path_lift, validate_path_geometry

SAMPLE_ROOT = "/public/home/iai806/abacus/agent-runs/20261009-varneb-G1-peaks-1438/prepared"


def prepare(repository, completed, output, *, job_reader=terminal_job_record):
    if output.exists():
        raise FileExistsError("refusing existing M refinement namespace")
    states = [job_reader(job) for job in ("28298794", "28380672")]
    if any(s["job_id"] != job or s["state"] != "COMPLETED" or s["exit_code"] != "0:0"
           for s, job in zip(states, ("28298794", "28380672"))):
        raise ValueError("both original M and sampling jobs must be confirmed completed")
    with tempfile.TemporaryDirectory(prefix="hfo2-M-sample-audit-") as folder:
        analysis = analyze(repository, completed, Path(folder) / "audit.json")
    result = analysis["channels"]["PO_to_M"]
    if result["ordinary_residual_passed"]:
        raise ValueError("inserted band already passes; no refinement required")
    m = load_manifest(completed)
    original, raw, _ = read_evaluated_observation(repository / BASE / SOURCES["PO_to_M"]["folder"])
    inserted = sorted((p for p in m["points"] if p["channel"] == "PO_to_M"), key=lambda p: p["fraction"])
    images, sources = [], []
    for i, atoms in enumerate(original):
        images.append(atoms.copy())
        sources.append(raw["raw_image_evaluations"][i]["raw_source"])
        if i == 3:
            for p in inserted:
                images.append(read(completed / p["geometry"], format="vasp"))
                sources.append(f"{SAMPLE_ROOT}/calculations/{p['index']:02d}/scf_000000")
    if len(images) != 12 or len(inserted) != 3:
        raise ValueError("registered twelve-image reconstruction required")
    lift = validate_periodic_path_lift(images)
    geometry = validate_path_geometry(images, minimum_distance=1.6, maximum_deformation=.25)
    output.mkdir(parents=True)
    write(output / "seed.traj", images)
    write(output / "initial.vasp", images[0], format="vasp", direct=True, sort=False)
    write(output / "final.vasp", images[-1], format="vasp", direct=True, sort=False)
    save(output / "factory_parameters.json", {"source_directory": sources[0], "seed_static_directories": sources})
    report = {"purpose": "G1_M_twelve_image_sampling_refinement", "source_jobs": states,
              "source_summary_sha256": analysis["summary_sha256"], "prepare_driver_sha256": sha256(Path(__file__)),
              "n_total_images": 12, "n_active_images": 10, "n_fixed_endpoints": 2,
              "restart_fmax_eV_A": result["replayed_ordinary_fmax_eV_A"], "ordinary_target_eV_A": .10,
              "step_cap": 20, "wall_cap_hours": 8, "maximum_moving_geometry_SCFs": 200,
              "pressure_GPa": 0., "climb": False, "optimizer": "FIRE", "maxstep_A": .02, "spring_eV_A2": .2,
              "all_twelve_initial_E_F_stress_reused": True, "new_DFT_calls_for_preparation": 0,
              "physical_inputs_changed": False, "old_optimizer_state_reused": False,
              "new_independent_channel": False, "G2_or_holdout_submission": False,
              "periodic_lift": lift, "geometry": geometry,
              "seed_file_sha256": {n: sha256(output/n) for n in ("seed.traj", "initial.vasp", "final.vasp", "factory_parameters.json")},
              "limits": "geometry/cache continuation with fresh optimizer, not acceleration or TS certification; audit cap termination before any further segment"}
    save(output / "manifest.json", report)
    return report


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repository", type=Path, default=Path("."))
    p.add_argument("--completed", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(prepare(a.repository, a.completed, a.output)["source_jobs"]))
