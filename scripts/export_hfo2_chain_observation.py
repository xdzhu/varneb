"""Freeze one fully evaluated chain observation without touching a live run.

Every image is matched to its exact ordered periodic SCF. Missing/partial
SCFs, a broken coordinate lift, or a mismatch with the optimizer log fail
closed. This exports evidence, not a restart or a converged-path claim.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.audit_hfo2_static_replica import sha256
from scripts.prepare_hfo2_lifted_restart import exact_cached_source
from vcneb import VCNEB, validate_path_geometry, validate_periodic_path_lift


def optimizer_row(path, step):
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) == 5 and fields[0].endswith(":") and fields[1] == str(step):
            rows.append({"optimizer": fields[0][:-1], "step": step, "time_CST": fields[2],
                         "highest_energy_eV_cell": float(fields[3]), "fmax_eV_A": float(fields[4])})
    if len(rows) != 1:
        raise ValueError("one completed optimizer log row for the requested snapshot required")
    return rows[0]


def replay(images):
    """Reproduce ordinary production forces using numeric-only calculators."""
    validate_periodic_path_lift(images)
    geometry = validate_path_geometry(images, minimum_distance=1.6, maximum_deformation=.25)
    chain = VCNEB(images, k=.2, climb=False, pressure=0.)
    forces = chain.get_forces()
    diagnostics = chain.path_diagnostics()
    return chain, forces, diagnostics, geometry


def export(workdir, step, output, source_job_id):
    if output.exists():
        raise FileExistsError("refusing existing observation namespace")
    if not isinstance(step, int) or isinstance(step, bool) or step < 0:
        raise ValueError("a nonnegative snapshot step required")
    snapshot = workdir / "snapshots" / f"step_{step:04d}"
    paths = sorted(snapshot.glob("POSCAR_*"))
    if len(paths) not in (9, 10, 12) or [p.name for p in paths] != [f"POSCAR_{i:02d}" for i in range(len(paths))]:
        raise ValueError("a complete nine/ten/twelve-image snapshot required")
    before = {p.name: sha256(p) for p in paths}
    row = optimizer_row(workdir / "vcneb.opt.log", step)
    images, evaluations = [], []
    for index, path in enumerate(paths):
        image = read(path, format="vasp")
        if image.get_chemical_symbols() != ["Hf"] * 4 + ["O"] * 8 or not image.pbc.all():
            raise ValueError("ordered periodic Hf4O8 images required")
        source, raw = exact_cached_source(image, workdir / f"image_{index:04d}")
        image.calc = SinglePointCalculator(image, **raw)
        images.append(image)
        evaluations.append({
            "image_index": index, "raw_source": str(source),
            "raw_log_sha256": sha256(source / "OUT.ABACUS/running_scf.log"),
            "input_sha256": {name: sha256(source / name) for name in (*CONTRACT, "STRU")},
            "snapshot_POSCAR_sha256": before[path.name],
            "ordered_periodic_geometry_matched": True,
            "new_SCF_for_export": False,
            "energy_eV_cell": float(raw["energy"]),
            "forces_eV_A": np.asarray(raw["forces"]).tolist(),
            "stress_ASE_voigt_eV_A3": np.asarray(raw["stress"]).tolist(),
        })
    chain, forces, diagnostics, geometry = replay(images)
    fmax = float(np.linalg.norm(forces, axis=1).max())
    if abs(fmax - row["fmax_eV_A"]) > 7e-7:
        raise ValueError("exact-SCF replay differs from the production optimizer force")
    if abs(float(chain.enthalpies.max()) - row["highest_energy_eV_cell"]) > 7e-7:
        raise ValueError("exact-SCF replay differs from the production energy")
    if any(sha256(path) != before[path.name] for path in paths):
        raise ValueError("snapshot changed while being inspected; export a completed older step")
    report = {
        "format_version": 1, "status": "complete_observation_not_final_result",
        "source_job_id": source_job_id, "source_workdir": str(workdir),
        "snapshot_step": step, "optimizer_log_row": row,
        "export_script_sha256": sha256(Path(__file__)),
        "n_total_images": len(images), "n_active_images": len(images) - 2,
        "formula_units": 4, "pressure_GPa": 0., "climb": False,
        "ordered_symbols": images[0].get_chemical_symbols(),
        "k_eV_A2": .2, "fmax_target_eV_A": .10, "cell_scale_A": chain.cell_scale,
        "replayed_fmax_eV_A": fmax, "new_DFT_calls": 0,
        "extended_reaction_coordinate_A": chain.reaction_coordinate().tolist(),
        "raw_image_evaluations": evaluations, "path_diagnostics": diagnostics,
        "initial_geometry_gate": geometry,
        "limitations": ["fixed ordinary-production contract only; no physical-input changes",
                        "a converged NEB residual does not certify a stationary first-order saddle",
                        "discrete peak energy is not a sampling-converged barrier",
                        "live run/source/raw SCFs remain unmodified"],
    }
    output.mkdir(parents=True, exist_ok=False)
    for path in paths:
        shutil.copyfile(path, output / path.name)
    write(output / "evaluated_chain.traj", images)
    report["evaluated_chain_sha256"] = sha256(output / "evaluated_chain.traj")
    (output / "observation.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--step", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-job-id", required=True)
    args = parser.parse_args()
    report = export(args.workdir, args.step, args.output, args.source_job_id)
    print(json.dumps({key: report[key] for key in
                      ("status", "source_job_id", "snapshot_step", "replayed_fmax_eV_A", "new_DFT_calls")}))


if __name__ == "__main__":
    main()
