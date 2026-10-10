"""Read-only, same-boundary audit of an evaluated Hf4O8 G2 chain.

Unlike the historical free-cell exporter, this replays the exact clamped
subspace and its registered coordinate scale. It never invokes a calculator.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write

from examples.hfo2_fixed_input_factory import (
    CONTRACT, read_fixed_hfo2_stru, same_ordered_geometry,
)
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.export_hfo2_chain_observation import optimizer_row
from vcneb import VCNEB, clamped_plane_vcneb_boundary
from vcneb import validate_path_geometry, validate_periodic_path_lift


CELL_SCALE = 5.12968067458423
PILOT_SCRIPT_SHA256 = "7bf764238c4cc857ca2b7e3e94642f0517e1f3b84f4e4c5b17833379ac1101ed"


def audited_source(image, directory, *, endpoint_record=None):
    """Require pinned raw data; never accept a nearby or permuted geometry."""
    if endpoint_record is not None:
        audit_path = directory / "seed_cache_audit.json"
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        if (audit.get("policy") != "identical_ordered_clamped_endpoint_hash_pinned"
                or audit.get("raw_source") != endpoint_record["directory"]
                or audit.get("input_sha256") != endpoint_record["input_sha256"]
                or audit.get("raw_log_sha256") != endpoint_record["raw_log_sha256"]):
            raise ValueError("fixed clamped endpoint cache provenance differs")
        candidates = [(Path(endpoint_record["directory"]), audit_path, audit)]
    else:
        candidates = []
        for path in sorted(directory.glob("scf_*/call_audit.json"), reverse=True):
            candidates.append((path.parent, path, json.loads(path.read_text(encoding="utf-8"))))
    for source, audit_path, audit in candidates:
        audit_hash = sha256(audit_path)
        hashes = audit.get("input_sha256", {})
        if (set(hashes) != {*CONTRACT, "STRU"}
                or any(hashes[n] != h for n, h in CONTRACT.items())
                or {n: sha256(source/n) for n in hashes} != hashes
                or sha256(source/"OUT.ABACUS/running_scf.log") != audit.get("raw_log_sha256")):
            raise ValueError("cached SCF physical inputs or pinned raw log changed")
        if not same_ordered_geometry(image, read_fixed_hfo2_stru(source/"STRU")):
            continue
        raw = audited_results(source)
        if "results" in audit and any(
                not np.allclose(raw[k], audit["results"][k], rtol=0., atol=1e-12)
                for k in ("energy", "forces", "stress")):
            raise ValueError("SCF audit and freshly parsed raw results differ")
        if sha256(audit_path) != audit_hash:
            raise ValueError("SCF audit changed during inspection")
        if ({n:sha256(source/n) for n in hashes} != hashes
                or sha256(source/"OUT.ABACUS/running_scf.log") != audit["raw_log_sha256"]):
            raise ValueError("SCF input/raw log changed during inspection")
        return source, raw, {"audit_path": str(audit_path), "audit_sha256": audit_hash,
                             "input_sha256": hashes, "raw_log_sha256": audit["raw_log_sha256"]}
    raise ValueError("no exact ordered, complete and hash-pinned SCF for this snapshot")


def load_boundary(workdir, images, production_script):
    """This bounded exporter accepts only the registered first G2 pilot."""
    if sha256(production_script) != PILOT_SCRIPT_SHA256:
        raise ValueError("unregistered production script: spring/driver contract not established")
    path = workdir/"vcneb_preflight.json"
    record = json.loads(path.read_text(encoding="utf-8"))
    expected = {"driver": "vcneb.material_runner", "n_images": 9,
                "external_pressure_gpa": 0., "cell_mode": "full", "cell_scale_A": CELL_SCALE,
                "fmax_target_eV_per_A": .10, "climbing_image_requested": False,
                "climb_after_steps": None, "align_cells": False,
                "cell_interpolation": "linear", "mode_subspace": None,
                "factory": "examples.hfo2_fixed_input_factory:make_clamped_seed_cached_factory",
                "calculator_validation": "runtime_instantiated"}
    if any(record.get(k) != v for k, v in expected.items()):
        raise ValueError("registered clamped G2 runtime contract differs")
    mechanical = record.get("mechanical_boundary") or {}
    if (mechanical.get("kind") != "clamped_plane" or mechanical.get("allow_tilt") is not True
            or mechanical.get("cell_dofs") != 3 or mechanical.get("fixed_cell_rows") != [0, 1]
            or mechanical.get("cell_scale_A") != CELL_SCALE):
        raise ValueError("clamped mechanical boundary or coordinate metric differs")
    reference = Path(mechanical["reference_file"])
    if sha256(reference) != mechanical["reference_file_sha256"]:
        raise ValueError("fixed substrate reference file changed")
    substrate = read(reference, format="vasp")
    if not np.allclose(substrate.cell.array, mechanical["reference_cell_A"], atol=1e-10, rtol=0):
        raise ValueError("fixed substrate reference matrix differs")
    boundary = clamped_plane_vcneb_boundary(12, substrate.cell.array, allow_tilt=True)
    boundary.validate_images(images)
    params = record["calculator_parameters"]
    if set(params) != {"source_directory", "seed_cache_records"}:
        raise ValueError("physical input overrides prohibited")
    caches = params["seed_cache_records"]
    if (len(caches) != 9 or caches[0] is None or caches[-1] is None
            or any(r is not None for r in caches[1:-1])):
        raise ValueError("two fixed clamped caches and seven active images required")
    return boundary, record, caches, {"runtime_preflight_sha256": sha256(path),
                                     "production_script_sha256": sha256(production_script)}


def export(workdir, step, output, source_job_id, *, production_script):
    if output.exists():
        raise FileExistsError("refusing existing observation namespace")
    if type(step) is not int or step < 0:
        raise ValueError("a nonnegative snapshot step required")
    paths = sorted((workdir/"snapshots"/f"step_{step:04d}").glob("POSCAR_*"))
    if [p.name for p in paths] != [f"POSCAR_{i:02d}" for i in range(9)]:
        raise ValueError("complete nine-image clamped snapshot required")
    before = {p.name: sha256(p) for p in paths}
    images = [read(p, format="vasp") for p in paths]
    if any(a.get_chemical_symbols() != ["Hf"]*4+["O"]*8 or not a.pbc.all() for a in images):
        raise ValueError("ordered periodic Hf4O8 images required")
    boundary, metadata, caches, provenance = load_boundary(workdir, images, production_script)
    row = optimizer_row(workdir/"vcneb.opt.log", step)
    validate_periodic_path_lift(images)
    geometry = validate_path_geometry(images, cell_scale=CELL_SCALE,
                                     minimum_distance=1.6, maximum_deformation=.25)
    evaluations = []
    for i, image in enumerate(images):
        source, raw, pinned = audited_source(image, workdir/f"image_{i:04d}", endpoint_record=caches[i])
        image.calc = SinglePointCalculator(image, **raw)
        evaluations.append({"image_index": i, "raw_source": str(source), **pinned,
                            "snapshot_POSCAR_sha256": before[paths[i].name],
                            "energy_eV_cell": float(raw["energy"]),
                            "forces_eV_A": np.asarray(raw["forces"]).tolist(),
                            "stress_ASE_voigt_eV_A3": np.asarray(raw["stress"]).tolist(),
                            "new_SCF_for_export": False, "ordered_periodic_geometry_matched": True})
    # material_runner passes MIC only to initial interpolation; run_vcneb
    # uses the supplied continuous lift with VCNEB's default mic=False.
    chain = VCNEB(images, cell_scale=CELL_SCALE, k=.2, climb=False, pressure=0.,
                  **boundary.vcneb_kwargs(images))
    forces = chain.get_forces()
    fmax = float(np.linalg.norm(forces, axis=1).max())
    if (abs(fmax-row["fmax_eV_A"]) > 7e-7
            or abs(float(chain.enthalpies.max())-row["highest_energy_eV_cell"]) > 7e-7):
        raise ValueError("same-boundary exact-SCF replay differs from production log")
    if (any(sha256(p) != before[p.name] for p in paths)
            or sha256(workdir/"vcneb_preflight.json") != provenance["runtime_preflight_sha256"]):
        raise ValueError("observation source changed during inspection")
    energies = chain.enthalpies
    relative = (energies-energies[0])*1000/4
    report = {"format_version": 1, "status": "complete_observation_not_final_result",
              "source_job_id": source_job_id, "source_workdir": str(workdir), "snapshot_step": step,
              "export_script_sha256": sha256(Path(__file__)), **provenance,
              "mechanical_boundary": metadata["mechanical_boundary"],
              "n_total_images": 9, "n_active_images": 7, "formula_units": 4,
              "pressure_GPa": 0., "climb": False, "k_eV_A2": .2, "cell_scale_A": CELL_SCALE,
              "fmax_target_eV_A": .10, "replayed_fmax_eV_A": fmax,
              "ordinary_residual_pass": bool(fmax <= .10), "new_DFT_calls": 0,
              "optimizer_log_row": row, "raw_image_evaluations": evaluations,
              "extended_reaction_coordinate_A": chain.reaction_coordinate().tolist(),
              "relative_enthalpy_meV_fu": relative.tolist(),
              "sampled_forward_barrier_meV_fu": float(relative.max()),
              "sampled_reverse_barrier_meV_fu": float((energies.max()-energies[-1])*1000/4),
              "path_diagnostics": chain.path_diagnostics(), "geometry_gate": geometry,
              "limitations": ["not a converged barrier when ordinary residual fails",
                              "discrete peaks do not certify sampling convergence or a stationary TS",
                              "clamped reaction stress is not a free-cell convergence criterion",
                              "raw SCFs, production source and run directories remain unmodified"]}
    output.mkdir(parents=True, exist_ok=False)
    for p in paths:
        shutil.copyfile(p, output/p.name)
    shutil.copyfile(workdir/"vcneb_preflight.json", output/"runtime_preflight.json")
    write(output/"evaluated_chain.traj", images)
    report["evaluated_chain_sha256"] = sha256(output/"evaluated_chain.traj")
    (output/"observation.json").write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--step", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-job-id", required=True)
    parser.add_argument("--production-script", type=Path, required=True)
    args = parser.parse_args()
    report = export(args.workdir, args.step, args.output, args.source_job_id,
                    production_script=args.production_script)
    print(json.dumps({k: report[k] for k in ("status", "replayed_fmax_eV_A", "ordinary_residual_pass")}))


if __name__ == "__main__":
    main()
