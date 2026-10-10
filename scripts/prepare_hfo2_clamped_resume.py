"""Prepare a geometry-only continuation of one audited G2 pilot (no DFT)."""
import argparse
import json
from pathlib import Path
import shutil

from ase.io import read

from scripts.audit_hfo2_static_replica import sha256
from scripts.export_hfo2_clamped_observation import CELL_SCALE, REGISTERED_FACTORIES
from vcneb import clamped_plane_vcneb_boundary, validate_periodic_path_lift


def prepare(observation, output):
    if output.exists():
        raise FileExistsError("fresh continuation namespace required")
    path = observation/"observation.json"
    r = json.loads(path.read_text(encoding="utf-8"))
    if (r["status"] != "complete_observation_not_final_result"
            or r["n_total_images"] != 9 or r["ordinary_residual_pass"]
            or r["pressure_GPa"] != 0 or r["climb"] or r["cell_scale_A"] != CELL_SCALE
            or r["production_script_sha256"] not in REGISTERED_FACTORIES
            or sha256(observation/"evaluated_chain.traj") != r["evaluated_chain_sha256"]):
        raise ValueError("an audited incomplete ordinary clamped G2 pilot required")
    preflight = observation/"runtime_preflight.json"
    if sha256(preflight) != r["runtime_preflight_sha256"]:
        raise ValueError("runtime preflight hash differs")
    metadata = json.loads(preflight.read_text(encoding="utf-8"))
    if metadata["mechanical_boundary"] != r["mechanical_boundary"]:
        raise ValueError("mechanical boundary differs from observed chain")
    substrate = Path(r["mechanical_boundary"]["reference_file"])
    if sha256(substrate) != r["mechanical_boundary"]["reference_file_sha256"]:
        raise ValueError("substrate changed")
    images = read(observation/"evaluated_chain.traj",index=":")
    if len(images)!=9 or any(a.get_chemical_symbols()!=["Hf"]*4+["O"]*8 for a in images):
        raise ValueError("nine ordered Hf4O8 frames required")
    clamped_plane_vcneb_boundary(12,read(substrate).cell.array,allow_tilt=True).validate_images(images)
    validate_periodic_path_lift(images)
    records = []
    for i,e in enumerate(r["raw_image_evaluations"]):
        if (e["image_index"] != i or sha256(observation/f"POSCAR_{i:02d}") != e["snapshot_POSCAR_sha256"]
                or sha256(Path(e["audit_path"])) != e["audit_sha256"]):
            raise ValueError("snapshot/raw audit provenance changed")
        records.append({"directory":e["raw_source"],"input_sha256":e["input_sha256"],
                        "raw_log_sha256":e["raw_log_sha256"]})
    if len(records)!=9:
        raise ValueError("all nine cache records required")
    parameters = {"source_directory":metadata["calculator_parameters"]["source_directory"],
                  "seed_cache_records":records}
    output.mkdir(parents=True,exist_ok=False)
    for source,target in (("POSCAR_00","initial.vasp"),("POSCAR_08","final.vasp"),
                          ("evaluated_chain.traj","seed.traj")):
        shutil.copyfile(observation/source,output/target)
    shutil.copyfile(substrate,output/"substrate.vasp")
    (output/"factory_parameters.json").write_text(json.dumps(parameters,indent=2)+"\n",encoding="utf-8")
    manifest = {"status":"audited_G2_geometry_continuation_prepared", "source_job_id":r["source_job_id"],
                "source_observation_sha256":sha256(path), "source_step":r["snapshot_step"],
                "cell_scale_A":CELL_SCALE,"pressure_GPa":0,"climb":False,"ordinary_fmax_eV_A":.10,
                "n_total_images":9,"n_fixed_endpoints":2,"n_active_images":7,"current_frame_exact_caches":9,
                "new_independent_chains":0,"new_DFT_calls":0,"physical_inputs_changed":False,
                "holdout_generated_or_read":False,"optimizer_state_restored":False,
                "optimizer_policy":"new FIRE state at identical audited step; no velocities restored",
                "files_sha256":{n:sha256(output/n) for n in
                    ("initial.vasp","final.vasp","seed.traj","substrate.vasp","factory_parameters.json")}}
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    return manifest


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observation",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(prepare(args.observation,args.output)))


if __name__=="__main__":
    main()
