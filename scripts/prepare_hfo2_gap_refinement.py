"""Insert three audited real SCFs into the historical guided chain, no DFT."""

import argparse
import json
from pathlib import Path

from ase.io import read, write

from examples.hfo2_fixed_input_factory import CONTRACT, same_ordered_geometry
from scripts.audit_hfo2_sparse_gap import audit_point
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.prepare_hfo2_reference_variants import write_clean_poscar
from vcneb import endpoint_structure_record, validate_path_geometry, minimum_image_path_lift, validate_periodic_path_lift


def prepare(source_root, chain, gap_root, output):
    if output.exists():
        raise FileExistsError("refusing existing refined-chain namespace")
    originals = [read(chain / f"POSCAR_{i:02d}", format="vasp") for i in range(7)]
    gaps = []
    for i in range(3):
        audit_point(gap_root, i)
        gaps.append(read(gap_root / "calculations" / f"{i:02d}" / "STRU", format="abacus"))
    images = originals[:2] + gaps + originals[2:]
    images, lift_audit = minimum_image_path_lift(images)
    sources = ([source_root / "00", source_root / "01"]
               + [gap_root / "calculations" / f"{i:02d}" for i in range(3)]
               + [source_root / f"{i:02d}" for i in range(2, 7)])
    audits = []
    for i, (atoms, source) in enumerate(zip(images, sources)):
        if not same_ordered_geometry(atoms, read(source / "STRU", format="abacus")):
            raise ValueError(f"seed/raw geometry mismatch at refined image{i}")
        if any(sha256(source / n) != h for n, h in CONTRACT.items()):
            raise ValueError("historical/gap electronic contract differs")
        raw = audited_results(source)
        audits.append({"index": i, "source": str(source), "energy_eV": float(raw["energy"]),
                       "structure": endpoint_structure_record(atoms),
                       "raw_log_sha256": sha256(source / "OUT.ABACUS/running_scf.log")})
    geometry = validate_path_geometry(images, minimum_distance=1.6, maximum_deformation=.25)
    output.mkdir(parents=True)
    write(output / "seed.traj", images)
    write_clean_poscar(output / "initial.vasp", images[0])
    write_clean_poscar(output / "final.vasp", images[-1])
    parameters = {"source_directory": str(source_root / "00"),
                  "seed_static_directories": [str(p) for p in sources]}
    (output / "factory_parameters.json").write_text(json.dumps(parameters, indent=2) + "\n", encoding="utf-8")
    report = {"purpose": "repair_undetected_internal_peak_sampling_before_claiming_low_barrier",
              "n_total_images": 10, "n_fixed_endpoints": 2, "n_active_images": 8,
              "old_images_retained": 7, "new_already_evaluated_images": 3,
              "first_optimization_step_cap": 10, "fmax_eV_A": .10, "climb": False,
              "optimizer": "FIRE", "maxstep_A": .02, "spring_eV_A2": .2,
              "pressure_GPa": 0, "physical_inputs_changed": False,
              "cached_seed_results": audits, "initial_geometry_gate": geometry,
              "periodic_lift_audit": lift_audit,
              "limitations": "initial sampled profile is not the optimized MEP; final interval resolution must be audited separately"}
    (output / "manifest.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    validate_periodic_path_lift(images)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source-root", "chain", "gap-root", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    report = prepare(args.source_root, args.chain, args.gap_root, args.output)
    print(json.dumps({k: report[k] for k in ("n_total_images", "n_active_images", "first_optimization_step_cap")}))


if __name__ == "__main__":
    main()
