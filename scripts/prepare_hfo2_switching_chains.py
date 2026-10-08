"""Stage two ordered switching candidates only after electronic Berry gate.

No DFT or job submission. Keep the explicit parent-defined atom identity;
do not auto-map distorted endpoints into the same permutation. These are
unconstrained linear starting bands, not a proof of two distinct final MEPs.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ase.io import read, write
import numpy as np

from examples.hfo2_fixed_input_factory import CONTRACT, same_ordered_geometry
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.hfo2_endpoint_polarization import LABELS
from scripts.prepare_hfo2_reference_variants import write_clean_poscar
from vcneb import endpoint_structure_record, interpolate_vcneb, validate_path_geometry


def ordered_seed(initial, final):
    if (initial.get_chemical_symbols() != ["Hf"]*4 + ["O"]*8
            or final.get_chemical_symbols() != initial.get_chemical_symbols()
            or not np.allclose(initial.cell.array, final.cell.array, atol=1e-10, rtol=0)):
        raise ValueError("same-cell ordered Hf4O8 inversion endpoints required")
    images = interpolate_vcneb(initial, final, 9, align_cells=False, align_translation=False,
                               mic=True, mapping=list(range(12)), cell_interpolation="linear",
                               minimum_distance=1.6, maximum_deformation=.25)
    # Preserve exact endpoint representation as well as its periodic geometry.
    images[0], images[-1] = initial.copy(), final.copy()
    for image in images:
        image.calc = None
    validate_path_geometry(images, minimum_distance=1.6, maximum_deformation=.25)
    return images


def prepare(variants, polar_root, output):
    if output.exists():
        raise FileExistsError("refusing existing switching-seed namespace")
    summary_path, manifest_path = polar_root / "summary.json", polar_root / "manifest.json"
    summary = json.loads(summary_path.read_text())
    if summary.get("status") != "passed" or summary.get("manifest_sha256") != sha256(manifest_path):
        raise ValueError("actual electronic inversion gate not passed or manifest changed")
    manifest = json.loads(manifest_path.read_text())
    if tuple(p["label"] for p in manifest["points"]) != LABELS:
        raise ValueError("Berry endpoints are not the declared three variants")
    atoms, raw_sources, endpoint_audits = [], [], []
    files = ("PO.vasp", "PO_minus_T_preserving_inversion.vasp", "PO_minus_T_reversing_inversion.vasp")
    for point, filename in zip(manifest["points"], files):
        source = Path(point["source_directory"])
        if (any(sha256(source / n) != h for n, h in CONTRACT.items())
                or any(sha256(source / n) != h for n, h in point["baseline_input_sha256"].items())
                or sha256(source / "OUT.ABACUS/running_scf.log") != point["baseline_log_sha256"]):
            raise ValueError("original endpoint physical files or baseline evidence changed")
        geometry = read(variants / filename)
        if not same_ordered_geometry(geometry, read(source / "STRU", format="abacus")):
            raise ValueError("ordered geometric variant and actual static source differ")
        raw = audited_results(source)
        if (np.max(np.linalg.norm(raw["forces"], axis=1)) > .03
                or np.max(np.abs(raw["stress"]))*1602.176634 > 2):
            raise ValueError("original endpoint force/stress contract fails")
        endpoint_audits.append({"label": point["label"], "source": str(source),
                                "energy_eV_cell": float(raw["energy"]),
                                "structure": endpoint_structure_record(geometry),
                                "raw_log_sha256": sha256(source / "OUT.ABACUS/running_scf.log")})
        atoms.append(geometry)
        raw_sources.append(source)
    seeds = [ordered_seed(atoms[0], a) for a in atoms[1:]]
    if all(same_ordered_geometry(a, b) for a, b in zip(*seeds)):
        raise ValueError("two purported ordered candidate chains coincide")
    output.mkdir(parents=True)
    reports = []
    for index, images in enumerate(seeds, start=1):
        label = LABELS[index]
        folder = output / label
        folder.mkdir()
        write(folder / "seed.traj", images)
        write_clean_poscar(folder / "initial.vasp", images[0])
        write_clean_poscar(folder / "final.vasp", images[-1])
        (folder / "factory_parameters.json").write_text(json.dumps({
            "source_directory": str(raw_sources[0]),
            "seed_static_directories": [str(raw_sources[0])] + [None]*7 + [str(raw_sources[index])]}, indent=2) + "\n")
        report = {"purpose": "common_POplus_to_ordered_inversion_variant_switching_candidate",
                  "variant": label, "n_total_images": 9, "n_fixed_endpoints": 2, "n_active_images": 7,
                  "fmax_eV_A": .10, "climb": False, "pressure_GPa": 0,
                  "optimizer": "FIRE", "maxstep_A": .02, "spring_eV_A2": .2,
                  "first_optimization_step_cap": 10, "physical_inputs_changed": False,
                  "endpoint_audits": [endpoint_audits[0], endpoint_audits[index]],
                  "polarization_summary_sha256": sha256(summary_path),
                  "polarization_manifest_sha256": sha256(manifest_path),
                  "initial_geometry_gate": validate_path_geometry(images, minimum_distance=1.6, maximum_deformation=.25),
                  "atom_mapping": "identity after explicit parent operation; no per-image/endpoint auto-remapping",
                  "cell_boundary": "initial cells equal; subsequent joint-cell optimization remains full at P=0",
                  "limitations": "distinct seeds/orderings only, not topological inequivalence or optimized barriers; continuous polarization branch still pending"}
        (folder / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")
        reports.append(report)
    return reports


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("variants", "polar-root", "output"):
        parser.add_argument(f"--{key}", type=Path, required=True)
    args = parser.parse_args()
    reports = prepare(args.variants, args.polar_root, args.output)
    print(json.dumps({"prepared_candidates": len(reports), "n_active_images_each": 7, "DFT_submitted": False}))


if __name__ == "__main__":
    main()
