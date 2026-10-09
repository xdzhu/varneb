"""Six bounded G1 statics selected from a frozen ordinary-converged band.

No automatic MEP restart, CI, geometry relaxation or electronic parameter
override. Interpolate the existing ordered fractional lift and cell, not MIC.
Hermite estimates only select samples; they are not DFT barrier evidence.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read, write

from examples.hfo2_fixed_input_factory import CONTRACT, FixedHfo2Calculator, same_ordered_geometry
from scripts.analyze_hfo2_channel_network import read_evaluated_observation
from scripts.audit_hfo2_static_replica import audited_results, sha256
from vcneb.core import cell_force_from_arrays, deformation_from_cell, fractional_force_from_arrays
from vcneb import validate_path_geometry, validate_periodic_path_lift

COMMAND = "mpirun -np 32 /public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus"
SOURCE_OBSERVATION_SHA256 = "2058e4045e57ddfa23c95d551593aa861e6ccca169b9e98d3f12bf112bd5b229"
SEGMENTS = ((2, 3), (5, 6))
FRACTIONS = (.25, .50, .75)


def segment_derivatives(left, right, reference_cell):
    dq = right.get_scaled_positions(wrap=False) - left.get_scaled_positions(wrap=False)
    dF = (deformation_from_cell(right.cell.array, reference_cell)
          - deformation_from_cell(left.cell.array, reference_cell))
    result = []
    for atoms in (left, right):
        force_q = fractional_force_from_arrays(atoms.get_forces(), atoms.cell.array)
        force_F = cell_force_from_arrays(atoms.get_stress(voigt=False), atoms, reference_cell, pressure=0.)
        result.append(-float(np.sum(force_q * dq) + np.sum(force_F * dF)))
    return result  # eV/cell per dimensionless segment fraction


def hermite_candidates(e0, e1, d0, d1):
    coefficients = [2*e0-2*e1+d0+d1, -3*e0+3*e1-2*d0-d1, d0, e0]
    roots = np.roots([3*coefficients[0], 2*coefficients[1], coefficients[2]])
    return [{"fraction": float(z.real), "energy_meV_fu": float(np.polyval(coefficients, z.real))}
            for z in roots if abs(z.imag) < 1e-8 and 0 < z.real < 1]


def interpolate_segment(left, right, fraction):
    if not np.isfinite(fraction) or not 0 < fraction < 1:
        raise ValueError("a finite interior segment fraction required")
    if left.get_chemical_symbols() != right.get_chemical_symbols() or not np.array_equal(left.pbc, right.pbc):
        raise ValueError("ordered periodic endpoint identities differ")
    q = ((1-fraction)*left.get_scaled_positions(wrap=False)
         + fraction*right.get_scaled_positions(wrap=False))
    atoms = left.copy()
    atoms.calc = None
    atoms.set_cell((1-fraction)*left.cell.array + fraction*right.cell.array, scale_atoms=False)
    atoms.set_scaled_positions(q)
    return atoms


def prepare(observation, output):
    if output.exists():
        raise FileExistsError("refusing an existing six-point namespace")
    images, raw, digest = read_evaluated_observation(observation)
    if (digest != SOURCE_OBSERVATION_SHA256 or raw["snapshot_step"] != 69
            or raw["source_job_id"] != "28319570" or len(images) != 9
            or raw["replayed_fmax_eV_A"] > .10):
        raise ValueError("the frozen ordinary-preserving step69 source required")
    validate_periodic_path_lift(images)
    reference = images[0].cell.array
    energy0 = images[0].get_potential_energy()
    segment_records, candidates = [], []
    for left, right in SEGMENTS:
        a, b = images[left], images[right]
        energies = [(x.get_potential_energy()-energy0)*250 for x in (a, b)]
        derivatives = np.asarray(segment_derivatives(a, b, reference))*250
        if not (derivatives[0] > 0 and derivatives[1] < 0):
            raise ValueError("the registered positive-to-negative derivative bracket is absent")
        segment_records.append({"segment": [left, right], "energy_meV_fu": energies,
                                "dE_dlambda_meV_fu": derivatives.tolist(),
                                "Hermite_candidates_not_DFT": hermite_candidates(*energies, *derivatives)})
        for fraction in FRACTIONS:
            atoms = interpolate_segment(a, b, fraction)
            validate_periodic_path_lift([a, atoms, b])
            check = validate_path_geometry([a, atoms, b], minimum_distance=1.6, maximum_deformation=.25)
            candidates.append((atoms, {"index": len(candidates), "segment": [left, right],
                                        "fraction": fraction, "geometry_check": check}))
    output.mkdir(parents=True)
    geometry_root = output / "geometries"
    geometry_root.mkdir()
    points = []
    for atoms, record in candidates:
        path = geometry_root / f"POSCAR_{record['index']:02d}.vasp"
        write(path, atoms, format="vasp", direct=True, vasp5=True, sort=False)
        record.update({"geometry": path.relative_to(output).as_posix(), "geometry_sha256": sha256(path)})
        points.append(record)
    manifest = {"schema_version": 1, "purpose": "G1_preserving_step69_six_sampling_statics",
                "driver_sha256": sha256(Path(__file__)), "source_observation_sha256": digest,
                "source_job_id": raw["source_job_id"], "snapshot_step": 69,
                "physical_contract": CONTRACT, "source_static_directory": raw["raw_image_evaluations"][3]["raw_source"],
                "common_PO_energy_eV_cell": float(energy0), "formula_units": 4, "pressure_GPa": 0.,
                "sampled_band_maximum_meV_fu": max((x.get_potential_energy()-energy0)*250 for x in images),
                "segments": segment_records, "points": points, "maximum_new_SCF_calls": 6,
                "ordinary_fmax_eV_A": .10, "climb": False,
                "interpolation": "linear existing unwrapped fractional lift and cell; no MIC, atom remapping, alignment or relaxation",
                "limitations": "linear reconstruction statics only; neither an optimized continuous MEP nor stationary TS or certified barrier error",
                "new_DFT_calls_for_preparation": 0, "automatic_restart_or_G2_submission": False}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    return manifest


def load_manifest(root):
    m = json.loads((root / "manifest.json").read_text())
    if (m["purpose"] != "G1_preserving_step69_six_sampling_statics"
            or m["driver_sha256"] != sha256(Path(__file__)) or m["physical_contract"] != CONTRACT
            or len(m["points"]) != 6 or m["maximum_new_SCF_calls"] != 6
            or m["source_observation_sha256"] != SOURCE_OBSERVATION_SHA256
            or m["source_job_id"] != "28319570" or m["snapshot_step"] != 69
            or m["formula_units"] != 4 or m["pressure_GPa"] != 0
            or m["ordinary_fmax_eV_A"] != .10 or m["climb"]
            or abs(m["common_PO_energy_eV_cell"] + 9783.249675811956) > 1e-10):
        raise ValueError("sampling manifest/source/physical contract changed")
    for i, (point, expected) in enumerate(zip(m["points"],
                                             ((s, f) for s in SEGMENTS for f in FRACTIONS))):
        segment, fraction = expected
        if (point["index"] != i or point["segment"] != list(segment) or point["fraction"] != fraction
                or point["geometry"] != f"geometries/POSCAR_{i:02d}.vasp"):
            raise ValueError("declared sampling recipe changed")
    return m


def run_point(root, index):
    manifest = load_manifest(root)
    if isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < 6:
        raise ValueError("declared six-point index required")
    point = manifest["points"][index]
    path = root / point["geometry"]
    target = root / "calculations" / f"{index:02d}"
    if target.exists():
        raise FileExistsError("refusing repeat evaluation or existing calculation directory")
    if sha256(path) != point["geometry_sha256"]:
        raise ValueError("prepared geometry changed")
    atoms = read(path, format="vasp")
    calc = FixedHfo2Calculator(source=manifest["source_static_directory"], command=COMMAND, directory=target)
    atoms.calc = calc
    energy = atoms.get_potential_energy()
    directory = target / "scf_000000"
    raw = audited_results(directory)
    report = {"index": index, "segment": point["segment"], "fraction": point["fraction"],
              "manifest_sha256": sha256(root / "manifest.json"), "geometry_sha256": point["geometry_sha256"],
              "raw_log_sha256": sha256(directory / "OUT.ABACUS/running_scf.log"),
              "input_STRU_sha256": sha256(directory / "STRU"),
              "physical_input_sha256": {name: sha256(directory / name) for name in CONTRACT},
              "relative_energy_meV_fu": float((energy-manifest["common_PO_energy_eV_cell"])*250),
              "results": {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in raw.items()},
              "status": "audited_static_not_MEP_or_TS", "new_SCF_calls": 1}
    (target / "point_result.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    return report


def summarize(root):
    output = root / "summary.json"
    if output.exists():
        raise FileExistsError("refusing an existing summary")
    m = load_manifest(root)
    reports = []
    for i, point in enumerate(m["points"]):
        directory = root / "calculations" / f"{i:02d}"
        r = json.loads((directory / "point_result.json").read_text())
        rawdir = directory / "scf_000000"
        actual = audited_results(rawdir)
        if (r["index"] != i or r["segment"] != point["segment"] or r["fraction"] != point["fraction"]
                or r["manifest_sha256"] != sha256(root / "manifest.json")
                or r["physical_input_sha256"] != CONTRACT
                or any(sha256(rawdir / n) != h for n, h in CONTRACT.items())
                or sha256(rawdir / "OUT.ABACUS/running_scf.log") != r["raw_log_sha256"]
                or sha256(rawdir / "STRU") != r["input_STRU_sha256"]
                or sha256(root / point["geometry"]) != r["geometry_sha256"]
                or r["geometry_sha256"] != point["geometry_sha256"]
                or not same_ordered_geometry(read(root / point["geometry"], format="vasp"),
                                             read(rawdir / "STRU", format="abacus"))
                or abs(float(actual["energy"])-r["results"]["energy"]) > 1e-10
                or not np.allclose(actual["forces"], r["results"]["forces"], atol=1e-12, rtol=0)
                or not np.allclose(actual["stress"], r["results"]["stress"], atol=1e-12, rtol=0)
                or abs(r["relative_energy_meV_fu"] - (actual["energy"]-m["common_PO_energy_eV_cell"])*250) > 1e-10):
            raise ValueError("static result/raw provenance changed")
        reports.append(r)
    maximum = max(r["relative_energy_meV_fu"] for r in reports)
    summary = {"status": "six_sampling_statics_audited_not_relaxed_MEP", "points": reports,
               "highest_new_sample_meV_fu": maximum, "old_sampled_band_maximum_meV_fu": m["sampled_band_maximum_meV_fu"],
               "difference_meV_fu": maximum-m["sampled_band_maximum_meV_fu"], "new_SCF_calls": 6,
               "automatic_restart_or_G2_submission": False, "limitations": m["limitations"]}
    output.write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare", "run-point", "summarize"))
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--observation", type=Path)
    p.add_argument("--index", type=int)
    args = p.parse_args()
    if args.action == "prepare":
        if args.observation is None:
            p.error("prepare needs --observation")
        result = prepare(args.observation, args.root)
        print(json.dumps({"prepared_points": len(result["points"]), "new_DFT_calls": 0}))
    elif args.action == "run-point":
        result = run_point(args.root, args.index)
        print(json.dumps({k: result[k] for k in ("index", "relative_energy_meV_fu", "new_SCF_calls")}))
    else:
        result = summarize(args.root)
        print(json.dumps({k: result[k] for k in ("highest_new_sample_meV_fu", "difference_meV_fu", "new_SCF_calls")}))


if __name__ == "__main__":
    main()
