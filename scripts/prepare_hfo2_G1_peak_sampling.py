"""Five preregistered G1 peak samples; no band restart or G2 submission.

Reuse the frozen T/PO and PO/M ordered lifts, electronic bytes and calculator.
The Hermite screen selects points; it does not certify a DFT barrier or TS.
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
from scripts.prepare_hfo2_sampling_bridge import COMMAND, interpolate_segment, segment_derivatives, hermite_candidates
from vcneb import validate_path_geometry, validate_periodic_path_lift

BASE = Path("benchmarks/hfo2_channels/20261008")
SOURCES = {
    "T_to_PO": {"folder": "converged_gap/gap_converged_step06", "job": "28300425", "step": 6,
                "sha256": "b0d42fa6797986b038366f6856ba0943d2e2f4930356fd958ac638d1a2bef844", "images": 10},
    "PO_to_M": {"folder": "morning_update_20261009_0850/PO_M_step39", "job": "28298794", "step": 39,
                "sha256": "b1351c6cc87bbbcace3896808ea6a57c3bf66a4b913fbbac04da48a61a503526", "images": 9},
}
RECIPE = (("T_to_PO", 2, 3, .99), ("T_to_PO", 3, 4, .02),
          ("PO_to_M", 3, 4, .25), ("PO_to_M", 3, 4, .50), ("PO_to_M", 3, 4, .75))
PO_ENERGY = -9783.249675811956


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def prepare(repository, output):
    if output.exists():
        raise FileExistsError("refusing existing five-point namespace")
    bands, sources = {}, {}
    for name, spec in SOURCES.items():
        images, raw, digest = read_evaluated_observation(repository / BASE / spec["folder"])
        po = images[-1] if name == "T_to_PO" else images[0]
        if (digest != spec["sha256"] or raw["source_job_id"] != spec["job"]
                or raw["snapshot_step"] != spec["step"] or len(images) != spec["images"]
                or raw["replayed_fmax_eV_A"] > .10 or abs(po.get_potential_energy()-PO_ENERGY) > 1e-10):
            raise ValueError("frozen ordinary source/common PO reference changed")
        validate_periodic_path_lift(images)
        bands[name] = images
        sources[name] = {**spec, "source_static_directory": raw["raw_image_evaluations"][3]["raw_source"],
                         "sampled_maximum_relative_PO_meV_fu": max((a.get_potential_energy()-PO_ENERGY)*250 for a in images),
                         "source_fmax_eV_A": raw["replayed_fmax_eV_A"]}
    candidates, screens = [], {}
    for i, (name, left, right, fraction) in enumerate(RECIPE):
        images = bands[name]
        a, b = images[left], images[right]
        key = f"{name}:{left}->{right}"
        if key not in screens:
            energies = [(x.get_potential_energy()-PO_ENERGY)*250 for x in (a, b)]
            derivatives = np.asarray(segment_derivatives(a, b, images[0].cell.array))*250
            if not (derivatives[0] > 0 and derivatives[1] < 0):
                raise ValueError("registered derivative bracket absent")
            screens[key] = {"energies_relative_PO_meV_fu": energies, "dE_dlambda_meV_fu": derivatives.tolist(),
                            "Hermite_candidates_not_DFT": hermite_candidates(*energies, *derivatives)}
        atoms = interpolate_segment(a, b, fraction)
        validate_periodic_path_lift([a, atoms, b])
        check = validate_path_geometry([a, atoms, b], minimum_distance=1.6, maximum_deformation=.25)
        candidates.append((atoms, {"index": i, "channel": name, "segment": [left, right],
                                   "fraction": fraction, "geometry_check": check}))
    (output / "geometries").mkdir(parents=True)
    points = []
    for atoms, point in candidates:
        path = output / "geometries" / f"POSCAR_{point['index']:02d}.vasp"
        write(path, atoms, format="vasp", direct=True, vasp5=True, sort=False)
        points.append({**point, "geometry": path.relative_to(output).as_posix(), "geometry_sha256": sha256(path)})
    manifest = {"purpose": "G1_T_PO_M_five_peak_statics", "driver_sha256": sha256(Path(__file__)),
                "physical_contract": CONTRACT, "sources": sources, "screens_not_DFT": screens, "points": points,
                "common_PO_energy_eV_cell": PO_ENERGY, "formula_units": 4, "pressure_GPa": 0.,
                "ordinary_fmax_eV_A": .10, "climb": False, "maximum_new_SCF_calls": 5,
                "G1_peak_sampling_cap": 14, "prior_preserving_SCF_calls": 6,
                "remaining_reversing_sampling_cap": 3, "new_DFT_calls_for_preparation": 0,
                "interpolation": "existing unwrapped fractional lift and cell; no MIC/remapping/relaxation",
                "limitations": "linear reconstruction statics, not relaxed continuous MEP or stationary TS or certified error bound",
                "automatic_restart_or_G2_submission": False}
    save(output / "manifest.json", manifest)
    return manifest


def load_manifest(root):
    m = json.loads((root / "manifest.json").read_text())
    if (m["purpose"] != "G1_T_PO_M_five_peak_statics" or m["driver_sha256"] != sha256(Path(__file__))
            or m["physical_contract"] != CONTRACT or m["common_PO_energy_eV_cell"] != PO_ENERGY
            or m["formula_units"] != 4 or m["pressure_GPa"] != 0 or m["ordinary_fmax_eV_A"] != .10
            or m["climb"] or m["maximum_new_SCF_calls"] != 5 or len(m["points"]) != 5
            or m["G1_peak_sampling_cap"] != 14 or m["prior_preserving_SCF_calls"] != 6
            or m["remaining_reversing_sampling_cap"] != 3 or m["automatic_restart_or_G2_submission"]):
        raise ValueError("five-point physical contract changed")
    for name, spec in SOURCES.items():
        if any(m["sources"][name][k] != v for k, v in spec.items()):
            raise ValueError("frozen source changed")
    for i, (point, (name, left, right, fraction)) in enumerate(zip(m["points"], RECIPE)):
        if (point["index"] != i or point["channel"] != name or point["segment"] != [left, right]
                or point["fraction"] != fraction or point["geometry"] != f"geometries/POSCAR_{i:02d}.vasp"):
            raise ValueError("five-point recipe changed")
    return m


def run_point(root, index):
    m = load_manifest(root)
    if isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < 5:
        raise ValueError("registered five-point index required")
    point = m["points"][index]
    path, target = root / point["geometry"], root / "calculations" / f"{index:02d}"
    if target.exists():
        raise FileExistsError("refusing duplicate SCF namespace")
    if sha256(path) != point["geometry_sha256"]:
        raise ValueError("prepared geometry changed")
    atoms = read(path, format="vasp")
    atoms.calc = FixedHfo2Calculator(source=m["sources"][point["channel"]]["source_static_directory"],
                                    command=COMMAND, directory=target)
    energy = atoms.get_potential_energy()
    directory = target / "scf_000000"
    raw = audited_results(directory)
    r = {"index": index, "channel": point["channel"], "segment": point["segment"], "fraction": point["fraction"],
         "manifest_sha256": sha256(root / "manifest.json"), "geometry_sha256": point["geometry_sha256"],
         "raw_log_sha256": sha256(directory / "OUT.ABACUS/running_scf.log"),
         "input_STRU_sha256": sha256(directory / "STRU"),
         "physical_input_sha256": {n: sha256(directory / n) for n in CONTRACT},
         "relative_energy_meV_fu": float((energy-PO_ENERGY)*250),
         "results": {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in raw.items()},
         "new_SCF_calls": 1, "status": "audited_static_not_MEP_or_TS"}
    save(target / "point_result.json", r)
    return r


def summarize(root):
    if (root / "summary.json").exists():
        raise FileExistsError("refusing existing summary")
    m = load_manifest(root)
    reports = []
    for point in m["points"]:
        target = root / "calculations" / f"{point['index']:02d}"
        r = json.loads((target / "point_result.json").read_text())
        rawdir = target / "scf_000000"
        actual = audited_results(rawdir)
        if (any(r[k] != point[k] for k in ("index", "channel", "segment", "fraction", "geometry_sha256"))
                or r["manifest_sha256"] != sha256(root / "manifest.json")
                or r["physical_input_sha256"] != CONTRACT or any(sha256(rawdir/n) != h for n, h in CONTRACT.items())
                or r["raw_log_sha256"] != sha256(rawdir / "OUT.ABACUS/running_scf.log")
                or r["input_STRU_sha256"] != sha256(rawdir / "STRU")
                or sha256(root / point["geometry"]) != point["geometry_sha256"]
                or not same_ordered_geometry(read(root / point["geometry"], format="vasp"), read(rawdir / "STRU", format="abacus"))
                or any(not np.allclose(actual[k], r["results"][k], rtol=0, atol=1e-12) for k in ("energy", "forces", "stress"))
                or abs(r["relative_energy_meV_fu"]-(actual["energy"]-PO_ENERGY)*250) > 1e-10):
            raise ValueError("raw static provenance changed")
        reports.append(r)
    channels = {}
    for name, source in m["sources"].items():
        old = source["sampled_maximum_relative_PO_meV_fu"]
        highest = max(r["relative_energy_meV_fu"] for r in reports if r["channel"] == name)
        channels[name] = {"old_sampled_maximum_meV_fu": old, "highest_new_sample_meV_fu": highest,
                          "highest_union_sample_meV_fu": max(old, highest), "new_minus_old_meV_fu": highest-old}
    report = {"status": "five_peak_statics_audited_not_relaxed_MEP", "points": reports, "channels": channels,
              "new_SCF_calls": 5, "automatic_restart_or_G2_submission": False, "limitations": m["limitations"]}
    save(root / "summary.json", report)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare", "run-point", "summarize"))
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--repository", type=Path, default=Path("."))
    p.add_argument("--index", type=int)
    a = p.parse_args()
    result = (prepare(a.repository, a.root) if a.action == "prepare" else
              run_point(a.root, a.index) if a.action == "run-point" else summarize(a.root))
    print(json.dumps({k: v for k, v in result.items() if k in ("purpose", "status", "index", "relative_energy_meV_fu", "channels")}))


if __name__ == "__main__":
    main()
