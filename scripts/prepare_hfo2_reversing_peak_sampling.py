"""Two finite G1 reversing-peak statics; never submit/restart a band."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write

from examples.hfo2_fixed_input_factory import CONTRACT, FixedHfo2Calculator
from scripts.analyze_hfo2_channel_network import read_evaluated_observation
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.plot_hfo2_sampling_bridge import fixed_writer_geometry_matches
from scripts.prepare_hfo2_G1_peak_sampling import PO_ENERGY, save
from scripts.prepare_hfo2_sampling_bridge import COMMAND, hermite_candidates, interpolate_segment, segment_derivatives
from vcneb import VCNEB, validate_path_geometry, validate_periodic_path_lift

SOURCE = "benchmarks/hfo2_channels/20261008/terminal_G1_update_20261009_17/reversing_step45"
SOURCE_SHA = "23c1f7e3799ab0f1b388bd2791418c400f26eccf5679a3353c05dff326ef3f08"
RECIPE = ((3, 4, .99), (4, 5, .01))
LIMITATIONS = "two linear-reconstruction statics, not a continuous MEP, stationary TS, error bound or independent prediction"


def source_band(repository):
    images, raw, digest = read_evaluated_observation(repository / SOURCE)
    if (digest != SOURCE_SHA or raw["source_job_id"] != "28319571" or raw["snapshot_step"] != 45
            or len(images) != 9 or raw["replayed_fmax_eV_A"] > .10
            or any(abs(a.get_potential_energy()-PO_ENERGY) > 1e-10 for a in (images[0], images[-1]))):
        raise ValueError("frozen reversing ordinary source/common PO reference changed")
    validate_periodic_path_lift(images)
    return images, raw


def prepare(repository, root):
    if root.exists():
        raise FileExistsError("refusing existing reversing-peak namespace")
    images, raw = source_band(repository)
    candidates, screens = [], []
    for left, right, fraction in RECIPE:
        a, b = images[left], images[right]
        energies = [(x.get_potential_energy()-PO_ENERGY)*250 for x in (a, b)]
        derivatives = np.asarray(segment_derivatives(a, b, images[0].cell.array))*250
        if not (derivatives[0] > 0 > derivatives[1]):
            raise ValueError("registered physical derivative bracket absent")
        screens.append({"segment": [left, right], "energies_meV_fu": energies,
                        "dE_dlambda_meV_fu": derivatives.tolist(),
                        "Hermite_candidates_not_DFT": hermite_candidates(*energies, *derivatives)})
        atoms = interpolate_segment(a, b, fraction)
        validate_periodic_path_lift([a, atoms, b])
        check = validate_path_geometry([a, atoms, b], minimum_distance=1.6, maximum_deformation=.25)
        candidates.append((atoms, check))
    (root / "geometries").mkdir(parents=True)
    points = []
    for i, ((atoms, check), (left, right, fraction)) in enumerate(zip(candidates, RECIPE)):
        path = root / "geometries" / f"POSCAR_{i:02d}.vasp"
        write(path, atoms, format="vasp", direct=True, vasp5=True, sort=False)
        points.append({"index": i, "segment": [left, right], "fraction": fraction,
                       "geometry": path.relative_to(root).as_posix(), "geometry_sha256": sha256(path),
                       "geometry_check": check})
    manifest = {"purpose": "G1_reversing_two_peak_statics", "driver_sha256": sha256(Path(__file__)),
                "source": SOURCE, "source_sha256": SOURCE_SHA, "source_job": "28319571", "source_step": 45,
                "source_static_directory": raw["raw_image_evaluations"][4]["raw_source"],
                "source_fmax_eV_A": raw["replayed_fmax_eV_A"], "physical_contract": CONTRACT,
                "points": points, "screens_not_DFT": screens, "common_PO_energy_eV_cell": PO_ENERGY,
                "formula_units": 4, "pressure_GPa": 0., "ordinary_fmax_eV_A": .10, "climb": False,
                "maximum_new_SCF_calls": 2, "prior_peak_SCF_calls": 11, "G1_peak_sampling_cap": 14,
                "unused_cap_not_automatically_submitted": 1, "new_DFT_calls_for_preparation": 0,
                "interpolation": "existing unwrapped fractional lift and cell; no MIC/remapping/relaxation",
                "limitations": LIMITATIONS, "automatic_restart_or_G2_submission": False}
    save(root / "manifest.json", manifest)
    return manifest


def load_manifest(root):
    m = json.loads((root / "manifest.json").read_text())
    expected = {"purpose": "G1_reversing_two_peak_statics", "driver_sha256": sha256(Path(__file__)),
                "source": SOURCE, "source_sha256": SOURCE_SHA, "source_job": "28319571", "source_step": 45,
                "physical_contract": CONTRACT, "common_PO_energy_eV_cell": PO_ENERGY, "formula_units": 4,
                "pressure_GPa": 0., "ordinary_fmax_eV_A": .10, "climb": False, "maximum_new_SCF_calls": 2,
                "prior_peak_SCF_calls": 11, "G1_peak_sampling_cap": 14, "unused_cap_not_automatically_submitted": 1,
                "new_DFT_calls_for_preparation": 0, "limitations": LIMITATIONS,
                "automatic_restart_or_G2_submission": False}
    if any(m[k] != v for k, v in expected.items()) or len(m["points"]) != 2:
        raise ValueError("reversing-peak contract changed")
    for i, (p, (left, right, fraction)) in enumerate(zip(m["points"], RECIPE)):
        if (p["index"] != i or p["segment"] != [left, right] or p["fraction"] != fraction
                or p["geometry"] != f"geometries/POSCAR_{i:02d}.vasp"):
            raise ValueError("reversing-peak recipe changed")
    return m


def run_point(root, index):
    m = load_manifest(root)
    if isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < 2:
        raise ValueError("registered two-point index required")
    p = m["points"][index]
    target, geometry = root / "calculations" / f"{index:02d}", root / p["geometry"]
    if target.exists():
        raise FileExistsError("refusing duplicate SCF namespace")
    if sha256(geometry) != p["geometry_sha256"]:
        raise ValueError("prepared geometry changed")
    atoms = read(geometry, format="vasp")
    atoms.calc = FixedHfo2Calculator(source=m["source_static_directory"], command=COMMAND, directory=target)
    energy = atoms.get_potential_energy()
    rawdir = target / "scf_000000"
    actual = audited_results(rawdir)
    r = {**p, "manifest_sha256": sha256(root / "manifest.json"),
         "physical_input_sha256": {n: sha256(rawdir / n) for n in CONTRACT},
         "raw_log_sha256": sha256(rawdir / "OUT.ABACUS/running_scf.log"), "input_STRU_sha256": sha256(rawdir / "STRU"),
         "relative_energy_meV_fu": float((energy-PO_ENERGY)*250),
         "results": {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in actual.items()},
         "new_SCF_calls": 1, "status": "audited_static_not_MEP_or_TS"}
    save(target / "point_result.json", r)
    return r


def audit_points(root, m, *, require_all_bytes=False):
    probe = (root / "mpi_affinity.txt").read_text().splitlines()
    matches = [re.fullmatch(r"rank=(\d+) host=\S+", line) for line in probe]
    if any(x is None for x in matches) or sorted(int(x.group(1)) for x in matches) != list(range(32)):
        raise ValueError("unique actual 32-rank probe required")
    reports, atoms_list, elapsed = [], [], []
    for p in m["points"]:
        target = root / "calculations" / f"{p['index']:02d}"
        r = json.loads((target / "point_result.json").read_text())
        d = target / "scf_000000"
        actual = audited_results(d)
        call = json.loads((d / "call_audit.json").read_text())
        hashes = json.loads((d / "input_sha256.json").read_text())
        atoms = read(root / p["geometry"], format="vasp")
        byte_names = CONTRACT if require_all_bytes else ("INPUT", "KPT")
        if (any(r[k] != p[k] for k in p) or r["new_SCF_calls"] != 1
                or r["manifest_sha256"] != sha256(root / "manifest.json")
                or r["physical_input_sha256"] != CONTRACT or {n: hashes[n] for n in CONTRACT} != CONTRACT
                or any(sha256(d/n) != CONTRACT[n] for n in byte_names) or call["input_sha256"] != hashes
                or sha256(d / "STRU") != r["input_STRU_sha256"] or hashes["STRU"] != r["input_STRU_sha256"]
                or sha256(d / "OUT.ABACUS/running_scf.log") != r["raw_log_sha256"]
                or call["raw_log_sha256"] != r["raw_log_sha256"]
                or sha256(root / p["geometry"]) != p["geometry_sha256"]
                or not fixed_writer_geometry_matches(d / "STRU", atoms)
                or any(not np.allclose(actual[k], report["results"][k], rtol=0, atol=1e-12)
                       for k in ("energy", "forces", "stress") for report in (r, call))
                or abs(r["relative_energy_meV_fu"]-(actual["energy"]-PO_ENERGY)*250) > 1e-10
                or not np.isfinite(call["elapsed_seconds"]) or call["elapsed_seconds"] <= 0):
            raise ValueError("actual reversing sampling provenance changed")
        atoms.calc = SinglePointCalculator(atoms, **actual)
        reports.append(r)
        atoms_list.append(atoms)
        elapsed.append(call["elapsed_seconds"])
    return reports, atoms_list, elapsed


def summarize(root):
    if (root / "summary.json").exists():
        raise FileExistsError("refusing existing reversing summary")
    m = load_manifest(root)
    reports, _, elapsed = audit_points(root, m, require_all_bytes=True)
    r = {"status": "two_actual_reversing_statics_audited", "points": reports, "new_SCF_calls": 2,
         "SCF_seconds": elapsed, "physical_bytes_checked_on_HF": True,
         "automatic_restart_or_G2_submission": False, "limitations": LIMITATIONS}
    save(root / "summary.json", r)
    return r


def analyze(repository, root, output):
    if output.exists():
        raise FileExistsError("refusing existing reversing analysis")
    m = load_manifest(root)
    summary = json.loads((root / "summary.json").read_text())
    reports, new_atoms, elapsed = audit_points(root, m)
    if (summary["points"] != reports or summary["new_SCF_calls"] != 2
            or not summary["physical_bytes_checked_on_HF"] or summary["automatic_restart_or_G2_submission"]):
        raise ValueError("complete actual reversing summary required")
    original, _ = source_band(repository)
    inserted = dict(zip((left for left, _, _ in RECIPE), new_atoms))
    added = [a for i, old in enumerate(original) for a in ([old, inserted[i]] if i in inserted else [old])]
    validate_periodic_path_lift(added)
    band = VCNEB(added, pressure=0., k=.2, climb=False)
    fmax = float(np.linalg.norm(band.get_forces().reshape(-1, 3), axis=1).max())
    old = max((a.get_potential_energy()-PO_ENERGY)*250 for a in original)
    highest = max(r["relative_energy_meV_fu"] for r in reports)
    result = {"status": "two_actual_reversing_peak_checks_not_continuous_MEP_or_TS",
              "old_sampled_maximum_meV_fu": old, "highest_new_sample_meV_fu": highest,
              "highest_union_sample_meV_fu": max(old, highest), "new_minus_old_meV_fu": highest-old,
              "forward_sampled_barrier_meV_fu": max(old, highest), "reverse_sampled_barrier_meV_fu": max(old, highest),
              "endpoint_difference_meV_fu": 0., "inserted_total_images": len(added), "moving_images": len(added)-2,
              "replayed_ordinary_fmax_eV_A": fmax, "ordinary_residual_passed": fmax <= .10,
              "optimizer_steps": 0, "additional_SCF_calls": 0, "new_SCF_calls": 2, "total_G1_peak_SCF_calls": 13,
              "unique_MPI_ranks": 32, "SCF_seconds": elapsed, "total_SCF_seconds": sum(elapsed),
              "SCF_core_hours": sum(elapsed)*32/3600,
              "manifest_sha256": sha256(root / "manifest.json"), "summary_sha256": sha256(root / "summary.json"),
              "source_sha256": SOURCE_SHA, "analysis_driver_sha256": sha256(Path(__file__)),
              "new_DFT_calls_for_analysis": 0, "automatic_restart_or_G2_submission": False, "limitations": LIMITATIONS,
              "physical_bytes_checked_on_HF": "all six before SCF and summary; offline INPUT/KPT and raw records; repeated UPF/orbitals omitted"}
    save(output, result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare", "run-point", "summarize", "analyze"))
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--repository", type=Path, default=Path("."))
    p.add_argument("--index", type=int, choices=(0, 1))
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    result = (prepare(a.repository, a.root) if a.action == "prepare" else run_point(a.root, a.index)
              if a.action == "run-point" else summarize(a.root) if a.action == "summarize"
              else analyze(a.repository, a.root, a.output))
    print(json.dumps({k: v for k, v in result.items() if k in ("purpose", "status", "index", "relative_energy_meV_fu", "replayed_ordinary_fmax_eV_A", "highest_union_sample_meV_fu")}))


if __name__ == "__main__":
    main()
