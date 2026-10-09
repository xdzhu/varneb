"""E044: fourteen existing switching interiors, output-only SCF + native Berry.

Never alter the stationary-endpoint adapter, NEB factory, physical contract,
or source runs. The observable quadratures are not energies for barriers.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

import numpy as np

from examples.hfo2_fixed_input_factory import CONTRACT, same_ordered_geometry
from scripts import hfo2_endpoint_polarization as native
from scripts.analyze_hfo2_channel_network import read_evaluated_observation
from scripts.audit_hfo2_static_replica import INPUT_FILES, audited_results, sha256
from vcneb.polarization import modular_difference, quantum_lattice, sampled_reduced_branch

CASE = Path("benchmarks/hfo2_channels/20261008")
CHANNELS = (
    ("preserving", "switching_converged_update_20261009_1255/preserving_step69",
     "2058e4045e57ddfa23c95d551593aa861e6ccca169b9e98d3f12bf112bd5b229", "28319570", 69),
    ("reversing", "terminal_G1_update_20261009_17/reversing_step45",
     "23c1f7e3799ab0f1b388bd2791418c400f26eccf5679a3353c05dff326ef3f08", "28319571", 45),
)
TOLERANCES = {"SCF_energy_eV_cell": 1e-5, "SCF_force_eV_A": 1e-4,
              "SCF_stress_kbar": .02, "sampled_gap_min_eV": .1,
              "longitudinal_P_convergence_C_m2": .01}
AXIS_TRANSVERSE_MAX = 1e-3
LIMITATIONS = ("R3 component and sampled longitudinal quadrature only; conditional sampled "
               "branch lift, not full-vector polarization, spontaneous P, continuous-path "
               "winding/insulation, TS, stability or complete G1 certification")
# Accepted production driver before the post-processing-only numpy-bool fix.
# Historical physics, point recipe, native auditor and all other hashes still
# have to match. This is not permission for arbitrary old/mutated drivers.
HISTORICAL_DRIVER_SHA256 = "97ad87eab742adabf3ae96297efffdf619d1420966271d7481978df00ed83b82"
NATIVE_ADAPTER_LF_SHA256 = "39ae7c33b604a021a2722a15144ddf25a3fb47e4c37f0920d29e5bfb37699187"
NATIVE_ADAPTER_HISTORICAL_CRLF_SHA256 = "7147021eee5421cc283a15938c0d0f75aee236c488910f07bcaccc867ebcda53"


def save(path, value):
    def native_scalar(item):
        if isinstance(item, np.generic):
            return item.item()
        raise TypeError(f"not a JSON-native scalar: {type(item).__name__}")
    path.write_text(json.dumps(value, indent=2, default=native_scalar, allow_nan=False) + "\n", encoding="utf-8")


def blueprint(repository):
    root = repository / CASE
    endpoints = root / "polarization_summary.json"
    evidence = json.loads(endpoints.read_text())
    if evidence["status"] != "passed" or evidence["switching_path_branch_selected"]:
        raise ValueError("original modular endpoint audit required, not a selected path branch")
    channels, points = [], []
    for name, relative, digest, job, step in CHANNELS:
        images, raw, actual = read_evaluated_observation(root / relative)
        if (actual != digest or raw["source_job_id"] != job or raw["snapshot_step"] != step
                or len(images) != 9 or raw["replayed_fmax_eV_A"] > .10):
            raise ValueError("frozen ordinary switching source changed")
        cells = np.asarray([a.cell.array for a in images])
        transverse = np.linalg.norm(cells[:, 2, :2], axis=1) / np.linalg.norm(cells[:, 2], axis=1)
        if np.any(cells[:, 2, 2] <= 0) or np.max(transverse) > AXIS_TRANSVERSE_MAX:
            raise ValueError("R3 no longer tracks the declared common positive z direction")
        channels.append({"name": name, "observation": relative, "observation_sha256": digest,
                         "job": job, "step": step, "n_total_images": 9,
                         "ordinary_fmax_eV_A": raw["replayed_fmax_eV_A"],
                         "max_R3_transverse_fraction": float(transverse.max()),
                         "endpoint_labels": ["PO_plus", f"PO_minus_T_{name}"]})
        for i in range(1, 8):
            point = raw["raw_image_evaluations"][i]
            points.append({"index": len(points), "label": f"{name}_{i:02d}", "channel": name,
                           "image_index": i, "source_directory": point["raw_source"],
                           "baseline_log_sha256": point["raw_log_sha256"],
                           "baseline_results": {"energy": point["energy_eV_cell"],
                               "forces": point["forces_eV_A"], "stress": point["stress_ASE_voigt_eV_A3"]},
                           "baseline_input_sha256": point["input_sha256"],
                           "cell_A": images[i].cell.array.tolist(),
                           "positions_A": images[i].positions.tolist(),
                           "quantum_lattice_C_m2": quantum_lattice(images[i].cell.array).tolist()})
    return {"purpose": "E044_fourteen_existing_switching_interiors_native_Berry",
            "physical_contract": CONTRACT, "channels": channels, "points": points,
            "endpoint_summary_sha256": sha256(endpoints),
            "maximum_new_SCF_calls": 14, "maximum_new_NSCF_calls": 42,
            "driver_sha256": sha256(Path(__file__)),
            "native_adapter_sha256": sha256(Path(native.__file__)),
            "polarization_module_sha256": sha256(repository / "vcneb/polarization.py"),
            "new_DFT_calls_for_preparation": 0, "automatic_restart_or_G2_submission": False,
            "NSCF_energies_used_for_barriers": False, "SCF_physical_settings_changed": False,
            "SCF_delta": {"out_chg": "1", "out_bandgap": "1"},
            "NSCF_meshes": [[2, 2, n] for n in native.GRIDS], "gdir": 3,
            "abacus_commit": "f7cb1d3", "tolerances": TOLERANCES,
            "axis_transverse_max": AXIS_TRANSVERSE_MAX, "limitations": LIMITATIONS}


def replay_endpoints(repository, root):
    """Fresh read-only HF replay of all cached raw endpoint property evidence."""
    original = repository / CASE
    summary = json.loads((original / "polarization_summary.json").read_text())
    replay = root / "endpoint_replay"
    replay.mkdir()
    shutil.copyfile(original / "polarization_manifest_r10.json", replay / "manifest.json")
    for i, label in enumerate(native.LABELS):
        expected_scf = summary["SCF_output_only_audits"][label]
        scf = native.audit_scf(replay, i, Path(expected_scf["directory"]))
        if scf != {k: v for k, v in expected_scf.items() if k not in ("directory", "reused_completed_SCF")}:
            raise ValueError("cached endpoint SCF evidence changed")
        for expected in summary["native_Berry_audits"][label]:
            actual = native.audit_nscf(replay, i, expected["nz"], Path(expected["directory"]), scf)
            if actual != {k: v for k, v in expected.items() if k not in ("directory", "reused_completed_NSCF")}:
                raise ValueError("cached endpoint NSCF evidence changed")
    shutil.copyfile(original / "polarization_summary.json", root / "endpoint_summary.json")
    return {"status": "all_three_cached_endpoints_freshly_replayed", "new_DFT_calls": 0,
            "endpoint_summary_sha256": sha256(root / "endpoint_summary.json"),
            "cached_endpoint_manifest_sha256": sha256(replay / "manifest.json")}


def prepare(repository, root):
    if root.exists():
        raise FileExistsError("refusing existing path-property namespace")
    plan = blueprint(repository)
    # Validate every source before generating any executable staged input.
    for point in plan["points"]:
        source = Path(point["source_directory"])
        native.check_inputs(source, point["baseline_input_sha256"])
        if (any(point["baseline_input_sha256"][n] != h for n, h in CONTRACT.items())
                or sha256(source / "OUT.ABACUS/running_scf.log") != point["baseline_log_sha256"]):
            raise ValueError("original physical bytes or SCF log changed")
        atoms = native.read(source / "STRU", format="abacus")
        from ase import Atoms
        frozen = Atoms(["Hf"]*4 + ["O"]*8, positions=point["positions_A"], cell=point["cell_A"], pbc=True)
        if not same_ordered_geometry(atoms, frozen):
            raise ValueError("original ordered periodic geometry changed")
        results = audited_results(source)
        if any(not np.allclose(results[k], v, rtol=0, atol=1e-10)
               for k, v in point["baseline_results"].items()):
            raise ValueError("original E/F/stress differs from frozen observation")
    root.mkdir(parents=True)
    plan["endpoint_replay"] = replay_endpoints(repository, root)
    for point in plan["points"]:
        source = Path(point["source_directory"])
        base = (source / "INPUT").read_bytes()
        if {"out_chg", "out_bandgap"} & native.input_values(base.decode()).keys():
            raise ValueError("unexpected output flags in original source")
        point["nscf_inputs"] = []
        for stage in ("scf", *(f"nscf_22{n}" for n in native.GRIDS)):
            folder = root / "points" / point["label"] / stage
            folder.mkdir(parents=True)
            for name in INPUT_FILES:
                shutil.copyfile(source / name, folder / name)
            if stage == "scf":
                (folder / "INPUT").write_bytes(base + native.OUTPUT_ADDITION)
                point["scf_input_sha256"] = {n: sha256(folder / n) for n in INPUT_FILES}
            else:
                nz = int(stage[-1])
                (folder / "INPUT").write_bytes(native.berry_input(base))
                (folder / "KPT").write_bytes(f"K_POINTS\n0\nGamma\n2 2 {nz} 0 0 0\n".encode())
                point["nscf_inputs"].append({"nz": nz, "input_sha256": {n: sha256(folder / n) for n in INPUT_FILES}})
        point.update(occupied_bands=48, cached_scf_directory=None, cached_nscf_directories={})
    save(root / "manifest.json", plan)
    load_manifest(repository, root)
    return plan


def validate_blueprint_metadata(manifest, expected):
    for key, value in expected.items():
        if key == "driver_sha256":
            if manifest[key] not in (value, HISTORICAL_DRIVER_SHA256):
                raise ValueError("path-property contract changed: driver_sha256")
        elif key == "native_adapter_sha256":
            historical_line_endings = (value == NATIVE_ADAPTER_LF_SHA256
                                       and manifest[key] == NATIVE_ADAPTER_HISTORICAL_CRLF_SHA256)
            if manifest[key] != value and not historical_line_endings:
                raise ValueError("path-property contract changed: native_adapter_sha256")
        elif key != "points" and manifest[key] != value:
            raise ValueError(f"path-property contract changed: {key}")


def load_manifest(repository, root):
    manifest = json.loads((root / "manifest.json").read_text())
    expected = blueprint(repository)
    validate_blueprint_metadata(manifest, expected)
    if len(manifest["points"]) != 14 or sha256(root / "endpoint_summary.json") != expected["endpoint_summary_sha256"]:
        raise ValueError("fourteen-point recipe or cached endpoints changed")
    for point, target in zip(manifest["points"], expected["points"]):
        if any(point[k] != v for k, v in target.items()) or point["occupied_bands"] != 48 or point["cached_scf_directory"] is not None or point["cached_nscf_directories"]:
            raise ValueError("path-point/source/charge-cache recipe changed")
        # Independently define, rather than trust hashes of, allowed deltas.
        original_input = Path(point["source_directory"]) / "INPUT"
        if original_input.exists():
            base = original_input.read_bytes()
        else:
            staged = root / "points" / point["label"] / "scf/INPUT"
            body = staged.read_bytes()
            if not body.endswith(native.OUTPUT_ADDITION):
                raise ValueError("output-only SCF delta changed")
            base = body[:-len(native.OUTPUT_ADDITION)]
        import hashlib
        if hashlib.sha256(base).hexdigest() != CONTRACT["INPUT"]:
            raise ValueError("baseline INPUT changed")
        stages = [("scf", point["scf_input_sha256"], base + native.OUTPUT_ADDITION)]
        if [x["nz"] for x in point["nscf_inputs"]] != list(native.GRIDS):
            raise ValueError("fixed longitudinal quadrature recipe changed")
        stages += [(f"nscf_22{x['nz']}", x["input_sha256"], native.berry_input(base)) for x in point["nscf_inputs"]]
        for stage, hashes, body in stages:
            folder = root / "points" / point["label"] / stage
            native.check_inputs(folder, hashes)
            if (folder / "INPUT").read_bytes() != body or hashes["STRU"] != point["baseline_input_sha256"]["STRU"]:
                raise ValueError("output recipe or geometry changed")
            if any(hashes[n] != h for n, h in CONTRACT.items() if n not in ("INPUT", "KPT")):
                raise ValueError("basis/pseudopotential changed")
            mesh = 2 if stage == "scf" else int(stage[-1])
            if (folder / "KPT").read_bytes() != f"K_POINTS\n0\nGamma\n2 2 {mesh} 0 0 0\n".encode():
                raise ValueError("SCF/observable mesh contract changed")
    return manifest


def run(repository, root, index):
    load_manifest(repository, root)
    if isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < 14:
        raise ValueError("registered fourteen-point index required")
    native.run_endpoint(root, index)  # Generic executor; endpoint preparation stays strict.


def summarize(repository, root):
    manifest = load_manifest(repository, root)
    endpoints = json.loads((root / "endpoint_summary.json").read_text())
    point_audits = {}
    for i, point in enumerate(manifest["points"]):
        destination = root / "calculations" / point["label"]
        scf = native.audit_scf(root, i, destination / "scf")
        expected_scf = json.loads((destination / "scf_audit.json").read_text())
        if scf != {k: v for k, v in expected_scf.items() if k not in ("directory", "reused_completed_SCF")}:
            raise ValueError("completed SCF evidence changed")
        audits = [native.audit_nscf(root, i, n, destination / f"nscf_22{n}", scf) for n in native.GRIDS]
        for actual in audits:
            stored = json.loads((destination / f"berry_22{actual['nz']}.json").read_text())
            if actual != {k: v for k, v in stored.items() if k not in ("directory", "reused_completed_NSCF")}:
                raise ValueError("completed Berry evidence changed")
        point_audits[point["label"]] = {"SCF": scf, "NSCF": audits}
    channels = []
    for channel in manifest["channels"]:
        name = channel["name"]
        audits = [endpoints["native_Berry_audits"][channel["endpoint_labels"][0]]]
        audits += [point_audits[f"{name}_{i:02d}"]["NSCF"] for i in range(1, 8)]
        audits += [endpoints["native_Berry_audits"][channel["endpoint_labels"][1]]]
        changes = [abs(modular_difference(a[-1]["value_modern_SI_C_m2"], a[-2]["value_modern_SI_C_m2"],
                                        a[-1]["reported_modulus_modern_SI_C_m2"])) for a in audits]
        converged = bool(max(changes) <= TOLERANCES["longitudinal_P_convergence_C_m2"])
        lift = sampled_reduced_branch([a[-1]["value_modern_SI_C_m2"] for a in audits],
            [a[-1]["physical_quantum_C_m2"] for a in audits],
            [a[-1]["reported_modulus_modern_SI_C_m2"] for a in audits], changes)
        channels.append({"name": name, "native_Berry_audits_in_source_order": audits,
                         "longitudinal_224_to_228_modular_change_C_m2": changes,
                         "longitudinal_gate_passed": converged, "sampled_branch": lift})
    result = {"status": "bounded_path_property_measurements_complete",
              "manifest_sha256": sha256(root / "manifest.json"), "points": point_audits,
              "channels": channels, "new_SCF_calls": 14, "new_NSCF_calls": 42,
              "new_endpoints_peak_samples_or_channels": 0,
              "all_longitudinal_gates_passed": all(c["longitudinal_gate_passed"] for c in channels),
              "all_sampled_lifts_unique": all(c["sampled_branch"]["status"] == "conditional_sampled_lift" for c in channels),
              "automatic_restart_or_G2_submission": False, "full_G1_certified": False,
              "limitations": LIMITATIONS, "NSCF_energies_used_for_barriers": False}
    save(root / "summary.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("blueprint", "prepare", "run", "summarize"))
    parser.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int)
    args = parser.parse_args()
    if args.action == "blueprint":
        if args.root.exists():
            raise FileExistsError("refusing existing blueprint output")
        save(args.root, blueprint(args.repository))
    elif args.action == "prepare":
        print(json.dumps({"staged_existing_interiors": len(prepare(args.repository, args.root)["points"])}))
    elif args.action == "run":
        run(args.repository, args.root, args.index)
    else:
        result = summarize(args.repository, args.root)
        print(json.dumps({k: result[k] for k in ("status", "all_longitudinal_gates_passed", "all_sampled_lifts_unique", "full_G1_certified")}))


if __name__ == "__main__":
    main()
