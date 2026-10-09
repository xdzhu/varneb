"""Offline replay of E044's raw exports; no charge cubes/UPFs redistributed.

Full charge and basis-byte checks are the recorded HF audits, not claimed
from their omitted local files. Native raw logs and both eigenvalue formats
are re-parsed locally, and output-only SCF reproducibility is recomputed.
"""
from __future__ import annotations

import argparse
import json
import hashlib
from pathlib import Path
import re
import shutil

import numpy as np

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts import hfo2_endpoint_polarization as native
from scripts import hfo2_switching_path_polarization as flow
from scripts.audit_hfo2_static_replica import audited_results, sha256
from vcneb.polarization import (modular_difference, parse_abacus_berry, quantum_lattice,
                               sampled_band_gap, sampled_nscf_band_gap, sampled_reduced_branch)


def native_elapsed_seconds(body):
    rows = re.findall(r"Total\s+Time\s*:\s*(\d+)\s+h\s+(\d+)\s+mins\s+(\d+)\s+secs", body)
    if len(rows) != 1:
        raise ValueError("one complete native h/min/s runtime required")
    h, m, s = map(int, rows[0])
    if m >= 60 or s >= 60:
        raise ValueError("invalid native time fields")
    return h*3600 + m*60 + s


def export(source, output):
    """Copy only reviewed raw observables, not charge/basis caches or sources."""
    if output.exists():
        raise FileExistsError("refusing existing exported property namespace")
    output.mkdir(parents=True)
    for name in ("manifest.json", "summary.json", "endpoint_summary.json"):
        shutil.copyfile(source/name, output/name)
    m = json.loads((source/"manifest.json").read_text())
    for point in m["points"]:
        for stage in ("scf", *(f"nscf_22{n}" for n in native.GRIDS)):
            old = source/"calculations"/point["label"]/stage
            new = output/"calculations"/point["label"]/stage
            (new/"OUT.ABACUS").mkdir(parents=True)
            log = "running_scf.log" if stage == "scf" else "running_nscf.log"
            bands = "istate.info" if stage == "scf" else "BANDS_1.dat"
            for name in ("INPUT", "KPT", "STRU", "abacus.out", "abacus.err",
                         "OUT.ABACUS/"+log, "OUT.ABACUS/"+bands):
                shutil.copyfile(old/name, new/name)
    files = {p.relative_to(output).as_posix(): sha256(p) for p in sorted(output.rglob("*")) if p.is_file()}
    flow.save(output/"file_inventory.json", {"source_directory": str(source), "file_sha256": files,
        "omitted": "SCF charge cubes, UPFs/orbitals and repeated staged templates; verified in actual HF native audits",
        "new_DFT_calls": 0})
    return files


def audit(repository, root):
    inventory = json.loads((root/"file_inventory.json").read_text())
    for relative, expected_hash in inventory["file_sha256"].items():
        path = (root/relative).resolve()
        if root.resolve() not in path.parents or sha256(path) != expected_hash:
            raise ValueError("raw export inventory/path changed")
    m = json.loads((root/"manifest.json").read_text())
    summary = json.loads((root/"summary.json").read_text())
    expected = flow.blueprint(repository)
    flow.validate_blueprint_metadata(m, expected)
    if (sha256(root/"manifest.json") != summary["manifest_sha256"] or len(m["points"]) != 14
            or summary["new_SCF_calls"] != 14 or summary["new_NSCF_calls"] != 42
            or summary["full_G1_certified"] or summary["NSCF_energies_used_for_barriers"]
            or sha256(root/"endpoint_summary.json") != expected["endpoint_summary_sha256"]):
        raise ValueError("complete bounded measurements and unselected endpoints required")
    endpoints = json.loads((root/"endpoint_summary.json").read_text())
    elapsed, point_records = [], {}
    for point, frozen in zip(m["points"], expected["points"]):
        if any(point[k] != v for k, v in frozen.items()):
            raise ValueError("original terminal source changed")
        target = root/"calculations"/point["label"]
        records = summary["points"][point["label"]]
        scf = target/"scf"
        log = scf/"OUT.ABACUS/running_scf.log"
        if sha256(log) != records["SCF"]["log_sha256"]:
            raise ValueError("raw SCF log changed")
        native.check_runtime_version(log.read_text())
        raw = audited_results(scf)
        errors = {"energy_eV_cell": abs(float(raw["energy"]-point["baseline_results"]["energy"])),
            "force_eV_A": float(np.max(np.abs(raw["forces"]-point["baseline_results"]["forces"]))),
            "stress_kbar": float(np.max(np.abs(raw["stress"]-point["baseline_results"]["stress"]))) * 1602.176634}
        bands = sampled_band_gap((scf/"OUT.ABACUS/istate.info").read_text(), 48)
        if errors != records["SCF"]["SCF_output_only_errors"] or bands != records["SCF"]["bands"]:
            raise ValueError("raw SCF E/F/stress/band replay changed")
        if any(x > flow.TOLERANCES["SCF_"+k] for k, x in errors.items()):
            raise ValueError("output-only SCF reproduction failed")
        for stage in ("scf", *(f"nscf_22{n}" for n in native.GRIDS)):
            work = target/stage
            inputs = point["scf_input_sha256"] if stage == "scf" else next(x["input_sha256"] for x in point["nscf_inputs"] if x["nz"] == int(stage[-1]))
            for name in ("INPUT", "KPT", "STRU"):
                if sha256(work/name) != inputs[name]:
                    raise ValueError("exported input or geometry changed")
            if inputs["STRU"] != point["baseline_input_sha256"]["STRU"] or any(inputs[n] != h for n, h in CONTRACT.items() if n not in ("INPUT", "KPT")):
                raise ValueError("reported HF basis/geometry contract changed")
            body = (work/"INPUT").read_bytes()
            if stage == "scf":
                if not body.endswith(native.OUTPUT_ADDITION) or hashlib.sha256(body[:-len(native.OUTPUT_ADDITION)]).hexdigest() != CONTRACT["INPUT"]:
                    raise ValueError("SCF output-only delta changed")
                base = body[:-len(native.OUTPUT_ADDITION)]
                body_log = log.read_text()
            else:
                if body != native.berry_input(base):
                    raise ValueError("NSCF observable delta changed")
                nz = int(stage[-1])
                record = next(a for a in records["NSCF"] if a["nz"] == nz)
                file = work/"OUT.ABACUS/running_nscf.log"
                body_log = file.read_text()
                native.check_runtime_version(body_log)
                if sha256(file) != record["log_sha256"] or re.findall(r"\bDSIZE\s*=\s*(\d+)", body_log) != ["32"]:
                    raise ValueError("raw native NSCF log/ranks changed")
                berry = parse_abacus_berry(body_log)
                if any(record[k] != v for k, v in berry.items()):
                    raise ValueError("raw native polarization differs from summary")
                gap = sampled_nscf_band_gap((work/"OUT.ABACUS/BANDS_1.dat").read_text(), 48)
                if gap != record["bands"] or gap["n_kpoints"] != 4*nz or gap["sampled_indirect_gap_eV"] <= .1:
                    raise ValueError("raw NSCF band/grid/gap changed")
                q = float(np.linalg.norm(quantum_lattice(point["cell_A"])[2]))
                direction = np.array(point["cell_A"])[2]
                if (abs(q-record["physical_quantum_C_m2"]) > 1e-12
                        or abs(berry["reported_modulus_C_m2"]-2*q*native.NATIVE_SI_FACTOR) > 1.5e-7
                        or not np.allclose(berry["cartesian_projection_C_m2"], berry["value_C_m2"]*direction/np.linalg.norm(direction), rtol=0, atol=3e-7)
                        or record["value_modern_SI_C_m2"] != berry["value_C_m2"]/native.NATIVE_SI_FACTOR):
                    raise ValueError("actual-cell quantum, unit or vector projection changed")
            elapsed.append(native_elapsed_seconds(body_log))
        point_records[point["label"]] = records
    channel_records = []
    for original in summary["channels"]:
        name = original["name"]
        arrays = [endpoints["native_Berry_audits"]["PO_plus"]]
        arrays += [point_records[f"{name}_{i:02d}"]["NSCF"] for i in range(1, 8)]
        arrays += [endpoints["native_Berry_audits"][f"PO_minus_T_{name}"]]
        changes = [abs(modular_difference(a[-1]["value_modern_SI_C_m2"], a[-2]["value_modern_SI_C_m2"], a[-1]["reported_modulus_modern_SI_C_m2"])) for a in arrays]
        lift = sampled_reduced_branch([a[-1]["value_modern_SI_C_m2"] for a in arrays],
            [a[-1]["physical_quantum_C_m2"] for a in arrays],
            [a[-1]["reported_modulus_modern_SI_C_m2"] for a in arrays], changes)
        if (arrays != original["native_Berry_audits_in_source_order"]
                or changes != original["longitudinal_224_to_228_modular_change_C_m2"]
                or lift != original["sampled_branch"]):
            raise ValueError("sampled quadrature or branch replay changed")
        channel_records.append({"name": name, "maximum_224_to_228_change_C_m2": max(changes),
                                "sampled_branch": lift})
    return {"status": "raw_export_replayed", "new_DFT_calls": 0, "raw_SCF_logs": 14,
            "raw_NSCF_logs": 42, "raw_band_tables": 56, "native_integer_reported_DFT_seconds": sum(elapsed),
            "native_integer_reported_DFT_core_hours": sum(elapsed)*32/3600, "channels": channel_records,
            "runtime_resolution": "native h/min/s integer reports, not full-precision subprocess or allocation timing",
            "physical_files_and_charge_bytes_checked_on_HF_not_from_omitted_exports": True,
            "full_G1_certified": False, "limitations": flow.LIMITATIONS}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("export", "audit"))
    p.add_argument("--repository", type=Path, default=Path(__file__).resolve().parents[4])
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError("refusing previous offline audit")
    if a.action == "export":
        print(json.dumps({"exported_files": len(export(a.root, a.output)), "new_DFT_calls": 0}))
    else:
        result = audit(a.repository, a.root)
        flow.save(a.output, result)
        print(json.dumps({k: result[k] for k in ("status", "raw_SCF_logs", "raw_NSCF_logs", "native_integer_reported_DFT_core_hours")}))
