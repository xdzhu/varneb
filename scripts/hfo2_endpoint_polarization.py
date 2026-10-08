"""Bounded native LCAO Berry gate for PO+ and two inversion-mapped PO-.

SCF differs from the immutable historical contract only by output switches.
NSCF grids are observable quadrature on the SAME charge, never NEB energies.
Run endpoints only inside a checked allocation; no source output is changed.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess

from ase.io import read
import numpy as np

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.audit_hfo2_static_replica import INPUT_FILES, audited_results, sha256
from vcneb.polarization import (modular_difference, parse_abacus_berry, quantum_lattice,
                               sampled_band_gap, sampled_nscf_band_gap)


LABELS = ("PO_plus", "PO_minus_T_preserving", "PO_minus_T_reversing")
GRIDS = (2, 4, 8)  # transverse 2x2 unchanged; longitudinal convergence only
OUTPUT_ADDITION = b"\nout_chg 1\nout_bandgap 1\n"
ABI_SOURCE = "https://github.com/deepmodeling/abacus-develop/tree/f7cb1d3/examples/berryphase/lcao_PbTiO3"


def input_values(text):
    result = {}
    for line in text.splitlines():
        line = re.split(r"#|//", line, maxsplit=1)[0].strip()
        if not line or line == "INPUT_PARAMETERS":
            continue
        key, value = line.split(maxsplit=1)
        if key in result:
            raise ValueError(f"duplicate INPUT keyword {key}")
        result[key] = value
    return result


def berry_input(base):
    values = input_values(base.decode())
    values.update(calculation="nscf", init_chg="file", symmetry="-1", berry_phase="1", gdir="3",
                  cal_force="0", cal_stress="0", out_stru="0", out_band="1")
    # No changed cutoff, functional, spin, orbital integrals or solver.
    return ("INPUT_PARAMETERS\n" + "".join(f"{k} {v}\n" for k, v in values.items())).encode()


def prepare(sources, root, charge_caches=None):
    if root.exists():
        raise FileExistsError("refusing existing polarization namespace")
    if set(sources) != set(LABELS):
        raise ValueError("exactly the three preregistered ordered endpoints required")
    charge_caches = charge_caches or {}
    if set(charge_caches) - set(LABELS):
        raise ValueError("charge cache labels must match declared endpoints")
    records = []
    geometries = []
    for label in LABELS:
        source = Path(sources[label])
        if any(sha256(source / n) != h for n, h in CONTRACT.items()):
            raise ValueError(f"{label}: historical electronic contract changed")
        base = (source / "INPUT").read_bytes()
        if {"out_chg", "out_bandgap"} & input_values(base.decode()).keys():
            raise ValueError("unexpected output flags in original source")
        atoms = read(source / "STRU", format="abacus")
        if atoms.get_chemical_symbols() != ["Hf"] * 4 + ["O"] * 8:
            raise ValueError("ordered Hf4O8 required")
        if geometries and not np.allclose(atoms.cell.array, geometries[0].cell.array, atol=1e-10, rtol=0):
            raise ValueError("polarization comparison requires common cell")
        geometries.append(atoms)
        raw = audited_results(source)
        if (np.max(np.linalg.norm(raw["forces"], axis=1)) > .03
                or np.max(np.abs(raw["stress"])) * 1602.176634 > 2):
            raise ValueError("source is not an accepted stationary PO endpoint")
        valences = [float(re.search(r'z_valence\s*=\s*"\s*([^\"]+)"',
                                   (source / f"{s}.upf").read_text()).group(1)) for s in ("Hf", "O")]
        nelec = 4 * valences[0] + 8 * valences[1]
        if valences != [12., 6.] or nelec != 96:
            raise ValueError("audited nonmagnetic 96-electron pseudopotential contract differs")
        stage = root / "points" / label
        stage.mkdir(parents=True)
        scf = stage / "scf"
        scf.mkdir()
        for name in INPUT_FILES:
            shutil.copyfile(source / name, scf / name)
        (scf / "INPUT").write_bytes(base + OUTPUT_ADDITION)
        nscf_inputs = []
        for nz in GRIDS:
            folder = stage / f"nscf_22{nz}"
            folder.mkdir()
            for name in INPUT_FILES:
                shutil.copyfile(source / name, folder / name)
            (folder / "INPUT").write_bytes(berry_input(base))
            (folder / "KPT").write_bytes(f"K_POINTS\n0\nGamma\n2 2 {nz} 0 0 0\n".encode())
            nscf_inputs.append({"nz": nz, "input_sha256": {n: sha256(folder / n) for n in INPUT_FILES}})
        records.append({"label": label, "source_directory": str(source),
                        "cached_scf_directory": charge_caches.get(label),
                        "baseline_log_sha256": sha256(source / "OUT.ABACUS/running_scf.log"),
                        "baseline_results": {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in raw.items()},
                        "baseline_input_sha256": {n: sha256(source / n) for n in INPUT_FILES},
                        "scf_input_sha256": {n: sha256(scf / n) for n in INPUT_FILES},
                        "nscf_inputs": nscf_inputs, "occupied_bands": 48,
                        "cell_A": atoms.cell.array.tolist(), "quantum_lattice_C_m2": quantum_lattice(atoms.cell.array).tolist()})
    record = {"schema_version": 1, "purpose": "HfO2_three_endpoint_native_Berry_inversion_gate",
              "source_example": ABI_SOURCE, "abacus_commit": "f7cb1d3", "gdir": 3,
              "SCF_delta": {"out_chg": "1", "out_bandgap": "1"},
              "SCF_physical_settings_changed": False,
              "NSCF_delta": {"calculation": "nscf", "init_chg": "file", "symmetry": "-1",
                             "berry_phase": "1", "gdir": "3", "cal_force": "0", "cal_stress": "0", "out_stru": "0", "out_band": "1"},
              "NSCF_meshes": [[2, 2, n] for n in GRIDS], "NSCF_energies_used_for_barriers": False,
              "tolerances": {"SCF_energy_eV_cell": 1e-5, "SCF_force_eV_A": 1e-4, "SCF_stress_kbar": .02,
                             "sampled_gap_min_eV": .1, "longitudinal_P_convergence_C_m2": .01,
                             "inversion_modular_residual_C_m2": .01},
              "limitations": "R3 component only; longitudinal quadrature only; no automatic spontaneous-P/path-branch selection, no full-BZ insulating certificate",
              "points": records}
    (root / "manifest.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    for index, point in enumerate(records):
        if point["cached_scf_directory"]:
            cache = Path(point["cached_scf_directory"])
            audit = audit_scf(root, index, cache)
            point["cached_scf_audit"] = audit
    (root / "manifest.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return record


def check_inputs(folder, expected):
    if any(sha256(folder / n) != h for n, h in expected.items()):
        raise ValueError(f"staged input mutated: {folder}")


def audit_scf(root, index, work):
    manifest = json.loads((root / "manifest.json").read_text())
    point = manifest["points"][index]
    check_inputs(work, point["scf_input_sha256"])
    raw = audited_results(work)
    baseline = point["baseline_results"]
    errors = {"energy_eV_cell": abs(float(raw["energy"] - baseline["energy"])),
              "force_eV_A": float(np.max(np.abs(raw["forces"] - baseline["forces"]))),
              "stress_kbar": float(np.max(np.abs(raw["stress"] - baseline["stress"]))) * 1602.176634}
    if any(v > manifest["tolerances"]["SCF_" + k] for k, v in errors.items()):
        raise ValueError(f"output-only SCF reproducibility gate failed: {errors}")
    bands = sampled_band_gap((work / "OUT.ABACUS/istate.info").read_text(), point["occupied_bands"])
    if bands["sampled_indirect_gap_eV"] <= manifest["tolerances"]["sampled_gap_min_eV"]:
        raise ValueError("SCF sampled gap below predeclared gate")
    charges = list((work / "OUT.ABACUS").glob("SPIN*_CHG.cube"))
    if len(charges) != 1 or charges[0].stat().st_size == 0:
        raise ValueError("expected one completed nonmagnetic charge cube")
    return {"label": point["label"], "SCF_output_only_errors": errors, "bands": bands,
            "log_sha256": sha256(work / "OUT.ABACUS/running_scf.log"),
            "charge_sha256": {p.name: sha256(p) for p in charges}}


def audit_nscf(root, index, nz, work, scf_audit):
    manifest = json.loads((root / "manifest.json").read_text())
    point = manifest["points"][index]
    expected = next(x for x in point["nscf_inputs"] if x["nz"] == nz)
    check_inputs(work, expected["input_sha256"])
    log = work / "OUT.ABACUS/running_nscf.log"
    body = log.read_text()
    if re.findall(r"\bDSIZE\s*=\s*(\d+)", body) != ["32"] or not re.search(r"Total\s+Time\s*:", body):
        raise ValueError("incomplete or non-32-rank NSCF")
    for name, expected_hash in scf_audit["charge_sha256"].items():
        if sha256(work / "OUT.ABACUS" / name) != expected_hash:
            raise ValueError("fixed SCF charge changed during NSCF")
    berry = parse_abacus_berry(body)
    if berry["direction"] != 3:
        raise ValueError("unexpected Berry direction")
    bands = sampled_nscf_band_gap((work / "OUT.ABACUS/BANDS_1.dat").read_text(), point["occupied_bands"])
    if bands["n_kpoints"] != 4 * nz or bands["sampled_indirect_gap_eV"] <= manifest["tolerances"]["sampled_gap_min_eV"]:
        raise ValueError("NSCF grid/gap gate failed")
    q = float(np.linalg.norm(np.array(point["quantum_lattice_C_m2"])[2]))
    direction = np.array(point["cell_A"])[2]
    projection = berry["value_C_m2"] * direction / np.linalg.norm(direction)
    if not np.allclose(projection, berry["cartesian_projection_C_m2"], atol=3e-7, rtol=0):
        raise ValueError("Berry projection does not follow the actual cell vector")
    ratio = berry["reported_modulus_C_m2"] / q
    # This exact ABACUS source uses older SI constants and modulus=2 for
    # nspin=1 with even ionic valences. Do not silently divide printed P by 2.
    if abs(ratio - 2) > .001:
        raise ValueError("reported spin-paired modulus inconsistent with source/cell")
    return {"label": point["label"], "nz": nz, **berry, "bands": bands,
            "physical_quantum_C_m2": q, "native_modulus_over_physical_quantum": ratio,
            "log_sha256": sha256(log), "input_sha256": expected["input_sha256"]}


def run_endpoint(root, index):
    import os
    if (os.environ.get("RUN_DFT") != "1" or os.environ.get("SLURM_JOB_PARTITION") != "hfacnormal01"
            or os.environ.get("SLURM_JOB_NUM_NODES") != "1" or os.environ.get("SLURM_NTASKS") != "1"
            or os.environ.get("SLURM_CPUS_PER_TASK") != "32"):
        raise RuntimeError("explicit single-node32 allocation on hfacnormal01 required")
    manifest = json.loads((root / "manifest.json").read_text())
    point = manifest["points"][index]
    source = Path(point["source_directory"])
    check_inputs(source, point["baseline_input_sha256"])
    if sha256(source / "OUT.ABACUS/running_scf.log") != point["baseline_log_sha256"]:
        raise ValueError("baseline log changed")
    destination = root / "calculations" / point["label"]
    destination.mkdir(parents=True, exist_ok=False)
    argv = ["mpirun", "-np", "32", "/public/home/iai806/apprepo/abacus/v3.10.0LTS-intelmpi2025/app/bin/abacus"]

    def execute(stage, hashes, charge_source=None):
        work = destination / stage
        shutil.copytree(root / "points" / point["label"] / stage, work)
        check_inputs(work, hashes)
        if charge_source is not None:
            (work / "OUT.ABACUS").mkdir()
            for charge in (charge_source / "OUT.ABACUS").glob("SPIN*_CHG.cube"):
                shutil.copyfile(charge, work / "OUT.ABACUS" / charge.name)
        with (work / "abacus.out").open("w") as out, (work / "abacus.err").open("w") as err:
            subprocess.run(argv, cwd=work, stdout=out, stderr=err, check=True)
        return work

    scf = (Path(point["cached_scf_directory"]) if point.get("cached_scf_directory")
           else execute("scf", point["scf_input_sha256"]))
    scf_audit = audit_scf(root, index, scf)
    if point.get("cached_scf_audit") and scf_audit != point["cached_scf_audit"]:
        raise ValueError("cached SCF evidence changed after preflight")
    scf_audit["directory"] = str(scf)
    scf_audit["reused_completed_SCF"] = bool(point.get("cached_scf_directory"))
    (destination / "scf_audit.json").write_text(json.dumps(scf_audit, indent=2) + "\n")
    for inputs in point["nscf_inputs"]:
        nz = inputs["nz"]
        work = execute(f"nscf_22{nz}", inputs["input_sha256"], scf)
        audit = audit_nscf(root, index, nz, work, scf_audit)
        (destination / f"berry_22{nz}.json").write_text(json.dumps(audit, indent=2) + "\n")
        print(json.dumps({"label": point["label"], "nz": nz, "P_native_C_m2": audit["value_C_m2"]}), flush=True)


def summarize(root):
    manifest = json.loads((root / "manifest.json").read_text())
    audits = {label: [json.loads((root / "calculations" / label / f"berry_22{n}.json").read_text()) for n in GRIDS]
              for label in LABELS}
    convergence = {}
    for label, records in audits.items():
        mods = [r["reported_modulus_C_m2"] for r in records]
        if np.ptp(mods) > 1e-7:
            raise ValueError("modulus changes at fixed cell")
        convergence[label] = abs(modular_difference(records[-1]["value_C_m2"], records[-2]["value_C_m2"], mods[-1]))
    plus = audits["PO_plus"][-1]
    residuals = {label: abs(modular_difference(audits[label][-1]["value_C_m2"],
                                              -plus["value_C_m2"], plus["reported_modulus_C_m2"])) for label in LABELS[1:]}
    if any(abs(audits[label][-1]["reported_modulus_C_m2"] - plus["reported_modulus_C_m2"]) > 1e-7 for label in LABELS):
        raise ValueError("inversion endpoints have different polarization moduli")
    tol = manifest["tolerances"]
    passed = (all(x <= tol["longitudinal_P_convergence_C_m2"] for x in convergence.values())
              and all(x <= tol["inversion_modular_residual_C_m2"] for x in residuals.values()))
    # Modular opposition alone cannot distinguish a polar state from a
    # self-inverse 0/half-quantum class. Report, never erase, this ambiguity.
    self_inverse = abs(modular_difference(plus["value_C_m2"], -plus["value_C_m2"], plus["reported_modulus_C_m2"]))
    record = {"status": "passed" if passed and self_inverse > 2*tol["inversion_modular_residual_C_m2"] else "review_required",
              "manifest_sha256": sha256(root / "manifest.json"), "native_Berry_audits": audits,
              "SCF_output_only_audits": {l: json.loads((root / "calculations" / l / "scf_audit.json").read_text()) for l in LABELS},
              "longitudinal_224_to_228_modular_change_C_m2": convergence,
              "inversion_modular_residual_C_m2": residuals, "plus_self_inverse_distance_C_m2": self_inverse,
              "limitations": manifest["limitations"], "spontaneous_polarization_C_m2": None,
              "switching_path_branch_selected": False}
    (root / "summary.json").write_text(json.dumps(record, indent=2) + "\n")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "summarize"))
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--sources", type=Path)
    parser.add_argument("--charge-caches", type=Path)
    parser.add_argument("--index", type=int, choices=range(3))
    args = parser.parse_args()
    if args.action == "prepare":
        caches = json.loads(args.charge_caches.read_text()) if args.charge_caches else None
        report = prepare(json.loads(args.sources.read_text()), args.root, caches)
        print(json.dumps({"staged_endpoints": len(report["points"]), "NSCF_meshes": report["NSCF_meshes"]}))
    elif args.action == "run":
        if args.index is None:
            parser.error("--index required for run")
        run_endpoint(args.root, args.index)
    else:
        print(json.dumps({k: v for k, v in summarize(args.root).items() if k == "status" or "change" in k or "residual" in k}))


if __name__ == "__main__":
    main()
