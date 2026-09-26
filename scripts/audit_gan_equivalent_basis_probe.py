"""Audit attempted equivalent-cell VASP statics without inventing unrun results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read

from scripts.audit_gan_joint_curvature import completed_case
from scripts.audit_gan_ts_newton_probe import failed_bravais_case
from scripts.prepare_gan_equivalent_basis_probe import (
    TRANSFORMS, physically_equivalent,
)
from scripts.prepare_gan_ts_newton_probe import sha256


def audit(work_root: Path, previous_work: Path, previous_manifest: Path,
          previous_audit: Path) -> dict:
    manifest_path = work_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    old = json.loads(previous_manifest.read_text(encoding="utf-8"))
    previous = json.loads(previous_audit.read_text(encoding="utf-8"))
    source_case = previous_work / "cases/newton_full"
    if (manifest["purpose"] != "GaN_full_step_equivalent_basis_VASP_Bravais_diagnostic_only"
            or manifest["status"] != "inputs_finalized_no_DFT"
            or manifest["evaluation_encut_eV"] != 1000
            or manifest["kmesh"] != [8, 8, 6]
            or [case["name"] for case in manifest["cases"]] != list(TRANSFORMS)
            or manifest["source_sha256"]["previous_manifest"] != sha256(previous_manifest)
            or manifest["source_sha256"]["previous_audit"] != sha256(previous_audit)
            or manifest["source_sha256"]["full_failed_OUTCAR"] != sha256(source_case / "OUTCAR")
            or manifest["source_sha256"]["preparer"] != sha256(
                Path(__file__).with_name("prepare_gan_equivalent_basis_probe.py")
            ) or previous["cases"][1]["outcar_sha256"] != sha256(source_case / "OUTCAR")):
        raise ValueError("equivalent-basis stage and failed-source provenance mismatch")
    original = read(source_case / "POSCAR", format="vasp")
    if (manifest["original_full_step_POSCAR_sha256"] != sha256(source_case / "POSCAR")
            or old["cases"][1]["POSCAR_sha256"] != sha256(source_case / "POSCAR")):
        raise ValueError("failed full-step geometry changed")
    results = []
    for case in manifest["cases"]:
        directory = work_root / "cases" / case["name"]
        matrix = np.asarray(case["unimodular_matrix"], dtype=int)
        if (not np.array_equal(matrix, np.asarray(TRANSFORMS[case["name"]]))
                or not physically_equivalent(original, read(directory / "POSCAR", format="vasp"),
                                             matrix)):
            raise ValueError(f"candidate is not an approved equivalent basis: {case['name']}")
        for filename, expected in case["input_sha256"].items():
            if sha256(directory / filename) != expected:
                raise ValueError(f"staged input changed: {directory / filename}")
        if (case["input_sha256"]["KPOINTS"] != old["cases"][1]["input_sha256"]["KPOINTS"]
                or case["input_sha256"]["POTCAR"] != old["cases"][1]["input_sha256"]["POTCAR"]
                or case["input_sha256"]["INCAR"] != old["cases"][1]["input_sha256"]["INCAR"]):
            raise ValueError(f"electronic contract changed: {case['name']}")
        if not (directory / "OUTCAR").exists():
            if (directory / "vasp.stdout").exists():
                raise ValueError(f"stdout without OUTCAR: {case['name']}")
            results.append({"case": case["name"], "status": "not_evaluated",
                            "input_sha256": case["input_sha256"]})
            continue
        text = (directory / "OUTCAR").read_text(encoding="utf-8", errors="replace")
        if "Inconsistent Bravais lattice types found for crystalline and" in text:
            results.append(failed_bravais_case(directory, case))
        else:
            atoms, result = completed_case(directory, case)
            if not physically_equivalent(original, atoms, matrix):
                raise ValueError(f"VASP output geometry changed: {case['name']}")
            results.append({**result, "status": "complete_converged_static"})
    return {
        "status": "GaN_equivalent_basis_diagnostic_audited_not_TS_certificate",
        "results": results,
        "source_sha256": {
            "manifest": sha256(manifest_path),
            "previous_manifest": sha256(previous_manifest),
            "previous_audit": sha256(previous_audit),
            "original_failed_POSCAR": sha256(source_case / "POSCAR"),
            "auditor": sha256(Path(__file__)),
        },
        "interpretation_boundary": (
            "Only attempted equivalent representations are evidence. Failure of two "
            "signed permutations does not prove that no equivalent basis can work; "
            "the untouched third case is not a failed calculation."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--previous-work", type=Path, required=True)
    parser.add_argument("--previous-manifest", type=Path, required=True)
    parser.add_argument("--previous-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    report = audit(args.work_root, args.previous_work,
                   args.previous_manifest, args.previous_audit)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"],
                      "results": [[case["case"], case["status"]]
                                  for case in report["results"]]}))


if __name__ == "__main__":
    main()
