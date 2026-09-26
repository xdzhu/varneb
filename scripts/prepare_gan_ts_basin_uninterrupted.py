"""Stage a long, uninterrupted GaN plus-branch relaxation from an audited seed.

The ten-step pilot proved that the original seed parses and starts correctly.
Restarting from its tenth CONTCAR triggers VASP's inconsistent-Bravais check;
this test keeps all physical inputs fixed and changes only the NSW cap.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from scripts.prepare_gan_ts_newton_probe import sha256


def prepare(pilot_root: Path, pilot_audit_path: Path,
            failed_root: Path, work_root: Path) -> dict:
    if work_root.exists():
        raise FileExistsError(work_root)
    pilot_manifest_path = pilot_root / "manifest.json"
    pilot_manifest = json.loads(pilot_manifest_path.read_text(encoding="utf-8"))
    audit = json.loads(pilot_audit_path.read_text(encoding="utf-8"))
    if (pilot_manifest.get("purpose")
            != "GaN_1000eV_signed_basin_native_VASP_10step_pilot_not_endpoint_certificate"
            or audit.get("status")
            != "GaN_native_basin_10step_pilots_raw_audited_not_phase_or_endpoint_certificates"
            or audit.get("source_sha256", {}).get("manifest") != sha256(pilot_manifest_path)):
        raise ValueError("source ten-step pilot is not raw audited")
    name = "negative_mode_plus"
    source_record = next(case for case in pilot_manifest["cases"] if case["name"] == name)
    audit_record = next(case for case in audit["cases"] if case["case"] == name)
    source = pilot_root / "cases" / name
    if (not audit_record["contcar_is_last_evaluated_geometry"]
            or sha256(source / "OUTCAR") != audit_record["outcar_sha256"]
            or any(sha256(source / filename) != expected
                   for filename, expected in source_record["input_sha256"].items())):
        raise ValueError("source plus pilot output/input changed")
    failed = failed_root / "cases" / name / "OUTCAR"
    if "Inconsistent Bravais lattice types found" not in failed.read_text(
            encoding="utf-8", errors="replace"):
        raise ValueError("intermediate restart has not demonstrated the Bravais failure")
    incar = (source / "INCAR").read_bytes()
    if b"\r" in incar or incar.count(b"NSW = 10\n") != 1:
        raise ValueError("pilot INCAR is not the expected LF-only ten-step input")
    destination = work_root / "case"
    destination.mkdir(parents=True)
    for filename in ("POSCAR", "KPOINTS", "POTCAR"):
        shutil.copy2(source / filename, destination / filename)
    (destination / "INCAR").write_bytes(incar.replace(b"NSW = 10\n", b"NSW = 100\n"))
    hashes = {filename: sha256(destination / filename)
              for filename in ("POSCAR", "INCAR", "KPOINTS", "POTCAR")}
    (destination / "sha256.inputs.json").write_text(
        json.dumps(hashes, indent=2) + "\n", encoding="utf-8",
    )
    result = {
        "purpose": "GaN_plus_basin_100step_uninterrupted_from_audited_seed_not_endpoint_certificate",
        "status": "inputs_finalized_no_DFT", "case": name,
        "unchanged": ["source POSCAR", "POTCAR", "KPOINTS Gamma 8x8x6",
                      "ENCUT 1000 eV", "ISYM -1", "SYMPREC 1e-4",
                      "PSTRESS 457 kbar", "IBRION 2", "ISIF 3", "POTIM 0.25"],
        "changed": ["NSW cap 10 -> 100 from the original signed seed"],
        "input_sha256": hashes,
        "source_sha256": {"pilot_manifest": sha256(pilot_manifest_path),
                          "pilot_audit": sha256(pilot_audit_path),
                          "pilot_OUTCAR": audit_record["outcar_sha256"],
                          "failed_restart_OUTCAR": sha256(failed),
                          "preparer": sha256(Path(__file__))},
        "limitations": [
            "The first ten steps must match the audited pilot within numerical tolerance.",
            "A successful process exit or 100 steps does not certify force/stress convergence or phase identity.",
        ],
    }
    (work_root / "manifest.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-root", type=Path, required=True)
    parser.add_argument("--pilot-audit", type=Path, required=True)
    parser.add_argument("--failed-root", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.pilot_root, args.pilot_audit, args.failed_root, args.work_root)
    print(json.dumps({"status": result["status"], "changed": result["changed"]}))


if __name__ == "__main__":
    main()
