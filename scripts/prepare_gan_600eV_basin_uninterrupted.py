"""Stage two uninterrupted 600-eV GaN basin relaxations from original seeds.

The ten-step pilots used a constant plane-wave basis while changing the cell.
Their zero-byte WAVECAR files make VASP's ISTART=2 consistent continuation
impossible.  The same-geometry restart failed the 1-meV/cell energy gate on
one side, and restoring the original FFT grid did not fix it.  Rerunning from
the original signed seeds avoids silently splicing two basis contracts.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from scripts.prepare_gan_600eV_basin_restart_canary import CASES, sha256
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256


PURPOSE = "GaN_45p7_600eV_two_sided_uninterrupted_100step_basin_test"


def long_incar(source: bytes) -> bytes:
    if (b"\r" in source or source.count(b" NSW = 10\n") != 1
            or source.count(b" LWAVE = .FALSE.\n") != 1
            or source.count(b" ENCUT = 600.000000\n") != 1
            or source.count(b" PSTRESS = 457.0\n") != 1
            or source.count(b" ISIF = 3\n") != 1
            or source.count(b" IBRION = 2\n") != 1):
        raise ValueError("source pilot INCAR is not unchanged 600-eV/45.7-GPa contract")
    return (source.replace(b" NSW = 10\n", b" NSW = 100\n")
                  .replace(b" LWAVE = .FALSE.\n", b" LWAVE = .TRUE.\n"))


def prepare(pilot: Path, pilot_audit_path: Path,
            restart_audit_path: Path, fft_audit_path: Path,
            output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    pilot_manifest_path = pilot / "manifest.json"
    pilot_manifest = json.loads(pilot_manifest_path.read_text(encoding="utf-8"))
    pilot_audit = json.loads(pilot_audit_path.read_text(encoding="utf-8"))
    restart_audit = json.loads(restart_audit_path.read_text(encoding="utf-8"))
    fft_audit = json.loads(fft_audit_path.read_text(encoding="utf-8"))
    if (pilot_manifest.get("purpose")
            != "GaN_45p7_600eV_signed_native_VASP_basin_10step_pilot"
            or pilot_audit.get("status")
            != "GaN_600eV_basin_10step_pilots_raw_audited_not_phase_certificates"
            or pilot_audit.get("source_sha256", {}).get("manifest")
            != sha256(pilot_manifest_path)
            or restart_audit.get("status")
            != "GaN_600eV_basin_restart_canaries_raw_audited_not_future_cell_guarantee"
            or restart_audit.get("all_reproduce_within_1meV_per_cell") is not False
            or fft_audit.get("status")
            != "GaN_600eV_original_FFT_grid_one_step_raw_audited"
            or fft_audit.get("original_grid_reproduces_pilot_within_1meV_per_cell")
            is not False
            or fft_audit.get("source_sha256", {}).get("automatic_restart_audit")
            != sha256(restart_audit_path)
            or fft_audit.get("source_sha256", {}).get("pilot_OUTCAR")
            != sha256(pilot / "cases/grid_um_vz/OUTCAR")
            or any(record.get("pressure_GPa") != 45.7
                   for record in (pilot_manifest, pilot_audit, restart_audit, fft_audit))):
        raise ValueError("signed 600-eV pilot and failed-restart evidence differ")
    pilot_cases = {record["name"]: record for record in pilot_manifest["cases"]}
    audited_cases = {record["case"]: record for record in pilot_audit["cases"]}
    if set(pilot_cases) != set(CASES) or set(audited_cases) != set(CASES):
        raise ValueError("two signed pilot branches are required")
    checked = []
    for name in CASES:
        source = pilot / "cases" / name
        record, audited = pilot_cases[name], audited_cases[name]
        hashes = record["input_sha256"]
        if (sha256(source / "OUTCAR") != audited["outcar_sha256"]
                or audited["input_sha256"] != hashes
                or any(sha256(source / key) != digest for key, digest in hashes.items())
                or hashes["KPOINTS"] != PRODUCTION_INPUT_SHA256["KPOINTS"]
                or hashes["POTCAR"] != PRODUCTION_INPUT_SHA256["POTCAR"]):
            raise ValueError(f"signed pilot inputs/output changed: {name}")
        checked.append((name, source, record, audited))
    output.mkdir(parents=True)
    (output / "cases").mkdir()
    cases = []
    for name, source, record, audited in checked:
        destination = output / "cases" / name
        destination.mkdir()
        for key in ("POSCAR", "KPOINTS", "POTCAR"):
            shutil.copy2(source / key, destination / key)
        (destination / "INCAR").write_bytes(long_incar((source / "INCAR").read_bytes()))
        hashes = {key: sha256(destination / key)
                  for key in ("POSCAR", "INCAR", "KPOINTS", "POTCAR")}
        (destination / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8"
        )
        cases.append({
            "name": name,
            "source_pilot_OUTCAR_sha256": audited["outcar_sha256"],
            "source_pilot_input_sha256": record["input_sha256"],
            "input_sha256": hashes,
            "first_ten_step_comparator_H_eV_per_cell": (
                audited["first_evaluated_enthalpy_eV_per_cell"],
                audited["last_evaluated_enthalpy_eV_per_cell"]
            ),
        })
    summary = {
        "status": "inputs_finalized_no_DFT",
        "purpose": PURPOSE,
        "pressure_GPa": 45.7,
        "PSTRESS_kbar": 457.0,
        "encut_eV": 600,
        "ionic_settings": {"IBRION": 2, "ISIF": 3, "NSW": 100,
                           "EDIFFG_eV_per_A": -0.02, "POTIM": 0.25},
        "output_only_change": "LWAVE .FALSE. -> .TRUE. for possible ISTART=2 continuation",
        "n_cases": 2,
        "cases": cases,
        "source_sha256": {
            "pilot_manifest": sha256(pilot_manifest_path),
            "pilot_audit": sha256(pilot_audit_path),
            "automatic_restart_audit": sha256(restart_audit_path),
            "original_FFT_grid_audit": sha256(fft_audit_path),
            "preparer": sha256(Path(__file__)),
        },
        "limitations": [
            "Compare the first ten ionic enthalpies with the original pilot before using a long branch.",
            "A 100-step cap and process exit do not certify force/stress convergence or phase identity.",
            "LWAVE is output-only; VASP electronic, PAW, k-point, pressure and optimizer parameters remain unchanged.",
        ],
    }
    (output / "manifest.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("pilot", "pilot-audit", "restart-audit", "fft-audit", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    summary = prepare(args.pilot, args.pilot_audit, args.restart_audit,
                      args.fft_audit, args.output)
    print(json.dumps({"status": summary["status"], "n_cases": summary["n_cases"]}))


if __name__ == "__main__":
    main()
