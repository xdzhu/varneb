"""Recheck GaN's negative restart on its original pilot FFT grids at 600 eV.

The audited negative branch auto-switched coarse/fine grids on restart.
This isolates that numerical difference without changing ENCUT, geometry,
k-point mesh, pseudopotentials, pressure, or ionic algorithm.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from scripts.prepare_gan_600eV_basin_restart_canary import sha256


CASE = "grid_um_vz"
PURPOSE = "GaN_45p7_600eV_negative_basin_original_FFT_grid_one_step_canary"
FFT_TAGS = ("NGX = 24", "NGY = 24", "NGZ = 40",
            "NGXF = 48", "NGYF = 48", "NGZF = 80")


def grid_locked_incar(raw: bytes) -> bytes:
    if (b"\r" in raw or not raw.endswith(b"\n")
            or raw.count(b" NSW = 1\n") != 1
            or raw.count(b" ENCUT = 600.000000\n") != 1
            or raw.count(b" PSTRESS = 457.0\n") != 1
            or any(f" {tag.split(' =')[0]} =".encode() in raw for tag in FFT_TAGS)):
        raise ValueError("source is not the same-contract automatic-grid restart canary")
    return raw + b"".join(f" {tag}\n".encode("ascii") for tag in FFT_TAGS)


def prepare(canary_root: Path, canary_audit_path: Path,
            pilot_root: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    canary_manifest_path = canary_root / "manifest.json"
    canary_manifest = json.loads(canary_manifest_path.read_text(encoding="utf-8"))
    canary_audit = json.loads(canary_audit_path.read_text(encoding="utf-8"))
    negative = next(item for item in canary_audit["cases"] if item["case"] == CASE)
    source = canary_root / "cases" / CASE
    pilot_outcar = (pilot_root / "cases" / CASE / "OUTCAR").read_text(
        encoding="utf-8", errors="replace"
    )
    if (canary_manifest.get("purpose")
            != "GaN_45p7_600eV_native_VASP_basin_restart_one_step_canary"
            or canary_audit.get("source_sha256", {}).get("manifest")
            != sha256(canary_manifest_path)
            or negative["reproduces_pilot_enthalpy_within_1meV_per_cell"]
            or abs(negative["enthalpy_difference_from_pilot_last_eV_per_cell"]) < 0.005
            or sha256(source / "OUTCAR") != negative["outcar_sha256"]
            or "dimension x,y,z NGX =    24 NGY =   24 NGZ =   40" not in pilot_outcar
            or "dimension x,y,z NGXF=    48 NGYF=   48 NGZF=   80" not in pilot_outcar):
        raise ValueError("source does not show the measured negative-grid mismatch")
    source_outcar = (source / "OUTCAR").read_text(encoding="utf-8", errors="replace")
    if ("dimension x,y,z NGX =    24 NGY =   24 NGZ =   36" not in source_outcar
            or "dimension x,y,z NGXF=    48 NGYF=   48 NGZF=   72" not in source_outcar):
        raise ValueError("automatic restart did not have the observed smaller FFT grids")
    record = next(item for item in canary_manifest["cases"] if item["name"] == CASE)
    if any(sha256(source / key) != digest for key, digest in record["input_sha256"].items()):
        raise ValueError("automatic restart input hash changed")
    directory = output / "case"
    directory.mkdir(parents=True)
    for filename in ("POSCAR", "KPOINTS", "POTCAR"):
        shutil.copy2(source / filename, directory / filename)
    (directory / "INCAR").write_bytes(grid_locked_incar((source / "INCAR").read_bytes()))
    hashes = {key: sha256(directory / key)
              for key in ("POSCAR", "INCAR", "KPOINTS", "POTCAR")}
    (directory / "sha256.inputs.json").write_text(
        json.dumps(hashes, indent=2) + "\n", encoding="utf-8"
    )
    manifest = {
        "status": "inputs_finalized_no_DFT",
        "purpose": PURPOSE,
        "case": CASE,
        "pressure_GPa": 45.7,
        "PSTRESS_kbar": 457.0,
        "encut_eV": 600,
        "original_pilot_coarse_FFT_grid": [24, 24, 40],
        "original_pilot_fine_FFT_grid": [48, 48, 80],
        "automatic_restart_coarse_FFT_grid": [24, 24, 36],
        "automatic_restart_fine_FFT_grid": [48, 48, 72],
        "source_automatic_restart_enthalpy_eV_per_cell": negative["enthalpy_eV_per_cell"],
        "source_pilot_last_enthalpy_eV_per_cell": (
            negative["enthalpy_eV_per_cell"]
            - negative["enthalpy_difference_from_pilot_last_eV_per_cell"]
        ),
        "input_sha256": hashes,
        "source_sha256": {
            "automatic_restart_manifest": sha256(canary_manifest_path),
            "automatic_restart_audit": sha256(canary_audit_path),
            "automatic_restart_OUTCAR": sha256(source / "OUTCAR"),
            "pilot_OUTCAR": sha256(pilot_root / "cases" / CASE / "OUTCAR"),
            "preparer": sha256(Path(__file__)),
        },
        "limitations": [
            "The FFT mesh is a numerical representation; grid locking is not a higher ENCUT.",
            "One-step agreement does not ensure future variable-cell FFT adequacy.",
        ],
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("canary-root", "canary-audit", "pilot-root", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    manifest = prepare(args.canary_root, args.canary_audit,
                       args.pilot_root, args.output)
    print(json.dumps({"status": manifest["status"], "purpose": manifest["purpose"]}))


if __name__ == "__main__":
    main()
