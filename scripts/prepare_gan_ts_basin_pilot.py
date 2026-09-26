"""Stage short pressure-consistent native VASP relaxations from signed GaN probes.

These ten-step pilots check descent and phase identity. They do not alter the
original VCNEB chain or assert completed B4/B1 basin connections.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil

from ase.io import read

from scripts.prepare_gan_ts_newton_probe import same_geometry, sha256


def relaxation_incar(static_text: str) -> str:
    values = {
        "IBRION": "2", "ISIF": "3", "NSW": "10",
        "PSTRESS": "457.0", "EDIFFG": "-0.02", "POTIM": "0.25",
    }
    required = {
        "ENCUT": "1000.000000", "EDIFF": "1.00e-07", "ISYM": "-1",
        "SYMPREC": "1.00e-04", "IBRION": "-1", "ISIF": "2", "NSW": "0",
    }
    found = {}
    output = []
    for line in static_text.splitlines():
        match = re.fullmatch(r"\s*([A-Za-z_]+)\s*=\s*(.*?)\s*", line)
        if match:
            key = match.group(1).upper()
            if key in found:
                raise ValueError(f"duplicate static INCAR key: {key}")
            found[key] = match.group(2)
            if key in values:
                line = f" {key} = {values.pop(key)}"
        output.append(line)
    if any(found.get(key) != value for key, value in required.items()):
        raise ValueError("source static VASP INCAR contract changed")
    if "PSTRESS" in found or "EDIFFG" in found or "POTIM" in found:
        raise ValueError("static source unexpectedly has relaxation or pressure tags")
    output.extend(f" {key} = {value}" for key, value in values.items())
    return "\n".join(output) + "\n"


def prepare(signed_root: Path, signed_audit: Path, work_root: Path) -> dict:
    if work_root.exists():
        raise FileExistsError(work_root)
    signed_manifest_path = signed_root / "manifest.json"
    signed_manifest = json.loads(signed_manifest_path.read_text(encoding="utf-8"))
    audit = json.loads(signed_audit.read_text(encoding="utf-8"))
    if (signed_manifest.get("purpose") != "GaN_near_TS_1000eV_signed_unstable_mode_static_downhill_probe_not_basin_link"
            or signed_manifest.get("status") != "inputs_finalized_no_DFT"
            or audit.get("status") != "GaN_signed_unstable_mode_two_statics_audited_not_basin_links"
            or not audit.get("both_sides_downhill")
            or audit.get("source_sha256", {}).get("manifest") != sha256(signed_manifest_path)
            or signed_manifest.get("pressure_GPa") != 45.7
            or signed_manifest.get("encut_eV") != 1000):
        raise ValueError("signed negative-mode statics are not an audited downhill pair")
    audited = {item["case"]: item for item in audit["cases"]}
    if set(audited) != {"negative_mode_plus", "negative_mode_minus"}:
        raise ValueError("signed downhill audit does not cover both cases")
    work_root.mkdir(parents=True)
    cases = []
    for source_record in signed_manifest["cases"]:
        name = source_record["name"]
        source_dir = signed_root / "cases" / name
        source_outcar = source_dir / "OUTCAR"
        if sha256(source_outcar) != audited[name]["outcar_sha256"]:
            raise ValueError(f"signed source VASP output changed: {name}")
        if not same_geometry(read(source_dir / "POSCAR", format="vasp"),
                             read(source_outcar), tolerance=2e-5):
            raise ValueError(f"signed static structure/output differs: {name}")
        static_text = (source_dir / "INCAR").read_text(encoding="utf-8")
        new_incar = relaxation_incar(static_text)
        directory = work_root / "cases" / name
        directory.mkdir(parents=True)
        for filename in ("POSCAR", "KPOINTS", "POTCAR"):
            shutil.copy2(source_dir / filename, directory / filename)
        # VASP 6.3.2 on hf rejected CRLF here with IERR=5 while reading IBRION.
        # Explicit bytes prevent Windows text-mode newline conversion.
        (directory / "INCAR").write_bytes(new_incar.encode("ascii"))
        hashes = {filename: sha256(directory / filename)
                  for filename in ("POSCAR", "INCAR", "KPOINTS", "POTCAR")}
        (directory / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8",
        )
        cases.append({
            "name": name, "sign": source_record["sign"],
            "input_sha256": hashes,
            "source_OUTCAR_sha256": audited[name]["outcar_sha256"],
            "source_static_enthalpy_eV_per_cell": audited[name]["enthalpy_eV_per_cell"],
        })
    result = {
        "purpose": "GaN_1000eV_signed_basin_native_VASP_10step_pilot_not_endpoint_certificate",
        "status": "inputs_finalized_no_DFT", "pressure_GPa": 45.7,
        "PSTRESS_kbar": 457.0, "encut_eV": 1000,
        "ionic_settings": {"IBRION": 2, "ISIF": 3, "NSW": 10,
                            "EDIFFG_eV_per_A": -0.02, "POTIM": 0.25},
        "cases": cases,
        "source_sha256": {
            "signed_manifest": sha256(signed_manifest_path),
            "signed_audit": sha256(signed_audit),
            "preparer": sha256(Path(__file__)),
        },
        "limitations": [
            "The ten-step cap is a pilot, not an endpoint convergence criterion.",
            "VASP PSTRESS includes P*V in its relaxation energy; do not add P*V a second time.",
            "Final phase identity and comparable static E+P*V require independent review.",
        ],
    }
    (work_root / "manifest.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("signed-root", "signed-audit", "work-root"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.signed_root, args.signed_audit, args.work_root)
    print(json.dumps({"status": result["status"], "n_cases": len(result["cases"])}))


if __name__ == "__main__":
    main()
