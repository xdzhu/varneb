"""Stage two same-600-eV GaN pressure-relaxation basin-link pilots.

The two signed seeds are the raw-audited, downhill q_u=+/-0.02 Angstrom
points of the local joint-mode frozen cut.  Only ionic/cell-relaxation tags
are added; the original VASP electronic contract is copied byte-for-byte.
The ten-step native VASP runs are diagnostics, not saddle or phase proofs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil

from ase.io import read

from scripts.prepare_gan_600eV_ts_hessian import geometry_preflight, same_geometry
from scripts.prepare_gan_600eV_ts_newton_canary import PRODUCTION_INPUT_SHA256


CASES = ("grid_um_vz", "grid_up_vz")
EXPECTED_Q = {"grid_um_vz": -0.02, "grid_up_vz": 0.02}
REPLACEMENTS = {"IBRION": "2", "ISIF": "3", "NSW": "10"}
NEW_TAGS = ("PSTRESS = 457.0", "EDIFFG = -0.02", "POTIM = 0.25")
STATIC_TAGS = {
    "ENCUT": "600.000000", "EDIFF": "1.00e-07", "ISYM": "-1",
    "SYMPREC": "1.00e-04", "IBRION": "-1", "ISIF": "2", "NSW": "0",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relaxation_incar(raw: bytes) -> bytes:
    """Change only relaxation controls; keep every electronic tag unchanged."""
    if b"\r" in raw or not raw.endswith(b"\n"):
        raise ValueError("source VASP INCAR must be LF-only")
    lines = raw.decode("ascii").splitlines()
    found = {}
    for index, line in enumerate(lines):
        match = re.fullmatch(r"\s*([A-Za-z_]+)\s*=\s*(.*?)\s*", line)
        if match is None:
            continue
        key, value = match.group(1).upper(), match.group(2)
        if key in found:
            raise ValueError(f"duplicate VASP INCAR key: {key}")
        found[key] = value
        if key in REPLACEMENTS:
            lines[index] = f" {key} = {REPLACEMENTS[key]}"
    if (any(found.get(key) != value for key, value in STATIC_TAGS.items())
            or any(key in found for key in ("PSTRESS", "EDIFFG", "POTIM"))):
        raise ValueError("source is not the unchanged 600-eV static INCAR")
    lines.extend(f" {tag}" for tag in NEW_TAGS)
    return ("\n".join(lines) + "\n").encode("ascii")


def prepare(source_root: Path, audit_path: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(output)
    source_manifest_path = source_root / "manifest.json"
    source_manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if (source_manifest.get("status") != "inputs_finalized_no_DFT"
            or source_manifest.get("purpose") != "GaN_45p7_600eV_local_joint_mode_2D_frozen_pilot"
            or source_manifest.get("pressure_GPa") != 45.7
            or source_manifest.get("input_contract_sha256") != PRODUCTION_INPUT_SHA256
            or audit.get("status") != "GaN_600eV_local_2D_frozen_enthalpy_pilot_raw_audited"
            or audit.get("pressure_GPa") != 45.7
            or audit.get("source_sha256", {}).get("manifest") != sha256(source_manifest_path)
            or audit.get("predeclared_gate", {}).get("negative_u_positive_v_at_both_scales") is not True):
        raise ValueError("local GaN 600-eV source/audit contract failed")
    manifest_cases = {record["name"]: record for record in source_manifest["cases"]}
    audited_cases = {record["case"]: record for record in audit["cases"]}
    if len(manifest_cases) != 12 or len(audited_cases) != 12:
        raise ValueError("incomplete local joint-mode static grid")
    prepared = []
    for name in CASES:
        source_case = source_root / "cases" / name
        record = manifest_cases[name]
        measured = audited_cases[name]
        input_hashes = record["input_sha256"]
        if (float(record["q_u_A"]) != EXPECTED_Q[name]
                or float(record["q_v_A"]) != 0.0
                or float(measured["q_u_A"]) != EXPECTED_Q[name]
                or float(measured["q_v_A"]) != 0.0
                or float(measured["delta_enthalpy_meV_per_GaN"]) >= 0.0
                or measured["input_sha256"] != input_hashes
                or any(input_hashes.get(key) != digest
                       for key, digest in PRODUCTION_INPUT_SHA256.items())
                or any(sha256(source_case / key) != digest
                       for key, digest in input_hashes.items())
                or sha256(source_case / "OUTCAR") != measured["outcar_sha256"]):
            raise ValueError(f"signed downhill seed or 600-eV raw source differs: {name}")
        atoms = read(source_case / "POSCAR", format="vasp")
        final = read(source_case / "OUTCAR")
        metrics = geometry_preflight(atoms)
        if (atoms.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]
                or not same_geometry(atoms, final)
                or metrics["empirical_near_symmetry_warning"]):
            raise ValueError(f"signed seed geometry or Bravais-risk gate failed: {name}")
        incar = relaxation_incar((source_case / "INCAR").read_bytes())
        prepared.append((name, source_case, incar, measured, metrics))
    output.mkdir(parents=True)
    (output / "cases").mkdir()
    result_cases = []
    for name, source_case, incar, measured, metrics in prepared:
        case = output / "cases" / name
        case.mkdir()
        for filename in ("POSCAR", "KPOINTS", "POTCAR"):
            shutil.copy2(source_case / filename, case / filename)
        (case / "INCAR").write_bytes(incar)
        hashes = {filename: sha256(case / filename)
                  for filename in ("POSCAR", "INCAR", "KPOINTS", "POTCAR")}
        if (hashes["KPOINTS"] != PRODUCTION_INPUT_SHA256["KPOINTS"]
                or hashes["POTCAR"] != PRODUCTION_INPUT_SHA256["POTCAR"]):
            raise ValueError("600-eV electronic inputs changed during staging")
        (case / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8"
        )
        result_cases.append({
            "name": name,
            "q_u_A": EXPECTED_Q[name],
            "source_static_delta_enthalpy_meV_per_GaN": measured["delta_enthalpy_meV_per_GaN"],
            "source_OUTCAR_sha256": measured["outcar_sha256"],
            "input_sha256": hashes,
            "geometry_preflight": metrics,
        })
    result = {
        "status": "inputs_finalized_no_DFT",
        "purpose": "GaN_45p7_600eV_signed_native_VASP_basin_10step_pilot",
        "claim_limit": "Ten-step two-sided descent pilot; not a TS or basin certificate",
        "pressure_GPa": 45.7,
        "PSTRESS_kbar": 457.0,
        "encut_eV": 600,
        "ionic_settings": {"IBRION": 2, "ISIF": 3, "NSW": 10, "EDIFFG_eV_per_A": -0.02,
                           "POTIM": 0.25},
        "n_cases": 2,
        "cases": result_cases,
        "source_sha256": {
            "static_grid_manifest": sha256(source_manifest_path),
            "static_grid_audit": sha256(audit_path),
            "preparer": sha256(Path(__file__)),
        },
        "limitations": [
            "A ten-step cap is not a convergence or phase-identity criterion.",
            "VASP PSTRESS adds P*V during relaxation; do not add P*V again to its printed enthalpy.",
            "Native VASP cell relaxation may expose a Bravais classifier failure at a later step.",
            "The same-600-eV strain energy/stress derivative discrepancy remains unresolved.",
        ],
    }
    (output / "manifest.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.source_root, args.audit, args.output)
    print(json.dumps({"status": result["status"], "n_cases": result["n_cases"]}))


if __name__ == "__main__":
    main()
