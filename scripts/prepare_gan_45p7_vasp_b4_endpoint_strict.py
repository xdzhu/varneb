"""Continue the unchanged GaN/B4 electronic contract with a stricter ion stop.

The first native cell relaxation stopped at its initial geometry under
EDIFFG=-0.02 eV/A although the raw-stress residual was 2.912 kbar. This
second, isolated continuation changes only EDIFFG and the step cap.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
from ase.io import read


EXPECTED_SOURCE_SHA256 = {
    "manifest.json": "05c85802129f8f18aa207505d47789951309b5a171dce387867165da39ad1b82",
    "case/INCAR": "27865542559e4c88230e04291a2f6c6468b3bc1d09e91e8d4404be8bf2feab5f",
    "case/KPOINTS": "b5215f3e7608c27f49d45e8129d3edca2cfaa3ba88b8c389d4b52604715712ee",
    "case/POSCAR": "a59eb4fa87d0dea12cdb5f4856d65f3e4ecce11ffae2add4db2d22e0d10aeb47",
    "case/POTCAR": "f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168",
    "case/OUTCAR": "a1581054b55783598899907799ac9581d7fa8c066ff71fcb00c8dbd65753bdc8",
    "case/CONTCAR": "c6d64143848485f34b2f73c873346720ff3ca31631c98f99d923f3f846069370",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(source_root: Path, work_root: Path) -> dict:
    source_root = Path(source_root)
    work_root = Path(work_root)
    if work_root.exists():
        raise FileExistsError(work_root)
    for name, expected in EXPECTED_SOURCE_SHA256.items():
        if sha256(source_root / name) != expected:
            raise ValueError(f"first native relaxation source changed: {name}")
    old = json.loads((source_root / "manifest.json").read_text(encoding="utf-8"))
    if (old.get("purpose")
            != "GaN_45p7_600eV_original_B4_endpoint_cell_relax_not_path_replacement"
            or old.get("electronic_contract", {}).get("ENCUT_eV") != 600):
        raise ValueError("first relaxation was not the audited 600-eV B4 case")
    raw_text = (source_root / "case/OUTCAR").read_text(
        encoding="utf-8", errors="replace"
    )
    if ("aborting loop because EDIFF is reached" not in raw_text
            or "reached required accuracy - stopping structural energy minimisation"
               not in raw_text
            or "General timing and accounting informations for this job"
               not in raw_text):
        raise ValueError("first relaxation was not cleanly completed")
    initial = read(source_root / "case/POSCAR", format="vasp")
    final = read(source_root / "case/CONTCAR", format="vasp")
    if (np.max(np.abs(initial.cell.array - final.cell.array)) > 1e-8
            or np.max(np.abs(initial.positions - final.positions)) > 1e-8):
        raise ValueError("first relaxation unexpectedly changed B4 geometry")
    incar = (source_root / "case/INCAR").read_bytes()
    for old_bytes, new_bytes in (
        (b" NSW = 100\n", b" NSW = 30\n"),
        (b" EDIFFG = -0.02\n", b" EDIFFG = -0.005\n"),
    ):
        if incar.count(old_bytes) != 1:
            raise ValueError(f"missing unique first-run tag: {old_bytes!r}")
        incar = incar.replace(old_bytes, new_bytes)
    case = work_root / "case"
    case.mkdir(parents=True)
    shutil.copy2(source_root / "case/CONTCAR", case / "POSCAR")
    for name in ("KPOINTS", "POTCAR"):
        shutil.copy2(source_root / "case" / name, case / name)
    (case / "INCAR").write_bytes(incar)
    hashes = {name: sha256(case / name)
              for name in ("INCAR", "KPOINTS", "POSCAR", "POTCAR")}
    (case / "sha256.inputs.json").write_text(
        json.dumps(hashes, indent=2) + "\n", encoding="utf-8"
    )
    manifest = {
        "status": "inputs_finalized_no_DFT",
        "purpose": "GaN_45p7_600eV_original_B4_endpoint_strict_ionic_continuation",
        "pressure_GPa": 45.7, "PSTRESS_kbar": 457.0,
        "electronic_contract": old["electronic_contract"],
        "changed_from_first_relaxation": ["EDIFFG -0.02 -> -0.005 eV/A",
                                          "NSW 100 -> 30"],
        "input_sha256": hashes,
        "source_sha256": {**EXPECTED_SOURCE_SHA256,
                          "preparer": sha256(Path(__file__))},
        "limitations": [
            "No production VCNEB endpoint is replaced by preparing or running this case.",
            "The final raw pressure residual must be audited; the VASP stop is not a 2-kbar certificate.",
            "Physical electronic settings and external pressure are unchanged.",
        ],
    }
    (work_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    args = parser.parse_args()
    report = prepare(args.source_root, args.work_root)
    print(json.dumps({"status": report["status"],
                      "purpose": report["purpose"],
                      "input_sha256": report["input_sha256"]}))


if __name__ == "__main__":
    main()
