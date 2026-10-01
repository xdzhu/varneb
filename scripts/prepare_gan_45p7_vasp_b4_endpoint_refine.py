"""Stage an isolated 600-eV cell relaxation of the original GaN B4 endpoint.

Only ionic/cell optimization tags differ from the completed production
endpoint static. This script never edits the original result tree or VCNEB.
Run it on hf, where the licensed POTCAR remains in the source static folder.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


SOURCE_HASHES = {
    "INCAR": "83ba34d4b14ff4ea641bd74dec26e462ea1cc24b575db79a07138ccdfd552772",
    "KPOINTS": "b5215f3e7608c27f49d45e8129d3edca2cfaa3ba88b8c389d4b52604715712ee",
    "POSCAR": "a59eb4fa87d0dea12cdb5f4856d65f3e4ecce11ffae2add4db2d22e0d10aeb47",
    "POTCAR": "f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168",
    "OUTCAR": "3c6795acf28dad6d00dd33a6f99cd755bf0b0cf8c9c4490bcc1fb3de3005dc87",
}
IONIC_CHANGES = {
    b" IBRION = -1\n": b" IBRION = 2\n",
    b" ISIF = 2\n": b" ISIF = 3\n",
    b" NSW = 0\n": b" NSW = 100\n",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(source_dir: Path, work_root: Path) -> dict:
    source_dir = Path(source_dir)
    work_root = Path(work_root)
    if work_root.exists():
        raise FileExistsError(work_root)
    for filename, expected in SOURCE_HASHES.items():
        if sha256(source_dir / filename) != expected:
            raise ValueError(f"original B4 source changed: {filename}")
    incar = (source_dir / "INCAR").read_bytes()
    if b"\r" in incar or b"PSTRESS" in incar or b"EDIFFG" in incar:
        raise ValueError("original B4 static INCAR is not the audited source")
    for old, new in IONIC_CHANGES.items():
        if incar.count(old) != 1:
            raise ValueError(f"missing unique static tag: {old!r}")
        incar = incar.replace(old, new)
    incar += b" EDIFFG = -0.02\n POTIM = 0.25\n PSTRESS = 457.0\n"
    case = work_root / "case"
    case.mkdir(parents=True)
    for filename in ("POSCAR", "KPOINTS", "POTCAR"):
        shutil.copy2(source_dir / filename, case / filename)
    (case / "INCAR").write_bytes(incar)
    input_hashes = {name: sha256(case / name)
                    for name in ("INCAR", "KPOINTS", "POSCAR", "POTCAR")}
    (case / "sha256.inputs.json").write_text(
        json.dumps(input_hashes, indent=2) + "\n", encoding="utf-8"
    )
    manifest = {
        "status": "inputs_finalized_no_DFT",
        "purpose": "GaN_45p7_600eV_original_B4_endpoint_cell_relax_not_path_replacement",
        "pressure_GPa": 45.7,
        "PSTRESS_kbar": 457.0,
        "electronic_contract": {
            "ENCUT_eV": 600, "PBE": True, "Ga_d_plus_N_PAW": True,
            "KPOINTS": "Gamma 8 8 6", "EDIFF_eV": 1e-7,
            "ISYM": -1, "SYMPREC": 1e-4,
        },
        "changed_from_static": [
            "IBRION -1 -> 2", "ISIF 2 -> 3", "NSW 0 -> 100",
            "EDIFFG -0.02 eV/A", "POTIM 0.25", "PSTRESS 457.0 kbar",
        ],
        "input_sha256": input_hashes,
        "source_sha256": {**SOURCE_HASHES, "preparer": sha256(Path(__file__))},
        "limitations": [
            "This is a targeted check of the original B4 endpoint's 2-kbar stress miss, not a replacement VCNEB path.",
            "VASP's ionic stop alone will not certify the 2-kbar stress gate; audit the final raw force/stress.",
            "If geometry changes, evaluate endpoints and peak under one contract before making a barrier claim.",
        ],
    }
    (work_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    args = parser.parse_args()
    report = prepare(args.source_dir, args.work_root)
    print(json.dumps({"status": report["status"],
                      "purpose": report["purpose"],
                      "input_sha256": report["input_sha256"]}))


if __name__ == "__main__":
    main()
