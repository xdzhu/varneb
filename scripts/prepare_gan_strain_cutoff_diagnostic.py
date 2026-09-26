"""Stage seven isolated VASP cutoff-sensitivity statics for GaN strain.

Only ENCUT changes from the audited 600 eV image-15 contract. This is a
numerical *diagnostic*, not a replacement of the production VCNEB path or an
adjustment of the unrelated BTO/ABACUS 100 Ry orbital contract.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil

import numpy as np
from ase.io import read

def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def same_geometry(a, b) -> bool:
    delta = a.get_scaled_positions(wrap=False) - b.get_scaled_positions(wrap=False)
    delta -= np.rint(delta)
    return (a.get_chemical_symbols() == b.get_chemical_symbols()
            and np.allclose(a.cell.array, b.cell.array, rtol=0, atol=1e-7)
            and np.allclose(delta, 0, rtol=0, atol=1e-7))


def stage(source_work: Path, work_root: Path, encut_eV: int = 800) -> dict:
    if work_root.exists():
        raise FileExistsError(work_root)
    if encut_eV <= 600 or encut_eV > 1200:
        raise ValueError("diagnostic cutoff must be above 600 and at most 1200 eV")
    source_manifest_path = source_work / "manifest.json"
    original = json.loads(source_manifest_path.read_text(encoding="utf-8"))
    if (original.get("status") != "inputs_finalized_no_DFT"
            or original.get("step_A") != 0.02
            or original.get("coordinate_count") != 18
            or len(original.get("cases", [])) != 36):
        raise ValueError("source must be the complete 0.02 Å GaN joint-curvature stage")
    template = source_work / "cases/axis12_plus"
    initial_incar = (template / "INCAR").read_text(encoding="utf-8")
    revised_incar, replacement_count = re.subn(
        r"(?mi)^(\s*ENCUT\s*=\s*)600\.000000(\s*)$",
        lambda match: f"{match.group(1)}{encut_eV}.000000{match.group(2)}",
        initial_incar,
    )
    if replacement_count != 1:
        raise ValueError("cannot make a single isolated ENCUT substitution")
    center = read(source_work / "center_POSCAR", format="vasp")
    selected = ["center"] + [f"axis{axis:02d}_{side}" for axis in (12, 13, 14)
                             for side in ("plus", "minus")]
    work_root.mkdir(parents=True)
    (work_root / "cases").mkdir()
    (work_root / "logs").mkdir()
    records = []
    for name in selected:
        source_poscar = (source_work / "center_POSCAR" if name == "center"
                         else source_work / "cases" / name / "POSCAR")
        target = work_root / "cases" / name
        target.mkdir()
        shutil.copy2(source_poscar, target / "POSCAR")
        if (not same_geometry(read(target / "POSCAR", format="vasp"),
                              read(source_poscar, format="vasp"))
                or len(read(target / "POSCAR", format="vasp")) != len(center)):
            raise ValueError(f"geometry changed while staging {name}")
        (target / "INCAR").write_text(revised_incar, encoding="utf-8", newline="\n")
        for filename in ("KPOINTS", "POTCAR"):
            shutil.copy2(template / filename, target / filename)
        hashes = {filename: sha256(target / filename)
                  for filename in ("POSCAR", "INCAR", "KPOINTS", "POTCAR")}
        (target / "sha256.inputs.json").write_text(
            json.dumps(hashes, indent=2) + "\n", encoding="utf-8",
        )
        records.append({"name": name, "POSCAR_sha256": hashes["POSCAR"],
                        "input_sha256": hashes})
    manifest = {
        "purpose": "GaN_image15_VASP_strain_energy_stress_cutoff_sensitivity_only",
        "status": "inputs_finalized_no_DFT",
        "encut_eV": encut_eV,
        "baseline_encut_eV": 600,
        "step_A": 0.02,
        "pressure_GPa": original["pressure_GPa"],
        "cell_scale_A": original["cell_scale_A"],
        "source_600_work_root": str(source_work),
        "source_600_manifest_sha256": sha256(source_manifest_path),
        "source_600_trajectory_sha256": original["source_sha256"]["trajectory"],
        "source_600_summary_sha256": original["source_sha256"]["summary"],
        "source_static_OUTCAR_sha256": original["source_static_sha256"]["OUTCAR"],
        "preparer_sha256": sha256(Path(__file__)),
        "cases": records,
        "limitations": "Isolated local cutoff diagnostic; does not update production VCNEB or certify a TS.",
    }
    (work_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8",
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-work", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--encut-eV", type=int, default=800)
    args = parser.parse_args()
    result = stage(args.source_work, args.work_root, args.encut_eV)
    print(json.dumps({"status": result["status"], "n_cases": len(result["cases"]),
                      "encut_eV": result["encut_eV"]}))


if __name__ == "__main__":
    main()
