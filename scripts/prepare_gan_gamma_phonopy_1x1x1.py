"""Stage GaN endpoint Gamma-only Phonopy/VASP finite displacements.

The 1x1x1 matrix means *no expansion* of the four-atom VCNEB cell. Staging
is non-DFT; independent VASP statics run only under the Slurm array. Both
step sizes are retained so force-constant sensitivity can be audited before
any physical soft-mode claim.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read, write
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms


POTCAR_SHA256 = "f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168"
SOURCE_SUBDIRS = {"B4": "hf_static_initial/00", "B1": "hf_static_final/28"}
DISTANCES_A = (0.01, 0.02)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_static_source(source: Path) -> tuple[Atoms, dict[str, str]]:
    required = ("CONTCAR", "INCAR", "KPOINTS", "POTCAR", "OUTCAR")
    if any(not (source / name).is_file() or (source / name).stat().st_size == 0 for name in required):
        raise ValueError(f"incomplete VASP endpoint static source: {source}")
    if sha256(source / "POTCAR") != POTCAR_SHA256:
        raise ValueError("POTCAR differs from audited Ga_d+N production input")
    incar = {}
    for line in (source / "INCAR").read_text(encoding="utf-8").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            incar[key.strip().upper()] = value.strip().split("#", 1)[0].strip()
    for key, expected in {"ENCUT": 600, "EDIFF": 1e-7, "SYMPREC": 1e-4,
                          "ISYM": -1, "ISMEAR": 0, "SIGMA": 0.05,
                          "IBRION": -1, "NSW": 0}.items():
        if key not in incar or not np.isclose(float(incar[key]), expected, rtol=1e-9, atol=1e-12):
            raise ValueError(f"production VASP contract mismatch: {key}")
    if incar.get("GGA") != "PE" or incar.get("PREC", "").lower() != "accurate":
        raise ValueError("production PBE/PREC contract mismatch")
    kpoints = (source / "KPOINTS").read_text(encoding="utf-8").splitlines()
    if (len(kpoints) < 5 or kpoints[2].strip().lower() != "gamma"
            or kpoints[3].split() != ["8", "8", "6"]
            or kpoints[4].split() != ["0", "0", "0"]):
        raise ValueError("production electronic 8x8x6 Gamma k mesh mismatch")
    if "General timing and accounting informations" not in (source / "OUTCAR").read_text(
        encoding="utf-8", errors="replace"
    ):
        raise ValueError("endpoint static OUTCAR is not complete")
    atoms = read(source / "CONTCAR", format="vasp")
    if (atoms.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]
            or atoms.get_volume() <= 0 or not np.isfinite(atoms.positions).all()):
        raise ValueError("endpoint is not the aligned four-atom GaN path cell")
    return atoms, {name: sha256(source / name) for name in required}


def make_phonopy(atoms: Atoms, distance_A: float) -> Phonopy:
    unit = PhonopyAtoms(
        symbols=atoms.get_chemical_symbols(),
        cell=np.asarray(atoms.cell.array),
        scaled_positions=atoms.get_scaled_positions(wrap=False),
    )
    phonon = Phonopy(
        unit, supercell_matrix=np.eye(3, dtype=int),
        primitive_matrix=np.eye(3), symprec=1e-4, is_symmetry=True,
    )
    phonon.generate_displacements(
        distance=distance_A, is_plusminus=True, is_diagonal=False,
    )
    if not np.array_equal(phonon.supercell_matrix, np.eye(3, dtype=int)):
        raise ValueError("GaN Gamma preparation unexpectedly expanded the cell")
    return phonon


def stage(case_root: Path, work_root: Path) -> dict:
    if work_root.exists():
        raise FileExistsError(f"refusing to replace an existing work root: {work_root}")
    sources = {}
    for phase, relative in SOURCE_SUBDIRS.items():
        source = case_root / relative
        atoms, hashes = validate_static_source(source)
        sources[phase] = {"directory": str(source), "atoms": atoms, "sha256": hashes}
    if sources["B4"]["sha256"]["POTCAR"] != sources["B1"]["sha256"]["POTCAR"]:
        raise ValueError("endpoint POTCAR identity differs")
    work_root.mkdir(parents=True)
    (work_root / "cases").mkdir()
    (work_root / "logs").mkdir()
    records = []
    for phase in ("B4", "B1"):
        source = sources[phase]
        atoms = source["atoms"]
        for distance in DISTANCES_A:
            phonon = make_phonopy(atoms, distance)
            family_dir = work_root / f"{phase}_d{distance:.2f}"
            family_dir.mkdir()
            phonon.save(filename=str(family_dir / "phonopy_disp.yaml"))
            displacements = phonon.dataset["first_atoms"]
            for local_index, displaced in enumerate(phonon.supercells_with_displacements):
                if displaced is None:
                    raise ValueError("Phonopy generated an empty displacement")
                if len(displaced) != len(atoms):
                    raise ValueError("1x1x1 Phonopy changed the endpoint atom count")
                name = f"{phase}_d{distance:.2f}_{local_index:03d}"
                target = work_root / "cases" / name
                target.mkdir()
                for filename in ("INCAR", "KPOINTS", "POTCAR"):
                    shutil.copy2(Path(source["directory"]) / filename, target / filename)
                shifted = Atoms(
                    symbols=displaced.symbols,
                    cell=displaced.cell,
                    scaled_positions=displaced.scaled_positions,
                    pbc=True,
                )
                write(target / "POSCAR", shifted, format="vasp", direct=True, sort=False, vasp5=True)
                reread = read(target / "POSCAR", format="vasp")
                expected = np.asarray(displacements[local_index]["displacement"], dtype=float)
                atom_index = int(displacements[local_index]["number"])
                difference = reread.positions - atoms.positions
                difference -= np.rint(difference @ np.linalg.inv(atoms.cell.array)) @ atoms.cell.array
                if (reread.get_chemical_symbols() != atoms.get_chemical_symbols()
                        or not np.allclose(reread.cell.array, atoms.cell.array, atol=1e-8)
                        or not np.allclose(difference[atom_index], expected, atol=1e-7)
                        or np.max(np.abs(np.delete(difference, atom_index, axis=0))) > 1e-7):
                    raise ValueError(f"Phonopy/POSCAR displaced geometry mismatch: {name}")
                hashes = {filename: sha256(target / filename)
                          for filename in ("POSCAR", "INCAR", "KPOINTS", "POTCAR")}
                (target / "sha256.inputs.json").write_text(
                    json.dumps(hashes, indent=2) + "\n", encoding="utf-8"
                )
                records.append({
                    "name": name, "phase": phase, "distance_A": distance,
                    "local_index": local_index, "atom_index": atom_index,
                    "displacement_A": expected.tolist(), "input_sha256": hashes,
                })
    manifest = {
        "purpose": "GaN_1x1x1_endpoint_Gamma_force_constants_two_step_sizes",
        "pressure_GPa": 45.7,
        "supercell_matrix": np.eye(3, dtype=int).tolist(),
        "phonon_symprec": 1e-4,
        "force_calculator": "VASP 6.3.2 PBE/Ga_d+N PAW ENCUT600 Gamma8x8x6 EDIFF1e-7 ISYM=-1 SYMPREC1e-4",
        "source": {phase: {key: val for key, val in item.items() if key != "atoms"}
                   for phase, item in sources.items()},
        "cases": records,
        "limitations": "fixed-cell atomic Gamma modes only; no finite-q dispersion, NAC, or variable-cell TS certification",
    }
    (work_root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-root", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    args = parser.parse_args()
    manifest = stage(args.case_root, args.work_root)
    print(json.dumps({"status": "staged_no_DFT", "n_cases": len(manifest["cases"]),
                      "work_root": str(args.work_root)}))


if __name__ == "__main__":
    main()
