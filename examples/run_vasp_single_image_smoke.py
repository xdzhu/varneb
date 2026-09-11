"""Run one VASP energy/force/stress smoke calculation through ASE."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys

import numpy as np
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb.calculator import validate_image_calculators
from vcneb.vasp import attach_vasp_calculators, default_vasp_command


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--structure", required=True)
    parser.add_argument("--template-dir", required=True)
    parser.add_argument("--workdir", required=True)
    parser.add_argument("--vasp-bin", default=os.environ.get("VASP_BIN", "vasp_std"))
    parser.add_argument("--ncores", type=int, default=int(os.environ.get("NP", "40")))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    structure = Path(args.structure).resolve()
    template_dir = Path(args.template_dir).resolve()
    workdir = Path(args.workdir).resolve()
    atoms = read(structure)
    attach_vasp_calculators(
        [atoms],
        source_dir=template_dir,
        workdir=workdir,
        command=default_vasp_command(args.ncores, args.vasp_bin),
        overrides={"xc": "PBE", "pp": "PBE"},
    )
    report = validate_image_calculators(
        [atoms],
        require_stress=True,
        require_variable_cell=True,
        require_directory=True,
        require_unique_directories=True,
    )[0]
    energy = float(atoms.get_potential_energy())
    forces = np.asarray(atoms.get_forces(), dtype=float)
    stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
    max_force = float(np.linalg.norm(forces, axis=1).max())
    write(workdir / "POSCAR.final", atoms, format="vasp", direct=True, vasp5=True)
    result = {
        "calculator": report.to_dict(),
        "structure": str(structure),
        "template_dir": str(template_dir),
        "workdir": str(workdir),
        "natoms": len(atoms),
        "energy_eV": energy,
        "max_force_eV_per_A": max_force,
        "stress_eV_per_A3": stress.tolist(),
        "stress_finite": bool(np.all(np.isfinite(stress))),
        "volume_A3": float(atoms.get_volume()),
    }
    (workdir / "smoke_summary.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
