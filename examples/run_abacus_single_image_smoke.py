"""Run one physical ABACUS energy/force/stress smoke calculation.

This is deliberately separate from a VC-NEB production run.  It checks the
calculator contract, generated INPUT/KPT/STRU files, and the stress channel on
one mapped HfO2 endpoint before any multi-image job is launched.  Its defaults
follow the production 100-Ry, full 10-au DZP policy.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb.abacus import attach_abacus_calculators, make_ase_abacus_factory
from vcneb.calculator import validate_image_calculators


DEFAULT_PSEUDO_DIR = "/home/zhuxd/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Pseudopotential"
DEFAULT_BASIS_DIR = "/home/zhuxd/abacus/PSEUDO/ABACUS-orbitals/Dojo-NC-FR/Orb-DZP-10au"
DEFAULT_COMMAND = "/home/zhuxd/Software/abacus/INSTALL/3.10.0-LTS/bin/abacus"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--structure",
        default=str(ROOT / "validation" / "hfo2_t_to_po" / "image_00" / "POSCAR"),
    )
    parser.add_argument("--workdir", default=str(ROOT / "run_VCNEB" / "abacus_single_image_smoke"))
    parser.add_argument("--pseudo-dir", default=DEFAULT_PSEUDO_DIR)
    parser.add_argument("--basis-dir", default=DEFAULT_BASIS_DIR)
    parser.add_argument("--command", default=DEFAULT_COMMAND)
    parser.add_argument("--ecutwfc", type=float, default=100.0)
    parser.add_argument("--scf-thr", type=float, default=1e-8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    structure = Path(args.structure).resolve()
    workdir = Path(args.workdir).resolve()
    atoms = read(structure)

    parameters = {
        "calculation": "scf",
        "basis_type": "lcao",
        "dft_functional": "pbe",
        "ecutwfc": args.ecutwfc,
        "scf_thr": args.scf_thr,
        "scf_nmax": 150,
        "mixing_type": "pulay",
        "mixing_beta": 0.3,
        "kpts": [2, 2, 2],
        "pp": {"Hf": "Hf.upf", "O": "O.upf"},
        "basis": {
            "Hf": "Hf_gga_10au_100Ry_4s2p2d1f.orb",
            "O": "O_gga_10au_100Ry_2s2p1d.orb",
        },
        "pseudo_dir": args.pseudo_dir,
        "basis_dir": args.basis_dir,
        "cal_force": 1,
        "cal_stress": 1,
        "out_stru": 1,
    }
    factory = make_ase_abacus_factory(parameters=parameters, command=args.command)
    attach_abacus_calculators([atoms], workdir=workdir, factory=factory)
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
    result = {
        "calculator": report.to_dict(),
        "structure": str(structure),
        "workdir": str(workdir),
        "natoms": len(atoms),
        "energy_eV": energy,
        "max_force_eV_per_A": float(np.linalg.norm(forces, axis=1).max()),
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
