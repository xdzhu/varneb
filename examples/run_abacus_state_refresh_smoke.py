"""Check that the ABACUS ASE calculator refreshes after q/cell changes."""

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--structure",
        default=str(ROOT / "validation" / "batio3_cubic_to_tetragonal" / "hf_endpoint_relax_cubic" / "CONTCAR"),
    )
    parser.add_argument("--workdir", default=str(ROOT / "outputs" / "abacus_state_refresh_smoke"))
    parser.add_argument("--pseudo-dir", required=True)
    parser.add_argument("--basis-dir", required=True)
    parser.add_argument("--command", required=True)
    parser.add_argument("--ecutwfc", type=float, default=60.0)
    parser.add_argument("--scf-thr", type=float, default=1e-6)
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
        "scf_nmax": 80,
        "mixing_type": "pulay",
        "mixing_beta": 0.3,
        "kpts": [1, 1, 1],
        "pp": {"Ba": "Ba.upf", "Ti": "Ti.upf", "O": "O.upf"},
        "basis": {
            "Ba": "Ba_gga_10au_100Ry_4s2p1d.orb",
            "Ti": "Ti_gga_8au_100Ry_4s2p2d1f.orb",
            "O": "O_gga_7au_100Ry_2s2p1d.orb",
        },
        "pseudo_dir": args.pseudo_dir,
        "basis_dir": args.basis_dir,
        "cal_force": 1,
        "cal_stress": 1,
        "out_stru": 1,
    }
    factory = make_ase_abacus_factory(parameters=parameters, command=args.command)
    attach_abacus_calculators([atoms], workdir=workdir, factory=factory)
    validate_image_calculators(
        [atoms],
        require_stress=True,
        require_variable_cell=True,
        require_directory=True,
        require_unique_directories=True,
    )

    energy0 = float(atoms.get_potential_energy())
    forces0 = np.asarray(atoms.get_forces(), dtype=float)
    stress0 = np.asarray(atoms.get_stress(voigt=False), dtype=float)
    scaled = atoms.get_scaled_positions(wrap=False)
    scaled[1, 2] += 0.01
    atoms.set_scaled_positions(scaled)
    cell = atoms.cell.array.copy()
    cell[2, 2] *= 1.01
    atoms.set_cell(cell, scale_atoms=False)
    energy1 = float(atoms.get_potential_energy())
    forces1 = np.asarray(atoms.get_forces(), dtype=float)
    stress1 = np.asarray(atoms.get_stress(voigt=False), dtype=float)
    changed = bool(
        abs(energy1 - energy0) > 1e-8
        or np.max(np.abs(forces1 - forces0)) > 1e-8
        or np.max(np.abs(stress1 - stress0)) > 1e-10
    )
    result = {
        "structure": str(structure),
        "workdir": str(workdir),
        "energy_before_eV": energy0,
        "energy_after_eV": energy1,
        "energy_delta_eV": energy1 - energy0,
        "max_force_before_eV_per_A": float(np.linalg.norm(forces0, axis=1).max()),
        "max_force_after_eV_per_A": float(np.linalg.norm(forces1, axis=1).max()),
        "max_stress_before_eV_per_A3": float(np.max(np.abs(stress0))),
        "max_stress_after_eV_per_A3": float(np.max(np.abs(stress1))),
        "state_refresh_changed_results": changed,
    }
    if not changed:
        raise RuntimeError("ABACUS calculator results did not change after q/cell perturbation")
    output = workdir / "state_refresh_summary.json"
    output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
