"""Run a tiny real LAMMPS or CP2K static stress/force smoke.

This is deliberately not a material benchmark.  It verifies the backend
launch, isolated image directory, and ASE energy/force/stress contract before a
user selects a potential or basis set for a real VCNEB case.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ase import Atoms

from vcneb import (
    attach_image_calculators,
    inspect_calculator,
    make_ase_cp2k_factory,
    make_ase_lammps_factory,
    make_ase_espresso_factory,
    validate_qe_pseudopotentials,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("backend", choices=("lammps", "cp2k", "qe"))
    parser.add_argument("--workdir", default="backend_smoke")
    parser.add_argument("--command", default=None)
    parser.add_argument("--pseudo-dir", default=None, help="QE UPF directory")
    return parser.parse_args()


def run_lammps(workdir: Path, command: str | None) -> dict:
    # A two-atom Lennard-Jones cell exercises periodic stress and non-zero forces
    # without introducing a material-specific potential claim.
    images = [Atoms("Ar2", positions=[[0, 0, 0], [3.8, 0, 0]], cell=[12, 12, 12], pbc=True)]
    factory = make_ase_lammps_factory(
        parameters={
            "units": "metal",
            "atom_style": "atomic",
            "specorder": ["Ar"],
            "masses": ["1 39.948"],
            "pair_style": "lj/cut 10.0",
            "pair_coeff": ["* * 0.0103 3.4 10.0"],
        },
        command=command,
    )
    attach_image_calculators(images, workdir=workdir, factory=factory)
    return _evaluate(images[0])


def run_cp2k(workdir: Path, command: str | None) -> dict:
    # The HF CP2K module supplies the GTH-PBE basis/potential files.  This is a
    # one-molecule shell smoke, not a converged solid-state calculation.
    images = [Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]], cell=[8, 8, 8], pbc=True)]
    factory = make_ase_cp2k_factory(
        parameters={
            "basis_set": "DZVP-MOLOPT-SR-GTH",
            "pseudo_potential": "GTH-PBE",
            "cutoff": 300,
            "xc": "PBE",
            "charge": 0,
            "max_scf": 200,
            "inp": """
&FORCE_EVAL
  &DFT
    &SCF
      EPS_SCF 1.0E-6
      SCF_GUESS ATOMIC
      &OT
        MINIMIZER DIIS
        PRECONDITIONER FULL_SINGLE_INVERSE
      &END OT
    &END SCF
  &END DFT
&END FORCE_EVAL
""",
            "debug": False,
        },
        command=command,
    )
    attach_image_calculators(images, workdir=workdir, factory=factory)
    return _evaluate(images[0])


def run_qe(workdir: Path, command: str | None, pseudo_dir: str | None) -> dict:
    if not pseudo_dir:
        raise ValueError("QE smoke requires --pseudo-dir containing an explicit PBE UPF")
    pseudo_root = Path(pseudo_dir).resolve()
    pseudo_name = "H.SG15.PBE.UPF"
    validation = validate_qe_pseudopotentials(pseudo_root, {"H": pseudo_name})
    images = [Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.74]], cell=[8, 8, 8], pbc=True)]
    factory = make_ase_espresso_factory(
        parameters={
            "pseudopotentials": {"H": pseudo_name},
            "input_data": {
                "control": {"prefix": "varneb_qe_smoke", "outdir": str(workdir / "qe_tmp")},
                "system": {"ecutwfc": 30, "ecutrho": 240, "occupations": "fixed"},
                "electrons": {"conv_thr": 1.0e-8},
            },
            "kpts": (1, 1, 1),
        },
        command=command or "pw.x",
        pseudo_dir=pseudo_root,
    )
    attach_image_calculators(images, workdir=workdir, factory=factory)
    report = _evaluate(images[0])
    report["pseudopotential_validation"] = validation
    return report


def _evaluate(image: Atoms) -> dict:
    report = inspect_calculator(image.calc, require_stress=True, require_variable_cell=True)
    if not report.ok:
        raise RuntimeError("backend preflight failed: " + "; ".join(report.issues))
    energy = float(image.get_potential_energy())
    forces = image.get_forces()
    stress = image.get_stress()
    return {
        "calculator": report.to_dict(),
        "energy_eV": energy,
        "max_force_eV_per_A": float(abs(forces).max()),
        "stress_eV_per_A3": [float(value) for value in stress],
    }


def main() -> None:
    args = parse_args()
    workdir = Path(args.workdir).resolve()
    if args.backend == "lammps":
        report = run_lammps(workdir, args.command)
    elif args.backend == "cp2k":
        report = run_cp2k(workdir, args.command)
    else:
        report = run_qe(workdir, args.command, args.pseudo_dir)
    output = workdir / f"{args.backend}_smoke.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
