"""Regression smoke for the ABACUS-native endpoint input/output bridge."""

from __future__ import annotations

from argparse import Namespace
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.relax_abacus_native import (
    _last_force_stress,
    _parse_stru,
    _summary,
    _write_input,
    _write_kpt,
    _write_stru,
)


def main() -> None:
    source = ROOT / "validation" / "hfo2_t_to_po" / "T_HfO2_12.vasp"
    atoms = read(source)
    with TemporaryDirectory(prefix="vcneb-native-driver-") as tmp:
        root = Path(tmp)
        pp_dir = root / "pp"
        basis_dir = root / "basis"
        pp_dir.mkdir()
        basis_dir.mkdir()
        for name in ("Hf.upf", "O.upf"):
            (pp_dir / name).write_text("placeholder\n", encoding="utf-8")
        for name in ("Hf.orb", "O.orb"):
            (basis_dir / name).write_text("placeholder\n", encoding="utf-8")
        for name in ("Hf.upf", "O.upf", "Hf.orb", "O.orb"):
            (root / name).write_text("placeholder\n", encoding="utf-8")
        _write_stru(
            atoms,
            root / "STRU",
            {"Hf": "Hf.upf", "O": "O.upf"},
            {"Hf": "Hf.orb", "O": "O.orb"},
            pp_dir,
            basis_dir,
        )
        # A post-processing/recovery run may point pseudo/orbital directories
        # at the endpoint workdir itself; this must not raise SameFileError.
        _write_stru(
            atoms,
            root / "STRU_same_dir",
            {"Hf": "Hf.upf", "O": "O.upf"},
            {"Hf": "Hf.orb", "O": "O.orb"},
            root,
            root,
        )
        parsed = _parse_stru(root / "STRU")
        assert len(parsed) == 12
        assert parsed.get_chemical_formula() == "Hf4O8"
        assert np.allclose(parsed.cell.array, atoms.cell.array)
        assert parsed.get_chemical_symbols() == atoms.get_chemical_symbols()
        args = Namespace(
            structure=str(source),
            ecutwfc=100.0,
            scf_thr=1e-8,
            scf_nmax=150,
            kpts=[2, 2, 2],
            relax_nmax=200,
            force_thr=0.02,
            stress_thr=0.1,
            mixing_beta=0.3,
            relax_method="bfgs",
        )
        _write_input(root / "INPUT", args, 2)
        _write_kpt(root / "KPT", [2, 2, 2])
        text = (root / "INPUT").read_text(encoding="utf-8")
        assert "calculation       cell-relax" in text
        assert "relax_method      bfgs" in text
        assert "stress_thr        0.1" in text
        # Also accept the format emitted by ABACUS STRU_ION_D: species
        # labels, an explicit magnetic-moment line, and ``m 1 1 1`` flags.
        restart = root / "STRU_ION_D"
        restart.write_text(
            """ATOMIC_SPECIES
Hf 178.49 Hf.upf
O 15.999 O.upf

LATTICE_CONSTANT
1.8897261258

LATTICE_VECTORS
5 0 0
0 5 0
0 0 5

ATOMIC_POSITIONS
Direct

Hf #label
0.0 #magnetism
1 #number of atoms
0 0 0 m 1 1 1

O #label
0.0 #magnetism
1 #number of atoms
0.5 0.5 0.5 m 1 1 1
""",
            encoding="utf-8",
        )
        restart_parsed = _parse_stru(restart)
        assert len(restart_parsed) == 2
        assert restart_parsed.get_chemical_symbols() == ["Hf", "O"]
        out = root / "OUT.ABACUS"
        out.mkdir()
        (out / "running_cell-relax.log").write_text(
            "Largest gradient in force is 0.0123 eV/A\n"
            "TOTAL-PRESSURE: -0.40 KBAR\n"
            "Lattice relaxation is converged\n",
            encoding="utf-8",
        )
        summary = _summary(root, args, 0, atoms)
        assert summary["converged"] is False  # no STRU_ION_D in this fixture
        assert summary["last_largest_gradient_eV_per_A"] == 0.0123
        assert summary["last_total_pressure_kbar"] == -0.4
        force, stress = _last_force_stress(
            "TOTAL-FORCE (eV/Angstrom)\n"
            " Hf1 0.003 0.004 0.000\n O1 0.000 0.000 0.010\n"
            "TOTAL-STRESS (KBAR)\n1 2 3\n2 4 5\n3 5 6\n",
            2,
        )
        assert abs(force - 0.01) < 1e-12
        assert stress == 6.0
    print("native_abacus_driver_regression=ok")


if __name__ == "__main__":
    main()
