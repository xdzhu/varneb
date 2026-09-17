"""Prepare a Qian-compatible four-atom GaN B4-to-B1 VASP case.

The B4 endpoint is copied from an explicit, user-licensed VASP template.
The B1 endpoint contains the same two formula units in a tetragonal
representation: its basal angle is 90 degrees and c/a is sqrt(2).  This is
the compatible cell convention needed for the B4 (hexagonal) to B1
(rocksalt) variable-cell transformation; both endpoints therefore have the
same atom count and Ga/N ordering.

This script only creates input structures and copies the supplied licensed
POTCAR.  It writes a dedicated static VASP input (600 eV, 8x8x6 Gamma mesh)
for this cross-backend reproduction.  It deliberately performs no DFT
calculation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys

import numpy as np
from ase import Atoms
from ase.io import read, write

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb.provenance import endpoint_structure_record


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--b4-poscar", required=True, help="four-atom wurtzite GaN POSCAR/CONTCAR")
    parser.add_argument("--source-dir", required=True, help="directory with licensed INCAR, KPOINTS and POTCAR")
    parser.add_argument("--destination", required=True, help="new case directory; must not already exist")
    parser.add_argument(
        "--b1-cubic-a",
        type=float,
        default=4.263,
        help="initial conventional rocksalt lattice constant in Angstrom",
    )
    return parser.parse_args()


def b1_two_formula_unit_cell(a_cubic: float) -> Atoms:
    """Return B1 GaN as two formula units with gamma=90 degrees.

    The cell vectors are (a/2,a/2,0), (-a/2,a/2,0), and (0,0,a).
    It has half the conventional B1 volume and therefore exactly two GaN
    formula units, matching the conventional four-atom wurtzite cell.
    """

    if a_cubic <= 0.0:
        raise ValueError("--b1-cubic-a must be positive")
    cell = np.array(
        [
            [0.5 * a_cubic, 0.5 * a_cubic, 0.0],
            [-0.5 * a_cubic, 0.5 * a_cubic, 0.0],
            [0.0, 0.0, a_cubic],
        ]
    )
    # Ga first, then N: preserve species block ordering expected by POTCAR.
    scaled = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
            [0.5, 0.5, 0.0],
            [0.0, 0.0, 0.5],
        ]
    )
    return Atoms(["Ga", "Ga", "N", "N"], scaled_positions=scaled, cell=cell, pbc=True)


def main() -> int:
    args = parse_args()
    b4_path = Path(args.b4_poscar).resolve()
    source = Path(args.source_dir).resolve()
    destination = Path(args.destination).resolve()
    required = [b4_path, source / "POTCAR"]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing required input(s): " + ", ".join(missing))
    if destination.exists():
        raise FileExistsError(f"destination already exists: {destination}")

    b4 = read(b4_path)
    if b4.get_chemical_symbols() != ["Ga", "Ga", "N", "N"]:
        raise ValueError("B4 POSCAR must contain ordered Ga Ga N N")
    b1 = b1_two_formula_unit_cell(args.b1_cubic_a)
    template = destination / "template"
    template.mkdir(parents=True)
    write(template / "B4.POSCAR", b4, format="vasp", direct=True, vasp5=True)
    write(template / "B1.POSCAR", b1, format="vasp", direct=True, vasp5=True)
    # Qian et al. used a dense 8x8x6 mesh and a 75-Ry QE ultrasoft setup.
    # This VASP PAW reproduction instead uses 600 eV; it is intentionally
    # recorded as a cross-backend comparison rather than a bitwise duplicate.
    (template / "INCAR").write_text(
        "PREC = Accurate\nENCUT = 600\nEDIFF = 1E-7\nISMEAR = 0\nSIGMA = 0.05\n"
        "IBRION = -1\nNSW = 0\nISIF = 2\nISYM = 0\nLCHARG = .FALSE.\nLWAVE = .FALSE.\n",
        encoding="utf-8",
    )
    (template / "KPOINTS").write_text(
        "GaN B4-B1 8x8x6 Gamma mesh\n0\nGamma\n8 8 6\n0 0 0\n", encoding="utf-8"
    )
    shutil.copy2(source / "POTCAR", template / "POTCAR")

    payload = {
        "status": "prepared_no_dft",
        "case": "GaN B4(wurtzite) to B1(rocksalt)",
        "formula_units": 2,
        "atom_count": len(b4),
        "b1_cubic_a_A": args.b1_cubic_a,
        "endpoint_structures": {"B4": endpoint_structure_record(b4), "B1": endpoint_structure_record(b1)},
        "source": {"b4_poscar": str(b4_path), "licensed_vasp_source": str(source)},
    }
    (destination / "case_setup.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
