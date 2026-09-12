"""Generate a six-component strain finite-difference report.

The six variables are three diagonal deformation components and three
engineering-style symmetric shears.  The latter perturb ``F_ij`` and ``F_ji``
by the same scalar, so the analytic comparison is ``G_ij + G_ji``.  The
report is calculator-independent and uses the analytic VCNEB test potential.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
from ase import Atoms
from ase.units import GPa

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.run_toy_vcneb import ToyPhaseTransition
from vcneb.core import cell_force, cell_from_deformation


def make_atoms(reference_cell: np.ndarray, q: np.ndarray, deform: np.ndarray) -> Atoms:
    atoms = Atoms(
        "Ar",
        scaled_positions=[q],
        cell=cell_from_deformation(deform, reference_cell),
        pbc=True,
    )
    atoms.calc = ToyPhaseTransition(reference_cell)
    return atoms


def enthalpy(reference_cell: np.ndarray, q: np.ndarray, deform: np.ndarray, pressure: float) -> float:
    atoms = make_atoms(reference_cell, q, deform)
    return float(atoms.get_potential_energy() + pressure * atoms.get_volume())


def main() -> None:
    reference_cell = np.array(
        [[4.7, 0.2, 0.0], [0.4, 5.1, 0.3], [0.1, 0.2, 5.4]],
        dtype=float,
    )
    q = np.array([0.42, 0.57, 0.48], dtype=float)
    deform = np.array(
        [[1.08, 0.02, 0.01], [0.00, 0.97, 0.03], [0.01, 0.00, 1.03]],
        dtype=float,
    )
    epsilons = [1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
    pressures = [0.0, 0.5 * GPa]
    results = []
    for pressure in pressures:
        atoms = make_atoms(reference_cell, q, deform)
        analytic = cell_force(atoms, reference_cell, pressure=pressure)
        for epsilon in epsilons:
            diagonal_errors = []
            shear_errors = []
            for index in range(3):
                direction = np.zeros((3, 3))
                direction[index, index] = 1.0
                plus = enthalpy(reference_cell, q, deform + epsilon * direction, pressure)
                minus = enthalpy(reference_cell, q, deform - epsilon * direction, pressure)
                numeric = -(plus - minus) / (2.0 * epsilon)
                diagonal_errors.append(abs(numeric - analytic[index, index]))
            for row, col in [(0, 1), (0, 2), (1, 2)]:
                direction = np.zeros((3, 3))
                direction[row, col] = 1.0
                direction[col, row] = 1.0
                plus = enthalpy(reference_cell, q, deform + epsilon * direction, pressure)
                minus = enthalpy(reference_cell, q, deform - epsilon * direction, pressure)
                numeric = -(plus - minus) / (2.0 * epsilon)
                expected = analytic[row, col] + analytic[col, row]
                shear_errors.append(abs(numeric - expected))
            results.append(
                {
                    "pressure_GPa": pressure / GPa,
                    "epsilon": epsilon,
                    "max_diagonal_error_eV": max(diagonal_errors),
                    "max_symmetric_shear_error_eV": max(shear_errors),
                    "max_six_component_error_eV": max(diagonal_errors + shear_errors),
                }
            )

    report = {
        "calculator": "analytic ToyPhaseTransition",
        "cell_convention": "h = h0 F.T",
        "shear_convention": "F_ij and F_ji are perturbed by the same scalar",
        "pressures_GPa": [pressure / GPa for pressure in pressures],
        "epsilons": epsilons,
        "results": results,
    }
    output = ROOT / "outputs" / "vcneb_finite_difference_report.json"
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    for pressure in pressures:
        subset = [item for item in results if item["pressure_GPa"] == pressure / GPa]
        best = min(subset, key=lambda item: item["epsilon"])
        print(
            f"pressure_GPa={pressure / GPa:.6g} epsilon={best['epsilon']:.1e} "
            f"max_six_component_error_eV={best['max_six_component_error_eV']:.3e}"
        )
    print(f"output={output.resolve()}")


if __name__ == "__main__":
    main()
