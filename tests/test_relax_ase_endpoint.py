from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS

import numpy as np

from scripts.relax_ase_endpoint import EndpointBFGS, _snapshot, _stress_residual_kbar


def test_endpoint_snapshot_records_finite_force_stress_and_enthalpy() -> None:
    atoms = Atoms(
        "Cu2",
        positions=[[0.0, 0.0, 0.0], [2.5, 0.0, 0.0]],
        cell=[5.0, 5.0, 5.0],
        pbc=True,
    )
    atoms.calc = EMT()
    optimizer = BFGS(atoms, logfile=None)
    snapshot = _snapshot(atoms, atoms, optimizer, pressure_gpa=1.0, fmax=0.10)
    assert snapshot["max_atomic_force_eV_per_A"] < 1.0e-8
    assert len(snapshot["stress_eV_per_A3"]) == 3
    assert snapshot["enthalpy_eV"] > snapshot["potential_energy_eV"]
    assert snapshot["max_abs_stress_residual_kbar"] >= 0.0


def test_endpoint_bfgs_requires_force_and_explicit_stress_convergence() -> None:
    atoms = Atoms("Cu", positions=[[0.0, 0.0, 0.0]], cell=[4.0, 4.0, 4.0], pbc=True)
    atoms.calc = EMT()
    optimizer = EndpointBFGS(
        atoms,
        endpoint_atoms=atoms,
        pressure_gpa=0.0,
        stress_kbar=1.0,
        logfile=None,
    )
    gradient = -atoms.get_forces().ravel()
    optimizer.fmax = 0.10
    assert _stress_residual_kbar(atoms, 0.0) > 1.0
    assert optimizer.gradient_converged(gradient) is False


def test_endpoint_bfgs_preserves_force_only_mode() -> None:
    atoms = Atoms("Cu", positions=[[0.0, 0.0, 0.0]], cell=[4.0, 4.0, 4.0], pbc=True)
    atoms.calc = EMT()
    optimizer = EndpointBFGS(
        atoms,
        endpoint_atoms=atoms,
        pressure_gpa=0.0,
        stress_kbar=None,
        logfile=None,
    )
    gradient = np.zeros(3)
    optimizer.fmax = 0.10
    assert bool(optimizer.gradient_converged(gradient)) is True
