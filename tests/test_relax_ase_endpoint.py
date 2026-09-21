from ase import Atoms
from ase.calculators.emt import EMT
from ase.optimize import BFGS

from scripts.relax_ase_endpoint import _snapshot


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
