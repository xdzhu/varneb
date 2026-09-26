import pytest

from ase import Atoms

from scripts.audit_gan_ts_basin_pilot import energy_triples, ga_n_coordination


def test_vasp_pressure_work_pair_is_read_once():
    log = """  free  energy   TOTEN  =       -24.00000000 eV
  free  energy   TOTEN  =       -22.09836893 eV
  enthalpy is  TOTEN    =       -11.52495148 eV   P V=       10.57341745
"""
    assert energy_triples(log) == pytest.approx([(-22.09836893, -11.52495148, 10.57341745)])


def test_vasp_pressure_work_mismatch_rejected():
    log = """free energy TOTEN = -22.0 eV
enthalpy is TOTEN = -10.0 eV P V= 10.0
"""
    with pytest.raises(ValueError, match="do not close"):
        energy_triples(log)


def test_ga_n_coordination_uses_periodic_neighbors():
    atoms = Atoms("Ga2N2", scaled_positions=[[0, 0, 0], [0.5, 0.5, 0],
                                               [0.25, 0.25, 0.25], [0.75, 0.75, 0.75]],
                  cell=[4, 4, 4], pbc=True)
    counts = ga_n_coordination(atoms, 2.4)
    assert len(counts) == 2
    assert all(value > 0 for value in counts)
