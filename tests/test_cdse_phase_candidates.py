from pathlib import Path

import numpy as np
import pytest
from ase.io import read

from scripts.prepare_cdse_phase_candidates import phase_candidate
from scripts.prepare_cdse_sheppard_suite import symmetry, mapped_chain
from scripts.verify_cdse_endpoint_statics import raw_stationarity

FIXTURES = Path(__file__).resolve().parent / "fixtures/cdse_phase"


@pytest.mark.parametrize("label,expected", [("rs",225),("wz",186)])
def test_actual_failed_endpoint_is_bounded_candidate_not_a_convergence_claim(label, expected):
    path = FIXTURES / f"{label}_relaxed.vasp"
    atoms = read(path)
    cell, positions = atoms.cell.array.copy(), atoms.positions.copy()
    candidate, audit = phase_candidate(atoms, label)
    assert all(record["number"] == expected for record in symmetry(candidate))
    assert audit["candidate_only_not_stationarity_verified"]
    assert audit["maximum_principal_stretch_change"] < 1e-3
    assert audit["maximum_motif_displacement_A"] < 1e-3
    assert np.array_equal(atoms.cell.array, cell)
    assert np.array_equal(atoms.positions, positions)


def test_large_distortion_cannot_be_relabeled_as_target_phase():
    path = FIXTURES / "rs_relaxed.vasp"
    atoms = read(path)
    cell = atoms.cell.array.copy()
    cell[2] *= 1.02
    atoms.set_cell(cell, scale_atoms=True)
    with pytest.raises(ValueError, match="too far"):
        phase_candidate(atoms, "rs")


def test_raw_static_gate_checks_stress_without_symmetry_projection():
    path = FIXTURES / "wz_relaxed.vasp"
    atoms = read(path)
    summary = {"status":"completed", "potential_energy_eV":-20.,
               "forces_eV_per_A":np.zeros((8,3)).tolist(), "stress_eV_per_A3_voigt":np.zeros(6).tolist()}
    assert raw_stationarity(atoms, summary)["additional_DFT_evaluations"] == 0
    summary["stress_eV_per_A3_voigt"][2] = .1
    with pytest.raises(ValueError, match="not stationary"):
        raw_stationarity(atoms, summary)


def test_static_gate_rejects_missing_and_nonfinite_results():
    path = FIXTURES / "rs_relaxed.vasp"
    atoms = read(path)
    with pytest.raises(KeyError):
        raw_stationarity(atoms, {"status":"completed"})
    summary = {"status":"completed", "potential_energy_eV":float("nan"),
               "forces_eV_per_A":np.zeros((8,3)).tolist(), "stress_eV_per_A3_voigt":np.zeros(6).tolist()}
    with pytest.raises(ValueError, match="finite"):
        raw_stationarity(atoms, summary)


def test_actual_relaxed_endpoint_origin_noise_cannot_change_periodic_winding():
    rs = phase_candidate(read(FIXTURES / "rs_relaxed.vasp"), "rs")[0]
    wz = phase_candidate(read(FIXTURES / "wz_relaxed.vasp"), "wz")[0]
    chain = mapped_chain(rs, wz, "atomic_mapping")
    minimum = min(np.min(a.get_all_distances(mic=True) + np.eye(8)*100) for a in chain)
    assert minimum > 2.5  # real failure with nearest-image reselection: 1.674 A
    origin_shifted = wz.copy()
    origin_shifted.positions += np.array([1e-12,-1e-12,0.])
    other = mapped_chain(rs, origin_shifted, "atomic_mapping")
    assert all(np.allclose(a.positions,b.positions,atol=2e-12,rtol=0) for a,b in zip(chain,other))
