"""No-DFT extrema, unit normalization and evidence-gate checks."""
from pathlib import Path
import runpy

import numpy as np
import pytest
from ase import Atoms, units
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write
from ase.io.trajectory import Trajectory
from vcneb import endpoint_structure_record

ROOT = Path(__file__).resolve().parents[1]
MODULE = runpy.run_path(str(ROOT / "scripts" / "analyze_gan_qian_path.py"))
SUITE = runpy.run_path(str(ROOT / "scripts" / "prepare_gan_qian_suite.py"))


def test_three_peaks_are_measured_from_each_local_start_not_global_zero():
    result = MODULE["peak_segments"]([0, .57, -.01, .56, 0, .57, -.01])
    assert result["peak_indices"] == [1, 3, 5]
    assert result["valley_indices"] == [2, 4]
    assert [segment["forward_barrier_eV_per_GaN"] for segment in result["segments"]] == pytest.approx([.57, .57, .57])


def test_small_endpoint_wiggle_is_reported_only_when_requested():
    values = [0, .2, .34, .1, -.005, -.004, -.005]
    assert MODULE["peak_segments"](values)["peak_indices"] == [2]
    assert MODULE["peak_segments"](values, .0001)["peak_indices"] == [2, 5]


def test_analysis_reader_preserves_saved_energy_and_rejects_partial_tail(tmp_path):
    atom = Atoms("Ba", cell=[4, 4, 4], pbc=True)
    atom.calc = SinglePointCalculator(atom, energy=2.)
    path = tmp_path / "chain.traj"
    write(path, [atom] * 6)
    assert MODULE["latest_evaluated_chain"](path, 3)[0].get_potential_energy() == 2.
    with Trajectory(str(path), "a") as trajectory:
        trajectory.write(atom)
    with pytest.raises(ValueError, match="partial"):
        MODULE["latest_evaluated_chain"](path, 3)


def synthetic_case():
    images = SUITE["body_diagonal_chain"](SUITE["conventional_gan"](4.3, .25), SUITE["conventional_gan"](4.05, .5))
    pressure = 45 * units.GPa
    values = np.interp(np.arange(29), [0, 5, 10, 15, 20, 25, 28], [0, .57, -.01, .56, 0, .57, -.01])
    enthalpies = values * 4
    for atoms, enthalpy in zip(images, enthalpies):
        atoms.calc = SinglePointCalculator(atoms, energy=enthalpy - pressure * atoms.get_volume())
    summary = {"status": "completed", "converged": True, "final_max_generalized_force_eV_per_A": .09,
        "calculator_parameters": {"encut": 600., "ediff": 1e-7, "isym": -1, "symprec": 1e-4,
            "nsw": 0, "ibrion": -1, "isif": 2, "xc": "PBE", "pp": "PBE", "gamma": True,
            "kpts": [8, 8, 6], "setups": {"Ga": "_d"}},
        "licensed_input_fingerprints": {"POTCAR": {"sha256": MODULE["POTCAR_SHA256"]}},
        "runtime_parameter_changes_allowed": False, "endpoint_evaluation_policy": "fixed_cached_once",
        "n_images": 29, "image_enthalpies_eV": enthalpies.tolist(),
        "path_diagnostics": {"pressure_eV_per_A3": pressure,
            "images": [{"is_climbing_image": False}] * 29},
        "endpoint_structures": {"initial": endpoint_structure_record(images[0]), "final": endpoint_structure_record(images[-1])}}
    return summary, images, {"status": "ok", "issues": []}, [{"status": "ok", "image_indices": list(range(1, 28))}]


def test_full_synthetic_b3_analysis_normalizes_four_formula_units():
    summary, images, audit, manifest = synthetic_case()
    result = MODULE["analyze"](summary, images, "b3_diagonal", audit, manifest)
    assert result["n_formula_units"] == 4
    assert result["barrier_eV_per_GaN"] == pytest.approx(.57)
    assert result["resolved_extrema"]["peak_indices"] == [5, 15, 25]
    assert [segment["forward_barrier_eV_per_GaN"] for segment in result["resolved_extrema"]["segments"]] == pytest.approx([.57] * 3)


@pytest.mark.parametrize("failure", ["pressure", "paw", "endpoint_worker", "ci", "unconverged", "snapshot_energy"])
def test_analysis_rejects_insufficient_or_changed_evidence(failure):
    summary, images, audit, manifest = synthetic_case()
    if failure == "pressure":
        summary["path_diagnostics"]["pressure_eV_per_A3"] = 45.7 * units.GPa
    elif failure == "paw":
        summary["licensed_input_fingerprints"]["POTCAR"]["sha256"] = "wrong"
    elif failure == "endpoint_worker":
        manifest[0]["image_indices"] = list(range(29))
    elif failure == "ci":
        summary["path_diagnostics"]["images"][15] = {"is_climbing_image": True}
    elif failure == "unconverged":
        summary["converged"] = False
    else:
        summary["image_enthalpies_eV"][15] += .01
    with pytest.raises(ValueError):
        MODULE["analyze"](summary, images, "b3_diagonal", audit, manifest)
