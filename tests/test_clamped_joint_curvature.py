"""Joint curvature probes use the same exact substrate as the path, without DFT."""

import json
from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms

from scripts.relax_clamped_ase_endpoint import load_seed
from scripts.verify_clamped_joint_coordinates import main as verify_main, verify
from vcneb import ActiveJointCurvatureCoordinates as JointCurvatureCoordinates, clamped_plane_vcneb_boundary
from vcneb.joint_curvature import JointCurvatureCoordinates as LegacyCoordinates, symmetric_strain_basis


def _reference(rotated=False):
    atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
    cell = atoms.cell.array.copy()
    cell[0] += [0.0, 0.10, 0.15]
    cell[1] += [0.05, 0.0, 0.08]
    cell[2] += [0.11, -0.07, 0.0]
    if rotated:
        rotation, _ = np.linalg.qr(np.random.default_rng(81).normal(size=(3, 3)))
        assert np.linalg.det(rotation) > 0
        cell = cell @ rotation.T
    atoms.set_cell(cell, scale_atoms=True)
    return atoms


@pytest.mark.parametrize("allow_tilt", [False, True])
@pytest.mark.parametrize("deformed", [False, True])
@pytest.mark.parametrize("rotated", [False, True])
def test_joint_clamped_gradient_matches_energy_work(allow_tilt, deformed, rotated):
    reference = _reference(rotated)
    boundary = clamped_plane_vcneb_boundary(len(reference), reference.cell.array,
                                            allow_tilt=allow_tilt)
    chart = JointCurvatureCoordinates.for_clamped_plane(reference, cell_scale_A=3.6,
                                                       boundary=boundary)
    assert chart.cell_dofs == boundary.cell_dofs
    assert chart.size == 3 * len(reference) + boundary.cell_dofs
    center = np.zeros(chart.size)
    if deformed:
        center = np.random.default_rng(314).normal(size=chart.size) * 0.015
    at_center = chart.displaced(center)
    at_center.calc = EMT()
    pressure = 0.03  # eV/A^3, compressive-positive, analytic verification only
    gradient = chart.enthalpy_gradient(at_center, at_center.get_forces(),
                                      at_center.get_stress(voigt=False), pressure)
    h = 1e-5
    numerical = []
    for i in range(chart.size):
        pair = []
        for sign in (1, -1):
            delta = center.copy()
            delta[i] += sign * h
            probe = chart.displaced(delta)
            boundary.validate_images([probe])
            np.testing.assert_allclose(probe.cell[:2], reference.cell[:2], atol=1e-13)
            probe.calc = EMT()
            pair.append(probe.get_potential_energy() + pressure * probe.get_volume())
        numerical.append((pair[0] - pair[1]) / (2 * h))
    np.testing.assert_allclose(gradient, numerical, rtol=0, atol=3e-5)
    free = chart.translation_free_basis()
    np.testing.assert_allclose(free.T @ free, np.eye(chart.size - 3), atol=1e-12)
    translations = np.zeros((chart.size, 3))
    translations[:3 * len(reference)] = np.tile(np.eye(3), (len(reference), 1))
    np.testing.assert_allclose(translations.T @ free, 0, atol=1e-12)


def test_tilt_basis_is_not_replaced_with_symmetric_strain():
    reference = _reference(True)
    boundary = clamped_plane_vcneb_boundary(len(reference), reference.cell.array,
                                            allow_tilt=True)
    chart = JointCurvatureCoordinates.for_clamped_plane(reference, cell_scale_A=3.6,
                                                       boundary=boundary)
    assert max(np.linalg.norm(b - b.T) for b in chart.deformation_basis) > 0.5
    changes = np.array([reference.cell.array @ b.T for b in chart.deformation_basis])
    np.testing.assert_allclose(changes[:, :2], 0, atol=1e-14)
    assert np.linalg.matrix_rank(changes[:, 2]) == 3
    assert not chart.deformation_basis.flags.writeable
    default = JointCurvatureCoordinates(reference, 3.6)
    np.testing.assert_array_equal(default.deformation_basis, symmetric_strain_basis())
    assert default.size == 3 * len(reference) + 6


def test_default_extension_matches_legacy_geometry_and_gradients_exactly():
    reference = _reference(True)
    old = LegacyCoordinates(reference, 3.6)
    new = JointCurvatureCoordinates(reference, 3.6)
    for seed in range(4):
        delta = np.random.default_rng(seed).normal(size=new.size) * 0.01
        prior = old.displaced(delta)
        current = new.displaced(delta)
        np.testing.assert_array_equal(current.cell.array, prior.cell.array)
        np.testing.assert_array_equal(current.positions, prior.positions)
        current.calc = EMT()
        forces, stress = current.get_forces(), current.get_stress(voigt=False)
        np.testing.assert_array_equal(new.enthalpy_gradient(current, forces, stress, 0.03),
                                      old.enthalpy_gradient(prior, forces, stress, 0.03))
        np.testing.assert_array_equal(new.translation_free_basis(), old.translation_free_basis())


def test_incompatible_reference_and_evaluated_probe_are_rejected_unmodified():
    reference = _reference()
    boundary = clamped_plane_vcneb_boundary(len(reference), reference.cell.array,
                                            allow_tilt=True)
    wrong = reference.copy()
    wrong.cell[0, 0] += 0.001
    before = wrong.cell.array.copy()
    with pytest.raises(ValueError, match="substrate"):
        JointCurvatureCoordinates.for_clamped_plane(wrong, cell_scale_A=3.6,
                                                   boundary=boundary)
    chart = JointCurvatureCoordinates.for_clamped_plane(reference, cell_scale_A=3.6,
                                                       boundary=boundary)
    with pytest.raises(ValueError, match="substrate"):
        chart.enthalpy_gradient(wrong, np.zeros((len(wrong), 3)), np.zeros((3, 3)), 0.)
    np.testing.assert_array_equal(wrong.cell.array, before)


@pytest.mark.parametrize("basis", [np.ones((2, 3, 3)), np.zeros((1, 9)),
                                  np.full((1, 3, 3), np.nan)])
def test_invalid_deformation_basis_is_not_silently_normalized(basis):
    with pytest.raises(ValueError, match="deformation basis"):
        JointCurvatureCoordinates(_reference(), 3.6, deformation_basis=basis)


def test_explicit_fixed_cell_degenerates_to_atomic_chart_and_copies_basis():
    reference = _reference()
    chart = JointCurvatureCoordinates(reference, 3.6, deformation_basis=np.empty((0, 3, 3)))
    assert chart.cell_dofs == 0 and chart.size == 3 * len(reference)
    changed = chart.displaced(np.ones(chart.size) * 0.001)
    np.testing.assert_array_equal(changed.cell.array, reference.cell.array)
    basis = symmetric_strain_basis()
    generic = JointCurvatureCoordinates(reference, 3.6, deformation_basis=basis)
    basis[:] = 0
    np.testing.assert_array_equal(generic.deformation_basis, symmetric_strain_basis())


def test_extra_atomic_constraints_cannot_silently_change_active_variables():
    reference = _reference()
    reference.set_constraint(FixAtoms(indices=[0]))
    with pytest.raises(ValueError, match="explicit atomic subspace"):
        JointCurvatureCoordinates(reference, 3.6)


def test_actual_ten_hafnia_starters_have_same_path_and_probe_subspace():
    root = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds"
    manifest = json.loads((root / "clamped_seed_manifest.json").read_text())
    assert len(manifest["seeds"]) == 10
    for item in manifest["seeds"]:
        atoms, boundary, record = load_seed(root / item["manifest_file"])
        chart = JointCurvatureCoordinates.for_clamped_plane(
            atoms, cell_scale_A=record["cell_scale_A"], boundary=boundary,
        )
        assert chart.size == 39 and chart.translation_free_basis().shape == (39, 36)
        kwargs = boundary.vcneb_kwargs([atoms])
        cell_block = kwargs["mode_basis"][36:, 36:]
        np.testing.assert_allclose(chart.deformation_basis.reshape(3, 9).T, cell_block)
        before = atoms.positions.copy()
        for seed in range(3):
            delta = np.random.default_rng(seed).normal(size=chart.size) * 0.005
            probe = chart.displaced(delta)
            boundary.validate_images([probe])
            assert probe.calc is None
            assert probe.get_chemical_symbols() == atoms.get_chemical_symbols()
        np.testing.assert_array_equal(atoms.positions, before)


def test_verification_receipt_distinguishes_geometry_from_material_predictions():
    root = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds"
    report = verify(root)
    assert report["geometry_starters_checked"] == 10
    assert len(report["EMT_work"]) == 8 and report["EMT_evaluations"] == 232
    assert report["new_DFT_calls"] == 0
    assert not report["saddle_certified"]
    assert not report["material_joint_Hessian_measured"]
    assert not report["material_barrier_prediction_validated"]
    assert all(r["max_work_gradient_error_eV_A"] < report["work_tolerance_eV_A"]
               for r in report["EMT_work"])
    assert all(r["status"] == "geometry_seed_not_relaxed" for r in report["geometry"])
    assert all(len(value) == 64 for value in report["module_sha256"].values())
    json.dumps(report)


def test_verifier_refuses_existing_receipt_before_any_calculation(tmp_path, monkeypatch):
    report = tmp_path / "existing.json"
    report.write_text("original receipt")
    monkeypatch.setattr("sys.argv", ["verify", "--report", str(report)])
    with pytest.raises(FileExistsError, match="refusing existing"):
        verify_main()
    assert report.read_text() == "original receipt"
