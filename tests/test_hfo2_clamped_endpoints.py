"""Common-substrate preparation and endpoint physical gates use no DFT."""

import json
from pathlib import Path
import shutil

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read

from scripts.prepare_hfo2_clamped_endpoints import (AXES, PHASE_FILES, geometry_guard,
                                                   orient_long_axis_x, prepare)
from scripts.prepare_hfo2_reference_variants import write_clean_poscar
from scripts.audit_hfo2_static_replica import sha256
from scripts.relax_clamped_ase_endpoint import (ClampedEndpointBFGS, load_seed, relax, stress_report)
from vcneb import ClampedPlaneFilter, clamped_plane_vcneb_boundary


CASE = Path(__file__).resolve().parents[1] / "benchmarks/hfo2_channels/20261008"


def test_proper_rotation_preserves_geometry_identity_and_drops_calculator():
    atoms = read(CASE / "reference_variants/PO.vasp")
    atoms.calc = SinglePointCalculator(atoms, energy=123., forces=np.zeros((12, 3)), stress=np.zeros(6))
    before = atoms.positions.copy()
    rotated = orient_long_axis_x(atoms)
    assert rotated.calc is None
    assert np.linalg.det(np.eye(3)[list(AXES)]) == 1
    assert rotated.get_chemical_symbols() == atoms.get_chemical_symbols()
    np.testing.assert_array_equal(atoms.positions, before)
    np.testing.assert_allclose(rotated.get_all_distances(mic=True), atoms.get_all_distances(mic=True), atol=1e-12)
    assert rotated.get_volume() == pytest.approx(atoms.get_volume())
    np.testing.assert_allclose(rotated.get_scaled_positions(wrap=False),
                               atoms.get_scaled_positions(wrap=False)[:, AXES], atol=1e-12)


def test_actual_ten_seeds_have_one_substrate_per_condition_and_no_energy(tmp_path):
    output = tmp_path / "seeds"
    report = prepare(CASE, output)
    assert report["n_geometry_seeds"] == len(report["seeds"]) == 10
    assert report["new_DFT_calls"] == 0
    assert report["strain_conditions"] == [0., .01]
    assert report["holdout_strain_reserved_not_generated"] == .005
    assert report["orientation"]["atom_permutation"] == list(range(12))
    planes = {}
    for record in report["seeds"]:
        path = output / record["manifest_file"]
        assert sha256(path) == record["manifest_sha256"]
        atoms, boundary, manifest = load_seed(path)
        assert atoms.calc is None
        assert "energy_eV" not in manifest
        assert len(atoms) == 12 and record["minimum_distance_A"] > 1.6
        boundary.validate_images([atoms])
        strain = record["strain"]
        planes.setdefault(strain, atoms.cell.array[:2])
        np.testing.assert_allclose(atoms.cell[:2], planes[strain], atol=1e-14)
        free = orient_long_axis_x(read(CASE / PHASE_FILES[record["phase_label"]]))
        np.testing.assert_allclose(atoms.get_scaled_positions(wrap=False),
                                   free.get_scaled_positions(wrap=False), atol=1e-12)
        if strain == 0 and record["phase_label"] == "PO_plus":
            assert max(np.abs(manifest["free_phase_in_plane_length_changes"])) > .03
    np.testing.assert_allclose(planes[.01], planes[0.] * 1.01, atol=1e-14)
    with pytest.raises(FileExistsError, match="refusing existing"):
        prepare(CASE, output)


def test_changed_registered_source_refused_before_output_created(tmp_path):
    source = tmp_path / "case"
    (source / "reference_variants").mkdir(parents=True)
    for relative in (*PHASE_FILES.values(), "reference_variants/variant_manifest.json", "M_endpoint_audit.json"):
        shutil.copyfile(CASE / relative, source / relative)
    target = source / "reference_variants/PO.vasp"
    target.write_text(target.read_text() + "\n")
    with pytest.raises(ValueError, match="checksum changed"):
        prepare(source, tmp_path / "output")
    assert not (tmp_path / "output").exists()


@pytest.mark.parametrize("allow_tilt", [False, True])
@pytest.mark.parametrize("atomic_force,normal,shear,expected", [
    (0., 0., 0., True), (.04, 0., 0., False), (0., 2.1, 0., False), (0., 0., 2.1, None),
])
def test_bfgs_convergence_excludes_reactions_and_enforces_released_traction(
    allow_tilt, atomic_force, normal, shear, expected,
):
    atoms = Atoms("Ar", positions=[[1, 1, 1]], cell=np.diag([4., 4., 5.]), pbc=True)
    stress = np.array([[100., 10., shear], [10., -200., 0.], [shear, 0., normal]]) / 1602.176634
    atoms.calc = SinglePointCalculator(atoms, energy=0., forces=[[atomic_force, 0., 0.]], stress=stress)
    boundary = clamped_plane_vcneb_boundary(1, atoms.cell.array, allow_tilt=allow_tilt)
    target = ClampedPlaneFilter(atoms, boundary, cell_scale_A=5.)
    optimizer = ClampedEndpointBFGS(target, boundary=boundary, endpoint_atoms=atoms,
                                    pressure_gpa=0., stress_kbar=2., logfile=None)
    optimizer.fmax = .03
    if expected is None:
        expected = not allow_tilt
    assert optimizer.converged(target.get_forces()) == expected
    assert optimizer.gradient_converged(-target.get_forces().ravel()) == expected
    assert stress_report(atoms, boundary)["max_abs_full_stress_kbar_information_only"] == pytest.approx(200.)


def test_stress_pressure_sign_and_normal_only_shear_exclusion():
    atoms = Atoms("Ar", cell=np.diag([4., 4., 5.]), pbc=True)
    atoms.calc = SinglePointCalculator(atoms, energy=0., forces=np.zeros((1, 3)),
                                       stress=-np.eye(3) / 1602.176634 * 10)
    boundary = clamped_plane_vcneb_boundary(1, atoms.cell.array, allow_tilt=True)
    assert stress_report(atoms, boundary, pressure_gpa=1.)["open_traction_norm_kbar"] == pytest.approx(0.)


def test_complete_endpoint_bfgs_checkpoint_no_duplicate_evaluations(tmp_path):
    atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
    cell = atoms.cell.array.copy()
    cell[:2] *= 1.01
    cell[2, 2] *= 1.03
    atoms.set_cell(cell, scale_atoms=True)
    atoms.positions[0, 2] += .02
    boundary = clamped_plane_vcneb_boundary(len(atoms), cell, allow_tilt=True)
    manifest = {"minimum_distance_guard_A": 1.6, "cell_scale_A": 3.6, "pressure_gpa": 0.,
                "phase_label": "test_Cu", "strain": .01}
    calls = []

    class CountingEMT(EMT):
        def calculate(self, *args, **kwargs):
            calls.append(1)
            super().calculate(*args, **kwargs)

    def factory(image_index, image, directory):
        return CountingEMT(directory=str(directory))

    workdir = tmp_path / "endpoint"
    report = relax(atoms, boundary, manifest, workdir, factory=factory, steps=60, maxstep=.03)
    assert report["converged"] and report["status"] == "completed"
    assert report["open_traction_norm_kbar"] < 2 and report["max_atomic_force_eV_per_A"] < .03
    assert report["max_abs_full_stress_kbar_information_only"] > 2  # genuine clamping reaction
    assert report["substrate_plane_drift_A"] < 1e-13
    assert report["phase_and_variant_gate"] == "pending_independent_audit_no_symmetry_forcing"
    assert len(calls) == report["optimizer_steps"] + 1
    for image in read(workdir / "relax.traj", index=":"):
        boundary.validate_images([image])
    persisted = json.loads((workdir / "endpoint_relax_summary.json").read_text())
    assert persisted["endpoint"]["sha256"] == report["endpoint"]["sha256"]
    with pytest.raises(FileExistsError):
        relax(atoms, boundary, manifest, workdir, factory=factory)


def test_seed_checksum_and_namespace_escape_guards(tmp_path):
    atoms = bulk("Cu", "fcc", a=3.6, cubic=True)
    path = tmp_path / "POSCAR.seed"
    write_clean_poscar(path, atoms)
    manifest = {"format_version": 1, "status": "geometry_seed_not_relaxed", "seed_file": path.name,
                "seed_sha256": sha256(path), "ordered_species": atoms.get_chemical_symbols(),
                "reference_cell_A": atoms.cell.array.tolist(), "allow_tilt": True,
                "minimum_distance_guard_A": 1.6, "cell_scale_A": 3.6, "pressure_gpa": 0.}
    record = tmp_path / "endpoint_seed.json"
    record.write_text(json.dumps(manifest))
    load_seed(record)
    path.write_text(path.read_text() + "\n")
    with pytest.raises(ValueError, match="checksum changed"):
        load_seed(record)
    manifest["seed_file"] = "../another/POSCAR.seed"
    record.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="sibling"):
        load_seed(record)


def test_geometry_guard_rejects_overlap():
    atoms = Atoms("Ar2", positions=[[1, 1, 1], [1.1, 1, 1]], cell=[5, 5, 5], pbc=True)
    with pytest.raises(ValueError, match="unsafe minimum"):
        geometry_guard([atoms])
