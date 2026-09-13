"""Finite-difference checks for VC-NEB force transforms."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import types

import numpy as np
from ase import Atoms
from ase.calculators.calculator import Calculator, all_changes
from ase.io import write
from ase.io.trajectory import Trajectory

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.run_toy_vcneb import ToyPhaseTransition
from scripts.validate_vcneb_inputs import validate_abacus
from vcneb import (
    CalculatorCapabilityError,
    apply_chain_state,
    inspect_calculator,
    interpolate_vcneb,
    path_geometry_diagnostics,
    read_chain_trajectory,
    run_vcneb,
    validate_path_geometry,
    validate_image_calculators,
)
from vcneb import (
    Mode,
    build_direction_basis,
    build_mode_basis,
    direction_basis_conflicts,
    infer_atom_mapping,
    mode_guided_path,
    project_path_onto_modes,
    validate_atom_mapping,
)
from vcneb.abacus import make_ase_abacus_factory
from vcneb.core import VCNEB, cell_force, cell_from_deformation, deformation_from_cell, fractional_force
from vcneb.executor import ThreadedCalculatorExecutor
import vcneb.core as core_module


def make_atoms(reference_cell, q, deform):
    atoms = Atoms("Ar", scaled_positions=[q], cell=cell_from_deformation(deform, reference_cell), pbc=True)
    atoms.calc = ToyPhaseTransition(reference_cell)
    return atoms


class MetricCellCalculator(Calculator):
    """Rotation-invariant elastic model for non-diagonal cell tests."""

    implemented_properties = ["energy", "forces", "stress"]

    def __init__(self, reference_cell, *, stiffness=0.8):
        super().__init__()
        self.reference_cell = np.asarray(reference_cell, dtype=float)
        self.stiffness = float(stiffness)

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        deform = deformation_from_cell(atoms.cell.array, self.reference_cell)
        metric = deform @ deform.T
        delta = metric - np.eye(3)
        energy = 0.5 * self.stiffness * float(np.sum(delta * delta))
        grad_deform = 2.0 * self.stiffness * delta @ deform
        generalized_force = -grad_deform
        virial = generalized_force @ deform.T
        self.results["energy"] = energy
        self.results["forces"] = np.zeros((len(atoms), 3), dtype=float)
        self.results["stress"] = -virial / atoms.get_volume()


class LinearCoordinateCalculator(Calculator):
    """Energy linear in one fractional coordinate, with no cell force."""

    implemented_properties = ["energy", "forces", "stress"]

    def __init__(self, slope=0.7):
        super().__init__()
        self.slope = float(slope)

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        atoms = self.atoms
        q = atoms.get_scaled_positions(wrap=False)
        gradient = np.zeros_like(q)
        gradient[0, 0] = self.slope
        self.results["energy"] = float(self.slope * q[0, 0])
        self.results["forces"] = -gradient @ np.linalg.inv(atoms.cell.array.T)
        self.results["stress"] = np.zeros((3, 3), dtype=float)


class ZeroCalculator(Calculator):
    implemented_properties = ["energy", "forces", "stress"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        self.results["energy"] = 0.0
        self.results["forces"] = np.zeros((len(self.atoms), 3), dtype=float)
        self.results["stress"] = np.zeros((3, 3), dtype=float)


def check_cell_mask_regression() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = make_atoms(reference_cell, np.array([0.25, 0.5, 0.5]), np.eye(3))
    final_deform = np.eye(3)
    final_deform[0, 0] = 1.25
    final_deform[1, 1] = 1.30
    final = make_atoms(reference_cell, np.array([0.75, 0.5, 0.5]), final_deform)
    images = interpolate_vcneb(initial, final, n_images=3, align_cells=False)
    for image in images:
        image.calc = ToyPhaseTransition(reference_cell)

    mask = np.zeros((3, 3))
    mask[0, 0] = 1.0
    chain = VCNEB(images, k=1.0, climb=False, cell_mask=mask)

    before = deformation_from_cell(images[1].cell.array, reference_cell)
    x = chain.get_x()
    trial = x.copy()
    cell_offset = 3 * len(images[1])
    trial[cell_offset + 4] += 3.0 * chain.cell_scale
    chain.set_x(trial)
    after = deformation_from_cell(images[1].cell.array, reference_cell)
    inactive = np.ones((3, 3), dtype=bool)
    inactive[0, 0] = False
    if not np.allclose(after[inactive], before[inactive]):
        raise SystemExit("cell_mask failed to preserve inactive cell components")

    force_x = chain._compute_forces()
    cell_force_x = force_x[cell_offset : cell_offset + 9].reshape(3, 3)
    if np.max(np.abs(cell_force_x[inactive])) > 1e-12:
        raise SystemExit("cell_mask failed to remove inactive VC-NEB cell forces")


def check_atom_mask_regression() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = make_atoms(reference_cell, np.array([0.25, 0.5, 0.5]), np.eye(3))
    final = make_atoms(reference_cell, np.array([0.75, 0.8, 0.5]), np.eye(3))
    images = interpolate_vcneb(initial, final, n_images=3, align_cells=False)
    for image in images:
        image.calc = ToyPhaseTransition(reference_cell)

    mask = np.array([[1.0, 0.0, 1.0]])
    chain = VCNEB(images, k=1.0, climb=False, atom_mask=mask)
    before = chain.images[1].get_positions().copy()
    trial = chain.get_x()
    trial[1] += 2.0
    chain.set_x(trial)
    after = chain.images[1].get_positions()
    if not np.allclose(after[0, 1], before[0, 1]):
        raise SystemExit("atom_mask failed to preserve an inactive atomic component")
    force_x = chain._compute_forces()
    if abs(force_x[1]) > 1e-12:
        raise SystemExit("atom_mask failed to remove inactive atomic force")


def check_mode_guided_path() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = Atoms("Ar", scaled_positions=[[0.25, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final = Atoms("Ar", scaled_positions=[[0.75, 0.5, 0.5]], cell=reference_cell, pbc=True)
    mode = Mode([[0.0, 1.0, 0.0]])
    images = mode_guided_path(initial, final, 3, mode, amplitude=0.4, align_cells=False)
    if not np.allclose(images[0].positions, initial.positions) or not np.allclose(images[-1].positions, final.positions):
        raise SystemExit("mode_guided_path changed an endpoint")
    if abs(images[1].positions[0, 1] - 2.9) > 1e-12:
        raise SystemExit("mode_guided_path did not apply the normalized mode amplitude")
    projection = project_path_onto_modes(images, initial, mode)
    if abs(projection[1, 0] - 0.4) > 1e-12:
        raise SystemExit("project_path_onto_modes returned the wrong modal amplitude")
    cell_mode = Mode(np.zeros((1, 3)), cell=np.diag([0.2, 0.0, 0.0]))
    cell_images = mode_guided_path(
        initial,
        final,
        3,
        cell_mode,
        amplitude=0.5,
        normalize=False,
        align_cells=False,
    )
    if not np.allclose(cell_images[0].cell.array, initial.cell.array) or not np.allclose(
        cell_images[-1].cell.array, final.cell.array
    ):
        raise SystemExit("cell mode changed a mode-guided endpoint")
    cell_midpoint = deformation_from_cell(cell_images[1].cell.array, reference_cell)
    expected_cell_midpoint = np.eye(3) + 0.5 * cell_mode.cell / 5.0
    if not np.allclose(cell_midpoint, expected_cell_midpoint, rtol=0.0, atol=1e-12):
        raise SystemExit("mode_guided_path did not apply the cell mode component")


def check_cell_interpolation_strategies() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = make_atoms(reference_cell, np.array([0.25, 0.5, 0.5]), np.eye(3))
    final_deform = np.diag([1.6, 0.8, 1.2])
    final = make_atoms(reference_cell, np.array([0.75, 0.5, 0.5]), final_deform)

    linear = interpolate_vcneb(
        initial,
        final,
        n_images=3,
        align_cells=False,
        cell_interpolation="linear",
    )
    logarithmic = interpolate_vcneb(
        initial,
        final,
        n_images=3,
        align_cells=False,
        cell_interpolation="log_strain",
    )
    if not np.allclose(linear[-1].cell.array, final.cell.array) or not np.allclose(
        logarithmic[-1].cell.array, final.cell.array
    ):
        raise SystemExit("cell interpolation changed an endpoint")
    expected_log_midpoint = np.diag(np.sqrt(np.diag(final_deform)))
    log_midpoint = deformation_from_cell(logarithmic[1].cell.array, reference_cell)
    if not np.allclose(log_midpoint, expected_log_midpoint, rtol=0.0, atol=1e-12):
        raise SystemExit("log_strain interpolation returned the wrong geometric midpoint")
    if np.allclose(linear[1].cell.array, logarithmic[1].cell.array):
        raise SystemExit("linear and log_strain interpolation unexpectedly matched")

    callback_calls = []

    def custom_interpolator(parameter, deform0, deform1):
        callback_calls.append(parameter)
        return (1.0 - parameter) * deform0 + parameter * deform1

    custom = interpolate_vcneb(
        initial,
        final,
        n_images=4,
        align_cells=False,
        cell_interpolation=custom_interpolator,
    )
    if len(callback_calls) != 2 or not np.allclose(
        deformation_from_cell(custom[2].cell.array, reference_cell), 2.0 * final_deform / 3.0 + np.eye(3) / 3.0
    ):
        raise SystemExit("custom cell interpolation callback was not applied correctly")

    rotated = make_atoms(
        reference_cell,
        np.array([0.75, 0.5, 0.5]),
        np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]),
    )
    try:
        interpolate_vcneb(
            initial,
            rotated,
            n_images=3,
            align_cells=False,
            cell_interpolation="log_strain",
        )
    except ValueError:
        pass
    else:
        raise SystemExit("log_strain accepted an unaligned rigid rotation")


def check_atom_mapping() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = Atoms(
        "HHe",
        scaled_positions=[[0.10, 0.20, 0.30], [0.40, 0.50, 0.60]],
        cell=reference_cell,
        pbc=True,
    )
    final = Atoms(
        "HeH",
        scaled_positions=[[0.41, 0.50, 0.60], [0.10, 0.20, 0.30]],
        cell=reference_cell,
        pbc=True,
    )
    inferred = infer_atom_mapping(initial, final, mic=True)
    if inferred != [1, 0]:
        raise SystemExit(f"auto atom mapping returned {inferred}, expected [1, 0]")
    report = validate_atom_mapping(initial, final, "auto", mic=True)
    if report["mapping"] != [1, 0] or report["maximum_displacement_A"] > 0.1:
        raise SystemExit("atom mapping report did not preserve the expected permutation")
    images = interpolate_vcneb(
        initial,
        final,
        n_images=3,
        mapping="auto",
        align_cells=False,
        mic=True,
    )
    if images[-1].get_chemical_symbols() != initial.get_chemical_symbols():
        raise SystemExit("interpolate_vcneb did not apply the inferred atom permutation")
    try:
        validate_atom_mapping(initial, final, [0, 1])
    except ValueError:
        pass
    else:
        raise SystemExit("atom mapping accepted pairs with different elements")


def check_periodic_translation_alignment() -> None:
    cell = np.diag([4.0, 4.0, 4.0])
    initial = Atoms(
        "BaTiO3",
        scaled_positions=[
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
            [0.5, 0.5, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
        ],
        cell=cell,
        pbc=True,
    )
    final = Atoms(
        "BaTiO3",
        scaled_positions=[
            [0.5, 0.5, 0.06],
            [0.0, 0.0, 0.12],
            [0.0, 0.0, 0.54],
            [0.5, 0.0, 0.04],
            [0.0, 0.5, 0.04],
        ],
        cell=cell,
        pbc=True,
    )
    unaligned = interpolate_vcneb(
        initial,
        final,
        n_images=7,
        align_cells=False,
        mic=True,
        mapping="auto",
    )
    aligned = interpolate_vcneb(
        initial,
        final,
        n_images=7,
        align_cells=False,
        mic=True,
        mapping="auto",
        align_translation=True,
    )
    unaligned_minimum = min(
        float(np.min(image.get_all_distances(mic=True) + np.eye(len(image)) * 1e6))
        for image in unaligned
    )
    aligned_minimum = min(
        float(np.min(image.get_all_distances(mic=True) + np.eye(len(image)) * 1e6))
        for image in aligned
    )
    if unaligned_minimum > 1.1 or aligned_minimum < 1.6:
        raise SystemExit(
            "periodic translation alignment did not remove the endpoint-origin collision: "
            f"{unaligned_minimum}, {aligned_minimum}"
        )
    if not np.allclose(aligned[0].positions, initial.positions):
        raise SystemExit("periodic translation alignment changed the initial endpoint")
    metadata = aligned[0].info.get("vcneb_path_metadata", {})
    mapping = metadata.get("mapping", [])
    translation = np.asarray(metadata.get("translation_fractional_final_cell", []), dtype=float)
    if sorted(mapping) != list(range(len(initial))) or not metadata.get("align_translation"):
        raise SystemExit(f"periodic path metadata lost the chosen gauge: {metadata}")
    if translation.shape != (3,) or np.linalg.norm(translation) < 1e-8:
        raise SystemExit(f"periodic path metadata lost the chosen translation: {metadata}")
    print("periodic_translation_alignment_regression=ok")


def check_tangent_and_spring_components() -> None:
    cell = np.diag([5.0, 5.0, 5.0])
    monotonic = [
        Atoms("Ar", positions=[[5.0 * coordinate, 0.0, 0.0]], cell=cell, pbc=True)
        for coordinate in np.linspace(0.0, 1.0, 5)
    ]
    for image in monotonic:
        image.calc = LinearCoordinateCalculator()
    chain = VCNEB(monotonic, k=0.0, climb=False)
    force = chain._compute_forces()
    if np.max(np.abs(force)) > 1e-12:
        raise SystemExit(f"monotonic tangent did not remove parallel force: {force}")

    nonuniform = [
        Atoms("Ar", positions=[[5.0 * coordinate, 0.0, 0.0]], cell=cell, pbc=True)
        for coordinate in [0.0, 0.2, 0.7, 0.8, 1.0]
    ]
    for image in nonuniform:
        image.calc = ZeroCalculator()
    chain = VCNEB(nonuniform, k=1.0, climb=False)
    force = chain._compute_forces().reshape(3, 12)
    expected_x = np.array([1.5, -2.0, 0.5])
    if not np.allclose(force[:, 0], expected_x, rtol=0.0, atol=1e-12):
        raise SystemExit(f"spring force has the wrong extended-coordinate spacing: {force[:, 0]}")
    if np.max(np.abs(force[:, 1:])) > 1e-12:
        raise SystemExit("spring force leaked outside the path tangent")
    print("tangent_and_spring_regression=ok")


def check_energy_weighted_tangent_regimes() -> None:
    cell = np.diag([5.0, 5.0, 5.0])
    images = [
        Atoms("Ar", positions=[[0.0, 0.0, 0.0]], cell=cell, pbc=True)
        for _ in range(5)
    ]
    for image in images:
        image.calc = ZeroCalculator()
    chain = VCNEB(images, k=0.0, climb=False)
    image_x = []
    for x, y in [(0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (2.5, 1.0), (3.0, 2.0)]:
        vector = np.zeros(chain.image_ndofs)
        vector[:2] = [x, y]
        image_x.append(vector)

    tangent = chain._tangent(2, np.array([0.0, 1.0, 3.0, 2.0, 0.0]), image_x)
    expected = np.array([2.5, 1.0])
    expected /= np.linalg.norm(expected)
    if not np.allclose(tangent[:2], expected, rtol=0.0, atol=1e-12):
        raise SystemExit(f"energy-weighted peak tangent is wrong: {tangent[:2]}")

    tangent = chain._tangent(2, np.array([0.0, 3.0, 1.0, 2.0, 0.0]), image_x)
    if not np.allclose(tangent[:2], expected, rtol=0.0, atol=1e-12):
        raise SystemExit(f"energy-weighted valley tangent is wrong: {tangent[:2]}")

    monotonic = chain._tangent(2, np.array([0.0, 1.0, 2.0, 3.0, 4.0]), image_x)
    expected_monotonic = np.array([0.5, 1.0])
    expected_monotonic /= np.linalg.norm(expected_monotonic)
    if not np.allclose(monotonic[:2], expected_monotonic, rtol=0.0, atol=1e-12):
        raise SystemExit(f"monotonic tangent is wrong: {monotonic[:2]}")
    print("energy_weighted_tangent_regression=ok")


def check_nonorthogonal_mode_projection() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    reference = Atoms("Ar", scaled_positions=[[0.25, 0.5, 0.5]], cell=reference_cell, pbc=True)
    image = reference.copy()
    scaled = image.get_scaled_positions(wrap=False)
    scaled[0, :2] += [0.20, 0.10]
    image.set_scaled_positions(scaled)
    modes = [Mode([[1.0, 0.0, 0.0]]), Mode([[1.0, 1.0, 0.0]])]
    coefficients = project_path_onto_modes([reference, image], reference, modes)
    expected = np.array([0.50, 1.0 / np.sqrt(2.0)])
    if not np.allclose(coefficients[1], expected, rtol=0.0, atol=1e-12):
        raise SystemExit("non-orthogonal mode projection did not return least-squares coefficients")


def check_mask_validation() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = make_atoms(reference_cell, np.array([0.25, 0.5, 0.5]), np.eye(3))
    final = make_atoms(reference_cell, np.array([0.75, 0.5, 0.5]), np.eye(3))
    images = interpolate_vcneb(initial, final, n_images=3, align_cells=False)
    for image in images:
        image.calc = ToyPhaseTransition(reference_cell)
    for keyword, value in [("atom_mask", [[1.0, 0.5, 0.0]]), ("cell_mask", np.full((3, 3), 0.5))]:
        try:
            VCNEB(images, **{keyword: value})
        except ValueError:
            continue
        raise SystemExit(f"{keyword} accepted a non-binary mask")


def check_calculator_contract() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = make_atoms(reference_cell, np.array([0.25, 0.5, 0.5]), np.eye(3))
    final = make_atoms(reference_cell, np.array([0.75, 0.5, 0.5]), np.eye(3))
    images = interpolate_vcneb(initial, final, n_images=3, align_cells=False)
    for image in images:
        image.calc = ToyPhaseTransition(reference_cell)

    report = inspect_calculator(images[1].calc)
    if not report.ok or not report.has_stress or not report.variable_cell:
        raise SystemExit("calculator contract rejected a valid energy/force/stress calculator")
    if "stress" not in report.declared_properties:
        raise SystemExit("calculator contract did not expose declared stress capability")
    reports = validate_image_calculators(images)
    if len(reports) != len(images):
        raise SystemExit("calculator contract returned an incomplete image report")

    class MissingStress(Calculator):
        implemented_properties = ["energy", "forces"]

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            self.results["energy"] = 0.0
            self.results["forces"] = np.zeros((len(self.atoms), 3))

    images[1].calc = MissingStress()
    try:
        validate_image_calculators(images)
    except CalculatorCapabilityError as exc:
        message = str(exc)
        if "image 1" not in message or "stress" not in message:
            raise SystemExit("calculator contract error omitted image and stress context")
    else:
        raise SystemExit("calculator contract accepted a calculator without stress")

    try:
        run_vcneb(images, steps=0, logfile=None, trajectory=None)
    except CalculatorCapabilityError:
        pass
    else:
        raise SystemExit("run_vcneb did not preflight calculator capabilities")


def check_calculator_runtime_diagnostics() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])

    class FailingCalculator(Calculator):
        implemented_properties = ["energy", "forces", "stress"]

        def __init__(self):
            super().__init__(directory="runtime-image")
            self.command = "mpirun -np 4 fake-dft"

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            raise RuntimeError("SCF did not converge")

    initial = Atoms("Ar", scaled_positions=[[0.25, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final = initial.copy()
    final.set_scaled_positions([[0.75, 0.5, 0.5]])
    images = interpolate_vcneb(initial, final, n_images=3, align_cells=False)
    for image in images:
        image.calc = FailingCalculator()
    try:
        VCNEB(images, climb=False).get_forces()
    except RuntimeError as exc:
        message = str(exc)
        for token in ["image 0", "runtime-image", "fake-dft", "SCF did not converge"]:
            if token not in message:
                raise SystemExit(f"runtime calculator diagnostic omitted {token!r}")
    else:
        raise SystemExit("runtime calculator failure was not propagated")


def check_path_diagnostics() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = make_atoms(reference_cell, np.array([0.25, 0.5, 0.5]), np.eye(3))
    final_deform = np.eye(3)
    final_deform[0, 0] = 1.25
    final = make_atoms(reference_cell, np.array([0.75, 0.5, 0.5]), final_deform)
    images = interpolate_vcneb(initial, final, n_images=5, align_cells=False)
    for image in images:
        image.calc = ToyPhaseTransition(reference_cell)
    chain = VCNEB(images, k=0.5, climb=True)
    diagnostics = chain.path_diagnostics()
    if len(diagnostics["images"]) != 5 or diagnostics["highest_image_index"] is None:
        raise SystemExit("path diagnostics returned an incomplete image table")
    geometry = diagnostics.get("geometry", {})
    if not geometry.get("valid") or len(geometry.get("images", [])) != 5:
        raise SystemExit("path diagnostics omitted final-path geometry audit")
    interior = diagnostics["images"][2]
    required = {
        "volume_A3",
        "cell_lengths_A",
        "cell_angles_deg",
        "max_atom_force_eV_per_A",
        "max_stress_eV_per_A3",
        "max_cell_force_eV",
        "true_perpendicular_force_eV_per_A",
        "true_perpendicular_force_max_vector_eV_per_A",
        "spring_force_eV_per_A",
        "spring_force_max_vector_eV_per_A",
        "neb_residual_generalized_force_eV_per_A",
        "neb_residual_force_euclidean_eV_per_A",
    }
    if not required.issubset(interior):
        raise SystemExit("path diagnostics omitted physical or NEB force fields")
    if not np.isfinite(interior["max_stress_eV_per_A3"]):
        raise SystemExit("path diagnostics returned a non-finite stress summary")
    if not diagnostics["has_interior_barrier"] or diagnostics["interior_barrier_indices"] != [2]:
        raise SystemExit("path diagnostics failed to identify the analytic interior barrier")


def check_barrierless_path_diagnostics() -> None:
    """A monotonic band must not be reported as a validated transition state."""

    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = Atoms("Ar", scaled_positions=[[0.25, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final = Atoms("Ar", scaled_positions=[[0.75, 0.5, 0.5]], cell=reference_cell, pbc=True)
    images = interpolate_vcneb(initial, final, n_images=5, align_cells=False)
    for image in images:
        image.calc = LinearCoordinateCalculator()
    chain = VCNEB(images, k=0.2, climb=True)
    diagnostics = chain.path_diagnostics()
    saddle = chain.saddle_diagnostics()
    if diagnostics["has_interior_barrier"] or diagnostics["interior_peak_indices"]:
        raise SystemExit("monotonic path was incorrectly classified as an interior barrier")
    if saddle["is_local_peak"] or saddle["has_interior_barrier"]:
        raise SystemExit("monotonic path was incorrectly classified as a saddle")
    if not saddle["ci_warning"]:
        raise SystemExit("monotonic CI diagnostic omitted its warning")
    print("barrierless_path_diagnostics_regression=ok")


def check_path_geometry_preflight() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = Atoms(
        "Ar2",
        scaled_positions=[[0.10, 0.5, 0.5], [0.16, 0.5, 0.5]],
        cell=reference_cell,
        pbc=True,
    )
    final = initial.copy()
    final.set_scaled_positions([[0.20, 0.5, 0.5], [0.26, 0.5, 0.5]])
    images = interpolate_vcneb(initial, final, n_images=3, align_cells=False)
    report = path_geometry_diagnostics(images, minimum_distance=0.5)
    if report["valid"] or not report["issues"]:
        raise SystemExit("path geometry diagnostics missed an atom collision")
    try:
        interpolate_vcneb(
            initial,
            final,
            n_images=3,
            align_cells=False,
            minimum_distance=0.5,
        )
    except ValueError as exc:
        if "image" not in str(exc) or "distance" not in str(exc):
            raise SystemExit("path geometry preflight omitted image distance context")
    else:
        raise SystemExit("path geometry preflight accepted an overlapping path")
    valid = validate_path_geometry(images, minimum_distance=0.2)
    if not valid["valid"] or len(valid["images"]) != 3:
        raise SystemExit("path geometry preflight rejected a valid path")
    folded = [
        Atoms("Ar", scaled_positions=[[qx, 0.5, 0.5]], cell=reference_cell, pbc=True)
        for qx in (0.25, 0.45, 0.35, 0.60, 0.75)
    ]
    folded_report = path_geometry_diagnostics(folded, fold_cosine_threshold=0.0)
    if not folded_report["folded_junctions"] or 2 not in folded_report["folded_junctions"]:
        raise SystemExit("path geometry diagnostics missed a folded path")
    try:
        validate_path_geometry(folded, fold_cosine_threshold=0.0)
    except ValueError as exc:
        if "path fold" not in str(exc):
            raise SystemExit("fold preflight omitted path-fold context")
    else:
        raise SystemExit("folded path preflight unexpectedly succeeded")


def check_cell_validity_guards() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = make_atoms(reference_cell, np.array([0.25, 0.5, 0.5]), np.eye(3))
    final = make_atoms(reference_cell, np.array([0.75, 0.5, 0.5]), np.eye(3))
    invalid_final = final.copy()
    invalid_deform = np.eye(3)
    invalid_deform[0, 0] = 0.0
    invalid_final.set_cell(cell_from_deformation(invalid_deform, reference_cell), scale_atoms=False)
    try:
        interpolate_vcneb(initial, invalid_final, n_images=3, align_cells=False)
    except ValueError:
        pass
    else:
        raise SystemExit("interpolation accepted a singular endpoint cell")

    images = interpolate_vcneb(initial, final, n_images=3, align_cells=False)
    for image in images:
        image.calc = ToyPhaseTransition(reference_cell)
    chain = VCNEB(images, climb=False)
    trial = chain.get_x()
    cell_offset = 3 * len(images[1])
    trial[cell_offset] = -1.5 * chain.cell_scale
    try:
        chain.set_x(trial)
    except ValueError:
        pass
    else:
        raise SystemExit("VCNEB accepted an optimizer update with non-positive cell determinant")


def check_nonorthogonal_cell_force_regression() -> None:
    reference_cell = np.array(
        [[4.7, 0.2, 0.0], [0.4, 5.1, 0.3], [0.1, 0.2, 5.4]],
        dtype=float,
    )
    deform = np.array(
        [[1.08, 0.12, 0.03], [0.02, 0.97, 0.08], [0.04, 0.01, 1.03]],
        dtype=float,
    )
    q = np.array([0.42, 0.57, 0.48], dtype=float)

    def build(deformation: np.ndarray) -> Atoms:
        atoms = Atoms(
            "Ar",
            scaled_positions=[q],
            cell=cell_from_deformation(deformation, reference_cell),
            pbc=True,
        )
        atoms.calc = MetricCellCalculator(reference_cell)
        return atoms

    atoms = build(deform)
    analytic = cell_force(atoms, reference_cell)
    epsilon = 1e-6
    max_error = 0.0
    for row in range(3):
        for col in range(3):
            direction = np.zeros((3, 3))
            direction[row, col] = 1.0
            numeric = -(
                build(deform + epsilon * direction).get_potential_energy()
                - build(deform - epsilon * direction).get_potential_energy()
            ) / (2.0 * epsilon)
            max_error = max(max_error, abs(numeric - analytic[row, col]))
    if max_error > 1e-7:
        raise SystemExit(f"nonorthogonal cell-force regression failed: {max_error}")


def check_mode_subspace_constraint() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = Atoms("Ar", scaled_positions=[[0.25, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final = Atoms("Ar", scaled_positions=[[0.75, 0.5, 0.5]], cell=reference_cell, pbc=True)
    images = interpolate_vcneb(initial, final, n_images=3, align_cells=False)
    for image in images:
        image.calc = ToyPhaseTransition(reference_cell)

    basis = build_mode_basis(Mode([[1.0, 0.0, 0.0]]), initial)
    chain = VCNEB(images, k=1.0, climb=False, mode_basis=basis)
    before = chain.images[1].get_positions().copy()
    trial = chain.get_x()
    trial[1] += 1.75
    chain.set_x(trial)
    after = chain.images[1].get_positions()
    if not np.allclose(after[0, 1], before[0, 1], rtol=0.0, atol=1e-12):
        raise SystemExit("strict mode subspace did not remove an out-of-subspace update")

    force_x = chain._compute_forces()
    if abs(force_x[1]) > 1e-12:
        raise SystemExit("strict mode subspace leaked force outside the mode basis")

    incompatible = Atoms("Ar", scaled_positions=[[0.75, 0.6, 0.5]], cell=reference_cell, pbc=True)
    incompatible_images = interpolate_vcneb(initial, incompatible, n_images=3, align_cells=False)
    for image in incompatible_images:
        image.calc = ToyPhaseTransition(reference_cell)
    try:
        VCNEB(incompatible_images, climb=False, mode_basis=basis)
    except ValueError as error:
        if "cannot connect the endpoints" not in str(error):
            raise SystemExit("strict mode subspace reported the wrong endpoint error")
    else:
        raise SystemExit("strict mode subspace accepted incompatible endpoints")

    coupled_initial = Atoms("Ar", scaled_positions=[[0.25, 0.5, 0.5]], cell=reference_cell, pbc=True)
    coupled_final = Atoms("Ar", scaled_positions=[[0.75, 0.5, 0.5]], cell=reference_cell, pbc=True)
    coupled_cell = reference_cell.copy()
    coupled_cell[0, 0] *= 1.25
    coupled_final.set_cell(coupled_cell, scale_atoms=True)
    coupled_images = interpolate_vcneb(coupled_initial, coupled_final, n_images=7, align_cells=False)
    for image in coupled_images:
        image.calc = ToyPhaseTransition(reference_cell)
    coupled_mode = Mode([[1.0, 0.0, 0.0]], cell=np.diag([0.1, 0.0, 0.0]))
    coupled_basis = build_mode_basis(coupled_mode, coupled_initial, cell_scale=5.0)
    coupled_chain, _ = run_vcneb(
        coupled_images,
        k=0.15,
        climb=True,
        mode_basis=coupled_basis,
        optimizer="FIRE",
        fmax=0.01,
        steps=300,
        logfile=None,
        trajectory=None,
        snapshot_dir=None,
    )
    coupled_barrier, coupled_reaction = coupled_chain.barrier()
    coupled_saddle = coupled_chain.saddle_diagnostics()
    if abs(coupled_barrier - 0.25) > 5e-3 or abs(coupled_reaction) > 5e-6:
        raise SystemExit(
            f"coupled strict mode changed the known barrier: {coupled_barrier}, {coupled_reaction}"
        )
    if coupled_saddle["residual_generalized_force_eV_per_A"] > 0.01:
        raise SystemExit("coupled strict mode saddle residual exceeded tolerance")
    print("coupled_mode_constraint_regression=ok")


def check_projected_mode_constraint_and_direction_basis() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = Atoms("Ar", scaled_positions=[[0.25, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final = Atoms("Ar", scaled_positions=[[0.75, 0.5, 0.5]], cell=reference_cell, pbc=True)
    images = interpolate_vcneb(initial, final, n_images=3, align_cells=False)
    for image in images:
        image.calc = ToyPhaseTransition(reference_cell)

    basis = build_mode_basis(Mode([[0.0, 1.0, 0.0]]), initial)
    chain = VCNEB(
        images,
        k=1.0,
        climb=False,
        mode_basis=basis,
        constraint_mode="projected",
    )
    before = chain.images[1].get_positions().copy()
    trial = chain.get_x()
    trial[1] += 0.50
    chain.set_x(trial)
    after = chain.images[1].get_positions()
    if abs(after[0, 1] - before[0, 1] - 0.50) > 1e-12:
        raise SystemExit("projected mode constraint did not retain an allowed update")
    force_x = chain._compute_forces()
    if np.max(np.abs(np.delete(force_x, 1))) > 1e-12:
        raise SystemExit("projected mode constraint leaked force outside the mode basis")

    direction_basis = build_direction_basis(np.array([[0.0, 0.0, 1.0]]))
    if direction_basis.shape != (12, 1) or abs(direction_basis[2, 0] - 1.0) > 1e-12:
        raise SystemExit("direction basis has the wrong extended-coordinate layout")

    mixed_basis = build_direction_basis(np.array([[1.0, 1.0, 1.0]]))
    conflict = direction_basis_conflicts(
        mixed_basis,
        atom_mask=np.array([[1.0, 1.0, 0.0]]),
    )
    if not conflict["has_conflict"] or conflict["partially_clipped_columns"] != [0]:
        raise SystemExit("direction basis did not report a partially clipped component")
    fully_inactive = direction_basis_conflicts(
        mixed_basis,
        atom_mask=np.array([[0.0, 0.0, 0.0]]),
    )
    if fully_inactive["fully_inactive_columns"] != [0] or fully_inactive["rank_loss"] != 1:
        raise SystemExit("direction basis did not report a fully inactive component")


def check_climbing_image_saddle_diagnostics() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = Atoms("Ar", scaled_positions=[[0.25, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final = Atoms("Ar", scaled_positions=[[0.75, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final_cell = reference_cell.copy()
    final_cell[0, 0] *= 1.25
    final.set_cell(final_cell, scale_atoms=True)
    images = interpolate_vcneb(initial, final, n_images=7, align_cells=False)
    for index, image in enumerate(images[1:-1], start=1):
        scaled = image.get_scaled_positions(wrap=False)
        scaled[0, 1] += 0.04 * np.sin(index)
        image.set_scaled_positions(scaled)
        image.calc = ToyPhaseTransition(reference_cell)
    images[0].calc = ToyPhaseTransition(reference_cell)
    images[-1].calc = ToyPhaseTransition(reference_cell)

    chain, _ = run_vcneb(
        images,
        k=0.15,
        climb=True,
        optimizer="FIRE",
        fmax=0.01,
        steps=300,
        logfile=None,
        trajectory=None,
        snapshot_dir=None,
    )
    barrier, reaction = chain.barrier()
    diagnostics = chain.saddle_diagnostics()
    image_index = diagnostics["image_index"]
    if image_index is None:
        raise SystemExit("CI diagnostics did not identify an interior image")
    q = chain.images[image_index].get_scaled_positions(wrap=False)
    deform = deformation_from_cell(chain.images[image_index].cell.array, reference_cell)
    saddle_y = np.array([q[0, 0], deform[0, 0] - 1.0])
    if abs(barrier - 0.25) > 5e-3 or abs(reaction) > 5e-6:
        raise SystemExit(f"CI did not recover the analytic barrier: {barrier}, {reaction}")
    if np.linalg.norm(saddle_y - np.array([0.5, 0.125])) > 5e-3:
        raise SystemExit(f"CI saddle coordinate is inaccurate: {saddle_y}")
    if diagnostics["residual_generalized_force_eV_per_A"] > 0.01:
        raise SystemExit("CI saddle residual force exceeded the requested tolerance")
    curvature = diagnostics["tangent_curvature_eV_per_A2"]
    if curvature is None or curvature >= 0.0:
        raise SystemExit(f"CI saddle curvature is not negative: {curvature}")
    print("climbing_image_saddle_regression=ok")


def check_staged_climbing_protocol() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = Atoms("Ar", scaled_positions=[[0.25, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final = Atoms("Ar", scaled_positions=[[0.75, 0.5, 0.5]], cell=reference_cell, pbc=True)
    final_cell = reference_cell.copy()
    final_cell[0, 0] *= 1.25
    final.set_cell(final_cell, scale_atoms=True)
    images = interpolate_vcneb(initial, final, n_images=7, align_cells=False)
    for index, image in enumerate(images[1:-1], start=1):
        scaled = image.get_scaled_positions(wrap=False)
        scaled[0, 1] += 0.08 * np.sin(index)
        image.set_scaled_positions(scaled)
        image.calc = ToyPhaseTransition(reference_cell)
    images[0].calc = ToyPhaseTransition(reference_cell)
    images[-1].calc = ToyPhaseTransition(reference_cell)
    chain, optimizer = run_vcneb(
        images,
        k=0.15,
        climb=True,
        climb_after=2,
        optimizer="FIRE",
        fmax=0.01,
        steps=12,
        logfile=None,
        trajectory=None,
        snapshot_dir=None,
    )
    if not chain.climb or int(getattr(optimizer, "nsteps", 0)) < 2:
        raise SystemExit("staged climbing protocol did not enable CI after the requested steps")


def check_pressure_enthalpy_gradient() -> None:
    reference_cell = np.array(
        [
            [4.7, 0.2, 0.0],
            [0.4, 5.1, 0.3],
            [0.1, 0.2, 5.4],
        ],
        dtype=float,
    )
    q = np.array([0.42, 0.57, 0.48], dtype=float)
    deform = np.array(
        [
            [1.08, 0.02, 0.01],
            [0.00, 0.97, 0.03],
            [0.01, 0.00, 1.03],
        ],
        dtype=float,
    )
    pressure = 0.017
    atoms = Atoms("Ar", scaled_positions=[q], cell=cell_from_deformation(deform, reference_cell), pbc=True)
    atoms.calc = MetricCellCalculator(reference_cell)
    analytic_zero = cell_force(atoms, reference_cell)
    analytic = cell_force(atoms, reference_cell, pressure=pressure)
    eps = 1e-6
    max_zero_error = 0.0
    max_pressure_error = 0.0
    max_error = 0.0
    zero_errors = np.zeros((3, 3))
    pressure_errors = np.zeros((3, 3))
    volume = atoms.get_volume()
    pressure_term = -pressure * volume * np.linalg.inv(deform).T
    for row in range(3):
        for col in range(3):
            dd = np.zeros((3, 3))
            dd[row, col] = eps
            plus = Atoms("Ar", scaled_positions=[q], cell=cell_from_deformation(deform + dd, reference_cell), pbc=True)
            minus = Atoms("Ar", scaled_positions=[q], cell=cell_from_deformation(deform - dd, reference_cell), pbc=True)
            plus.calc = MetricCellCalculator(reference_cell)
            minus.calc = MetricCellCalculator(reference_cell)
            e_plus = plus.get_potential_energy()
            e_minus = minus.get_potential_energy()
            numeric_zero = -(e_plus - e_minus) / (2.0 * eps)
            numeric_volume = -(plus.get_volume() - minus.get_volume()) / (2.0 * eps)
            numeric = numeric_zero + numeric_volume * pressure
            zero_errors[row, col] = abs(numeric_zero - analytic_zero[row, col])
            pressure_errors[row, col] = abs(numeric_volume * pressure - pressure_term[row, col])
            max_zero_error = max(max_zero_error, zero_errors[row, col])
            max_pressure_error = max(max_pressure_error, pressure_errors[row, col])
            max_error = max(max_error, abs(numeric - analytic[row, col]))
    print(f"max_zero_pressure_cell_force_error={max_zero_error:.3e}")
    print(f"max_pressure_volume_force_error={max_pressure_error:.3e}")
    print(f"max_pressure_cell_force_error={max_error:.3e}")
    if max_error > 1e-7:
        print("zero_pressure_error_matrix=" + np.array2string(zero_errors, precision=6))
        print("pressure_volume_error_matrix=" + np.array2string(pressure_errors, precision=6))
        raise SystemExit("pressure enthalpy cell-force finite difference failed")


def check_trajectory_resume() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = make_atoms(reference_cell, np.array([0.25, 0.5, 0.5]), np.eye(3))
    final_deform = np.eye(3)
    final_deform[0, 0] = 1.25
    final = make_atoms(reference_cell, np.array([0.75, 0.5, 0.5]), final_deform)
    images = interpolate_vcneb(initial, final, n_images=4, align_cells=False)
    shifted = [image.copy() for image in images]
    for image in shifted[1:-1]:
        q = image.get_scaled_positions(wrap=False)
        q[:, 1] += 0.123
        image.set_scaled_positions(q)

    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "chain.traj"
        traj = Trajectory(str(path), "w")
        for image in images:
            traj.write(image)
        for image in shifted:
            traj.write(image)
        traj.write(shifted[0])
        traj.close()

        resumed = read_chain_trajectory(path, n_images=4)
        for got, expected in zip(resumed, shifted):
            if not np.allclose(got.positions, expected.positions):
                raise SystemExit("read_chain_trajectory did not ignore the partial tail")

        target = [image.copy() for image in images]
        calculators = [ToyPhaseTransition(reference_cell) for _ in target]
        for image, calc in zip(target, calculators):
            image.calc = calc
        apply_chain_state(target, resumed)
        for image, expected, calc in zip(target, shifted, calculators):
            if image.calc is not calc:
                raise SystemExit("apply_chain_state did not preserve calculators")
            if not np.allclose(image.positions, expected.positions):
                raise SystemExit("apply_chain_state did not copy resumed positions")


def check_snapshot_append_resume() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = make_atoms(reference_cell, np.array([0.25, 0.5, 0.5]), np.eye(3))
    final_deform = np.eye(3)
    final_deform[0, 0] = 1.25
    final = make_atoms(reference_cell, np.array([0.75, 0.5, 0.5]), final_deform)
    images = interpolate_vcneb(initial, final, n_images=4, align_cells=False)
    for image in images:
        image.calc = ToyPhaseTransition(reference_cell)

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        snapshot_dir = root / "snapshots"
        snapshot_dir.mkdir()
        existing = snapshot_dir / "chain_step_0000.traj"
        marker = b"existing snapshot marker"
        existing.write_bytes(marker)

        run_vcneb(
            images,
            k=0.15,
            climb=False,
            optimizer="FIRE",
            fmax=0.01,
            steps=1,
            logfile=None,
            trajectory=root / "chain.traj",
            trajectory_mode="a",
            snapshot_dir=snapshot_dir,
        )

        if existing.read_bytes() != marker:
            raise SystemExit("append-mode snapshots overwrote an existing snapshot")
        if not (snapshot_dir / "chain_step_0001.traj").exists():
            raise SystemExit("append-mode snapshots did not continue at the next step index")


def check_optimizer_api_shapes() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = make_atoms(reference_cell, np.array([0.25, 0.5, 0.5]), np.eye(3))
    final_deform = np.eye(3)
    final_deform[0, 0] = 1.25
    final = make_atoms(reference_cell, np.array([0.75, 0.5, 0.5]), final_deform)

    images = interpolate_vcneb(initial, final, n_images=4, align_cells=False)
    for image in images:
        image.calc = ToyPhaseTransition(reference_cell)
    chain = VCNEB(images, k=0.15, climb=False)
    if chain.get_positions().shape[0] != len(chain):
        raise SystemExit("VCNEB optimizer positions rows do not match __len__")
    if chain.get_masses().shape != (len(chain),):
        raise SystemExit("VCNEB get_masses must return one mass per optimizer row")

    for optimizer in ["BFGS", "LBFGS"]:
        trial_images = interpolate_vcneb(initial, final, n_images=4, align_cells=False)
        for image in trial_images:
            image.calc = ToyPhaseTransition(reference_cell)
        run_vcneb(
            trial_images,
            k=0.15,
            climb=False,
            optimizer=optimizer,
            optimizer_kwargs={"maxstep": 0.05},
            fmax=0.01,
            steps=1,
            logfile=None,
            trajectory=None,
            snapshot_dir=None,
        )


def check_set_x_cache_invalidation() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = make_atoms(reference_cell, np.array([0.25, 0.5, 0.5]), np.eye(3))
    final_deform = np.eye(3)
    final_deform[0, 0] = 1.25
    final = make_atoms(reference_cell, np.array([0.75, 0.5, 0.5]), final_deform)

    images = interpolate_vcneb(initial, final, n_images=3, align_cells=False)
    for image in images:
        image.calc = ToyPhaseTransition(reference_cell)
    chain = VCNEB(images, k=0.15, climb=False)

    cached_value = chain.get_value()
    x = chain.get_x()
    x[0] += 2.5
    chain.set_x(x)
    updated_value = chain.get_value()
    barrier, _ = chain.barrier()

    if abs(cached_value - updated_value) < 1e-6:
        raise SystemExit("set_x did not invalidate cached VC-NEB enthalpies")
    if abs(barrier - (updated_value - chain.enthalpies[0])) > 1e-12:
        raise SystemExit("barrier used stale enthalpies after set_x")


def check_parallel_endpoint_ownership() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    initial = make_atoms(reference_cell, np.array([0.25, 0.5, 0.5]), np.eye(3))
    final_deform = np.eye(3)
    final_deform[0, 0] = 1.25
    final = make_atoms(reference_cell, np.array([0.75, 0.5, 0.5]), final_deform)

    class FakeWorld:
        size = 4

        def __init__(self, rank):
            self.rank = rank

    original_world = core_module.world
    try:
        for rank in range(4):
            core_module.world = FakeWorld(rank)
            images = interpolate_vcneb(initial, final, n_images=5, align_cells=False)
            chain = VCNEB(images, k=0.15, climb=False, parallel=True)
            if chain._own_image(0) != (rank == 0) or chain._own_image(4) != (rank == 0):
                raise SystemExit("parallel endpoints must be owned only by rank 0")
            owned_interiors = [i for i in range(1, 4) if chain._own_image(i)]
            expected = [1 + rank] if rank < 3 else []
            if owned_interiors != expected:
                raise SystemExit("parallel interior image ownership changed unexpectedly")
    finally:
        core_module.world = original_world


def check_threaded_image_executor() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])

    def build_images():
        initial = make_atoms(reference_cell, np.array([0.25, 0.5, 0.5]), np.eye(3))
        final = make_atoms(reference_cell, np.array([0.75, 0.5, 0.5]), np.eye(3))
        images = interpolate_vcneb(initial, final, n_images=4, align_cells=False)
        for image in images:
            image.calc = ToyPhaseTransition(reference_cell)
        return images

    serial_chain = VCNEB(build_images(), k=0.15, climb=False)
    class RecordingExecutor(ThreadedCalculatorExecutor):
        def __init__(self):
            super().__init__(max_workers=2)
            self.batches = []

        def evaluate(self, images, *, indices=None):
            self.batches.append(list(range(len(images))) if indices is None else list(indices))
            return super().evaluate(images, indices=indices)

    recording_executor = RecordingExecutor()
    threaded_chain = VCNEB(
        build_images(),
        k=0.15,
        climb=False,
        image_executor=recording_executor,
    )
    serial_forces = serial_chain.get_forces()
    threaded_forces = threaded_chain.get_forces()
    if not np.allclose(serial_forces, threaded_forces, rtol=1e-11, atol=1e-11):
        raise SystemExit("threaded image executor changed VC-NEB forces")
    diagnostics = threaded_chain.path_diagnostics()
    if len(threaded_chain._last_evaluations or []) != 4:
        raise SystemExit("threaded image executor did not retain one result per image")
    if not recording_executor.batches or any(batch != [1, 2] for batch in recording_executor.batches):
        raise SystemExit(f"executor evaluated non-interior images: {recording_executor.batches}")
    endpoint_cache = dict(threaded_chain._endpoint_evaluations)
    batches_before_reset = len(recording_executor.batches)
    threaded_chain.set_x(threaded_chain.get_x())
    threaded_chain.get_forces()
    if len(recording_executor.batches) != batches_before_reset + 1 or recording_executor.batches[-1] != [1, 2]:
        raise SystemExit(f"executor batch did not stay interior-only: {recording_executor.batches}")
    if any(threaded_chain._endpoint_evaluations[index] is not endpoint_cache[index] for index in (0, 3)):
        raise SystemExit("fixed endpoint evaluations were not reused")
    if diagnostics["n_images"] != 4 or not all(np.isfinite(record["enthalpy_eV"]) for record in diagnostics["images"]):
        raise SystemExit("threaded image executor produced invalid diagnostics")

    class FailOnceCalculator(Calculator):
        implemented_properties = ["energy", "forces", "stress"]

        def __init__(self):
            super().__init__()
            self.calls = 0

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            self.calls += 1
            if self.calls == 1:
                raise RuntimeError("intentional one-time image failure")
            self.results["energy"] = 0.0
            self.results["forces"] = np.zeros((len(self.atoms), 3), dtype=float)
            self.results["stress"] = np.zeros((3, 3), dtype=float)

    retry_image = Atoms("Ar", positions=[[0.0, 0.0, 0.0]], cell=np.eye(3) * 5.0, pbc=True)
    retry_calc = FailOnceCalculator()
    retry_image.calc = retry_calc
    with tempfile.TemporaryDirectory() as tmp:
        manifest_path = Path(tmp) / "worker-manifest.jsonl"
        retry_executor = ThreadedCalculatorExecutor(
            max_workers=1, max_retries=1, manifest_path=manifest_path
        )
        retry_result = retry_executor.evaluate([retry_image])[0]
        if retry_result.energy != 0.0 or retry_executor.last_attempts.get(0) != 2:
            raise SystemExit("threaded image executor did not retry a failed image exactly once")
        records = [json.loads(line) for line in manifest_path.read_text(encoding="utf-8").splitlines()]
        if len(records) != 1 or records[0]["status"] != "ok" or records[0]["attempts"] != {"0": 2}:
            raise SystemExit("threaded image executor did not persist retry manifest")


def check_abacus_command_profile_factory() -> None:
    original_module = sys.modules.get("ase.calculators.abacus")
    fake_module = types.ModuleType("ase.calculators.abacus")
    calls = []

    class FakeProfile:
        def __init__(self, command):
            self.command = command

    class FakeAbacus:
        def __init__(self, directory=".", profile=None, **kwargs):
            self.directory = directory
            self.profile = profile
            self.parameters = kwargs
            calls.append(self)

    fake_module.Abacus = FakeAbacus
    fake_module.AbacusProfile = FakeProfile
    sys.modules["ase.calculators.abacus"] = fake_module
    try:
        reference_cell = np.diag([5.0, 5.0, 5.0])
        atoms = make_atoms(reference_cell, np.array([0.25, 0.5, 0.5]), np.eye(3))
        with tempfile.TemporaryDirectory() as tmp:
            factory = make_ase_abacus_factory(
                parameters={"calculation": "scf", "cal_force": 0, "cal_stress": 0},
                command="mpirun -np 4 abacus",
                pp={"Ar": "Ar.upf"},
                out_stru=0,
            )
            calc = factory(0, atoms, Path(tmp))

        if calc is not calls[-1]:
            raise SystemExit("ABACUS factory did not return the constructed calculator")
        if not isinstance(calc.profile, FakeProfile):
            raise SystemExit("ABACUS command was not converted to an AbacusProfile")
        if calc.profile.command != "mpirun -np 4 abacus":
            raise SystemExit("ABACUS command string was not preserved in the profile")
        if "command" in calc.parameters:
            raise SystemExit("ABACUS command leaked into calculator input parameters")
        for required in ["cal_force", "cal_stress", "out_stru"]:
            if calc.parameters.get(required) != 1:
                raise SystemExit(f"ABACUS factory did not enforce {required}=1")
        if calc.parameters.get("pp") != {"Ar": "Ar.upf"}:
            raise SystemExit("ABACUS factory did not preserve extra calculator kwargs")

        try:
            make_ase_abacus_factory(parameters={}, command="abacus", profile=object())
        except ValueError:
            pass
        else:
            raise SystemExit("ABACUS factory accepted both command and profile")
    finally:
        if original_module is None:
            sys.modules.pop("ase.calculators.abacus", None)
        else:
            sys.modules["ase.calculators.abacus"] = original_module


def check_abacus_validator_required_flags() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        template = root / "template"
        image_dir = root / "images" / "00"
        template.mkdir()
        image_dir.mkdir(parents=True)
        (template / "KPT").write_text("K_POINTS\n0\nGamma\n1 1 1 0 0 0\n", encoding="utf-8")
        (template / "STRU").write_text("ATOMIC_SPECIES\n", encoding="utf-8")
        (template / "INPUT").write_text(
            "calculation scf\nbasis_type lcao\ncal_force 1\ncal_stress 1\n",
            encoding="utf-8",
        )
        write(image_dir / "POSCAR.start", Atoms("Ar", scaled_positions=[[0.0, 0.0, 0.0]], cell=np.eye(3) * 5.0, pbc=True))

        report = validate_abacus(template, root / "images")
        if report["status"] != "failed" or not any("out_stru" in issue for issue in report["issues"]):
            raise SystemExit("ABACUS validator did not require out_stru for VC-NEB templates")

        with (template / "INPUT").open("a", encoding="utf-8") as handle:
            handle.write("out_stru 1\n")
        report = validate_abacus(template, root / "images")
        if report["status"] != "ok":
            raise SystemExit(f"ABACUS validator rejected a minimal valid template: {report['issues']}")

        for disabled_key in ["cal_force", "cal_stress", "out_stru"]:
            input_lines = {
                "calculation": "scf",
                "basis_type": "lcao",
                "cal_force": "1",
                "cal_stress": "1",
                "out_stru": "1",
            }
            input_lines[disabled_key] = "0"
            (template / "INPUT").write_text(
                "\n".join(f"{key} {value}" for key, value in input_lines.items()) + "\n",
                encoding="utf-8",
            )
            report = validate_abacus(template, root / "images")
            if report["status"] != "failed" or not any(disabled_key in issue for issue in report["issues"]):
                raise SystemExit(f"ABACUS validator accepted disabled {disabled_key}=0")


def main() -> None:
    reference_cell = np.array(
        [
            [4.7, 0.2, 0.0],
            [0.4, 5.1, 0.3],
            [0.1, 0.2, 5.4],
        ],
        dtype=float,
    )
    q = np.array([0.42, 0.57, 0.48], dtype=float)
    deform = np.diag([1.08, 0.97, 1.03])
    atoms = make_atoms(reference_cell, q, deform)
    assert np.allclose(deformation_from_cell(atoms.cell.array, reference_cell), deform)

    eps = 1e-6
    f_q = fractional_force(atoms)
    max_q_err = 0.0
    for axis in range(3):
        dq = np.zeros(3)
        dq[axis] = eps
        e_plus = make_atoms(reference_cell, q + dq, deform).get_potential_energy()
        e_minus = make_atoms(reference_cell, q - dq, deform).get_potential_energy()
        numeric_force = -(e_plus - e_minus) / (2.0 * eps)
        max_q_err = max(max_q_err, abs(numeric_force - f_q[0, axis]))

    f_cell = cell_force(atoms, reference_cell)
    max_cell_err = 0.0
    for row in range(3):
        for col in range(3):
            dd = np.zeros((3, 3))
            dd[row, col] = eps
            e_plus = make_atoms(reference_cell, q, deform + dd).get_potential_energy()
            e_minus = make_atoms(reference_cell, q, deform - dd).get_potential_energy()
            numeric_force = -(e_plus - e_minus) / (2.0 * eps)
            max_cell_err = max(max_cell_err, abs(numeric_force - f_cell[row, col]))

    print(f"max_fractional_force_error={max_q_err:.3e}")
    print(f"max_cell_force_error={max_cell_err:.3e}")
    if max_q_err > 1e-7 or max_cell_err > 1e-7:
        raise SystemExit(1)
    check_cell_mask_regression()
    print("cell_mask_regression=ok")
    check_atom_mask_regression()
    print("atom_mask_regression=ok")
    check_mode_guided_path()
    print("mode_guided_path_regression=ok")
    check_cell_interpolation_strategies()
    print("cell_interpolation_regression=ok")
    check_atom_mapping()
    print("atom_mapping_regression=ok")
    check_periodic_translation_alignment()
    check_tangent_and_spring_components()
    check_energy_weighted_tangent_regimes()
    check_nonorthogonal_mode_projection()
    print("nonorthogonal_mode_projection_regression=ok")
    check_calculator_contract()
    print("calculator_contract_regression=ok")
    check_calculator_runtime_diagnostics()
    print("calculator_runtime_diagnostics_regression=ok")
    check_path_diagnostics()
    print("path_diagnostics_regression=ok")
    check_barrierless_path_diagnostics()
    check_path_geometry_preflight()
    print("path_geometry_preflight_regression=ok")
    check_mask_validation()
    print("mask_validation_regression=ok")
    check_cell_validity_guards()
    print("cell_validity_guard_regression=ok")
    check_nonorthogonal_cell_force_regression()
    print("nonorthogonal_cell_force_regression=ok")
    check_mode_subspace_constraint()
    print("mode_subspace_constraint_regression=ok")
    check_projected_mode_constraint_and_direction_basis()
    print("projected_mode_constraint_regression=ok")
    check_climbing_image_saddle_diagnostics()
    check_staged_climbing_protocol()
    print("staged_climbing_protocol_regression=ok")
    check_pressure_enthalpy_gradient()
    print("pressure_enthalpy_gradient_regression=ok")
    check_trajectory_resume()
    print("trajectory_resume_regression=ok")
    check_snapshot_append_resume()
    print("snapshot_append_regression=ok")
    check_optimizer_api_shapes()
    print("optimizer_api_shapes_regression=ok")
    check_set_x_cache_invalidation()
    print("set_x_cache_invalidation_regression=ok")
    check_parallel_endpoint_ownership()
    print("parallel_endpoint_ownership_regression=ok")
    check_threaded_image_executor()
    print("threaded_image_executor_regression=ok")
    check_abacus_command_profile_factory()
    print("abacus_command_profile_regression=ok")
    check_abacus_validator_required_flags()
    print("abacus_validator_required_flags_regression=ok")


if __name__ == "__main__":
    main()
