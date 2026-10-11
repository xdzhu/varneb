"""Mode evaluations in the exact clamped-plane path/endpoint subspace.

This adapter reuses active joint coordinates, not a symmetric-strain mask.
It supplies the coordinate interface used by CalculatorModeEvaluator without
selecting modes, invoking a backend or releasing a prescribed substrate.
"""
from __future__ import annotations

from numbers import Real

import numpy as np
from ase import Atoms
from ase.stress import voigt_6_to_full_3x3_stress

from .active_curvature import ActiveJointCurvatureCoordinates
from .epitaxial_boundary import ClampedPlaneBoundary


class ClampedModeCoordinates:
    """Reference atomic u (Angstrom), then scaled permitted deformation (Angstrom).

    H=H0@(I+sum(y_cell/L)*D).T; s=s0+u@inv(H0). The Frobenius-orthonormal D
    are the same one/three directions as the declared endpoint/VCNEB boundary.
    With tilt allowed they are generally nonsymmetric. All atoms remain
    movable; three translations can be explicitly frozen as a coordinate gauge.
    There is no substrate-control coordinate eligible for orthogonal release.
    """

    coordinate_convention = "reference_cell_atomic_displacements_A_plus_scaled_clamped_deformation_A"

    def __init__(self, reference: Atoms, *, boundary: ClampedPlaneBoundary, cell_scale_A: float):
        if not isinstance(reference, Atoms) or not isinstance(boundary, ClampedPlaneBoundary):
            raise TypeError("ASE reference and explicit ClampedPlaneBoundary required")
        if (isinstance(cell_scale_A, (bool, np.bool_)) or not isinstance(cell_scale_A, Real)
                or not np.isfinite(cell_scale_A) or cell_scale_A <= 0):
            raise ValueError("cell_scale_A must be finite and positive in Angstrom")
        owned = reference.copy()
        owned.calc = None
        self._joint = ActiveJointCurvatureCoordinates.for_clamped_plane(
            owned, cell_scale_A=float(cell_scale_A), boundary=boundary,
        )
        self._boundary = boundary

    @property
    def coordinate_count(self) -> int:
        return self._joint.size

    @property
    def n_atoms(self) -> int:
        return len(self._joint.reference)

    @property
    def cell_dofs(self) -> int:
        return self._joint.cell_dofs

    @property
    def reference_atoms(self) -> Atoms:
        return self._joint.reference.copy()

    @property
    def reference_cell(self) -> np.ndarray:
        return self._joint.reference.cell.array.copy()

    @property
    def mechanical_contract(self) -> dict:
        """Numeric boundary/chart identity; no filename/optimizer identity."""
        return dict(kind="clamped_plane", reference_cell_A=self._boundary.reference_cell.tolist(),
            fixed_cell_rows=[0, 1], allow_tilt=self._boundary.allow_tilt,
            cell_dofs=self.cell_dofs, cell_scale_A=float(self._joint.cell_scale_A),
            deformation_basis=self._joint.deformation_basis.tolist(),
            raw_geometry_policy="reject_incompatible_no_projection",
            atomic_coordinate_unit="Angstrom", scaled_cell_coordinate_unit="Angstrom",
            gradient_unit="eV/Angstrom", prescribed_substrate_relaxable=False)

    def rigid_translation_directions(self) -> np.ndarray:
        directions = np.zeros((self.coordinate_count, 3))
        directions[:3*self.n_atoms] = np.tile(np.eye(3), (self.n_atoms, 1))
        return directions

    def to_atoms(self, coordinates) -> Atoms:
        """Construct the declared ordered/unwrapped geometry without calculation."""
        return self._joint.displaced(np.asarray(coordinates, dtype=float))

    def from_atoms(self, atoms: Atoms) -> np.ndarray:
        """Invert an already ordered continuous lift; never MIC or rematch atoms.

        The caller is responsible for frame, origin, atom order and periodic
        lift registration. Integer wraps are retained as different coordinates,
        not automatically removed as in a nearest-reference geometry matcher.
        """
        reference = self._joint.reference
        if (not isinstance(atoms, Atoms)
                or atoms.get_chemical_symbols() != reference.get_chemical_symbols()
                or not np.array_equal(atoms.pbc, reference.pbc)):
            raise ValueError("ordered atom/species/periodic identity differs from reference")
        self._joint._validate_candidate(atoms)
        deformation = np.linalg.solve(reference.cell.array, atoms.cell.array).T - np.eye(3)
        coefficients = np.einsum("aij,ij->a", self._joint.deformation_basis, deformation)
        recovered = np.einsum("a,aij->ij", coefficients, self._joint.deformation_basis)
        if not np.allclose(recovered, deformation, atol=1e-9, rtol=0):
            raise ValueError("cell contains deformation outside the declared clamped chart")
        u = (atoms.get_scaled_positions(wrap=False)-reference.get_scaled_positions(wrap=False)) @ reference.cell.array
        return np.concatenate((u.ravel(), coefficients*self._joint.cell_scale_A))

    def energy_and_gradient(self, coordinates, atoms: Atoms, *, pressure_eV_per_A3: float = 0.):
        """H=E+PV in eV and raw physical dH/dy in eV/Angstrom, not NEB forces."""
        expected = self.to_atoms(coordinates)
        if (not isinstance(atoms, Atoms)
                or atoms.get_chemical_symbols() != expected.get_chemical_symbols()
                or not np.array_equal(atoms.pbc, expected.pbc)
                or atoms.constraints
                or not np.allclose(atoms.cell.array, expected.cell.array, atol=1e-8, rtol=0)):
            raise ValueError("calculator structure disagrees with clamped coordinates")
        fractional_error = atoms.get_scaled_positions(wrap=False) - expected.get_scaled_positions(wrap=False)
        fractional_error[:, expected.pbc] -= np.rint(fractional_error[:, expected.pbc])
        if not np.allclose(fractional_error, 0., atol=1e-8, rtol=0):
            raise ValueError("calculator ordered positions disagree with clamped coordinates")
        pressure = float(pressure_eV_per_A3)
        if not np.isfinite(pressure):
            raise ValueError("pressure must be finite in eV/Angstrom^3")
        energy = float(atoms.get_potential_energy())
        forces = np.asarray(atoms.get_forces(apply_constraint=False), dtype=float)
        # ASE's public get_stress can average an asymmetric raw3x3 tensor
        # while converting to Voigt. Do not silently repair such backend data.
        raw_stress = np.asarray(atoms.calc.get_property("stress", atoms), dtype=float)
        stress = (voigt_6_to_full_3x3_stress(raw_stress) if raw_stress.shape == (6,)
                  else raw_stress)
        if (not np.isfinite(energy) or forces.shape != (self.n_atoms, 3)
                or stress.shape != (3, 3) or not np.isfinite(forces).all()
                or not np.isfinite(stress).all()
                or not np.allclose(stress, stress.T, atol=1e-8, rtol=0)):
            raise ValueError("calculator must supply finite energy/forces/symmetric ASE stress")
        gradient = self._joint.enthalpy_gradient(atoms, forces, stress, pressure)
        return energy + pressure*atoms.get_volume(), gradient
