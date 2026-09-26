"""Reference-cell atomic displacements and symmetric cell strain for ASE.

The independent coordinates are Cartesian atomic displacements in a fixed
reference cell (Angstrom), followed by six symmetric strain components in ASE
Voigt order (xx, yy, zz, yz, xz, xy). Atomic identity, origin, and periodic
image choices must be audited before constructing this object. No calculator
is launched by coordinate conversion alone.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np
from ase import Atoms


Array = np.ndarray


def _strain_matrix(voigt: Array) -> Array:
    xx, yy, zz, yz, xz, xy = voigt
    return np.array([[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]], dtype=float)


def _strain_voigt(matrix: Array) -> Array:
    return np.array([
        matrix[0, 0], matrix[1, 1], matrix[2, 2],
        matrix[1, 2], matrix[0, 2], matrix[0, 1],
    ])


class ReferenceCellCoordinates:
    """Map ``(u_1,...,u_N, eta_6)`` to periodic ASE structures.

    ASE stores lattice vectors as rows. For reference cell ``C0`` and symmetric
    ``A=I+eta``, this convention uses ``C=C0@A`` and
    ``r=(s0+u@inv(C0))@C``. Thus ``u`` is the reference-cell Cartesian
    displacement, not the current-cell Cartesian displacement. Rigid cell
    rotation is excluded from the six-strain chart.

    The reference atoms must already be reordered, translated and unwrapped
    into the desired atom gauge. Inverting another structure uses a nearest
    periodic image relative to this reference; it cannot discover a different
    atom permutation or resolve a half-cell ambiguity on its own.
    """

    def __init__(self, reference: Atoms):
        if not isinstance(reference, Atoms) or len(reference) < 1:
            raise ValueError("reference must be a nonempty ASE Atoms")
        cell = np.asarray(reference.cell.array, dtype=float)
        if not np.all(np.isfinite(cell)) or np.linalg.det(cell) <= 0.0:
            raise ValueError("reference cell must be finite and right-handed")
        self._reference = reference.copy()
        self._reference.calc = None
        self._cell = cell.copy()
        self._inverse_cell = np.linalg.inv(cell)
        self._scaled = reference.get_scaled_positions(wrap=False).copy()
        self._numbers = reference.get_atomic_numbers().copy()
        self._pbc = np.asarray(reference.pbc, dtype=bool).copy()

    @property
    def n_atoms(self) -> int:
        return len(self._reference)

    @property
    def coordinate_count(self) -> int:
        return 3 * self.n_atoms + 6

    @property
    def reference_cell(self) -> Array:
        return self._cell.copy()

    @property
    def reference_atoms(self) -> Atoms:
        return self._reference.copy()

    def rigid_translation_directions(self) -> Array:
        """Three Cartesian translations in the ``(u, eta)`` chart.

        These are useful as explicit frozen gauge directions during
        conditional mode relaxation. The caller must still check that its
        chosen mode axes are metric-orthogonal to them.
        """

        directions = np.zeros((self.coordinate_count, 3), dtype=float)
        directions[:3 * self.n_atoms] = np.tile(np.eye(3), (self.n_atoms, 1))
        return directions

    def _split(self, coordinates: Array | Sequence[float]) -> tuple[Array, Array, Array]:
        values = np.asarray(coordinates, dtype=float).reshape(-1)
        if values.shape != (self.coordinate_count,) or not np.all(np.isfinite(values)):
            raise ValueError(f"coordinates must contain {self.coordinate_count} finite values")
        displacement = values[:3 * self.n_atoms].reshape(self.n_atoms, 3)
        strain = _strain_matrix(values[-6:])
        deformation = np.eye(3) + strain
        if np.min(np.linalg.eigvalsh(deformation)) <= 1e-8:
            raise ValueError("symmetric deformation must be positive definite")
        return displacement, strain, deformation

    def to_atoms(self, coordinates: Array | Sequence[float]) -> Atoms:
        """Construct a calculator-free structure at the declared mode/strain point."""

        displacement, _, deformation = self._split(coordinates)
        atoms = self._reference.copy()
        atoms.set_cell(self._cell @ deformation, scale_atoms=False)
        atoms.set_scaled_positions(self._scaled + displacement @ self._inverse_cell)
        atoms.calc = None
        return atoms

    def remove_mass_weighted_translation(
        self, coordinates: Array | Sequence[float], masses_amu: Array | Sequence[float]
    ) -> tuple[Array, Array]:
        """Return a centered copy plus its removed Cartesian translation.

        Use this only when the mode provenance explicitly removes rigid
        translations. The resulting structure is a rigidly shifted gauge of
        the original and should not be compared by absolute positions.
        """

        displacement, _, _ = self._split(coordinates)
        masses = np.asarray(masses_amu, dtype=float).reshape(-1)
        if masses.shape != (self.n_atoms,) or not np.all(np.isfinite(masses)) or np.any(masses <= 0.0):
            raise ValueError("masses_amu must contain one positive mass per atom")
        translation = np.average(displacement, axis=0, weights=masses)
        centered = np.asarray(coordinates, dtype=float).reshape(-1).copy()
        centered[:3 * self.n_atoms] = (displacement - translation).reshape(-1)
        return centered, translation

    def from_atoms(self, atoms: Atoms, *, symmetry_tolerance: float = 1e-8) -> Array:
        """Invert a structure in the fixed reference atom and periodic gauge."""

        if len(atoms) != self.n_atoms or not np.array_equal(atoms.get_atomic_numbers(), self._numbers):
            raise ValueError("atom count/order/species differ from the reference")
        if not np.array_equal(np.asarray(atoms.pbc, dtype=bool), self._pbc):
            raise ValueError("periodic boundary flags differ from the reference")
        deformation = self._inverse_cell @ atoms.cell.array
        if not np.all(np.isfinite(deformation)) or not np.allclose(
            deformation, deformation.T, rtol=0.0, atol=symmetry_tolerance
        ):
            raise ValueError("cell contains an unsupported rotation or nonsymmetric deformation")
        deformation = 0.5 * (deformation + deformation.T)
        if np.min(np.linalg.eigvalsh(deformation)) <= 1e-8:
            raise ValueError("symmetric deformation must be positive definite")
        scaled_delta = atoms.get_scaled_positions(wrap=False) - self._scaled
        for axis in range(3):
            if self._pbc[axis]:
                scaled_delta[:, axis] -= np.rint(scaled_delta[:, axis])
                if np.any(np.abs(np.abs(scaled_delta[:, axis]) - 0.5) < 1e-8):
                    raise ValueError("half-cell atomic image is ambiguous; align the atom gauge first")
        displacement = scaled_delta @ self._cell
        return np.concatenate([displacement.reshape(-1), _strain_voigt(deformation - np.eye(3))])

    def energy_and_gradient(
        self,
        coordinates: Array | Sequence[float],
        atoms: Atoms,
        *,
        pressure_eV_per_A3: float = 0.0,
    ) -> tuple[float, Array]:
        """Return ``H=E+PV`` and its physical coordinate gradient.

        ``atoms`` must be the structure returned by :meth:`to_atoms` for these
        coordinates and must have an ASE calculator supplying energy, atomic
        forces and symmetric tensile-positive stress in eV/Angstrom^3. The
        pressure is positive for compression. ``dH/d(u,eta)`` is obtained from
        *raw* calculator forces/stress, never NEB spring or climbing forces.
        Strain off-diagonal coordinates are symmetric matrix entries, so their
        derivatives have the required factor of two.
        """

        displacement, _, deformation = self._split(coordinates)
        expected_cell = self._cell @ deformation
        expected_scaled = self._scaled + displacement @ self._inverse_cell
        if len(atoms) != self.n_atoms or not np.array_equal(atoms.get_atomic_numbers(), self._numbers):
            raise ValueError("calculator structure has incompatible atoms")
        if not np.allclose(atoms.cell.array, expected_cell, rtol=0.0, atol=1e-8):
            raise ValueError("calculator structure cell disagrees with coordinates")
        scaled_error = atoms.get_scaled_positions(wrap=False) - expected_scaled
        for axis in range(3):
            if self._pbc[axis]:
                scaled_error[:, axis] -= np.rint(scaled_error[:, axis])
        if not np.allclose(scaled_error, 0.0, rtol=0.0, atol=1e-8):
            raise ValueError("calculator structure positions disagree with coordinates")
        pressure = float(pressure_eV_per_A3)
        if not np.isfinite(pressure):
            raise ValueError("pressure_eV_per_A3 must be finite")
        energy = float(atoms.get_potential_energy())
        forces = np.asarray(atoms.get_forces(), dtype=float)
        stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
        if forces.shape != (self.n_atoms, 3) or stress.shape != (3, 3):
            raise ValueError("calculator must supply full atomic forces and 3x3 stress")
        if not np.isfinite(energy) or not np.all(np.isfinite(forces)) or not np.all(np.isfinite(stress)):
            raise ValueError("calculator returned non-finite energy/forces/stress")
        if not np.allclose(stress, stress.T, rtol=0.0, atol=1e-8):
            raise ValueError("calculator stress must be symmetric")
        volume = float(atoms.get_volume())
        gradient_atomic = -forces @ deformation.T
        stress_with_pressure = 0.5 * (stress + stress.T) + pressure * np.eye(3)
        strain_gradient = volume * np.linalg.solve(deformation.T, stress_with_pressure)
        strain_gradient = 0.5 * (strain_gradient + strain_gradient.T)
        gradient_voigt = _strain_voigt(strain_gradient)
        gradient_voigt[3:] *= 2.0
        gradient = np.concatenate([gradient_atomic.reshape(-1), gradient_voigt])
        return energy + pressure * volume, gradient


__all__ = ["ReferenceCellCoordinates"]
