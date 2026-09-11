"""Variable-cell NEB in an ASE-compatible extended coordinate space.

The implementation optimizes fractional atomic coordinates together with a
cell deformation gradient.  It deliberately stays calculator-agnostic: any ASE
calculator that can provide energy, forces, and stress can be used.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Iterator, Optional, Sequence

import numpy as np

from ase import Atoms
from ase.io import write
from ase.io.trajectory import Trajectory
from ase.optimize import BFGS, FIRE, LBFGS
from ase.parallel import world
from ase.units import GPa


Array = np.ndarray


@dataclass
class VCNEBState:
    """Extended coordinates for one image.

    ``q`` stores fractional coordinates without wrapping.  ``deform`` is the
    ASE/UnitCellFilter-style deformation gradient, with the current cell given
    by ``cell0 @ deform.T``.
    """

    q: Array
    deform: Array


def cell_matrix(atoms: Atoms) -> Array:
    return atoms.cell.array.copy()


def deformation_from_cell(cell: Array, reference_cell: Array) -> Array:
    """Return deformation gradient F for ASE row-vector cells.

    ASE stores the three cell vectors as rows and Cartesian positions are
    ``scaled @ cell``.  The UnitCellFilter convention is
    ``F = solve(cell0, cell).T`` and ``cell = cell0 @ F.T``.
    """

    return np.linalg.solve(reference_cell, cell).T


def cell_from_deformation(deform: Array, reference_cell: Array) -> Array:
    return reference_cell @ deform.T


def state_from_atoms(atoms: Atoms, reference_cell: Array) -> VCNEBState:
    return VCNEBState(
        q=atoms.get_scaled_positions(wrap=False).copy(),
        deform=deformation_from_cell(cell_matrix(atoms), reference_cell),
    )


def apply_state(
    atoms: Atoms,
    state: VCNEBState,
    reference_cell: Array,
    *,
    wrap_positions: bool = False,
) -> None:
    atoms.set_cell(cell_from_deformation(state.deform, reference_cell), scale_atoms=False)
    atoms.set_scaled_positions(state.q)
    if wrap_positions:
        atoms.wrap(eps=1e-12)


def polar_rotation(matrix: Array) -> Array:
    u, _, vt = np.linalg.svd(matrix)
    rot = u @ vt
    if np.linalg.det(rot) < 0.0:
        vt[-1] *= -1.0
        rot = u @ vt
    return rot


def remove_global_rotation(reference: Atoms, atoms: Atoms) -> None:
    """Rotate ``atoms`` into the gauge closest to ``reference``."""

    ref_cell = cell_matrix(reference)
    deform = deformation_from_cell(cell_matrix(atoms), ref_cell)
    rot = polar_rotation(deform)
    atoms.set_cell(cell_matrix(atoms) @ rot, scale_atoms=False)
    atoms.set_positions(atoms.positions @ rot)


def _fractional_delta(a: Array, b: Array, pbc: Sequence[bool], mic: bool) -> Array:
    delta = np.asarray(a) - np.asarray(b)
    if mic:
        pbc_arr = np.asarray(pbc, dtype=bool)
        delta[:, pbc_arr] -= np.rint(delta[:, pbc_arr])
    return delta


def interpolate_vcneb(
    initial: Atoms,
    final: Atoms,
    n_images: int,
    *,
    align_cells: bool = True,
    mic: bool = False,
    wrap_positions: bool = False,
) -> list[Atoms]:
    """Create an initial variable-cell band, including both endpoints."""

    if n_images < 2:
        raise ValueError("n_images must include endpoints and be at least 2")
    if len(initial) != len(final):
        raise ValueError("Initial and final structures have different atom counts")
    if initial.get_chemical_symbols() != final.get_chemical_symbols():
        raise ValueError("Initial and final structures must use the same atom order")

    first = initial.copy()
    last = final.copy()
    if align_cells:
        remove_global_rotation(first, last)

    reference_cell = cell_matrix(first)
    q0 = first.get_scaled_positions(wrap=False)
    q1 = last.get_scaled_positions(wrap=False)
    dq = _fractional_delta(q1, q0, first.pbc, mic)
    deform0 = deformation_from_cell(cell_matrix(first), reference_cell)
    deform1 = deformation_from_cell(cell_matrix(last), reference_cell)

    images = []
    for index in range(n_images):
        lam = index / (n_images - 1)
        image = first.copy()
        state = VCNEBState(
            q=q0 + lam * dq,
            deform=(1.0 - lam) * deform0 + lam * deform1,
        )
        apply_state(image, state, reference_cell, wrap_positions=wrap_positions)
        images.append(image)
    return images


def fractional_force(atoms: Atoms) -> Array:
    """Convert Cartesian forces to forces conjugate to fractional coordinates."""

    return np.asarray(atoms.get_forces()) @ cell_matrix(atoms).T


def cell_force(
    atoms: Atoms,
    reference_cell: Array,
    *,
    pressure: float = 0.0,
    mask: Optional[Array] = None,
) -> Array:
    """Force conjugate to the deformation gradient.

    This follows ASE's UnitCellFilter row-vector convention.  Stress and
    pressure are in eV/A^3; pressure is positive for compression.
    """

    stress = np.asarray(atoms.get_stress(voigt=False))
    volume = atoms.get_volume()
    virial = -volume * (stress + np.eye(3) * pressure)
    deform = deformation_from_cell(cell_matrix(atoms), reference_cell)
    force = np.linalg.solve(deform, virial.T).T
    if mask is not None:
        force = force * np.asarray(mask, dtype=float)
    return force


def _sum_array(value: Array) -> Array:
    summed = world.sum(value)
    return value if summed is None else summed


def _sum_scalar(value: float) -> float:
    try:
        return float(world.sum_scalar(value))
    except AttributeError:
        summed = world.sum(value)
        return float(value if summed is None else summed)


class VCNEB:
    """ASE optimizer target for variable-cell NEB.

    Coordinates used by the optimizer are length-like:

    - fractional coordinates are mapped through the reference cell,
    - deformation-gradient components are multiplied by ``cell_scale``.

    This keeps atomic and cell steps on comparable numerical footing and makes
    ``fmax`` roughly an eV/A criterion for both blocks.
    """

    def __init__(
        self,
        images: Sequence[Atoms],
        *,
        pressure: float = 0.0,
        k: float | Iterable[float] = 0.2,
        climb: bool = True,
        cell_scale: Optional[float] = None,
        atom_mask: Optional[Array] = None,
        cell_mask: Optional[Array] = None,
        mic: bool = False,
        wrap_positions: bool = False,
        parallel: bool = False,
        dynamic_relaxation: float = 1.0,
        dynamic_energy_scale: float = 0.5,
        log: Optional[Callable[[str], None]] = None,
    ) -> None:
        if len(images) < 2:
            raise ValueError("VCNEB needs at least two images")
        if any(len(image) != len(images[0]) for image in images):
            raise ValueError("All images must have the same atom count")
        symbols = images[0].get_chemical_symbols()
        if any(image.get_chemical_symbols() != symbols for image in images):
            raise ValueError("All images must have the same atom order")

        self.images = list(images)
        self.n_images = len(self.images)
        self.reference_cell = cell_matrix(self.images[0])
        self.pressure = float(pressure)
        self.climb = bool(climb)
        self.cell_scale = float(
            cell_scale if cell_scale is not None else abs(np.linalg.det(self.reference_cell)) ** (1.0 / 3.0)
        )
        if self.cell_scale <= 0.0:
            raise ValueError("cell_scale must be positive")
        self.atom_mask = None
        if atom_mask is not None:
            atom_mask_array = np.asarray(atom_mask, dtype=float)
            if atom_mask_array.size != 3 * self.n_atoms:
                raise ValueError("atom_mask must have shape (n_atoms, 3)")
            self.atom_mask = atom_mask_array.reshape(self.n_atoms, 3).copy()
        self.cell_mask = None if cell_mask is None else np.asarray(cell_mask, dtype=float).reshape(3, 3)
        self.mic = bool(mic)
        self.wrap_positions = bool(wrap_positions)
        self.parallel = bool(parallel)
        self.dynamic_relaxation = float(dynamic_relaxation)
        self.dynamic_energy_scale = float(dynamic_energy_scale)
        if not 0.0 <= self.dynamic_relaxation <= 1.0:
            raise ValueError("dynamic_relaxation must be between 0 and 1")
        if self.dynamic_energy_scale <= 0.0:
            raise ValueError("dynamic_energy_scale must be positive")
        self.log = log

        k_arr = np.atleast_1d(np.asarray(k, dtype=float))
        if k_arr.size == 1:
            self.k = np.full(self.n_images - 1, float(k_arr[0]))
        elif k_arr.size == self.n_images - 1:
            self.k = k_arr.copy()
        else:
            raise ValueError("k must be a scalar or have length n_images - 1")

        self._last_enthalpies: Optional[Array] = None
        self._last_forces_x: Optional[Array] = None
        self._owners = self._assign_image_owners()

    def __ase_optimizable__(self) -> "VCNEB":
        return self

    def __len__(self) -> int:
        return self.ndofs() // 3

    @property
    def n_atoms(self) -> int:
        return len(self.images[0])

    @property
    def image_ndofs(self) -> int:
        return 3 * self.n_atoms + 9

    def ndofs(self) -> int:
        return max(0, self.n_images - 2) * self.image_ndofs

    def _assign_image_owners(self) -> list[int]:
        owners = [0] * self.n_images
        if not self.parallel or world.size == 1:
            for image_index in range(1, self.n_images - 1):
                owners[image_index] = 0
            return owners
        for rank_index, image_index in enumerate(range(1, self.n_images - 1)):
            owners[image_index] = rank_index % world.size
        return owners

    def _own_image(self, image_index: int) -> bool:
        if image_index in (0, self.n_images - 1):
            return (not self.parallel) or world.size == 1 or world.rank == 0
        if not self.parallel or world.size == 1:
            return True
        return self._owners[image_index] == world.rank

    def _state(self, image_index: int) -> VCNEBState:
        return state_from_atoms(self.images[image_index], self.reference_cell)

    def _state_to_x(self, state: VCNEBState) -> Array:
        x_atoms = state.q @ self.reference_cell
        x_cell = self.cell_scale * (state.deform - np.eye(3))
        return np.concatenate([x_atoms.reshape(-1), x_cell.reshape(-1)])

    def _active_x_mask(self) -> Array:
        atom_mask = (
            np.ones((self.n_atoms, 3), dtype=float)
            if self.atom_mask is None
            else self.atom_mask
        ).reshape(-1)
        cell_mask = (
            np.ones((3, 3), dtype=float)
            if self.cell_mask is None
            else self.cell_mask
        ).reshape(-1)
        return np.concatenate([atom_mask, cell_mask])

    def _x_to_state(self, x: Array, image_index: int) -> VCNEBState:
        atom_size = 3 * self.n_atoms
        x_atoms = x[:atom_size].reshape(self.n_atoms, 3)
        x_cell = x[atom_size:].reshape(3, 3)
        q = x_atoms @ np.linalg.inv(self.reference_cell)
        if self.atom_mask is not None:
            current = state_from_atoms(self.images[image_index], self.reference_cell)
            current_x_atoms = current.q @ self.reference_cell
            q = (
                (x_atoms * self.atom_mask + current_x_atoms * (1.0 - self.atom_mask))
                @ np.linalg.inv(self.reference_cell)
            )
        deform = np.eye(3) + x_cell / self.cell_scale
        if self.cell_mask is not None:
            current = deformation_from_cell(cell_matrix(self.images[image_index]), self.reference_cell)
            deform = np.eye(3) + (deform - np.eye(3)) * self.cell_mask + (current - np.eye(3)) * (1.0 - self.cell_mask)
        return VCNEBState(q=q, deform=deform)

    def _force_to_x(self, force: VCNEBState) -> Array:
        f_atoms = force.q @ np.linalg.inv(self.reference_cell.T)
        f_cell = force.deform / self.cell_scale
        return np.concatenate([f_atoms.reshape(-1), f_cell.reshape(-1)])

    def _x_to_force(self, x_force: Array) -> VCNEBState:
        atom_size = 3 * self.n_atoms
        f_atoms_x = x_force[:atom_size].reshape(self.n_atoms, 3)
        f_cell_x = x_force[atom_size:].reshape(3, 3)
        return VCNEBState(
            q=f_atoms_x @ self.reference_cell.T,
            deform=f_cell_x * self.cell_scale,
        )

    def _image_x(self, image_index: int) -> Array:
        return self._state_to_x(self._state(image_index))

    def get_x(self) -> Array:
        if self.ndofs() == 0:
            return np.zeros(0)
        return np.concatenate([self._image_x(i) for i in range(1, self.n_images - 1)])

    def set_x(self, x: Array) -> None:
        x = np.asarray(x, dtype=float).reshape(-1)
        if x.size != self.ndofs():
            raise ValueError(f"Expected {self.ndofs()} coordinates, got {x.size}")
        for offset, image_index in enumerate(range(1, self.n_images - 1)):
            lo = offset * self.image_ndofs
            hi = lo + self.image_ndofs
            apply_state(
                self.images[image_index],
                self._x_to_state(x[lo:hi], image_index),
                self.reference_cell,
                wrap_positions=self.wrap_positions,
            )
        self._last_enthalpies = None
        self._last_forces_x = None

    def get_positions(self) -> Array:
        return self.get_x().reshape((-1, 3))

    def set_positions(self, positions: Array) -> None:
        self.set_x(np.asarray(positions, dtype=float).reshape(-1))

    def get_masses(self) -> Array:
        return np.ones(len(self))

    def iterimages(self) -> Iterator[Atoms]:
        yield from self.images

    def _enthalpy_and_force(self, image_index: int) -> tuple[float, VCNEBState]:
        enthalpy = 0.0
        f_q: Optional[Array] = None
        f_cell: Optional[Array] = None
        if self._own_image(image_index):
            atoms = self.images[image_index]
            energy = atoms.get_potential_energy()
            enthalpy = float(energy + self.pressure * atoms.get_volume())
            f_q = fractional_force(atoms)
            f_cell = cell_force(
                atoms,
                self.reference_cell,
                pressure=self.pressure,
                mask=self.cell_mask,
            )
        enthalpy = _sum_scalar(enthalpy)
        if f_q is None:
            f_q = np.zeros((self.n_atoms, 3), dtype=float)
        if f_cell is None:
            f_cell = np.zeros((3, 3), dtype=float)
        return enthalpy, VCNEBState(q=_sum_array(f_q), deform=_sum_array(f_cell))

    def _tangent(self, image_index: int, enthalpies: Array, image_x: list[Array]) -> Array:
        d_minus = image_x[image_index] - image_x[image_index - 1]
        d_plus = image_x[image_index + 1] - image_x[image_index]
        if enthalpies[image_index + 1] > enthalpies[image_index] > enthalpies[image_index - 1]:
            tangent = d_plus
        elif enthalpies[image_index + 1] < enthalpies[image_index] < enthalpies[image_index - 1]:
            tangent = d_minus
        else:
            tangent = (
                abs(enthalpies[image_index + 1] - enthalpies[image_index]) * d_plus
                + abs(enthalpies[image_index - 1] - enthalpies[image_index]) * d_minus
            )
        norm = np.linalg.norm(tangent)
        if norm < 1e-14:
            tangent = d_plus
            norm = np.linalg.norm(tangent)
        if norm < 1e-14:
            return np.zeros_like(tangent)
        return tangent / norm

    def _dynamic_weights(self, enthalpies: Array) -> Array:
        weights = np.ones(self.n_images)
        alpha = self.dynamic_relaxation
        if alpha >= 1.0 or self.n_images <= 2:
            return weights
        scale = max(abs(self.dynamic_energy_scale), 1e-12)
        emax = float(enthalpies.max())
        for image_index in range(1, self.n_images - 1):
            weights[image_index] = alpha + (1.0 - alpha) * np.exp(-(emax - enthalpies[image_index]) / scale)
        return weights

    def get_forces(self) -> Array:
        return self._compute_forces().reshape((-1, 3))

    def get_gradient(self) -> Array:
        return -self._compute_forces()

    def _compute_forces(self) -> Array:
        enthalpies = np.zeros(self.n_images)
        true_forces = []
        image_x = []
        image_x_active = []
        active_mask = self._active_x_mask()
        for image_index in range(self.n_images):
            h_i, f_i = self._enthalpy_and_force(image_index)
            enthalpies[image_index] = h_i
            true_forces.append(self._force_to_x(f_i))
            x_i = self._image_x(image_index)
            image_x.append(x_i)
            image_x_active.append(x_i * active_mask)

        self._last_enthalpies = enthalpies.copy()
        force_x = np.zeros(self.ndofs())
        if self.n_images <= 2:
            self._last_forces_x = force_x
            return force_x

        climbing_image = int(1 + np.argmax(enthalpies[1:-1]))
        weights = self._dynamic_weights(enthalpies)

        for offset, image_index in enumerate(range(1, self.n_images - 1)):
            tangent = self._tangent(image_index, enthalpies, image_x_active)
            true_force = true_forces[image_index] * active_mask
            true_parallel = np.dot(true_force, tangent) * tangent
            force_perp = true_force - true_parallel

            d_plus = np.linalg.norm(image_x_active[image_index + 1] - image_x_active[image_index])
            d_minus = np.linalg.norm(image_x_active[image_index] - image_x_active[image_index - 1])
            spring = self.k[image_index - 1] * (d_plus - d_minus) * tangent
            neb_force = force_perp + spring

            if self.climb and image_index == climbing_image:
                neb_force = true_force - 2.0 * true_parallel
            neb_force *= weights[image_index] * active_mask

            lo = offset * self.image_ndofs
            hi = lo + self.image_ndofs
            force_x[lo:hi] = neb_force

        self._last_forces_x = force_x.copy()
        if self.log is not None:
            max_force = 0.0 if force_x.size == 0 else float(np.linalg.norm(force_x.reshape(-1, 3), axis=1).max())
            self.log(f"max_force={max_force:.6g} enthalpies=" + " ".join(f"{e:.8f}" for e in enthalpies))
        return force_x

    def get_potential_energy(self, force_consistent: bool = False) -> float:
        return self.get_value()

    def get_value(self) -> float:
        if self._last_enthalpies is None:
            enthalpies = [self._enthalpy_and_force(i)[0] for i in range(self.n_images)]
            self._last_enthalpies = np.asarray(enthalpies, dtype=float)
        return float(self._last_enthalpies.max())

    def converged(self, gradient: Array, fmax: float) -> bool:
        max_force = self.gradient_norm(gradient)
        return bool(np.isfinite(max_force) and max_force < fmax)

    def gradient_norm(self, gradient: Array) -> float:
        if gradient.size == 0:
            return 0.0
        forces = (-gradient).reshape(-1, 3)
        return float(np.linalg.norm(forces, axis=1).max())

    @property
    def enthalpies(self) -> Array:
        if self._last_enthalpies is None:
            self.get_forces()
        assert self._last_enthalpies is not None
        return self._last_enthalpies.copy()

    def reaction_coordinate(self) -> Array:
        coords = [0.0]
        active_mask = self._active_x_mask()
        for image_index in range(1, self.n_images):
            dx = (self._image_x(image_index) - self._image_x(image_index - 1)) * active_mask
            coords.append(coords[-1] + float(np.linalg.norm(dx)))
        return np.asarray(coords)

    def barrier(self) -> tuple[float, float]:
        enthalpies = self.enthalpies
        return float(enthalpies.max() - enthalpies[0]), float(enthalpies[-1] - enthalpies[0])

    def write_chain(self, path: str | Path) -> None:
        write(str(path), self.images)

    def write_step_directory(self, directory: str | Path, step: int) -> None:
        root = Path(directory)
        root.mkdir(parents=True, exist_ok=True)
        write(str(root / f"chain_step_{step:04d}.traj"), self.images)
        step_dir = root / f"step_{step:04d}"
        step_dir.mkdir(exist_ok=True)
        for image_index, image in enumerate(self.images):
            write(str(step_dir / f"POSCAR_{image_index:02d}"), image, format="vasp", direct=True, vasp5=True)

    def plot_band(self, filename: str | Path) -> None:
        import matplotlib.pyplot as plt

        s = self.reaction_coordinate()
        h = self.enthalpies
        fig, ax = plt.subplots()
        ax.plot(s, h - h[0], marker="o")
        ax.set_xlabel("Reaction coordinate (A, extended)")
        ax.set_ylabel("Relative enthalpy (eV)")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(filename, dpi=200)
        plt.close(fig)


def read_chain_trajectory(
    path: str | Path,
    *,
    n_images: int,
    step: int = -1,
) -> list[Atoms]:
    """Read one complete VC-NEB chain from a flat ASE trajectory.

    ``run_vcneb`` stores each optimizer snapshot as ``n_images`` consecutive
    trajectory frames.  Interrupted writes can leave a partial tail; negative
    steps count over complete snapshots only, so ``step=-1`` returns the latest
    complete chain and ignores any partial trailing frames.
    """

    if n_images < 2:
        raise ValueError("n_images must include endpoints and be at least 2")
    frames = list(Trajectory(str(path), "r"))
    complete_steps = len(frames) // n_images
    if complete_steps == 0:
        raise ValueError(f"No complete {n_images}-image chain found in {path}")
    step_index = complete_steps + step if step < 0 else step
    if step_index < 0 or step_index >= complete_steps:
        raise IndexError(f"Trajectory {path} has {complete_steps} complete chain snapshots; got step {step}")
    chain = [frame.copy() for frame in frames[step_index * n_images : (step_index + 1) * n_images]]
    _validate_chain_compatibility(chain)
    return chain


def apply_chain_state(target: Sequence[Atoms], source: Sequence[Atoms]) -> None:
    """Copy cells and positions from ``source`` images while preserving calculators."""

    if len(target) != len(source):
        raise ValueError(f"Expected {len(target)} source images, got {len(source)}")
    _validate_chain_compatibility(source)
    target_symbols = target[0].get_chemical_symbols()
    for image_index, (dst, src) in enumerate(zip(target, source)):
        if dst.get_chemical_symbols() != target_symbols or src.get_chemical_symbols() != target_symbols:
            raise ValueError(f"Image {image_index} atom symbols/order are incompatible with the target chain")
        dst.set_cell(src.cell.array, scale_atoms=False)
        dst.set_positions(src.positions)
        dst.pbc = src.pbc


def _validate_chain_compatibility(images: Sequence[Atoms]) -> None:
    if not images:
        raise ValueError("Chain is empty")
    symbols = images[0].get_chemical_symbols()
    for image_index, image in enumerate(images):
        if image.get_chemical_symbols() != symbols:
            raise ValueError(f"Image {image_index} atom symbols/order differ from image 0")
        if image.cell.rank != 3 or image.get_volume() <= 0.0:
            raise ValueError(f"Image {image_index} has an invalid 3D cell")


def _make_optimizer(name: str, chain: VCNEB, logfile: str | Path | None):
    key = name.upper()
    if key == "FIRE":
        return FIRE(chain, logfile=logfile)
    if key == "LBFGS":
        return LBFGS(chain, logfile=logfile)
    if key == "BFGS":
        return BFGS(chain, logfile=logfile)
    raise ValueError(f"Unknown optimizer {name!r}; choose BFGS, LBFGS, or FIRE")


def _next_snapshot_step(directory: str | Path) -> int:
    root = Path(directory)
    if not root.exists():
        return 0
    next_step = 0
    for path in root.iterdir():
        name = path.name
        if name.startswith("chain_step_") and name.endswith(".traj"):
            token = name.removeprefix("chain_step_").removesuffix(".traj")
        elif name.startswith("step_"):
            token = name.removeprefix("step_")
        else:
            continue
        try:
            next_step = max(next_step, int(token) + 1)
        except ValueError:
            continue
    return next_step


def run_vcneb(
    images: Sequence[Atoms],
    *,
    pressure_gpa: float = 0.0,
    k: float | Iterable[float] = 0.2,
    climb: bool = True,
    cell_scale: Optional[float] = None,
    atom_mask: Optional[Array] = None,
    cell_mask: Optional[Array] = None,
    mic: bool = False,
    wrap_positions: bool = False,
    parallel: bool = False,
    optimizer: str = "FIRE",
    fmax: float = 0.05,
    steps: int = 300,
    logfile: str | Path | None = "vcneb-opt.log",
    trajectory: str | Path | None = "vcneb.traj",
    trajectory_mode: str = "w",
    snapshot_dir: str | Path | None = None,
    snapshot_start: Optional[int] = None,
    log: Optional[Callable[[str], None]] = None,
) -> tuple[VCNEB, object]:
    """Run a VC-NEB optimization with an ASE optimizer."""

    chain = VCNEB(
        images,
        pressure=pressure_gpa * GPa,
        k=k,
        climb=climb,
        cell_scale=cell_scale,
        atom_mask=atom_mask,
        cell_mask=cell_mask,
        mic=mic,
        wrap_positions=wrap_positions,
        parallel=parallel,
        log=log,
    )
    opt = _make_optimizer(optimizer, chain, logfile)

    if trajectory_mode not in {"w", "a"}:
        raise ValueError("trajectory_mode must be 'w' or 'a'")
    traj = Trajectory(str(trajectory), trajectory_mode) if trajectory else None
    if snapshot_start is None:
        snapshot_start = _next_snapshot_step(snapshot_dir) if trajectory_mode == "a" and snapshot_dir is not None else 0
    if snapshot_start < 0:
        raise ValueError("snapshot_start must be non-negative")
    step_counter = {"value": int(snapshot_start)}

    def save_snapshot() -> None:
        if traj is not None:
            for image in chain.images:
                traj.write(image)
        if snapshot_dir is not None:
            chain.write_step_directory(snapshot_dir, step_counter["value"])
        step_counter["value"] += 1

    opt.attach(save_snapshot, interval=1)
    opt.run(fmax=fmax, steps=steps)
    if traj is not None:
        traj.close()
    return chain, opt
