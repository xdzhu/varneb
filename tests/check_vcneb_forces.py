"""Finite-difference checks for VC-NEB force transforms."""

from __future__ import annotations

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
from vcneb import apply_chain_state, interpolate_vcneb, read_chain_trajectory, run_vcneb
from vcneb import Mode, build_direction_basis, build_mode_basis, mode_guided_path, project_path_onto_modes
from vcneb.abacus import make_ase_abacus_factory
from vcneb.core import VCNEB, cell_force, cell_from_deformation, deformation_from_cell, fractional_force
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


def check_nonorthogonal_mode_projection() -> None:
    reference_cell = np.diag([5.0, 5.0, 5.0])
    reference = Atoms("Ar", scaled_positions=[[0.25, 0.5, 0.5]], cell=reference_cell, pbc=True)
    image = reference.copy()
    scaled = image.get_scaled_positions(wrap=False)
    scaled[0, :2] += [0.20, 0.10]
    image.set_scaled_positions(scaled)
    modes = [Mode([[1.0, 0.0, 0.0]]), Mode([[1.0, 1.0, 0.0]])]
    coefficients = project_path_onto_modes([reference, image], reference, modes)
    expected = np.array([0.50, 0.50])
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
    check_nonorthogonal_mode_projection()
    print("nonorthogonal_mode_projection_regression=ok")
    check_mask_validation()
    print("mask_validation_regression=ok")
    check_mode_subspace_constraint()
    print("mode_subspace_constraint_regression=ok")
    check_projected_mode_constraint_and_direction_basis()
    print("projected_mode_constraint_regression=ok")
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
    check_abacus_command_profile_factory()
    print("abacus_command_profile_regression=ok")
    check_abacus_validator_required_flags()
    print("abacus_validator_required_flags_regression=ok")


if __name__ == "__main__":
    main()
