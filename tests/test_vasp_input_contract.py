"""Input-source prevention regressions from both real GaN Bravais failures."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.calculator import CalculationFailed
from ase.calculators.calculator import Calculator
from ase.calculators.vasp import Vasp
from ase.io import read, write

from vcneb.calculator import classify_calculator_failure
from vcneb.executor import ThreadedCalculatorExecutor
from ase.calculators.singlepoint import SinglePointCalculator
from vcneb.vasp import (
    REQUIRED_VCNEB_STATIC_PARAMETERS, ExplicitPotcarVasp,
    cached_vasp_static_endpoint_calculator, collect_vasp_params,
    prepare_vasp_static_parameters, validate_vasp_static_parameters,
    vasp_input_fingerprints,
)
from vcneb.provenance import endpoint_structure_record
from vcneb.vasp_contract import (
    POSCAR_LATTICE_SIGNIFICANT_DIGITS, VaspInputContractError,
    canonical_parameters, poscar_lattice_roundtrip, rewrite_poscar_lattice_exact,
    validate_vasp_image_geometry,
)

FIXTURES = Path(__file__).parent / "fixtures" / "vasp_bravais"


def source_inputs(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    (source / "INCAR").write_text("ENCUT=600\nISYM=2\nSYMPREC=1e-5\nNSW=40\nIBRION=2\nISIF=3\n")
    (source / "KPOINTS").write_text("Automatic\n0\nGamma\n8 8 6\n0 0 0\n")
    (source / "POTCAR").write_text(" TITEL = PAW_PBE Ga_d 06Sep2000\n TITEL = PAW_PBE N 08Apr2002\n")
    return source


def calculator(tmp_path, atoms):
    source = source_inputs(tmp_path)
    parameters, potcar = prepare_vasp_static_parameters(source)
    directory = tmp_path / "image"
    directory.mkdir()
    calc = ExplicitPotcarVasp(source_potcar=potcar, directory=directory,
                             minimum_distance=1.4, **parameters)
    calc.lock_input_contract(atoms, source)
    return calc, source


@pytest.mark.parametrize("name", ["gan_image15.vasp", "gan_image24.vasp"])
def test_real_failed_cells_are_valid_and_not_symmetrized(name):
    atoms = read(FIXTURES / name)
    positions, cell = atoms.positions.copy(), atoms.cell.array.copy()
    report = validate_vasp_image_geometry(atoms, minimum_distance=1.4)
    assert 30 < report["volume_A3"] < 40
    assert report["minimum_distance_A"] > 1.4
    np.testing.assert_array_equal(atoms.positions, positions)
    np.testing.assert_array_equal(atoms.cell.array, cell)


def test_source_relaxation_settings_become_one_explicit_static_policy(tmp_path):
    params, _ = prepare_vasp_static_parameters(source_inputs(tmp_path))
    assert {key: params[key] for key in REQUIRED_VCNEB_STATIC_PARAMETERS} == REQUIRED_VCNEB_STATIC_PARAMETERS
    assert params["isym"] == -1 and params["symprec"] == 1e-4
    np.testing.assert_array_equal(params["kpts"], [8, 8, 6])


@pytest.mark.parametrize("value", [None, 0, -1, float("nan"), float("inf")])
def test_nonfinite_or_implicit_symprec_is_rejected(value):
    with pytest.raises(VaspInputContractError, match="SYMPREC"):
        validate_vasp_static_parameters(dict(REQUIRED_VCNEB_STATIC_PARAMETERS, symprec=value))


def test_noninteger_static_setting_is_not_truncated():
    with pytest.raises(VaspInputContractError, match="NSW"):
        validate_vasp_static_parameters(dict(REQUIRED_VCNEB_STATIC_PARAMETERS, nsw=0.5))


@pytest.mark.parametrize("change", ["parameters", "INCAR", "KPOINTS", "POTCAR", "symbols", "magmoms", "charges"])
def test_changed_contract_stops_before_input_generation(tmp_path, monkeypatch, change):
    atoms = read(FIXTURES / "gan_image24.vasp")
    calc, source = calculator(tmp_path, atoms)
    calls = []
    monkeypatch.setattr(Vasp, "write_input", lambda *args, **kwargs: calls.append(1))
    if change == "parameters":
        calc.set(symprec=1e-8)
    elif change == "symbols":
        atoms[0].symbol = "Al"
    elif change == "magmoms":
        atoms.set_initial_magnetic_moments([1, 0, 0, 0])
    elif change == "charges":
        atoms.set_initial_charges([1, 0, 0, 0])
    else:
        with (source / change).open("a") as handle:
            handle.write("\nchanged\n")
    with pytest.raises(VaspInputContractError):
        calc.write_input(atoms)
    assert not calls


@pytest.mark.parametrize("kind", ["nan", "left_handed", "singular", "ill_conditioned", "collision", "nonperiodic"])
def test_invalid_updated_image_stops_before_launch(kind):
    atoms = Atoms("GaN", positions=[[0, 0, 0], [1, 1, 1]], cell=[4, 4, 4], pbc=True)
    if kind == "nan":
        atoms.positions[0, 0] = np.nan
    elif kind == "left_handed":
        atoms.cell[0, 0] = -4
    elif kind == "singular":
        atoms.cell[0, 0] = 0
    elif kind == "ill_conditioned":
        atoms.cell[0, 0] = 1e-9
    elif kind == "collision":
        atoms.positions[1] = atoms.positions[0]
    else:
        atoms.pbc = False
    with pytest.raises(VaspInputContractError):
        validate_vasp_image_geometry(atoms, minimum_distance=1.4)


def test_serialized_poscar_is_checked_and_provenance_written(tmp_path, monkeypatch):
    atoms = read(FIXTURES / "gan_image24.vasp")
    calc, _ = calculator(tmp_path, atoms)

    def writer(self, atoms, **kwargs):
        self.sort = np.arange(len(atoms))
        self.spinpol = False
        directory = Path(self.directory)
        write(directory / "POSCAR", atoms, format="vasp")
        self.write_incar(atoms, directory=self.directory)
        self.write_kpoints(atoms, directory=self.directory)
        (directory / "POTCAR").write_text("replacement")

    monkeypatch.setattr(Vasp, "write_input", writer)
    calc.write_input(atoms)
    report = json.loads((Path(calc.directory) / "vasp_input_contract.json").read_text())
    assert report["effective_symmetry"] == {"isym": -1, "symprec": 1e-4}
    assert report["structure"]["sha256"] == endpoint_structure_record(atoms)["sha256"]
    assert set(report["input_sha256"]) == {"POSCAR", "INCAR", "KPOINTS", "POTCAR"}
    assert "not_bravais_certification" in report["validation_scope"]
    assert report["poscar_lattice_serialization"] == {
        "policy": "exact_binary64_roundtrip",
        "significant_decimal_digits": POSCAR_LATTICE_SIGNIFICANT_DIGITS,
    }


def test_poscar_lattice_policy_preserves_real_manager_bits(tmp_path):
    # This is image 7 immediately before ASE's old fixed-decimal writer
    # changed its Bravais classification on hf VASP 6.3.2.
    cell = np.array([
        [1.5961020861014259, -2.7648985555405456, -1.990128040149855e-9],
        [1.5961020777947577, 2.7648985507446855, 8.506919139724213e-10],
        [-1.902038257670883e-9, 3.220630587909376e-9, 4.480852661673487],
    ])
    atoms = Atoms("Ga2N2", cell=cell, scaled_positions=np.zeros((4, 3)), pbc=True)
    path = tmp_path / "POSCAR"
    write(path, atoms, format="vasp", direct=True, vasp5=True)
    old = read(path, format="vasp").cell.array
    assert not np.array_equal(old, cell)
    assert np.max(np.abs(old - cell)) < 1e-12
    rewrite_poscar_lattice_exact(path, cell)
    new = read(path, format="vasp").cell.array
    assert np.array_equal(new, cell)
    assert np.array_equal(poscar_lattice_roundtrip(cell), cell)


def test_bravais_diagnostic_takes_precedence_over_mpi_wrapper(tmp_path):
    log = tmp_path / "vasp.out"
    log.write_text("Inconsistent Bravais lattice types found for crystalline and reciprocal lattice\nMPI_ABORT\n")
    assert classify_calculator_failure(RuntimeError("mpirun returned error 1"), diagnostic_paths=[log]) == "vasp_bravais_lattice_inconsistency"
    assert classify_calculator_failure(RuntimeError("mpirun: VASP SCF did not converge")) == "scf_nonconvergence"


@pytest.mark.parametrize("converged, valid", [(False, True), (True, False)])
def test_unconverged_or_nonfinite_vasp_results_are_never_accepted(tmp_path, monkeypatch, converged, valid):
    atoms = read(FIXTURES / "gan_image24.vasp")
    calc, _ = calculator(tmp_path, atoms)

    def fake_calculate(self, atoms, *args):
        self.atoms = atoms.copy()
        self.converged = converged
        self.results = {"energy": -1.0 if valid else np.nan,
                        "forces": np.zeros((len(atoms), 3)), "stress": np.zeros(6)}

    monkeypatch.setattr(Vasp, "calculate", fake_calculate)
    with pytest.raises(CalculationFailed):
        calc.calculate(atoms)
    assert not calc.results


def test_endpoint_cache_cannot_hide_effective_parameter_change(tmp_path):
    source = source_inputs(tmp_path)
    atoms = read(FIXTURES / "gan_image24.vasp")
    params, _ = prepare_vasp_static_parameters(source)
    summary = {
        "status": "completed", "execution_mode": "fixed_initial_endpoint_static_scf",
        "evaluated_image_index": 0, "n_images": 29,
        "endpoint_structures": {"initial": endpoint_structure_record(atoms)},
        "licensed_input_fingerprints": vasp_input_fingerprints(source),
        "calculator_parameters": dict(params, symprec=1e-8),
        "potential_energy_eV": -1.0, "forces_eV_per_A": np.zeros((4, 3)).tolist(),
        "stress_eV_per_A3_voigt": [0.0] * 6,
    }
    path = tmp_path / "summary.json"
    path.write_text(canonical_parameters(summary))
    with pytest.raises(ValueError, match="effective calculator parameters"):
        cached_vasp_static_endpoint_calculator(path, atoms, endpoint="initial", n_images=29,
            source_dir=source, directory=tmp_path / "00", calculator_parameters=params)


def test_exact_state_cache_hit_cannot_hide_parameter_drift(tmp_path):
    atoms = read(FIXTURES / "gan_image24.vasp")
    calc, _ = calculator(tmp_path, atoms)
    atoms.calc = SinglePointCalculator(atoms, energy=-1.0, forces=np.zeros((4, 3)), stress=np.zeros(6))
    executor = ThreadedCalculatorExecutor(1, cache_dir=tmp_path / "cache", cache_namespace="same-run")
    executor.evaluate([atoms], indices=[24])
    atoms.calc = calc
    calc.set(symprec=1e-8)
    with pytest.raises(VaspInputContractError, match="parameters changed"):
        executor.evaluate([atoms], indices=[24])


def test_static_validator_serializes_real_kpoint_array(tmp_path):
    source = source_inputs(tmp_path)
    root = Path(__file__).resolve().parents[1]
    result = subprocess.run([sys.executable, str(root / "scripts" / "validate_vasp_vcneb_static.py"),
                             "--source", str(source), "--output", str(tmp_path / "report.json")],
                            capture_output=True, text=True, check=True)
    report = json.loads(result.stdout)
    assert report["effective_image_parameters"]["kpts"] == [8, 8, 6]
    assert report["effective_image_parameters"]["symprec"] == 1e-4
    assert "not_bravais_certification" in report["validation_scope"]


def test_gan_recovery_gates_whole_saved_chain_before_endpoint_statics():
    root = Path(__file__).resolve().parents[1]
    runner = (root / "cluster" / "cu17_gan_b4_b1_vcneb_resume_isym_minus1.sh").read_text()
    assert "vasp_symprec=${VASP_SYMPREC:-1e-4}" in runner
    assert runner.index("scripts/probe_vasp_input_contract.py") < runner.index("for endpoint in initial final")
    assert '--trajectory "${old_workdir}/vcneb.traj" --n-images 29' in runner
    assert '--ncores 40 --vasp-bin "${vasp_bin}" --run' in runner


def test_gan_core_gates_do_not_confuse_openmp_threads_with_cpu_affinity():
    root = Path(__file__).resolve().parents[1]
    for name in ("cu17_gan_b4_b1_vcneb_serial.sh", "cu17_gan_b4_b1_vcneb_resume_isym_minus1.sh"):
        runner = (root / "cluster" / name).read_text()
        assert "visible_cores=$(env -u OMP_NUM_THREADS -u OMP_THREAD_LIMIT nproc)" in runner
        assert '[[ "${visible_cores}" == 40 ]]' in runner
        assert "export OMP_NUM_THREADS=1" in runner


def test_potcar_species_order_is_a_prelaunch_contract(tmp_path, monkeypatch):
    atoms = read(FIXTURES / "gan_image24.vasp")
    source = source_inputs(tmp_path)
    (source / "POTCAR").write_text(" TITEL = PAW_PBE N 08Apr2002\n TITEL = PAW_PBE Ga_d 06Sep2000\n")
    params, potcar = prepare_vasp_static_parameters(source)
    calc = ExplicitPotcarVasp(source_potcar=potcar, directory=tmp_path / "image", **params)
    calc.lock_input_contract(atoms, source)
    calls = []
    monkeypatch.setattr(Vasp, "write_input", lambda *args, **kwargs: calls.append(1))
    with pytest.raises(VaspInputContractError, match="dataset elements/order"):
        calc.write_input(atoms)
    assert not calls


def test_serial_failure_preserves_success_and_does_not_start_remaining_images(tmp_path):
    calls = []

    class Fake(Calculator):
        implemented_properties = ["energy", "forces", "stress"]

        def __init__(self, index):
            super().__init__()
            self.index = index

        def calculate(self, atoms, *args, **kwargs):
            super().calculate(atoms, *args, **kwargs)
            calls.append(self.index)
            if self.index == 2:
                raise RuntimeError("injected calculator failure")
            self.results = {"energy": -1.0, "forces": np.zeros((len(atoms), 3)), "stress": np.zeros(6)}

    images = [Atoms("Ga", cell=[4, 4, 4], pbc=True, calculator=Fake(index)) for index in [1, 2, 3]]
    executor = ThreadedCalculatorExecutor(1, cache_dir=tmp_path / "cache", cache_namespace="serial",
                                         manifest_path=tmp_path / "manifest.jsonl")
    with pytest.raises(RuntimeError, match="image 2"):
        executor.evaluate(images, indices=[1, 2, 3])
    assert calls == [1, 2]
    assert executor._load_cached(1, images[0]) is not None
    assert executor._load_cached(2, images[1]) is None
    assert executor._load_cached(3, images[2]) is None
    assert json.loads((tmp_path / "manifest.jsonl").read_text())["not_evaluated"] == [3]
