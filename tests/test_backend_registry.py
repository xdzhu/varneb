"""Regression tests for the public backend registry and factories."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

from ase import Atoms
import pytest

from vcneb.backends import (
    attach_image_calculators,
    backend_capability_matrix,
    get_backend_spec,
    make_ase_cp2k_factory,
    make_ase_abinit_factory,
    make_ase_calculator_factory,
    make_ase_lammps_factory,
)
from vcneb.abacus import _minimal_abacus_results, _read_vcneb_results
from vcneb.calculator import classify_calculator_failure, inspect_calculator
from vcneb.config import RunConfig
from vcneb.optimizer_registry import get_optimizer_spec, optimizer_capability_matrix


ROOT = Path(__file__).resolve().parents[1]


def test_backend_matrix_has_all_supported_adapters() -> None:
    names = [row["name"] for row in backend_capability_matrix()]
    assert names == ["abacus", "vasp", "qe", "lammps", "cp2k", "abinit"]
    assert get_backend_spec("LAMMPS").variable_cell


def test_optimizer_registry_is_independent_from_backend_registry() -> None:
    assert get_optimizer_spec("split-fire").name == "SplitFIRE"
    rows = optimizer_capability_matrix()
    assert rows
    assert {row["backend_independent"] for row in rows} == {"true"}


def test_run_config_accepts_cross_backend_optimizer_combinations(tmp_path) -> None:
    common = {
        "initial": tmp_path / "initial.vasp",
        "final": tmp_path / "final.vasp",
        "workdir": tmp_path / "run",
    }
    vasp_split_fire = RunConfig(backend="vasp", optimizer="SplitFIRE", **common)
    abacus_bfgs = RunConfig(backend="abacus", optimizer="BFGS", **common)
    assert (vasp_split_fire.backend, vasp_split_fire.optimizer) == ("vasp", "SplitFIRE")
    assert (abacus_bfgs.backend, abacus_bfgs.optimizer) == ("abacus", "BFGS")


def test_abacus_reader_ignores_optional_eigenvalue_parser(monkeypatch, tmp_path) -> None:
    abacus_io = pytest.importorskip("ase.io.abacus")

    class Chunk:
        energy = -1.25
        free_energy = -1.25
        forces_sort = [[0.0, 0.0, 0.0]] * 4
        stress = [[0.0, 0.0, 0.0]] * 3
        magmom = None
        dipole = None

    output = tmp_path / "OUT.ABACUS"
    output.mkdir()
    (output / "running_scf.log").write_text("minimal", encoding="utf-8")
    monkeypatch.setattr(abacus_io, "_get_abacus_chunks", lambda *args, **kwargs: [Chunk()])
    result = _read_vcneb_results(tmp_path, output_suffix="ABACUS", calculation="scf")
    assert result["energy"] == -1.25
    assert result["forces"] == [[0.0, 0.0, 0.0]] * 4
    assert result["stress"] == [[0.0, 0.0, 0.0]] * 3


def test_abacus_reader_rejects_singleton_mpi_fallback(tmp_path) -> None:
    output = tmp_path / "OUT.ABACUS"
    output.mkdir()
    (output / "running_scf.log").write_text("apparently valid", encoding="utf-8")
    (tmp_path / "abacus.out").write_text(
        "MPI startup(): PMI server not found\n", encoding="utf-8"
    )
    with pytest.raises(RuntimeError, match="PMI server not found"):
        _read_vcneb_results(tmp_path, output_suffix="ABACUS", calculation="scf")


def test_abacus_minimal_contract_parser_handles_usable_damaged_log(tmp_path) -> None:
    output = tmp_path / "running_scf.log"
    output.write_text(
        """TOTAL-FORCE (eV/Angstrom)
Ga1 0.1 0.2 0.3
N1 -0.1 -0.2 -0.3

TOTAL-STRESS (KBAR)
x 1.0 2.0 3.0
y 4.0 5.0 6.0
z 7.0 8.0 9.0

final etot is -1.25 eV
""",
        encoding="utf-8",
    )
    result = _minimal_abacus_results(output)
    assert result["forces"].shape == (2, 3)
    assert result["stress"].shape == (6,)
    assert result["energy"] == -1.25


def test_abacus_minimal_contract_parser_accepts_final_etot_marker(tmp_path) -> None:
    output = tmp_path / "running_scf.log"
    output.write_text(
        "TOTAL-FORCE (eV/Angstrom)\n"
        "Ga1 0.1 0.2 0.3\n"
        "N1 -0.1 -0.2 -0.3\n\n"
        "TOTAL-STRESS (KBAR)\n"
        "x 1.0 2.0 3.0\n"
        "y 4.0 5.0 6.0\n"
        "z 7.0 8.0 9.0\n\n"
        "!FINAL_ETOT_IS -4726.5466304132096411 eV\n",
        encoding="utf-8",
    )
    result = _minimal_abacus_results(output)
    assert result["energy"] == pytest.approx(-4726.54663041321)


def test_lammps_factory_is_explicit_and_isolated(tmp_path) -> None:
    factory = make_ase_lammps_factory(
        parameters={"pair_style": "zero 10.0", "pair_coeff": ["* *"]},
        command="lmp_mpi",
    )
    atoms = Atoms("Ar", cell=[10, 10, 10], pbc=True)
    image_dir = tmp_path / "image_0001"
    image_dir.mkdir()
    calculator = factory(1, atoms, image_dir)
    assert calculator.parameters["command"] == "lmp_mpi"
    assert Path(calculator.parameters["tmp_dir"]) == image_dir
    report = inspect_calculator(calculator)
    assert report.directory == str(image_dir)
    assert report.command == "lmp_mpi -echo log -screen none -log /dev/stdout"
    calculator.clean()


def test_cp2k_factory_keeps_output_inside_image(monkeypatch, tmp_path) -> None:
    import ase.calculators.cp2k as cp2k

    captured = {}

    class FakeCP2K:
        def __init__(self, **kwargs):
            captured.update(kwargs)

    monkeypatch.setattr(cp2k, "CP2K", FakeCP2K)
    factory = make_ase_cp2k_factory(parameters={"cutoff": 300}, command="cp2k_shell")
    image_dir = tmp_path / "image_0002"
    image_dir.mkdir()
    factory(2, Atoms("H", cell=[5, 5, 5], pbc=True), image_dir)
    assert captured["label"].endswith("cp2k")
    assert captured["command"] == "cp2k_shell"
    assert captured["stress_tensor"] is True


def test_cp2k_factory_converts_explicit_rydberg_cutoff(monkeypatch, tmp_path) -> None:
    import ase.calculators.cp2k as cp2k
    from ase.units import Rydberg

    captured = {}

    class FakeCP2K:
        def __init__(self, **kwargs):
            captured.update(kwargs)

    monkeypatch.setattr(cp2k, "CP2K", FakeCP2K)
    factory = make_ase_cp2k_factory(parameters={"cutoff_ry": 400}, command="cp2k_shell")
    image_dir = tmp_path / "image_0003"
    image_dir.mkdir()
    factory(3, Atoms("H", cell=[5, 5, 5], pbc=True), image_dir)
    assert captured["cutoff"] == pytest.approx(400 * Rydberg)


def test_cp2k_factory_rejects_ambiguous_cutoff_units(monkeypatch) -> None:
    import ase.calculators.cp2k as cp2k

    monkeypatch.setattr(cp2k, "CP2K", lambda **kwargs: None)
    with pytest.raises(ValueError, match="both cutoff and cutoff_ry"):
        make_ase_cp2k_factory(parameters={"cutoff": 400, "cutoff_ry": 400})


def test_abinit_factory_uses_profile_and_image_directory(monkeypatch, tmp_path) -> None:
    import ase.calculators.abinit as abinit

    captured = {}

    class FakeProfile:
        def __init__(self, command, *, pp_paths=None):
            captured["command"] = command
            captured["pp_paths"] = pp_paths

    class FakeAbinit:
        def __init__(self, **kwargs):
            captured["kwargs"] = kwargs

    monkeypatch.setattr(abinit, "AbinitProfile", FakeProfile)
    monkeypatch.setattr(abinit, "Abinit", FakeAbinit)
    (tmp_path / "pseudo").mkdir()
    factory = make_ase_abinit_factory(
        parameters={"ecut": 20, "toldfe": 1.0e-5},
        command="abinit",
        pp_paths=tmp_path / "pseudo",
    )
    image_dir = tmp_path / "image_0001"
    image_dir.mkdir()
    factory(1, Atoms("H", cell=[5, 5, 5], pbc=True), image_dir)
    assert captured["command"].strip("'").endswith("varneb_abinit_runner.sh")
    assert captured["pp_paths"] == [str(tmp_path / "pseudo")]
    assert captured["kwargs"]["directory"] == str(image_dir)
    assert captured["kwargs"]["nsym"] == 1


def test_abinit_factory_rejects_missing_pseudopotential_contract(tmp_path) -> None:
    with pytest.raises(ValueError, match="require explicit pp_paths"):
        make_ase_abinit_factory(
            parameters={"pps": "hgh"}, command="abinit"
        )
    with pytest.raises(FileNotFoundError, match="pseudopotential directories"):
        make_ase_abinit_factory(
            parameters={"pps": "hgh"},
            command="abinit",
            pp_paths=tmp_path / "missing-pseudo",
        )


@pytest.mark.parametrize(
    ("message", "category"),
    [
        ("ABINIT chkorthsy: inconsistent lattice vectors; MPI_ABORT", "abinit_symmetry_failure"),
        ("ABINIT requires pp_paths: pseudopotential directory missing", "pseudopotential_contract"),
        ("ABACUS result parsing failed: no parseable final stress", "output_contract_violation"),
    ],
)
def test_backend_contract_failures_are_not_misclassified_as_mpi(message, category) -> None:
    assert classify_calculator_failure(RuntimeError(message)) == category


def test_generic_ase_factory_injects_private_directory_and_command(tmp_path) -> None:
    captured = {}

    class FakeCalculator:
        def __init__(self, *, directory, command, cutoff):
            captured.update(directory=directory, command=command, cutoff=cutoff)

    factory = make_ase_calculator_factory(
        FakeCalculator,
        parameters={"cutoff": 400},
        command="srun fake-code",
    )
    image_dir = tmp_path / "image_0003"
    image_dir.mkdir()
    factory(3, Atoms("H", cell=[5, 5, 5], pbc=True), image_dir)
    assert captured == {
        "directory": str(image_dir),
        "command": "srun fake-code",
        "cutoff": 400,
    }


def test_attach_image_calculators_uses_one_directory_per_image(tmp_path) -> None:
    seen = []

    def factory(index, image, directory):
        seen.append((index, directory))
        return object()

    images = [Atoms("H", cell=[4, 4, 4], pbc=True) for _ in range(3)]
    attach_image_calculators(images, workdir=tmp_path, factory=factory)
    assert [index for index, _ in seen] == [0, 1, 2]
    assert all(directory.name == f"image_{index:04d}" for index, directory in seen)
    assert all((directory / "structure.start.vasp").is_file() for _, directory in seen)


def test_cli_backend_and_init_commands(tmp_path) -> None:
    result = subprocess.run(
        [sys.executable, "-m", "vcneb", "backends", "--json"],
        cwd=ROOT, capture_output=True, text=True, check=True,
    )
    assert {row["name"] for row in json.loads(result.stdout)} >= {"lammps", "cp2k"}
    target = tmp_path / "case" / "varneb.json"
    result = subprocess.run(
        [sys.executable, "-m", "vcneb", "init", str(target)],
        cwd=ROOT, capture_output=True, text=True, check=True,
    )
    assert json.loads(target.read_text(encoding="utf-8"))["fmax_ev_per_angstrom"] == 0.10
    result = subprocess.run(
        [sys.executable, "-m", "vcneb", "validate-config", str(target)],
        cwd=ROOT, capture_output=True, text=True, check=True,
    )
    assert "valid VARNEB config" in result.stdout
    result = subprocess.run(
        [sys.executable, "-m", "vcneb", "doctor", "--backend", "cp2k", "--json"],
        cwd=ROOT, capture_output=True, text=True, check=True,
    )
    doctor = json.loads(result.stdout)[0]
    assert "cp2k_shell.psmp" in doctor["executable_candidates"]
