"""Regression tests for the public backend registry and factories."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import sys

from ase import Atoms

from vcneb.backends import (
    attach_image_calculators,
    backend_capability_matrix,
    get_backend_spec,
    make_ase_cp2k_factory,
    make_ase_lammps_factory,
)
from vcneb.calculator import inspect_calculator


ROOT = Path(__file__).resolve().parents[1]


def test_backend_matrix_has_all_supported_adapters() -> None:
    names = [row["name"] for row in backend_capability_matrix()]
    assert names == ["abacus", "vasp", "qe", "lammps", "cp2k"]
    assert get_backend_spec("LAMMPS").variable_cell


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
