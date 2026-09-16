"""Unit tests for the Quantum ESPRESSO VCNEB adapter boundary."""

from __future__ import annotations

from ase import Atoms
import pytest

from vcneb.qe import attach_qe_calculators, make_ase_espresso_factory, static_qe_input_data


def test_static_qe_input_data_enforces_scf_force_and_stress() -> None:
    data = static_qe_input_data({"system": {"ecutwfc": 100}, "control": {"prefix": "bto"}})
    assert data["system"] == {"ecutwfc": 100}
    assert data["control"] == {
        "prefix": "bto",
        "calculation": "scf",
        "tstress": True,
        "tprnfor": True,
    }


@pytest.mark.parametrize(
    "input_data",
    [
        {"control": {"calculation": "vc-relax"}},
        {"calculation": "relax"},
        {"control": {"tstress": False}},
        {"control": {"tprnfor": False}},
    ],
)
def test_static_qe_input_data_rejects_nonstatic_or_incomplete_images(input_data) -> None:
    with pytest.raises(ValueError):
        static_qe_input_data(input_data)


def test_qe_factory_and_image_directories(monkeypatch, tmp_path) -> None:
    import ase.calculators.espresso as espresso

    captured: list[dict] = []

    class FakeProfile:
        def __init__(self, command, pseudo_dir):
            self.command = command
            self.pseudo_dir = pseudo_dir

    class FakeEspresso:
        implemented_properties = ("energy", "forces", "stress")

        def __init__(self, *, profile, directory, **parameters):
            captured.append({"profile": profile, "directory": directory, "parameters": parameters})
            self.profile = profile
            self.directory = directory

    monkeypatch.setattr(espresso, "Espresso", FakeEspresso)
    monkeypatch.setattr(espresso, "EspressoProfile", FakeProfile)
    factory = make_ase_espresso_factory(
        parameters={"pseudopotentials": {"Ba": "Ba.upf"}, "input_data": {"system": {"ecutwfc": 100}}},
        command="srun pw.x",
        pseudo_dir="/pseudo",
    )
    images = [Atoms("Ba", cell=[4, 4, 4], pbc=True), Atoms("Ba", cell=[4, 4, 4], pbc=True)]
    attach_qe_calculators(images, workdir=tmp_path, factory=factory)

    assert [entry["directory"] for entry in captured] == [str(tmp_path / "00"), str(tmp_path / "01")]
    assert all((tmp_path / f"{index:02d}" / "POSCAR.start").is_file() for index in range(2))
    assert captured[0]["profile"].command == "srun pw.x"
    assert captured[0]["parameters"]["input_data"]["control"]["calculation"] == "scf"


def test_qe_factory_requires_one_explicit_launch_path() -> None:
    with pytest.raises(ValueError, match="requires either profile"):
        make_ase_espresso_factory(parameters={})
    with pytest.raises(ValueError, match="either command or profile"):
        make_ase_espresso_factory(parameters={}, command="pw.x", pseudo_dir="/pseudo", profile=object())
