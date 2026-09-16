"""Unit tests for the Quantum ESPRESSO VCNEB adapter boundary."""

from __future__ import annotations

from ase import Atoms
import pytest

from vcneb.qe import (
    attach_qe_calculators,
    make_ase_espresso_factory,
    load_approved_qe_pseudopotential_manifest,
    static_qe_input_data,
    validate_qe_pseudopotentials,
)


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


def test_qe_pseudopotential_gate_validates_element_pbe_and_sha256(tmp_path) -> None:
    upf = tmp_path / "Ba-pbe.UPF"
    upf.write_text('<UPF version="2.0.1" element="Ba" functional="PBE">\n', encoding="utf-8")
    report = validate_qe_pseudopotentials(tmp_path, {"Ba": upf.name})
    assert report[0]["species"] == "Ba"
    assert report[0]["functional_marker"] == "PBE"
    assert len(report[0]["sha256"]) == 64


def test_qe_pseudopotential_gate_rejects_missing_or_wrong_metadata(tmp_path) -> None:
    wrong = tmp_path / "wrong.UPF"
    wrong.write_text('<UPF element="Ti" functional="LDA">\n', encoding="utf-8")
    with pytest.raises(ValueError, match="does not match"):
        validate_qe_pseudopotentials(tmp_path, {"Ba": wrong.name})
    with pytest.raises(FileNotFoundError, match="missing"):
        validate_qe_pseudopotentials(tmp_path, {"Ba": "Ba.UPF"})


def test_approved_qe_manifest_pins_filenames_and_md5(tmp_path) -> None:
    manifest = tmp_path / "manifest.json"
    manifest.write_text(
        '{"approval_status":"approved","species":{"Ba":{"filename":"Ba.UPF","md5":"00000000000000000000000000000000"}}}',
        encoding="utf-8",
    )
    loaded = load_approved_qe_pseudopotential_manifest(manifest, {"Ba"})
    assert loaded["pseudopotentials"] == {"Ba": "Ba.UPF"}
    assert loaded["expected_md5"]["Ba"] == "00000000000000000000000000000000"
    manifest.write_text('{"approval_status":"candidate_requires_user_approval","species":{}}', encoding="utf-8")
    with pytest.raises(ValueError, match="approval_status"):
        load_approved_qe_pseudopotential_manifest(manifest, {"Ba"})


def test_qe_pseudopotential_gate_rejects_a_manifest_md5_mismatch(tmp_path) -> None:
    upf = tmp_path / "Ba.UPF"
    upf.write_text('<UPF element="Ba" functional="PBE">\n', encoding="utf-8")
    with pytest.raises(ValueError, match="MD5"):
        validate_qe_pseudopotentials(tmp_path, {"Ba": upf.name}, expected_md5={"Ba": "0" * 32})
