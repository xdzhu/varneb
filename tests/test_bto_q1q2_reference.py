"""Cross-platform provenance checks for the BTO example loader."""

from __future__ import annotations

import json

import pytest

from examples.bto_q1q2_reference import (
    _require_report_hash, sha256, validate_bto_gamma_source,
)


def test_report_hash_keys_work_after_windows_to_linux_copy(tmp_path) -> None:
    source = tmp_path / "cubic_CONTCAR"
    source.write_text("audited structure\n", encoding="utf-8")
    digest = sha256(source)
    for reported in (r"outputs\bto\cubic_CONTCAR", "outputs/bto/cubic_CONTCAR"):
        assert _require_report_hash({"input_sha256": {reported: digest}}, source) == digest
    with pytest.raises(ValueError, match="mismatched"):
        _require_report_hash({"input_sha256": {r"outputs\bto\cubic_CONTCAR": "bad"}}, source)


def test_bto_gamma_source_keeps_phonon_supercell_separate_from_kmesh(tmp_path) -> None:
    force_sets = tmp_path / "FORCE_SETS"
    force_sets.write_text("5\n6\n", encoding="utf-8")
    gamma_provenance = tmp_path / "gamma.json"
    eigen_provenance = tmp_path / "eigen.json"
    gamma = {
        "endpoint": "cubic BaTiO3",
        "supercell": [1, 1, 1],
        "n_signed_displacements": 6,
        "finite_difference_amplitude_A": 0.01,
        "calculator": {
            "code": "ABACUS", "functional": "PBE", "ecutwfc_Ry": 100,
            "basis": "full Ba/Ti/O DZP 10 au", "kpoints": [4, 4, 4],
        },
        "force_constants": {
            "raw_force_sets": force_sets.name,
            "phonopy_gamma_eigenpairs_provenance": eigen_provenance.name,
        },
    }
    eigen = {
        "force_sets_sha256": sha256(force_sets),
        "method": "Phonopy run_qpoints([[0, 0, 0]], with_eigenvectors=True)",
        "n_modes": 15,
    }
    gamma_provenance.write_text(json.dumps(gamma), encoding="utf-8")
    eigen_provenance.write_text(json.dumps(eigen), encoding="utf-8")
    hashes = validate_bto_gamma_source(gamma_provenance, force_sets, eigen_provenance)
    assert hashes["force_sets"] == sha256(force_sets)

    gamma["supercell"] = [2, 2, 2]
    gamma_provenance.write_text(json.dumps(gamma), encoding="utf-8")
    with pytest.raises(ValueError, match="1x1x1"):
        validate_bto_gamma_source(gamma_provenance, force_sets, eigen_provenance)

    gamma["supercell"] = [1, 1, 1]
    gamma["calculator"]["kpoints"] = [1, 1, 1]
    gamma_provenance.write_text(json.dumps(gamma), encoding="utf-8")
    with pytest.raises(ValueError, match="1x1x1"):
        validate_bto_gamma_source(gamma_provenance, force_sets, eigen_provenance)
