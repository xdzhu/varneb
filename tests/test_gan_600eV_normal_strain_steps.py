"""Guard the bounded 600-eV GaN normal-strain finite-difference experiment."""

import json
from pathlib import Path

import pytest

from scripts.prepare_gan_600eV_normal_strain_steps import (
    DEFAULT_AUDIT, DEFAULT_SOURCE, INPUT_SHA256, prepare, sha256,
)
import scripts.audit_gan_600eV_normal_strain_steps as step_audit


def test_preparation_keeps_original_electronic_inputs(tmp_path):
    if not (DEFAULT_SOURCE / "POTCAR").is_file():
        pytest.skip("licensed VASP POTCAR is not redistributed in a clean checkout")
    root = tmp_path / "normal-strain"
    report = prepare(DEFAULT_SOURCE, DEFAULT_AUDIT, root)
    assert report["n_cases"] == 12
    assert report["steps_A"] == [0.01, 0.005]
    assert report["axes"] == [12, 13, 14]
    assert len({item["name"] for item in report["cases"]}) == 12
    assert all(not item["empirical_near_symmetry_warning"]
               and item["minimum_distance_A"] > 1.4
               for item in report["cases"])
    for item in report["cases"]:
        case = root / "cases" / item["name"]
        assert all(sha256(case / key) == value
                   for key, value in item["input_sha256"].items())
        assert all(item["input_sha256"][key] == digest
                   for key, digest in INPUT_SHA256.items())
    with pytest.raises(FileExistsError):
        prepare(DEFAULT_SOURCE, DEFAULT_AUDIT, root)


def test_archived_normal_strain_outputs_keep_same_contract_and_nonzero_residual():
    report = step_audit.audit()
    archived = json.loads((step_audit.ROOT / "benchmarks/numerical_integrity/"
                           "gan_600eV_normal_strain_step_dependence_20261001.json")
                          .read_text(encoding="utf-8"))
    assert archived["source_sha256"]["auditor"] == sha256(Path(step_audit.__file__))
    assert archived["source_sha256"]["manifest"] == sha256(
        step_audit.DEFAULT_WORK / "manifest.json")
    assert archived["comparison_to_prior_0p02_A"] == report["comparison_to_prior_0p02_A"]
    assert report["n_new_statics"] == 12
    assert report["electronic_input_sha256"] == INPUT_SHA256
    assert not report["POTCAR_locally_reverified"]
    assert len(report["rows"]) == 6
    for row in report["rows"]:
        assert 0.015 < abs(row["secant_minus_simpson_eV_per_A"]) < 0.03
        assert row["kpoint_plane_wave_counts"]["n_kpoints"] == 384
        assert row["kpoint_plane_wave_counts"]["n_changed_minus_vs_plus"] > 300
        assert row["fft_grids_unchanged"]
    assert all(len(item["residual_eV_per_A"]) == 3
               for item in report["comparison_to_prior_0p02_A"])
