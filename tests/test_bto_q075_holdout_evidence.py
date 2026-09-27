"""The pre-registered BTO holdout's static branch canary is not a PES."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest


BUNDLE = (Path(__file__).resolve().parents[1] / "benchmarks" / "numerical_integrity" /
          "bto_q075_holdout_hf_2026-09-27")


def _load(name: str) -> tuple[dict, str]:
    payload = (BUNDLE / name).read_bytes()
    return json.loads(payload), hashlib.sha256(payload).hexdigest()


def test_holdout_static_branches_have_audited_sources() -> None:
    preflight, preflight_sha = _load(
        "bto_transverse_soft_conditional_q075_q015_preflight_2026-09-27.json")
    canary, canary_sha = _load("canary_starts_result.json")
    audit, _ = _load("audit-q075_q015-canary-27787471.json")
    assert preflight["q_parallel_q_transverse_sqrt_amu_A"] == [0.75, 0.15]
    assert preflight["phonon_supercell"] == [1, 1, 1]
    assert preflight["electronic_kpoints"] == [4, 4, 4]
    assert preflight["no_dft_launched"] is True
    assert [row["label"] for row in preflight["branch_starts"]] == [
        "frozen", "+Q_y", "-Q_y"]
    assert canary["status"] == "three_start_statics_complete_optimizer_not_run"
    assert canary["n_new_evaluations"] == 3
    assert canary["preflight_sha256"] == preflight_sha
    assert str(canary["slurm_job_id"]) == "27787471"
    assert audit["status"] == "all_three_raw_DFT_points_validated"
    assert audit["source_sha256"]["preflight"] == preflight_sha
    assert audit["source_sha256"]["canary_result"] == canary_sha
    assert [row["label"] for row in audit["points"]] == [
        "frozen", "+Q_y", "-Q_y"]
    assert audit["frozen_minus_plus_branch_meV_per_BTO"] == pytest.approx(
        12.846519838149106, abs=1e-8)
    assert abs(audit["signed_branch_energy_difference_meV_per_BTO"]) < 1e-6
