"""Regression gate for the pre-declared BTO internal conditional holdout."""

from __future__ import annotations

import json
from pathlib import Path
import shutil

import pytest

from scripts.analyze_bto_conditional_holdout import EVIDENCE, analyze


def test_holdout_prediction_pass_is_not_a_continuous_pes_certificate() -> None:
    report = analyze()
    assert report["single_holdout_energy_prediction_passed"] is True
    assert abs(report["low_prediction_error_meV_per_BTO"]) < 2.0
    assert report["observed_low_minus_zero_gap_meV_per_BTO"] == pytest.approx(43.8404470010)
    assert set(report["branches_from_unique_audited_cache_points"]) == {"Q_y=0", "+Q_y", "-Q_y"}
    assert report["n_audited_DFT_points"] == 77
    assert report["n_new_DFT_points"] == 74
    assert report["source_branch_terminals_explicit"] is False
    assert report["branch_terminal_mapping_replayed_from_audited_cache"] is True
    assert [report["branch_replay"][key]["final_evaluation_directory"].split("-")[1]
            for key in ("frozen", "+Q_y", "-Q_y")] == ["000018", "000047", "000076"]
    assert report["curvature_checked_at_holdout"] is False
    assert report["continuous_conditioned_PES_certified"] is False


def test_holdout_rejects_a_changed_result_file(tmp_path: Path) -> None:
    for source in EVIDENCE.glob("*.json"):
        if source.name == "holdout_analysis.json":
            continue
        shutil.copy2(source, tmp_path / source.name)
    path = tmp_path / "conditional_q075_q015_result.json"
    changed = json.loads(path.read_text(encoding="utf-8"))
    changed["energy_minus_c_eV_per_BTO"] -= 0.001
    path.write_text(json.dumps(changed), encoding="utf-8")
    with pytest.raises(ValueError, match="source hashes"):
        analyze(evidence=tmp_path)
