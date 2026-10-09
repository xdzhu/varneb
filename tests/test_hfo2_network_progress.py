"""Current measured progress must not overwrite historical evidence or certify G1."""

import copy
import json
from pathlib import Path

import pytest

pytest.importorskip("pymatgen")

import scripts.analyze_hfo2_network_progress as progress
from scripts.audit_hfo2_static_replica import sha256


ROOT = Path(__file__).resolve().parents[1]
OLD = ROOT / "benchmarks/hfo2_channels/20261008/network_update_20261009/network_specification.json"
NEW = ROOT / "benchmarks/hfo2_channels/20261008/morning_update_20261009_0850/network_specification.json"


def test_original_stage_is_still_one_pass_without_changing_its_physical_record():
    original = progress.analyze_case(OLD)
    result = progress.analyze_progress(OLD)
    assert result["channels"] == original["channels"]
    assert result["ordinary_residual_passed_sources"] == 1
    assert result["ordinary_residual_passed_channels"] == ["PO_to_T"]
    assert result["physical_analysis_source_sha256"] == original["analysis_source_sha256"]
    assert result["analysis_source_sha256"] == sha256(Path(progress.__file__))
    assert result["historical_caption_replaced"] in original["limitations"]
    assert result["historical_caption_replaced"] not in result["limitations"]
    assert not result["full_G1_passed"] and not result["TS_certified"]


def test_actual_morning_stage_recognises_both_measured_passes_not_final_gates():
    result = progress.analyze_progress(NEW)
    assert result["existing_image_records"] == 37
    assert result["ordinary_residual_passed_sources"] == 2
    assert result["ordinary_residual_passed_channels"] == ["PO_to_T", "PO_to_M"]
    assert result["unconverged_channels"] == ["PO_flip_T_pattern_preserving", "PO_flip_T_pattern_reversing"]
    m = result["channels"][1]
    assert m["source_job_id"] == "28298794" and m["snapshot_step"] == 39
    assert m["source_NEB_fmax_eV_A"] == pytest.approx(.09790448957151193)
    assert m["peak_energy_relative_PO_meV_fu"] == pytest.approx(71.58220716746655)
    assert not result["full_G1_passed"] and not result["TS_certified"]
    assert result["new_DFT_calls"] == 0 and not result["physical_parameters_changed"]
    terminal = json.loads((NEW.parent / "PO_M_terminal_summary.json").read_text())
    assert terminal["converged"] and terminal["termination"] == "force_threshold"
    assert terminal["saddle_diagnostics"]["true_tangential_force_eV_per_A"] == pytest.approx(-.3834950037218957)


def test_even_all_ordinary_passes_do_not_promote_unmeasured_scientific_gates(monkeypatch):
    original = progress.analyze_case(OLD)
    altered = copy.deepcopy(original)
    for candidate in altered["channels"]:
        candidate["ordinary_residual_passed"] = True
    before = copy.deepcopy(altered)
    monkeypatch.setattr(progress, "analyze_case", lambda _: copy.deepcopy(altered))
    result = progress.analyze_progress(OLD)
    assert result["ordinary_residual_passed_sources"] == 4
    assert result["unconverged_channels"] == []
    assert not result["full_G1_passed"] and not result["TS_certified"]
    assert altered == before


def test_unknown_historical_schema_fails_closed(monkeypatch):
    result = progress.analyze_case(OLD)
    result["limitations"] = []
    monkeypatch.setattr(progress, "analyze_case", lambda _: result)
    with pytest.raises(ValueError, match="historical caption"):
        progress.analyze_progress(OLD)


def test_existing_output_refused_before_any_input_read(tmp_path, monkeypatch):
    output = tmp_path / "keep.json"
    output.write_text("original bytes")
    monkeypatch.setattr("sys.argv", ["progress", "--specification", "missing", "--output", str(output)])
    with pytest.raises(FileExistsError):
        progress.main()
    assert output.read_text() == "original bytes"
