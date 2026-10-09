"""Actual frozen-source review plus protocol refusal tests; no new DFT."""
import copy
from pathlib import Path

import pytest

from scripts.audit_hfo2_G1_gate import audit, validate_topology
from scripts.plot_hfo2_G1_terminal_network import build_data

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def actual():
    return build_data(ROOT)


def test_actual_review_has_four_candidate_edges_not_a_TS_or_winding_certificate():
    r = audit(ROOT)
    assert r["candidate_G1_gate_passed"] and r["new_DFT_calls"] == 0
    assert r["existing_terminal_records"] == 40 and r["fresh_raw_terminal_records"] == 0
    assert r["sampling"]["peak_SCF_calls"] == 13
    assert r["sampling"]["preserving_union_peak_meV_fu"] == pytest.approx(38.50844808448528)
    assert not r["TS_global_escape_stability_spontaneous_P_or_distinct_winding_certified"]
    assert not r["G2_matrix_or_holdout_automatically_submitted"]
    for b in r["terminal_barriers"]:
        assert b["terminal_sampled_forward_meV_fu"]-b["terminal_sampled_reverse_meV_fu"] == pytest.approx(b["endpoint_difference_meV_fu"])


@pytest.mark.parametrize("mutation", ["missing", "unconverged", "unknown_source", "phase", "nonfinite"])
def test_candidate_protocol_refuses_changed_or_incomplete_terminal_evidence(actual, mutation):
    channels = copy.deepcopy(actual["channels"])
    if mutation == "missing":
        channels.pop()
    elif mutation == "unconverged":
        channels[0]["source_NEB_fmax_eV_A"] = .101
    elif mutation == "unknown_source":
        channels[0]["source_observation_sha256"] = "0"*64
    elif mutation == "phase":
        channels[0]["rows"][-1]["structure_audit"]["symmetry_sweep"][0]["symbol"] = "P1"
    else:
        channels[0]["rows"][1]["Q_pattern_x_A"] = float("nan")
    with pytest.raises(ValueError):
        validate_topology(channels)
