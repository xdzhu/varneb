"""The published launch counts are reproducible from archived sacct rows."""

from pathlib import Path

import pytest

from scripts.audit_acceleration_slurm_accounting import audit, parse_accounting


RAW = (Path(__file__).resolve().parents[1] / "benchmarks" / "convergence"
       / "hf_slurm_abacus_steps_20260928.psv")


def test_archived_slurm_steps_and_cutoff_boundaries():
    result = audit(RAW)
    assert result["source_sha256"]["sacct_raw_psv"] == (
        "2a6f31d279e301fbee2d4b920877ef5e43f771471b5eff2678f6d82cc8c7ce3a"
    )
    routes = result["routes"]
    assert routes["bto_global"]["completed_launches_by_cutoff_second"] == 205
    assert routes["bto_block001"]["completed_launches_by_cutoff_second"] == 44
    assert routes["hfo2_global"]["completed_launches_by_cutoff_second"] == 217
    assert routes["hfo2_staged_coarse"]["completed_launches_by_cutoff_second"] == 122
    assert routes["hfo2_staged_refine"]["completed_launches_by_cutoff_second"] == 12
    assert routes["bto_global"]["launched_in_same_cutoff_second"] == 4
    assert routes["hfo2_staged_coarse"]["launched_in_same_cutoff_second"] == 5
    assert result["comparison"]["bto_block001_saved_percent"] == pytest.approx(78.5365853659)
    assert result["comparison"]["hfo2_staged_saved_percent"] == pytest.approx(38.2488479263)


def test_missing_step_and_failed_step_are_rejected():
    raw = RAW.read_text(encoding="utf-8")
    rows = [line for line in raw.splitlines() if not line.startswith("27729028.10|")]
    with pytest.raises(ValueError, match="missing or repeated"):
        parse_accounting("\n".join(rows))
    with pytest.raises(ValueError, match="non-successful"):
        parse_accounting(raw.replace("27729028.10|abacus|COMPLETED|",
                                     "27729028.10|abacus|FAILED|", 1))
