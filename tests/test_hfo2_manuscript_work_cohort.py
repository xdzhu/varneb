"""The actual six-band Table 7 cohort, without author arrays or new DFT."""
from pathlib import Path

import pytest

from scripts.audit_hfo2_manuscript_work_cohort import audit_cohort, audit_manuscript


@pytest.fixture(scope="module")
def cohort():
    return audit_cohort()


def manuscript():
    return (Path(__file__).resolve().parents[1]/"paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md").read_text(encoding="utf-8")


def test_actual_native_counts_and_reuse_are_not_statistical_replicates(cohort):
    assert (cohort["bands"], cohort["frame_records"], cohort["distinct_raw_log_STRU_identities"]) == (6, 54, 50)
    assert sorted(len(row["occurrences"]) for row in cohort["reused_identities"]) == [2, 4]
    assert all(index == 0 for row in cohort["reused_identities"] for _, index in row["occurrences"])
    assert cohort["maximum_native_dual_work_defect_eV_cell"] < 1.2e-16
    assert cohort["new_DFT_calls"] == cohort["scheduler_mutations"] == 0
    assert cohort["holdout_or_G3_access"] is False
    assert cohort["full_native_runtime_or_assets_reaudited_locally"] is False
    assert cohort["physical_error_bound_or_TS_certified"] is False


def test_actual_manuscript_has_every_table_and_component_row(cohort):
    assert audit_manuscript(manuscript(), cohort)


@pytest.mark.parametrize("old,new", [
    ("six completed clamped training bands", "five completed clamped training bands"),
    ("54 frame records", "45 frame records"),
    ("50 distinct raw-log/geometry", "42 distinct raw-log/geometry"),
    ("| 97.227719 | 75.922579 |", "| 97.227719 | 75.922580 |"),
    ("67.611606; open-cell", "67.611605; open-cell"),
])
def test_missing_cohort_and_altered_data_fail_closed(cohort, old, new):
    text = manuscript()
    assert old in text
    with pytest.raises(ValueError, match="stale|differs|omit"):
        audit_manuscript(text.replace(old, new), cohort)
