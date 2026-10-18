"""Source-derived prose and preservation checks, not independent material tests."""
from collections import Counter
import hashlib, json, re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper/VARNEB_JCTC"
CASE = ROOT / "benchmarks/hfo2_channels/20261008/G3_sheet_consistency_E172_20261015"


def result_prose(text):
    return text.split("The subsequently completed T candidate sheet", 1)[1].split(
        "### E.5.", 1
    )[0]


def require_source_numbers(text, data):
    s = data["summary"]
    p = result_prose(text)
    numbers = [
        f"{s['total_radius_range_A'][1]:.5f}",
        *[f"{v:.5f}" for v in s["released_gradient_range_eV_A"]],
        f"{s['relative_to_released_origin_energy_range_meV_cell'][0]:.3f}",
        f"{s['largest_retained_only_energy_work_difference']['retained_only_energy_minus_work_meV_cell']:.3f}",
        f"{s['largest_retained_only_energy_work_difference']['released_endpoint_work_meV_cell']:.3f}",
        f"{s['largest_retained_only_energy_work_difference']['full_energy_minus_trapezoid_work_meV_cell']:.3f}",
        f"{s['largest_full_endpoint_quadrature_difference']['full_energy_minus_trapezoid_work_meV_cell']:.3f}",
        f"{abs(s['largest_retained_loop_circulation']['retained_endpoint_circulation_meV_cell']):.3f}",
        f"{abs(s['largest_full_loop_circulation']['full_endpoint_circulation_meV_cell']):.3f}",
    ]
    if any(p.count(n) != 1 for n in numbers):
        raise ValueError(
            "source-derived diagnostic number missing, changed or duplicated"
        )


def test_manuscript_numbers_match_all_observed_diagnostics_and_units():
    data = json.loads((CASE / "analysis.json").read_bytes())
    text = (PAPER / "MANUSCRIPT_DRAFT.md").read_text(encoding="utf-8")
    require_source_numbers(text, data)
    p = " ".join(result_prose(text).split())
    assert data["formula_units"] == 4
    for phrase in (
        "25 points required 118 new SCFs",
        "40 nearest grid edges",
        "16 elementary closed loops",
        "meV/cell",
        "different loops",
        "not DFT, interpolation or barrier-error bounds",
        "not an end-to-end MEP or certified stable conditional minima",
    ):
        assert phrase in p


@pytest.mark.parametrize("old", ["1.209", "1.094", "0.115", "1.039", "2.000", "0.578"])
def test_changed_or_duplicated_numbers_cannot_silently_enter_prose(old):
    data = json.loads((CASE / "analysis.json").read_bytes())
    text = (PAPER / "MANUSCRIPT_DRAFT.md").read_text(encoding="utf-8")
    for bad in (text.replace(old, "999.999"), text.replace(old, old + " " + old)):
        with pytest.raises(ValueError, match="source-derived"):
            require_source_numbers(bad, data)


def test_old_113_table_lines_nine_figures_citations_and_wrapper_are_retained():
    before = json.loads((CASE / "baseline_preservation.json").read_bytes())
    text = (PAPER / "MANUSCRIPT_DRAFT.md").read_text(encoding="utf-8")
    counts = Counter(
        hashlib.sha256(line.encode()).hexdigest()
        for line in text.splitlines()
        if line.startswith("|")
    )
    assert before["baseline_commit"] == "6472844bf4d7689a5106c1e3a5c6d930b03c75d0"
    assert before["original_table_line_count"] == 113
    assert not Counter(before["original_table_hash_counts"])-counts
    # Original119 plus shared-reference(16), nonlinear(7), atomic(6), statics(6).
    # The completed P8 comparison adds6main and25appendix lines, no old loss.
    assert sum(counts.values())==218
    assert 'Table A21. Prospective T local errors' in text and 'Table A18. All frozen T test arms' in text
    assert 'Table A17. Four frozen directional training statics' in text
    figures = re.findall(r"!\[[^\]]*\]\(([^)]+)\)", text)
    assert not Counter(before["original_figure_targets"])-Counter(figures)
    assert len(figures)==13 and len(set(figures))==13
    assert all(
        hashlib.sha256((PAPER / n).read_bytes()).hexdigest() == h
        for n, h in before["original_figure_sha256"].items()
    )
    assert (
        re.findall(r"https?://[^\s)]+", text) == before["original_external_citations"]
    )
    assert (
        hashlib.sha256(
            (PAPER / "varneb_JCTC_draft.tex").read_text().encode()
        ).hexdigest()
        == before["tex_wrapper_sha256"]
    )


def test_methods_ledger_and_interpretation_keep_finite_residual_and_budget_boundaries():
    main = (PAPER / "MANUSCRIPT_DRAFT.md").read_text(encoding="utf-8")
    methods = (PAPER / "METHODS_DRAFT.md").read_text(encoding="utf-8")
    ledger = (PAPER / "TERMINOLOGY_LEDGER.md").read_text(encoding="utf-8")
    assert (
        "| Released endpoint work |" in ledger
        and "| Closed-loop endpoint circulation |" in ledger
    )
    for phrase in (
        "Cauchy inequality bounds only",
        "does not bound",
        "exact-envelope assumption",
        "nonconservative DFT forces",
        "does not identify one mode as the sole cause",
    ):
        assert phrase in main
    for phrase in (
        "118newSCFs",
        "preserving30/81",
        "T29/81",
        "8",
        "zero endpoint gradients",
        "not energy error",
        "40nearest axial edges",
    ):
        assert phrase in methods
