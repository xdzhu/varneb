"""Figure data must preserve the audited GaN TS/endpoint protocol boundary."""

import json

import matplotlib.pyplot as plt
import pytest

from scripts.plot_gan_joint_mode_coupling import (
    DEFAULT_REPORT, draw_figure, load_evidence, source_rows,
)


def test_gan_coupling_figure_uses_pinned_audited_data_and_distinct_contracts():
    report, digest = load_evidence(DEFAULT_REPORT)
    rows = source_rows(report, digest)
    assert len(rows) == 14
    joint = [row for row in rows if row["panel"] == "a" and row["category"] == "joint"]
    assert len(joint) == 2
    assert all(row["value"] < 0 and row["electronic_contract"] == "VASP 1000 eV"
               for row in joint)
    optical = [row for row in rows if row["panel"] == "b"]
    assert len(optical) == 6
    assert all(row["electronic_contract"] == "VASP 600 eV endpoint Gamma basis"
               and row["unit"] == "% atomic optical direction"
               for row in optical)
    for phase in ("B4", "B1"):
        assert sum(row["value"] for row in optical if row["phase"] == phase) == pytest.approx(100.0)
    figure = draw_figure(rows)
    assert len(figure.axes) == 2
    assert figure.axes[0].get_position().y0 == figure.axes[1].get_position().y0
    assert figure.axes[0].get_position().y1 == figure.axes[1].get_position().y1
    plt.close(figure)


def test_gan_coupling_figure_rejects_modified_report(tmp_path):
    altered = json.loads(DEFAULT_REPORT.read_text(encoding="utf-8"))
    altered["joint_hessian_negative_eigenvalues_eV_per_A2"]["0p01"] = +4.0
    path = tmp_path / "altered.json"
    path.write_text(json.dumps(altered), encoding="utf-8")
    with pytest.raises(ValueError, match="pinned report"):
        load_evidence(path)
