from pathlib import Path
import json

import pytest

from scripts import prepare_hfo2_G1_peak_sampling as peak

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture()
def prepared(tmp_path):
    root = tmp_path / "five"
    return root, peak.prepare(ROOT, root)


def test_actual_sources_and_fixed_five_point_budget(prepared):
    root, m = prepared
    assert len(m["points"]) == 5 and len(list(root.glob("geometries/*.vasp"))) == 5
    assert [(p["channel"], *p["segment"], p["fraction"]) for p in m["points"]] == list(peak.RECIPE)
    assert m["sources"]["PO_to_M"]["sampled_maximum_relative_PO_meV_fu"] == pytest.approx(71.58220716746655)
    assert m["sources"]["T_to_PO"]["sampled_maximum_relative_PO_meV_fu"] == pytest.approx(115.21016416145358)
    assert m["screens_not_DFT"]["PO_to_M:3->4"]["Hermite_candidates_not_DFT"][0]["energy_meV_fu"] == pytest.approx(81.7188487782859)
    assert m["G1_peak_sampling_cap"] == 14
    assert m["new_DFT_calls_for_preparation"] == 0 and not m["automatic_restart_or_G2_submission"]
    assert peak.load_manifest(root) == m


@pytest.mark.parametrize("key,value", [("climb", True), ("pressure_GPa", 1), ("common_PO_energy_eV_cell", 0),
                                      ("formula_units", 12), ("ordinary_fmax_eV_A", .01), ("maximum_new_SCF_calls", 6)])
def test_tampered_contract_rejected(prepared, key, value):
    root, m = prepared
    m[key] = value
    peak.save(root / "manifest.json", m)
    with pytest.raises(ValueError, match="contract"):
        peak.load_manifest(root)


def test_recipe_and_source_tampering_rejected(prepared):
    root, m = prepared
    m["sources"]["T_to_PO"]["step"] = 7
    peak.save(root / "manifest.json", m)
    with pytest.raises(ValueError, match="source"):
        peak.load_manifest(root)
    m["sources"]["T_to_PO"]["step"] = 6
    m["points"][0]["geometry"] = "../other.vasp"
    peak.save(root / "manifest.json", m)
    with pytest.raises(ValueError, match="recipe"):
        peak.load_manifest(root)


def test_existing_and_changed_geometry_refused_before_DFT(prepared, monkeypatch):
    root, m = prepared
    monkeypatch.setattr(peak, "FixedHfo2Calculator", lambda **_: pytest.fail("DFT loaded"))
    with pytest.raises(FileExistsError):
        peak.prepare(Path("missing"), root)
    p = root / m["points"][0]["geometry"]
    p.write_bytes(p.read_bytes()+b"\n")
    with pytest.raises(ValueError, match="geometry changed"):
        peak.run_point(root, 0)
    (root / "calculations/00").mkdir(parents=True)
    with pytest.raises(FileExistsError):
        peak.run_point(root, 0)
    (root / "summary.json").write_text("retain")
    with pytest.raises(FileExistsError):
        peak.summarize(root)


@pytest.mark.parametrize("index", [True, -1, 5, None])
def test_out_of_budget_index_rejected(prepared, index):
    with pytest.raises(ValueError, match="index"):
        peak.run_point(prepared[0], index)
