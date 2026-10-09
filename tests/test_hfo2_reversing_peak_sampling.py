from pathlib import Path
import json
import shutil
import pytest

from scripts import prepare_hfo2_reversing_peak_sampling as peak

ROOT = Path(__file__).resolve().parents[1]
COMPLETED = ROOT / "benchmarks/hfo2_channels/20261008/reversing_peak_sampling_20261009/completed_HF/prepared"


@pytest.fixture()
def prepared(tmp_path):
    path = tmp_path / "two"
    return path, peak.prepare(ROOT, path)


def test_actual_reversing_source_two_point_cap(prepared):
    path, m = prepared
    assert len(m["points"]) == 2
    assert [(p["segment"][0], p["segment"][1], p["fraction"]) for p in m["points"]] == list(peak.RECIPE)
    assert m["source_fmax_eV_A"] == pytest.approx(.09402284566626806, abs=1e-12, rel=0)
    assert m["screens_not_DFT"][0]["Hermite_candidates_not_DFT"][0]["energy_meV_fu"] == pytest.approx(392.82995563637496)
    assert m["prior_peak_SCF_calls"] + m["maximum_new_SCF_calls"] == 13 < m["G1_peak_sampling_cap"]
    assert m["new_DFT_calls_for_preparation"] == 0 and not m["automatic_restart_or_G2_submission"]
    assert peak.load_manifest(path) == m


@pytest.mark.parametrize("key,value", [("pressure_GPa", 1), ("climb", True), ("maximum_new_SCF_calls", 3),
                                      ("source_step", 44), ("common_PO_energy_eV_cell", 0), ("ordinary_fmax_eV_A", .01)])
def test_contract_tamper_rejected(prepared, key, value):
    path, m = prepared
    m[key] = value
    peak.save(path / "manifest.json", m)
    with pytest.raises(ValueError, match="contract"):
        peak.load_manifest(path)


def test_recipe_tamper_rejected(prepared):
    path, m = prepared
    m["points"][0]["geometry"] = "../outside.vasp"
    peak.save(path / "manifest.json", m)
    with pytest.raises(ValueError, match="recipe"):
        peak.load_manifest(path)


@pytest.mark.parametrize("index", [True, -1, 2, None])
def test_out_of_budget_no_DFT(prepared, monkeypatch, index):
    monkeypatch.setattr(peak, "FixedHfo2Calculator", lambda **_: pytest.fail("DFT loaded"))
    with pytest.raises(ValueError, match="index"):
        peak.run_point(prepared[0], index)


def test_preserve_outputs_changed_geometry_rejected(prepared, monkeypatch):
    path, m = prepared
    monkeypatch.setattr(peak, "FixedHfo2Calculator", lambda **_: pytest.fail("DFT loaded"))
    with pytest.raises(FileExistsError):
        peak.prepare(ROOT, path)
    geometry = path / m["points"][0]["geometry"]
    geometry.write_bytes(geometry.read_bytes()+b"\n")
    with pytest.raises(ValueError, match="geometry changed"):
        peak.run_point(path, 0)
    (path / "calculations/00").mkdir(parents=True)
    with pytest.raises(FileExistsError):
        peak.run_point(path, 0)
    (path / "summary.json").write_text("keep")
    with pytest.raises(FileExistsError):
        peak.summarize(path)


def test_actual_two_SCFs_no_extra_optimizer_or_G2(tmp_path):
    report = peak.analyze(ROOT, COMPLETED, tmp_path / "analysis.json")
    assert report["new_SCF_calls"] == 2 and report["total_G1_peak_SCF_calls"] == 13
    assert report["highest_union_sample_meV_fu"] == pytest.approx(392.8229051971357, rel=0, abs=1e-8)
    assert report["highest_new_sample_meV_fu"] < report["old_sampled_maximum_meV_fu"]
    assert report["inserted_total_images"] == 11 and report["moving_images"] == 9
    assert report["ordinary_residual_passed"]
    assert report["replayed_ordinary_fmax_eV_A"] == pytest.approx(.09402284566626806, rel=0, abs=1e-12)
    assert report["optimizer_steps"] == report["additional_SCF_calls"] == report["new_DFT_calls_for_analysis"] == 0
    assert not report["automatic_restart_or_G2_submission"]


def test_delivered_prepared_geometries_preserve_declared_raw_bytes():
    case = COMPLETED.parents[1]
    for root in (case / "prepared_local", case / "prepared_HF_LF", COMPLETED):
        m = peak.load_manifest(root)
        for point in m["points"]:
            assert peak.sha256(root / point["geometry"]) == point["geometry_sha256"]


@pytest.mark.parametrize("kind", ["force", "mpi", "log"])
def test_actual_evidence_tamper_rejected(tmp_path, kind):
    root = tmp_path / "changed"
    shutil.copytree(COMPLETED, root)
    if kind == "force":
        p = root / "calculations/00/point_result.json"
        record = json.loads(p.read_text())
        record["results"]["forces"][0][0] += .01
        peak.save(p, record)
    elif kind == "mpi":
        (root / "mpi_affinity.txt").write_text("rank=0 host=node61\n" * 32)
    else:
        p = root / "calculations/00/scf_000000/OUT.ABACUS/running_scf.log"
        p.write_bytes(p.read_bytes()+b"\nchanged\n")
    with pytest.raises(ValueError):
        peak.analyze(ROOT, root, tmp_path / "rejected.json")
