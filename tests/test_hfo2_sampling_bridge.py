from pathlib import Path
import json

import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.io import read

from scripts import prepare_hfo2_sampling_bridge as bridge

ROOT = Path(__file__).resolve().parents[1]
OBSERVATION = ROOT / "benchmarks/hfo2_channels/20261008/switching_converged_update_20261009_1255/preserving_step69"


@pytest.fixture()
def prepared(tmp_path):
    output = tmp_path / "six"
    return output, bridge.prepare(OBSERVATION, output)


def test_actual_complete_band_selects_only_six_registered_statics(prepared):
    output, m = prepared
    assert m["maximum_new_SCF_calls"] == 6 and m["new_DFT_calls_for_preparation"] == 0
    assert [p["segment"] for p in m["points"]] == [[2, 3]]*3 + [[5, 6]]*3
    assert [p["fraction"] for p in m["points"]] == [.25, .5, .75]*2
    assert m["sampled_band_maximum_meV_fu"] == pytest.approx(32.806023216835456)
    for segment in m["segments"]:
        assert segment["dE_dlambda_meV_fu"][0] > 0
        assert segment["dE_dlambda_meV_fu"][1] < 0
        assert segment["Hermite_candidates_not_DFT"][0]["energy_meV_fu"] == pytest.approx(38.217947, abs=1e-5)
    assert not m["automatic_restart_or_G2_submission"]
    assert not m["climb"] and m["ordinary_fmax_eV_A"] == .10
    assert len(list(output.glob("geometries/*.vasp"))) == 6


def test_fractional_lift_and_cell_are_exact_linear_not_MIC(prepared):
    output, m = prepared
    a = read(OBSERVATION / "POSCAR_02", format="vasp")
    b = read(OBSERVATION / "POSCAR_03", format="vasp")
    # A common integer representative can lie outside the principal cell.
    shift = np.tile([1, -1, 2], (len(a), 1))
    for x in (a, b):
        x.set_scaled_positions(x.get_scaled_positions(wrap=False)+shift)
    mid = bridge.interpolate_segment(a, b, .5)
    np.testing.assert_allclose(mid.get_scaled_positions(wrap=False),
                               .5*(a.get_scaled_positions(wrap=False)+b.get_scaled_positions(wrap=False)), atol=1e-14)
    np.testing.assert_allclose(mid.cell.array, .5*(a.cell.array+b.cell.array), atol=1e-14)
    assert mid.calc is None
    assert mid.get_chemical_symbols() == a.get_chemical_symbols()


def test_segment_work_matches_energy_derivative_with_atomic_and_oblique_cell_change():
    a = bulk("Cu", "fcc", a=3.6, cubic=True)
    b = a.copy()
    b.set_cell(a.cell.array @ np.array([[1.01, .003, 0], [0, .99, .002], [.004, 0, 1.005]]).T,
               scale_atoms=True)
    q = b.get_scaled_positions(wrap=False)
    q[0] += [.012, -.004, .007]
    b.set_scaled_positions(q)
    for x in (a, b):
        x.calc = EMT()
    expected = bridge.segment_derivatives(a, b, a.cell.array)
    h = 1e-5
    def energy(f):
        if f in (0, 1):
            return (a if f == 0 else b).get_potential_energy()
        x = bridge.interpolate_segment(a, b, f)
        x.calc = EMT()
        return x.get_potential_energy()
    observed = [(-3*energy(0)+4*energy(h)-energy(2*h))/(2*h),
                (3*energy(1)-4*energy(1-h)+energy(1-2*h))/(2*h)]
    np.testing.assert_allclose(expected, observed, rtol=0, atol=1e-5)


@pytest.mark.parametrize("key,value", [("climb", True), ("pressure_GPa", 1.), ("formula_units", 12),
                                      ("common_PO_energy_eV_cell", 0.), ("ordinary_fmax_eV_A", .01)])
def test_tampered_physical_interpretation_rejected(prepared, key, value):
    root, m = prepared
    m[key] = value
    (root / "manifest.json").write_text(json.dumps(m))
    with pytest.raises(ValueError, match="contract"):
        bridge.load_manifest(root)


def test_changed_geometry_rejected_before_calculator(prepared, monkeypatch):
    root, m = prepared
    p = root / m["points"][0]["geometry"]
    p.write_bytes(p.read_bytes()+b"\n")
    monkeypatch.setattr(bridge, "FixedHfo2Calculator", lambda **_: pytest.fail("calculator loaded"))
    with pytest.raises(ValueError, match="geometry changed"):
        bridge.run_point(root, 0)


def test_existing_output_or_point_or_summary_refused_early(prepared, monkeypatch):
    root, _ = prepared
    with pytest.raises(FileExistsError):
        bridge.prepare(Path("missing"), root)
    (root / "calculations/00").mkdir(parents=True)
    monkeypatch.setattr(bridge, "FixedHfo2Calculator", lambda **_: pytest.fail("calculator loaded"))
    with pytest.raises(FileExistsError):
        bridge.run_point(root, 0)
    (root / "summary.json").write_text("retain")
    with pytest.raises(FileExistsError):
        bridge.summarize(root)
    assert (root / "summary.json").read_text() == "retain"


def test_wrong_recipe_cannot_retarget_geometry(prepared):
    root, m = prepared
    m["points"][0]["geometry"] = "../different.vasp"
    (root / "manifest.json").write_text(json.dumps(m))
    with pytest.raises(ValueError, match="recipe"):
        bridge.load_manifest(root)
