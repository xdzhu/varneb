"""Material boundary/dual-frame tests without importing author mode arrays."""
from pathlib import Path
import json
import numpy as np
import pytest

from scripts.analyze_hfo2_clamped_mode_gradient import analyze
from scripts.analyze_hfo2_clamped_reference import terminal_images, original_chart_images
from scripts.analyze_hfo2_clamped_external_work import SEEDS, SEEDS_SHA
from scripts.analyze_hfo2_clamped_reference import pinned_json
from vcneb.continuous_projection import project_reference_atomic_gradient
from scripts.analyze_hfo2_cmma_path_reference import MASS
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/"benchmarks/hfo2_channels/20261008"
AUDIT = "56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c"


def test_unknown_training_receipt_rejected_before_any_paths_are_opened():
    with pytest.raises(ValueError, match="before label access"):
        analyze(None, None, None, None, None, terminal_audit_sha256="0"*64,
                clamped_seed_manifest=None, clamped_seed_manifest_sha256=None)


def test_native_covector_transport_matches_direct_production_frame():
    images, _, residual, _ = terminal_images(CASE/"clamped_M_terminal_E068_20261011", AUDIT)
    seeds = pinned_json(SEEDS, SEEDS_SHA)
    t = read(CASE/"reference_variants/T.vasp")
    original = original_chart_images(images, seeds["orientation"], t)
    r = np.asarray(seeds["orientation"]["cartesian_rotation_rows"])
    rng = np.random.default_rng(81)
    basis = np.linalg.qr(rng.normal(size=(36, 4)))[0]
    rotated = np.einsum("nak,ba->nbk", basis.reshape(12, 3, 4), r).reshape(36, 4)
    for native, old in zip(images, original):
        f = native.get_forces(apply_constraint=False)
        a = project_reference_atomic_gradient(f @ r, old.cell.array, t.cell.array, basis, masses_amu=MASS)
        b = project_reference_atomic_gradient(f, native.cell.array, r @ t.cell.array @ r.T, rotated, masses_amu=MASS)
        for key in ("internal_gradient_norm", "retained_gradient_norm", "omitted_gradient_norm", "retained_squared_gradient_fraction"):
            assert a[key] == pytest.approx(b[key], abs=1e-13, rel=0)
        assert old.calc is None
    assert residual["fmax_eV_A"] == pytest.approx(.09921408814764762, abs=1e-12)


def test_actual_native_report_preserves_all_references_and_claim_limits():
    result = json.loads((CASE/"clamped_mode_gradient_E078_20261011/analysis.json").read_text())
    assert result["terminal_audit_sha256"] == AUDIT
    assert result["source_job_id"] == "28722320" and result["sampled_peak_image"] == 3
    assert result["new_DFT_calls"] == result["scheduler_mutations"] == 0
    assert result["G3_selection_or_holdout_access"] is False
    assert result["complete_six_physical_bytes_rechecked_locally"] is False
    assert result["physical_inputs_production_modes_geometry_or_optimizer_changed"] is False
    rows = result["images"]
    assert [r["image_index"] for r in rows] == list(range(9))
    expected = ["T_geometric_triplet", "T_two_lowest_complete_optical_doublets",
                "Cmma_frame_103", "Cmma_frame_106", "Cmma_frame_114", "Cmma_frame_127"]
    assert all(list(r["representations"]) == expected for r in rows)
    peak = result["sampled_peak"]["representations"]
    assert peak[expected[0]]["retained_squared_gradient_fraction"] == pytest.approx(.09028221034427668, abs=1e-12)
    assert peak[expected[1]]["retained_squared_gradient_fraction"] == pytest.approx(.32709626956356713, abs=1e-12)
    assert peak[expected[2]]["retained_squared_gradient_fraction"] == pytest.approx(.9350846598522286, abs=1e-12)
    for row in rows:
        for name, span in row["representations"].items():
            assert span["source_basis_rank"] == span["translation_free_basis_rank"] == (3 if name == expected[0] else 4)
            assert abs(span["norm_identity_defect"]) < 1e-12
            assert span["source_modes_or_forces_ASR_corrected"] is False
            assert span["translation_free_chart_restriction_applied"] is True
