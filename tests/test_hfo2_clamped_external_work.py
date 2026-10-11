"""External physical work retains reactions and declares its frozen chart."""
from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.calculators.singlepoint import SinglePointCalculator

from scripts.analyze_hfo2_clamped_external_work import analyze, image_work
from vcneb import clamped_plane_vcneb_boundary

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT/"benchmarks/hfo2_channels/20261008/clamped_M_terminal_E068_20261011"
SHA = "56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c"


@pytest.mark.parametrize("rotated", [False, True])
@pytest.mark.parametrize("strain", [0., .01])
def test_real_energy_derivative_of_each_declared_frozen_chart(rotated, strain):
    image = bulk("Cu", "fcc", a=3.6, cubic=True)
    cell = image.cell.array.copy()
    cell[2] += [.27, -.19, .11]
    if rotated:
        rotation, _ = np.linalg.qr(np.random.default_rng(71).normal(size=(3, 3)))
        cell = cell @ rotation
    image.set_cell(cell, scale_atoms=True)
    image.calc = EMT()
    boundary = clamped_plane_vcneb_boundary(len(image), cell, allow_tilt=True)
    result = image_work(image, boundary, strain=strain, cell_scale_A=3.6)
    h = 1e-5
    for key, jacobian_key in (
        ("fixed_open_row_partial_eV_cell_per_strain", "fixed_open_row_cell_jacobian_A_per_strain"),
        ("affine_open_row_partial_eV_cell_per_strain", "affine_cell_jacobian_A_per_strain"),
    ):
        values = []
        for sign in (1, -1):
            probe = image.copy()
            probe.set_cell(cell + sign*h*np.array(result[jacobian_key]), scale_atoms=True)
            probe.calc = EMT()
            values.append(probe.get_potential_energy())
        assert result[key] == pytest.approx((values[0]-values[1])/(2*h), abs=3e-5)


def test_open_traction_vanishes_but_substrate_reaction_still_does_work():
    image = bulk("Cu", "fcc", a=3.6, cubic=True)
    image.cell[2] += [.27, -.19, .11]
    boundary = clamped_plane_vcneb_boundary(len(image), image.cell.array, allow_tilt=True)
    projector = np.eye(3) - np.outer(boundary.normal, boundary.normal)
    stress = .04*projector
    image.calc = SinglePointCalculator(image, energy=0., forces=np.zeros((len(image), 3)), stress=stress)
    result = image_work(image, boundary, strain=0., cell_scale_A=3.6)
    np.testing.assert_allclose(boundary.open_traction(stress), 0., atol=1e-15)
    assert result["fixed_open_row_partial_eV_cell_per_strain"] == pytest.approx(2*.04*image.get_volume())
    assert result["affine_minus_fixed_open_row_eV_cell_per_strain"] == pytest.approx(0., abs=1e-13)


def test_tilted_nonstationary_open_row_has_chart_dependent_partial_work():
    image = bulk("Cu", "fcc", a=3.6, cubic=True)
    image.cell[2, 0] = .4
    stress = np.array([[.03, 0., .02], [0., .04, 0.], [.02, 0., .05]])
    image.calc = SinglePointCalculator(image, energy=0., forces=np.zeros((len(image), 3)), stress=stress)
    boundary = clamped_plane_vcneb_boundary(len(image), image.cell.array, allow_tilt=True)
    result = image_work(image, boundary, strain=0., cell_scale_A=3.6)
    expected = image.get_volume()*.02*.4/3.6
    assert result["affine_minus_fixed_open_row_eV_cell_per_strain"] == pytest.approx(expected)
    assert result["affine_open_row_partial_eV_cell_per_strain"] != result["fixed_open_row_partial_eV_cell_per_strain"]


def test_actual_native_M_trace_keeps_units_inputs_and_prediction_limit():
    result = analyze(CASE, SHA)
    assert result["native_portable_frames_reparsed"] == 9
    assert result["sampled_peak_image"] == 3
    assert result["strain_training_condition"] == 0.
    assert result["ordinary_fmax_eV_A"] == pytest.approx(.09921408814764762)
    assert result["new_DFT_calls"] == result["scheduler_mutations"] == 0
    assert not result["G3_selection_or_holdout_access"]
    for row in result["images"]:
        assert row["fixed_open_row_partial_meV_fu_per_unit_strain"] == pytest.approx(
            row["fixed_open_row_partial_eV_cell_per_strain"]*250)
    assert result["images"][0]["fixed_open_row_relative_PO_partial_meV_fu_per_unit_strain"] == 0.


def test_unadmitted_receipt_is_rejected_before_any_label_read(monkeypatch):
    import scripts.analyze_hfo2_clamped_external_work as module

    def forbidden_read(*args):
        raise AssertionError("an unadmitted label must not be opened")

    monkeypatch.setattr(module, "terminal_images", forbidden_read)
    with pytest.raises(ValueError):
        module.analyze(CASE, "0"*64)


def test_actual_reversing_work_remains_frozen_chart_not_stationary_barrier_slope():
    case = ROOT/"benchmarks/hfo2_channels/20261008/clamped_reversing_terminal_E083_20261011"
    result = analyze(case, "632e0e9b2b821b3bfc7d66b361c44c7db734988fb5d018ea2f31441fab14ed41")
    assert result["status"] == "native_training_band_frozen_external_work_not_stationary_barrier_response"
    assert result["new_DFT_calls"] == result["scheduler_mutations"] == 0
    assert result["G3_selection_or_holdout_access"] is False
    assert result["source_job_id"] == "28692776" and result["sampled_peak_image"] == 4
    peak = result["sampled_peak"]
    assert peak["fixed_open_row_relative_PO_partial_meV_fu_per_unit_strain"] == pytest.approx(201.23094619250415, abs=1e-9)
    assert peak["fixed_open_row_partial_meV_fu_per_unit_strain"] == pytest.approx(-27.183577479466106, abs=1e-9)
    # Almost untilted geometry explains tiny reparameterization work; this
    # does not certify that the sampled peak is a stationary TS branch.
    assert abs(peak["affine_minus_fixed_open_row_eV_cell_per_strain"]) < 2e-11
