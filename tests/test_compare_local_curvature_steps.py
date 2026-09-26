import pytest

from scripts.compare_local_curvature_steps import compare


def bto_report(step, lowest, *, refined="same"):
    return {
        "source_sha256": {
            "preflight": "same", "refined_result": refined, "branch_replay": "same",
        },
        "q1_q2_sqrt_amu_A": [0.6, 0.0], "axis_kind": "transverse_soft",
        "n_matched_signed_probes": 32, "step_sqrt_amu_A": step,
        "eigenvalues_eV_per_amu_A2": [lowest] + [0.1 + i for i in range(15)],
        "energy_hessian_diagonal_max_abs_difference_eV_per_amu_A2": 0.008,
        "energy_gradient_max_abs_difference_eV_per_sqrt_amu_A": 0.001,
    }


def test_bto_comparison_retains_numerical_resolution_warning():
    result = compare(bto_report(0.05, 0.005), bto_report(0.10, 0.004), "bto")
    assert result["negative_index_stable_by_cutoff"]["0.0"]
    assert result["smallest_absolute_curvature_below_observed_diagonal_mismatch"]
    assert result["status"].endswith("not_minimum_or_TS_certificate")


def test_bto_comparison_rejects_different_stationary_center():
    with pytest.raises(ValueError, match="physical source"):
        compare(bto_report(0.05, 0.005), bto_report(0.10, 0.004, refined="other"), "bto")


def test_bto_two_step_energy_derivative_extrapolation_is_diagnostic():
    first = bto_report(0.05, 0.005)
    second = bto_report(0.10, 0.004)
    for report in (first, second):
        step = report["step_sqrt_amu_A"]
        report["force_gradient_per_open_direction_eV_per_sqrt_amu_A"] = [0.4] * 16
        report["energy_gradient_per_open_direction_eV_per_sqrt_amu_A"] = [
            0.4 + 2.0 * step**2
        ] * 16
    diagnostic = compare(first, second, "bto")["energy_force_gradient_step_extrapolation"]
    assert diagnostic["max_abs_error_first_eV_per_sqrt_amu_A"] == pytest.approx(0.005)
    assert diagnostic["max_abs_error_second_eV_per_sqrt_amu_A"] == pytest.approx(0.02)
    assert diagnostic["max_abs_error_extrapolated_eV_per_sqrt_amu_A"] == pytest.approx(0.0)


def test_bto_two_step_gradient_requires_same_center():
    first = bto_report(0.05, 0.005)
    second = bto_report(0.10, 0.004)
    for report in (first, second):
        report["force_gradient_per_open_direction_eV_per_sqrt_amu_A"] = [0.4] * 16
        report["energy_gradient_per_open_direction_eV_per_sqrt_amu_A"] = [0.4] * 16
    second["force_gradient_per_open_direction_eV_per_sqrt_amu_A"][0] = 0.5
    with pytest.raises(ValueError, match="shared center"):
        compare(first, second, "bto")


def test_gan_comparison_checks_one_negative_mode_across_steps():
    base = {
        "source_sha256": {"center_OUTCAR": "same"},
        "pressure_GPa": 45.7, "encut_eV": 1000,
        "n_static_displacements": 36, "center_enthalpy_eV_per_cell": -11.5,
        "energy_hessian_diagonal_max_abs_difference_eV_per_A2": 0.2,
        "energy_gradient_max_abs_difference_eV_per_A": 0.01,
    }
    first = dict(base, step_A=0.02,
                 translation_free_eigenvalues_eV_per_A2=[-4.0] + [2.0 + i for i in range(14)])
    second = dict(base, step_A=0.01,
                  translation_free_eigenvalues_eV_per_A2=[-4.1] + [2.1 + i for i in range(14)])
    result = compare(first, second, "gan")
    assert all(result["negative_index_stable_by_cutoff"].values())
    assert result["lowest_eigenvalue_change"] == pytest.approx(0.1)
