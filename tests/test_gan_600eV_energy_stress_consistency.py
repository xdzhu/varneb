import json

import pytest

from scripts.audit_gan_600eV_energy_stress_consistency import (
    OUTPUT,
    basis_metadata,
    simpson_average,
)


def test_simpson_average_integrates_a_cubic_exactly():
    polynomial = lambda x: 1 + 2 * x + 3 * x**2 + 4 * x**3
    assert simpson_average(polynomial(-1), polynomial(0), polynomial(1)) == pytest.approx(2)


def test_basis_metadata_extracts_both_fft_grids_and_plane_waves(tmp_path):
    outcar = tmp_path / "OUTCAR"
    outcar.write_text(
        "dimension x,y,z NGX = 24 NGY = 24 NGZ = 40\n"
        "dimension x,y,z NGXF= 48 NGYF= 48 NGZF= 80\n"
        "maximum number of plane-waves: 1242\n",
        encoding="utf-8",
    )
    assert basis_metadata(outcar) == {
        "maximum_plane_waves": 1242,
        "coarse_fft_grid": [24, 24, 40],
        "fine_fft_grid": [48, 48, 80],
    }


def test_committed_600eV_audit_does_not_certify_ts_or_change_protocol():
    report = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert report["status"] == "GaN_600eV_raw_energy_stress_discrepancy_localized_to_volume_strain"
    assert report["pressure_GPa"] == 45.7
    assert report["step_A"] == 0.02
    assert "600 eV" in report["electronic_contract"]
    assert len(report["axes"]) == 18
    assert all(row["fft_grids_unchanged"] for row in report["axes"])
    residual = report["max_abs_residual_eV_per_A_by_group"]
    assert residual["atomic"] < 0.001
    assert residual["shear_strain"] < 0.0001
    assert residual["diagonal_strain"] > 0.02
    assert all(
        0.3 < row["observed_minus_integrated_pressure_GPa"] < 0.32
        for row in report["axes"]
        if row["group"] == "diagonal_strain"
    )
    assert "do not uniquely prove" in report["interpretation"]
