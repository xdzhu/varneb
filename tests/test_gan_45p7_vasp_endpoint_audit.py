"""Production GaN endpoint statics must not be confused with basin returns."""

import pytest

from scripts import audit_gan_45p7_vasp_endpoints as endpoint_audit


def test_archived_production_endpoints_match_raw_vasp_statics():
    report = endpoint_audit.audit()
    assert report["status"] == "GaN_45p7_VASP_production_endpoint_raw_statics_audited"
    assert report["total_images"] == 29
    assert report["contract"]["ENCUT_eV"] == 600
    assert report["contract"]["k_mesh"] == [8, 8, 6]
    initial, final = report["endpoints"]
    assert (initial["phase"], initial["image_index"]) == ("B4", 0)
    assert (final["phase"], final["image_index"]) == ("B1", 28)
    assert all(item["force_pass_0p02eV_per_A"] for item in (initial, final))
    assert initial["max_raw_stress_residual_kbar"] == pytest.approx(2.91154)
    assert not initial["stress_pass_2kbar"]
    assert final["max_raw_stress_residual_kbar"] == pytest.approx(0.23863)
    assert final["stress_pass_2kbar"]
    assert all(item["raw_minus_chain_energy_eV"] == 0 for item in (initial, final))
    assert all(item["max_POSCAR_vs_chain_position_difference_A"] == 0
               for item in (initial, final))


def test_endpoint_audit_rejects_source_hash_drift(monkeypatch):
    monkeypatch.setitem(endpoint_audit.EXPECTED_SOURCE_SHA256,
                        "initial_00_OUTCAR", "0" * 64)
    with pytest.raises(ValueError, match="source hash changed"):
        endpoint_audit.audit()
