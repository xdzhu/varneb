"""The peak raw-output evidence must not silently promote CP2K's tail."""

from __future__ import annotations

import json
from pathlib import Path

from scripts.audit_gan_image15_raw_outputs import audit_peak_outputs


ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "paper/VARNEB_CPC/evidence/gan_45p7_peak_raw_outputs_20260927/audit.json"


def test_frozen_image15_raw_output_audit() -> None:
    report = audit_peak_outputs()
    assert report == json.loads(FROZEN.read_text(encoding="utf-8"))
    for backend in ("qe", "abinit"):
        assert report["findings"][backend]["status"] == "raw_output_reproduces_final_peak_energy_force_stress"
        assert report["findings"][backend]["energy_difference_eV"] == 0.0
        assert report["findings"][backend]["max_absolute_force_difference_eV_per_A"] == 0.0
        assert report["findings"][backend]["max_absolute_stress_difference_eV_per_A3"] == 0.0
    cp2k = report["findings"]["cp2k"]
    assert cp2k["status"] == "last_visible_text_output_does_not_certify_final_peak"
    assert cp2k["last_program_end_precedes_last_energy"]
    assert 0.0003 < abs(cp2k["last_visible_energy_minus_chain_eV"]) < 0.0004
    assert cp2k["exact_geometry_cache_key"].startswith("image_0015_cc482fbb3276f4d8")
    assert cp2k["worker_cache_energy_difference_eV"] == 0.0
    assert cp2k["worker_cache_max_absolute_force_difference_eV_per_A"] == 0.0
    assert cp2k["worker_cache_max_absolute_stress_difference_eV_per_A3"] == 0.0
