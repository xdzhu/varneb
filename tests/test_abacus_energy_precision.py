"""ABACUS machine energy should not lose precision to the display line."""

from __future__ import annotations

import pytest

import numpy as np

from vcneb.abacus import _minimal_abacus_results, _prefer_machine_precision_energy


def test_prefers_matching_full_precision_marker(tmp_path) -> None:
    log = tmp_path / "running_scf.log"
    log.write_text(
        " final etot is -3737.9367135 eV\n"
        " !FINAL_ETOT_IS -3.7379367134697677102D+03 eV\n",
        encoding="utf-8",
    )
    assert _prefer_machine_precision_energy(log, -3737.9367135) == pytest.approx(
        -3737.9367134697677, abs=1e-12,
    )


def test_uses_last_complete_marker_and_keeps_legacy_without_one(tmp_path) -> None:
    log = tmp_path / "running_scf.log"
    log.write_text(" final etot is -10.0 eV\n", encoding="utf-8")
    assert _prefer_machine_precision_energy(log, -10.0) == -10.0
    log.write_text(" !FINAL_ETOT_IS -11.0 eV\n !FINAL_ETOT_IS -10.00000004 eV\n",
                   encoding="utf-8")
    assert _prefer_machine_precision_energy(log, -10.0) == pytest.approx(-10.00000004)


def test_large_marker_disagreement_fails_instead_of_changing_energy_semantics(tmp_path) -> None:
    log = tmp_path / "running_scf.log"
    log.write_text(" !FINAL_ETOT_IS -11.0 eV\n", encoding="utf-8")
    with pytest.raises(ValueError, match="disagrees"):
        _prefer_machine_precision_energy(log, -10.0)
    with pytest.raises(ValueError, match="non-finite"):
        _prefer_machine_precision_energy(log, float("nan"))


def test_minimal_abacus_fallback_converts_compressive_kbar_to_ase_stress(tmp_path) -> None:
    log = tmp_path / "running_scf.log"
    log.write_text(
        "TOTAL-FORCE (eV/Angstrom)\n"
        "Ba1  0.01 -0.02 0.03\n"
        "TOTAL-STRESS (KBAR)\n"
        "-1.602176634 0.0 -0.801088317\n"
        "0.0 -3.204353268 0.0\n"
        "-0.801088317 0.0 -4.806529902\n"
        "!FINAL_ETOT_IS -3.7379367171902390510D+03 eV\n",
        encoding="utf-8",
    )
    result = _minimal_abacus_results(log)
    assert result["energy"] == pytest.approx(-3737.936717190239)
    assert np.allclose(result["forces"], [[0.01, -0.02, 0.03]])
    assert np.allclose(result["stress"], [0.001, 0.002, 0.003, 0.0, 0.0005, 0.0], atol=1e-12)
