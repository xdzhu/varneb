"""Prevent the 600-eV GaN basin pilot from changing electronic settings."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

from scripts.prepare_gan_600eV_ts_basin_pilot import relaxation_incar


STATIC = b"""INCAR created by Atomic Simulation Environment
 ENCUT = 600.000000
 SIGMA = 0.050000
 EDIFF = 1.00e-07
 SYMPREC = 1.00e-04
 GGA = PE
 PREC = Accurate
 IBRION = -1
 ISIF = 2
 ISMEAR = 0
 ISYM = -1
 NSW = 0
 LCHARG = .FALSE.
 LWAVE = .FALSE.
"""


def test_only_ionic_pressure_tags_change():
    result = relaxation_incar(STATIC)
    assert b"ENCUT = 600.000000\n" in result
    assert b"EDIFF = 1.00e-07\n" in result
    assert b"SYMPREC = 1.00e-04\n" in result
    assert b"IBRION = 2\n" in result
    assert b"ISIF = 3\n" in result
    assert b"NSW = 10\n" in result
    assert b"PSTRESS = 457.0\n" in result
    assert b"EDIFFG = -0.02\n" in result
    assert b"POTIM = 0.25\n" in result
    assert b"\r" not in result


def test_archived_submitted_input_contract_matches_preparer():
    root = Path(__file__).resolve().parents[1]
    record = json.loads((
        root / "benchmarks/numerical_integrity/gan_600eV_ts_basin_pilot_inputs_20260928.json"
    ).read_text(encoding="utf-8"))
    assert record["purpose"] == "GaN_45p7_600eV_signed_native_VASP_basin_10step_pilot"
    assert record["pressure_GPa"] == 45.7
    assert record["encut_eV"] == 600
    assert record["PSTRESS_kbar"] == 457.0
    assert [(case["name"], case["q_u_A"]) for case in record["cases"]] == [
        ("grid_um_vz", -0.02), ("grid_up_vz", 0.02),
    ]
    digest = hashlib.sha256(relaxation_incar(STATIC)).hexdigest()
    assert all(case["input_sha256"]["INCAR"] == digest for case in record["cases"])


@pytest.mark.parametrize("bad", [
    STATIC.replace(b"600.000000", b"1000.000000"),
    STATIC.replace(b"SYMPREC = 1.00e-04", b"SYMPREC = 1.00e-12"),
    STATIC + b" PSTRESS = 457.0\n",
    STATIC.replace(b"\n", b"\r\n"),
])
def test_rejects_altered_static_contract(bad):
    with pytest.raises(ValueError):
        relaxation_incar(bad)
