"""Distinguish initial structure snapshots from actual final calculator inputs."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
from ase.io import read


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "paper" / "VARNEB_CPC" / "evidence"
SAVED = EVIDENCE / "gan_45p7_image15_geometry_sanity_20260927"
CHAIN = EVIDENCE / "gan_45p7_final_chains_20260927"
ACTUAL = EVIDENCE / "gan_45p7_final_peak_inputs_20260927"
HASHES = {
    "abacus": "f1050345c9750c8adb67bf3b82c4e75746e5d30dd95f75c7f00e295b27ed6e07",
    "vasp": "3e323ad457fb458106c1e67de5ff5f747e8229013be6b0ed5ef09f069f425865",
    "qe": "99ced52beab323983966e67bef15a76f29afda436d7841ccafe6d0334e389570",
    "abinit": "eee69c323aac34a66eda681442766d44dbf62a229a4661766cfa4fe7b43250d5",
    "cp2k": "4c23aa4c6de89b1466467136d6242fa1741a89f6e15293a2ac7a811515953798",
}


def test_three_initial_structure_snapshots_differ_from_final_image_15() -> None:
    for backend, digest in HASHES.items():
        saved_path = SAVED / f"{backend}.vasp"
        assert hashlib.sha256(saved_path.read_bytes()).hexdigest() == digest
        saved = read(saved_path, format="vasp")
        final = read(CHAIN / f"gan_{backend}_45p7_final_chain.traj", index="15")
        assert saved.get_chemical_formula() == final.get_chemical_formula() == "Ga2N2"
        delta_volume = saved.get_volume() - final.get_volume()
        if backend in {"abacus", "vasp"}:
            assert abs(delta_volume) < 1e-10
        else:
            assert 0.2 < delta_volume < 0.3


def test_actual_qe_abinit_cp2k_inputs_match_the_final_peak_geometry() -> None:
    cases = {
        "qe": ("qe.pwi", "espresso-in",
               "b558fff1f3b16379a009f905514fa92c141062f24bb2859d1b323e57b190c767"),
        "abinit": ("abinit.in", "abinit-in",
                   "b4d4388996655c4d13aae3c786dd3bd15f01c81ab13a23198a2d74e13b2069f3"),
        "cp2k": ("cp2k.inp", "cp2k-restart",
                 "935170d769aa368ddc8f1bc3089223953b275ea5b52bc70168b84285ab93cff1"),
    }
    for backend, (filename, fmt, digest) in cases.items():
        path = ACTUAL / filename
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
        actual = read(path, format=fmt)
        final = read(CHAIN / f"gan_{backend}_45p7_final_chain.traj", index="15")
        assert actual.get_chemical_symbols() == final.get_chemical_symbols()
        np.testing.assert_allclose(actual.cell.array, final.cell.array, atol=1e-9, rtol=0)
        np.testing.assert_allclose(actual.positions, final.positions, atol=1e-9, rtol=0)
