"""Do not mistake mutable backend workdirs for final-chain raw evidence."""

from __future__ import annotations

import hashlib
from pathlib import Path

from ase.io import read


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "paper" / "VARNEB_CPC" / "evidence"
RAW = EVIDENCE / "gan_45p7_peak_raw_input_20260927"
CHAIN = EVIDENCE / "gan_45p7_final_chains_20260927"
HASHES = {
    "abacus": "f1050345c9750c8adb67bf3b82c4e75746e5d30dd95f75c7f00e295b27ed6e07",
    "vasp": "3e323ad457fb458106c1e67de5ff5f747e8229013be6b0ed5ef09f069f425865",
    "qe": "99ced52beab323983966e67bef15a76f29afda436d7841ccafe6d0334e389570",
    "abinit": "eee69c323aac34a66eda681442766d44dbf62a229a4661766cfa4fe7b43250d5",
    "cp2k": "4c23aa4c6de89b1466467136d6242fa1741a89f6e15293a2ac7a811515953798",
}


def test_three_mutable_peak_directories_are_not_final_image_15() -> None:
    for backend, digest in HASHES.items():
        raw_path = RAW / f"{backend}.vasp"
        assert hashlib.sha256(raw_path.read_bytes()).hexdigest() == digest
        current = read(raw_path, format="vasp")
        final = read(CHAIN / f"gan_{backend}_45p7_final_chain.traj", index="15")
        assert current.get_chemical_formula() == final.get_chemical_formula() == "Ga2N2"
        delta_volume = current.get_volume() - final.get_volume()
        if backend in {"abacus", "vasp"}:
            assert abs(delta_volume) < 1e-10
        else:
            assert 0.2 < delta_volume < 0.3
