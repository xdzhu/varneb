from __future__ import annotations

import runpy
from pathlib import Path

from ase import Atoms
from ase.io import read, write
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "setup_batio3_vca_case.py"


def _potcar(path: Path, symbol: str, zval: float) -> None:
    path.write_text(f"   TITEL  = PAW_PBE {symbol}\n   POMASS = 1; ZVAL = {zval:.3f}\n", encoding="utf-8")


def test_setup_duplicates_only_the_ba_site_and_writes_vca_weights(tmp_path: Path) -> None:
    module = runpy.run_path(str(SCRIPT))
    atoms = Atoms("BaTiO3", scaled_positions=[[0, 0, 0]] * 5, cell=[4, 4, 4], pbc=True)
    initial, final = tmp_path / "T.vasp", tmp_path / "C.vasp"
    write(initial, atoms, format="vasp")
    write(final, atoms, format="vasp")
    incar, kpoints = tmp_path / "INCAR", tmp_path / "KPOINTS"
    incar.write_text("ENCUT = 600\n", encoding="utf-8")
    kpoints.write_text("Gamma\n", encoding="utf-8")
    potentials = {}
    for symbol, zval in (("Ba_sv", 10), ("Sr_sv", 10), ("Ti_sv", 12), ("O", 6)):
        path = tmp_path / f"POTCAR.{symbol}"
        _potcar(path, symbol, zval)
        potentials[symbol] = path

    result = module["setup"](
        initial_path=initial, final_path=final, destination=tmp_path / "case",
        incar=incar, kpoints=kpoints, ba_potcar=potentials["Ba_sv"],
        sr_potcar=potentials["Sr_sv"], ti_potcar=potentials["Ti_sv"], o_potcar=potentials["O"],
    )
    expanded = read(tmp_path / "case" / "initial" / "POSCAR.vca")
    assert expanded.get_chemical_symbols() == ["Ba", "Sr", "Ti", "O", "O", "O"]
    assert (expanded.positions[0] == expanded.positions[1]).all()
    assert "VCA = 0.50000000 0.50000000 1.00000000 1.00000000" in (tmp_path / "case" / "initial" / "INCAR").read_text()
    assert result["composition"] == "Ba0.5Sr0.5TiO3"
    assert result["physical_atom_count"] == 5
    assert result["expanded_vasp_atom_count"] == 6


def test_setup_rejects_nonisovalent_a_site_potentials(tmp_path: Path) -> None:
    module = runpy.run_path(str(SCRIPT))
    atoms = Atoms("BaTiO3", cell=[4, 4, 4], pbc=True)
    for name in ("T.vasp", "C.vasp"):
        write(tmp_path / name, atoms, format="vasp")
    (tmp_path / "INCAR").write_text("ENCUT=600\n")
    (tmp_path / "KPOINTS").write_text("Gamma\n")
    pots = []
    for index, (symbol, zval) in enumerate((("Ba_sv", 10), ("La", 11), ("Ti_sv", 12), ("O", 6))):
        path = tmp_path / f"P{index}"
        _potcar(path, symbol, zval)
        pots.append(path)
    with pytest.raises(ValueError, match="same ZVAL"):
        module["setup"](
            initial_path=tmp_path / "T.vasp", final_path=tmp_path / "C.vasp", destination=tmp_path / "case",
            incar=tmp_path / "INCAR", kpoints=tmp_path / "KPOINTS", ba_potcar=pots[0], sr_potcar=pots[1],
            ti_potcar=pots[2], o_potcar=pots[3],
        )
