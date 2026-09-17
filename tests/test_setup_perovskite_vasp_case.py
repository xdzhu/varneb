from __future__ import annotations

import runpy
from pathlib import Path

from ase import Atoms
from ase.io import read, write
import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "setup_perovskite_vasp_case.py"


def _potcar(path: Path, symbol: str, zval: float) -> None:
    path.write_text(f" TITEL = PAW_PBE {symbol}\n POMASS = 1; ZVAL = {zval:.3f}\n", encoding="latin-1")


def _templates(tmp_path: Path) -> tuple[Path, Path, Path, Path, dict[str, Path]]:
    atoms = Atoms("BaTiO3", scaled_positions=[[0, 0, 0]] * 5, cell=[4, 4, 4], pbc=True)
    initial, final = tmp_path / "initial.vasp", tmp_path / "final.vasp"
    write(initial, atoms, format="vasp")
    write(final, atoms, format="vasp")
    incar, kpoints = tmp_path / "INCAR", tmp_path / "KPOINTS"
    incar.write_text("ENCUT = 600\n", encoding="utf-8")
    kpoints.write_text("Gamma\n", encoding="utf-8")
    pots = {}
    for label, zval in (("Pb_d", 14), ("Ti_sv", 12), ("Zr_sv", 12), ("O", 6)):
        path = tmp_path / f"POTCAR.{label}"
        _potcar(path, label, zval)
        pots[label] = path
    return initial, final, incar, kpoints, pots


def test_setup_pto_replaces_template_a_site_without_vca(tmp_path: Path) -> None:
    module = runpy.run_path(str(SCRIPT))
    initial, final, incar, kpoints, pots = _templates(tmp_path)
    result = module["setup"](
        template_initial=initial, template_final=final, destination=tmp_path / "pto", incar=incar, kpoints=kpoints,
        a_potcar=pots["Pb_d"], b_potcar=pots["Ti_sv"], o_potcar=pots["O"],
    )
    assert read(tmp_path / "pto" / "initial" / "CONTCAR").get_chemical_symbols() == ["Pb", "Ti", "O", "O", "O"]
    assert "VCA" not in (tmp_path / "pto" / "initial" / "INCAR").read_text()
    assert result["composition"] == "PbTiO3"


def test_setup_pzt50_expands_only_b_site_and_requires_equal_zval(tmp_path: Path) -> None:
    module = runpy.run_path(str(SCRIPT))
    initial, final, incar, kpoints, pots = _templates(tmp_path)
    result = module["setup"](
        template_initial=initial, template_final=final, destination=tmp_path / "pzt", incar=incar, kpoints=kpoints,
        a_potcar=pots["Pb_d"], b_potcar=pots["Ti_sv"], o_potcar=pots["O"],
        vca_b_potcar=pots["Zr_sv"], vca_b_symbol="Zr",
    )
    expanded = read(tmp_path / "pzt" / "initial" / "POSCAR.vca")
    assert expanded.get_chemical_symbols() == ["Pb", "Ti", "Zr", "O", "O", "O"]
    assert (expanded.positions[1] == expanded.positions[2]).all()
    assert "VCA = 1.00000000 0.50000000 0.50000000 1.00000000" in (tmp_path / "pzt" / "initial" / "INCAR").read_text()
    assert result["composition"] == "Pb(Ti0.5Zr0.5)O3"

    _potcar(pots["Zr_sv"], "Zr_sv", 10)
    with pytest.raises(ValueError, match="equal ZVAL"):
        module["setup"](
            template_initial=initial, template_final=final, destination=tmp_path / "bad", incar=incar, kpoints=kpoints,
            a_potcar=pots["Pb_d"], b_potcar=pots["Ti_sv"], o_potcar=pots["O"],
            vca_b_potcar=pots["Zr_sv"], vca_b_symbol="Zr",
        )
