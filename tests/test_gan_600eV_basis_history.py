from pathlib import Path

import pytest

from scripts.audit_gan_600eV_basis_history import audit_pair, parse_outcar


def outcar(*, initial_volume: float, final_volume: float, pw: int,
           stress: tuple[float, float, float], pstress: float,
           energy: float, n_kpoints: int = 384) -> str:
    external = sum(stress) / 3.0 - pstress
    lines = [
        " total plane-waves  NPLWV =  8",
        " dimension x,y,z NGX = 2 NGY = 2 NGZ = 2",
        " dimension x,y,z NGXF= 4 NGYF= 4 NGZF= 4",
        f" volume of cell : {initial_volume:.4f}",
    ]
    lines.extend(
        f" k-point {index:3d} : 0.0000 0.0000 0.0000  plane waves: {pw}"
        for index in range(1, n_kpoints + 1)
    )
    lines.extend([
        f" volume of cell : {final_volume:.4f}",
        f" in kB {stress[0]:.5f} {stress[1]:.5f} {stress[2]:.5f} 0.0 0.0 0.0",
        f" external pressure = {external:.2f} kB  Pullay stress = {pstress:.2f} kB",
        f" free  energy   TOTEN  = {energy:.8f} eV",
    ])
    return "\n".join(lines) + "\n"


def test_outcar_parser_uses_raw_stress_not_pressure_shift():
    parsed = parse_outcar(outcar(
        initial_volume=40, final_volume=35, pw=100,
        stress=(456, 455, 457), pstress=457, energy=-20,
    ))
    assert parsed["coarse_fft_grid"] == [2, 2, 2]
    assert parsed["fine_fft_grid"] == [4, 4, 4]
    assert parsed["last_external_pressure_kbar"] == -1
    assert parsed["last_max_raw_stress_residual_from_45p7GPa_kbar"] == 2
    assert parsed["n_kpoints"] == 384


def test_outcar_parser_rejects_inconsistent_pressure():
    raw = outcar(initial_volume=40, final_volume=35, pw=100,
                 stress=(456, 455, 457), pstress=457, energy=-20)
    with pytest.raises(ValueError, match="identity fails"):
        parse_outcar(raw.replace("external pressure = -1.00", "external pressure = 456.00"))


def test_pair_audit_requires_same_geometry_and_electronics(tmp_path: Path):
    native = tmp_path / "native"
    static = tmp_path / "static"
    native.mkdir()
    static.mkdir()
    common = "ENCUT = 600\nSIGMA = 0.05\nEDIFF = 1e-7\nSYMPREC = 1e-4\nGGA = PE\nPREC = Accurate\nISMEAR = 0\nISYM = -1\n"
    (native / "INCAR").write_text(common + "PSTRESS = 457\n", encoding="utf-8")
    (static / "INCAR").write_text(common, encoding="utf-8")
    for name in ("KPOINTS", "POTCAR"):
        (native / name).write_text(name, encoding="utf-8")
        (static / name).write_text(name, encoding="utf-8")
    (native / "CONTCAR").write_text("terminal geometry", encoding="utf-8")
    (static / "POSCAR").write_text("terminal geometry", encoding="utf-8")
    (native / "OUTCAR").write_text(outcar(
        initial_volume=40, final_volume=35, pw=100,
        stress=(456, 455, 457), pstress=457, energy=-20,
    ), encoding="utf-8")
    (static / "OUTCAR").write_text(outcar(
        initial_volume=35, final_volume=35, pw=101,
        stress=(453, 453, 453), pstress=0, energy=-19.993,
    ), encoding="utf-8")
    result = audit_pair(native, static)
    assert result["kpoints_with_changed_plane_wave_count"] == 384
    assert result["static_minus_native_TOTEN_meV_per_cell"] == pytest.approx(7)
    (static / "INCAR").write_text(common.replace("ENCUT = 600", "ENCUT = 1000"), encoding="utf-8")
    with pytest.raises(ValueError, match="electronic settings"):
        audit_pair(native, static)
