"""Compare the VASP basis and raw stress of GaN basin relaxations and statics.

This is a read-only, case-specific audit. It does not infer a basis-set
convergence limit from one pair of runs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
COMMON_ELECTRONIC_KEYS = (
    "ENCUT", "SIGMA", "EDIFF", "SYMPREC", "GGA", "PREC", "ISMEAR", "ISYM"
)
PRESSURE_KBAR = 457.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _last(pattern: str, text: str, *, flags: int = re.MULTILINE) -> re.Match[str]:
    matches = list(re.finditer(pattern, text, flags))
    if not matches:
        raise ValueError(f"OUTCAR field missing: {pattern}")
    return matches[-1]


def _grid(text: str, prefix: str) -> list[int]:
    suffix = "F" if prefix.endswith("F") else ""
    match = _last(
        rf"^\s*dimension x,y,z NGX{suffix}\s*=\s*(\d+)\s+"
        rf"NGY{suffix}\s*=\s*(\d+)\s+"
        rf"NGZ{suffix}\s*=\s*(\d+)",
        text,
    )
    return [int(value) for value in match.groups()]


def parse_outcar(text: str) -> dict:
    kpoint_counts = [
        (int(match.group(1)), int(match.group(2)))
        for match in re.finditer(
            r"^\s*k-point\s+(\d+)\s+:.*?plane waves:\s*(\d+)\s*$",
            text,
            re.MULTILINE,
        )
    ]
    if not kpoint_counts or [number for number, _ in kpoint_counts] != list(
        range(1, len(kpoint_counts) + 1)
    ):
        raise ValueError("k-point plane-wave list is absent or incomplete")
    values = [count for _, count in kpoint_counts]
    stress = [
        float(value)
        for value in _last(
            rf"^\s*in kB\s+({FLOAT})\s+({FLOAT})\s+({FLOAT})\s+"
            rf"({FLOAT})\s+({FLOAT})\s+({FLOAT})\s*$",
            text,
        ).groups()
    ]
    pressure_line = _last(
        rf"external pressure\s*=\s*({FLOAT})\s*kB\s+"
        rf"Pullay stress\s*=\s*({FLOAT})\s*kB",
        text,
    )
    external_pressure, pstress = (float(value) for value in pressure_line.groups())
    mean_raw_pressure = sum(stress[:3]) / 3.0
    if abs(external_pressure - (mean_raw_pressure - pstress)) > 0.02:
        raise ValueError("VASP external-pressure/PSTRESS/raw-stress identity fails")
    coarse = _grid(text, "NGX")
    fine = _grid(text, "NGXF")
    nplwv = int(_last(r"total plane-waves\s+NPLWV\s*=\s*(\d+)", text).group(1))
    if nplwv != coarse[0] * coarse[1] * coarse[2]:
        raise ValueError("NPLWV is inconsistent with the reported coarse FFT grid")
    volumes = [
        float(match.group(1))
        for match in re.finditer(rf"volume of cell\s*:\s*({FLOAT})", text)
    ]
    if not volumes:
        raise ValueError("OUTCAR has no volume records")
    energy = float(_last(rf"free\s+energy\s+TOTEN\s*=\s*({FLOAT})\s+eV", text).group(1))
    return {
        "initial_volume_A3_reported": volumes[0],
        "last_volume_A3_reported": volumes[-1],
        "coarse_fft_grid": coarse,
        "fine_fft_grid": fine,
        "NPLWV_fft_grid_points": nplwv,
        "kpoint_plane_wave_counts": values,
        "n_kpoints": len(values),
        "first_kpoint_plane_waves": values[0],
        "sum_kpoint_plane_waves": sum(values),
        "last_TOTEN_eV_per_cell": energy,
        "last_raw_stress_kbar_vasp_order": stress,
        "last_max_raw_stress_residual_from_45p7GPa_kbar": max(
            [abs(value - PRESSURE_KBAR) for value in stress[:3]]
            + [abs(value) for value in stress[3:]]
        ),
        "last_external_pressure_kbar": external_pressure,
        "PSTRESS_kbar_reported": pstress,
    }


def parse_incar(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        active = line.split("#", 1)[0].strip()
        if "=" in active:
            key, value = active.split("=", 1)
            values[key.strip().upper()] = value.strip().upper()
    return values


def compact_outcar(record: dict) -> dict:
    counts = record["kpoint_plane_wave_counts"]
    summary = {key: value for key, value in record.items()
               if key != "kpoint_plane_wave_counts"}
    summary["min_kpoint_plane_waves"] = min(counts)
    summary["max_kpoint_plane_waves"] = max(counts)
    summary["kpoint_plane_wave_counts_sha256"] = hashlib.sha256(
        json.dumps(counts, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    return summary


def audit_pair(native_dir: Path, static_dir: Path) -> dict:
    files = {
        "native_OUTCAR": native_dir / "OUTCAR",
        "native_INCAR": native_dir / "INCAR",
        "native_KPOINTS": native_dir / "KPOINTS",
        "native_POTCAR": native_dir / "POTCAR",
        "native_CONTCAR": native_dir / "CONTCAR",
        "static_OUTCAR": static_dir / "OUTCAR",
        "static_INCAR": static_dir / "INCAR",
        "static_KPOINTS": static_dir / "KPOINTS",
        "static_POTCAR": static_dir / "POTCAR",
        "static_POSCAR": static_dir / "POSCAR",
    }
    hashes = {label: sha256(path) for label, path in files.items()}
    if (hashes["native_CONTCAR"] != hashes["static_POSCAR"]
            or hashes["native_KPOINTS"] != hashes["static_KPOINTS"]
            or hashes["native_POTCAR"] != hashes["static_POTCAR"]):
        raise ValueError("terminal geometry, KPOINTS, or POTCAR differs")
    native_input = parse_incar(files["native_INCAR"])
    static_input = parse_incar(files["static_INCAR"])
    if (any(native_input.get(key) != static_input.get(key)
            for key in COMMON_ELECTRONIC_KEYS)
            or float(native_input["ENCUT"]) != 600.0
            or float(native_input["PSTRESS"]) != PRESSURE_KBAR
            or "PSTRESS" in static_input):
        raise ValueError("electronic settings, 600-eV ENCUT, or pressure contract differs")
    native = parse_outcar(files["native_OUTCAR"].read_text(encoding="utf-8", errors="replace"))
    static = parse_outcar(files["static_OUTCAR"].read_text(encoding="utf-8", errors="replace"))
    if (native["n_kpoints"] != static["n_kpoints"]
            or native["n_kpoints"] != 384
            or abs(native["PSTRESS_kbar_reported"] - PRESSURE_KBAR) > 0.01
            or abs(static["PSTRESS_kbar_reported"]) > 0.01
            or abs(native["last_volume_A3_reported"] - static["initial_volume_A3_reported"]) > 0.01):
        raise ValueError("k-point count, printed pressure, or geometry volume differs")
    counts = list(zip(native["kpoint_plane_wave_counts"], static["kpoint_plane_wave_counts"]))
    return {
        "native": compact_outcar(native),
        "fresh_static": compact_outcar(static),
        "kpoints_with_changed_plane_wave_count": sum(a != b for a, b in counts),
        "static_minus_native_TOTEN_meV_per_cell": 1000.0 * (
            static["last_TOTEN_eV_per_cell"] - native["last_TOTEN_eV_per_cell"]
        ),
        "same_terminal_POSCAR_bytes": True,
        "same_electronic_INCAR_keys": list(COMMON_ELECTRONIC_KEYS),
        "same_KPOINTS_and_POTCAR_bytes": True,
        "source_sha256": hashes,
    }


def audit(root: Path) -> dict:
    cases = {}
    for label, native, static in (
        ("B1", "grid_um_vz", "negative"),
        ("B4", "grid_up_vz", "positive"),
    ):
        cases[label] = audit_pair(
            root / "outputs" / "gan_600eV_basin_uninterrupted_20260928" / "cases" / native,
            root / "outputs" / f"gan_600eV_basin_final_static_{static}_20260928" / "case",
        )
    return {
        "status": "raw_basis_history_difference_measured_not_cutoff_convergence",
        "material": "GaN",
        "pressure_GPa": 45.7,
        "ENCUT_eV_both": 600,
        "cases": cases,
        "interpretation": (
            "The native ISIF=3 run and fresh static use identical terminal geometry, "
            "electronic INCAR keys, KPOINTS and POTCAR, but their plane-wave "
            "membership differs. This is consistent with VASP's documented "
            "fixed-basis cell-relaxation and fresh-static basis reset. It is "
            "not by itself a full Pulay-stress or cutoff-convergence bound."
        ),
        "reference": "https://vasp.at/wiki/ISTART",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.root)
    if args.output.exists():
        raise FileExistsError(args.output)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({label: {
        "changed_kpoints": item["kpoints_with_changed_plane_wave_count"],
        "delta_E_meV_per_cell": item["static_minus_native_TOTEN_meV_per_cell"],
    } for label, item in result["cases"].items()}))


if __name__ == "__main__":
    main()
