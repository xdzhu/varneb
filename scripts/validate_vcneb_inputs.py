"""Static input checks for VC-NEB image directories.

This script intentionally does not launch VASP or ABACUS.  It catches common
setup mistakes before a calculator-backed VC-NEB run is submitted.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import sys

from ase.io import read


def _existing(path: Path, names: list[str]) -> Path | None:
    for name in names:
        candidate = path / name
        if candidate.exists():
            return candidate
    return None


def _image_dirs(root: Path) -> list[Path]:
    dirs = []
    for child in sorted(root.iterdir()):
        if not child.is_dir():
            continue
        if _existing(child, ["POSCAR", "POSCAR.start", "STRU.start.vasp", "CONTCAR"]) is not None:
            dirs.append(child)
    return dirs


def _read_image(path: Path):
    image_file = _existing(path, ["POSCAR", "POSCAR.start", "STRU.start.vasp", "CONTCAR"])
    if image_file is None:
        raise FileNotFoundError(f"No POSCAR-like image file in {path}")
    return image_file, read(image_file)


def _parse_abacus_input(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    for raw in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        fields = line.split()
        if len(fields) >= 2:
            values[fields[0].lower()] = " ".join(fields[1:])
    return values


def _abacus_flag_enabled(value: str) -> bool:
    first = value.split()[0].strip().lower() if value.split() else ""
    if first in {"true", "t", ".true.", "yes", "on"}:
        return True
    if first in {"false", "f", ".false.", "no", "off"}:
        return False
    try:
        return float(first) != 0.0
    except ValueError:
        return bool(first)


def _abacus_species_files(stru: Path) -> tuple[list[str], list[str]]:
    pp_files: list[str] = []
    orbital_files: list[str] = []
    if not stru.exists():
        return pp_files, orbital_files
    lines = stru.read_text(encoding="utf-8", errors="ignore").splitlines()
    section = None
    for raw in lines:
        line = raw.strip()
        upper = line.upper()
        if not line:
            continue
        if upper in {"ATOMIC_SPECIES", "NUMERICAL_ORBITAL"}:
            section = upper
            continue
        if re.match(r"^[A-Z_]+$", upper):
            section = None
        if section == "ATOMIC_SPECIES":
            fields = line.split()
            if len(fields) >= 3:
                pp_files.append(fields[2])
        elif section == "NUMERICAL_ORBITAL":
            fields = line.split()
            if fields:
                orbital_files.append(fields[0])
    return pp_files, orbital_files


def _common_structure_report(image_dirs: list[Path]) -> tuple[list[str], list[str], list[int]]:
    issues: list[str] = []
    image_files: list[str] = []
    atom_counts: list[int] = []
    reference_symbols: list[str] | None = None
    for image_dir in image_dirs:
        image_file, atoms = _read_image(image_dir)
        image_files.append(str(image_file))
        symbols = atoms.get_chemical_symbols()
        atom_counts.append(len(symbols))
        if reference_symbols is None:
            reference_symbols = symbols
        elif symbols != reference_symbols:
            issues.append(f"{image_file}: chemical symbols/order differs from first image")
        if atoms.cell.rank != 3 or atoms.get_volume() <= 0.0:
            issues.append(f"{image_file}: invalid 3D periodic cell")
    return issues, image_files, atom_counts


def validate_vasp(template_dir: Path, image_root: Path) -> dict:
    issues: list[str] = []
    for name in ["INCAR", "KPOINTS", "POTCAR"]:
        path = template_dir / name
        if not path.exists():
            issues.append(f"Missing VASP template file: {path}")
        elif path.stat().st_size == 0:
            issues.append(f"Empty VASP template file: {path}")
    image_dirs = _image_dirs(image_root)
    if not image_dirs:
        issues.append(f"No image directories with POSCAR-like files under {image_root}")
    else:
        structure_issues, image_files, atom_counts = _common_structure_report(image_dirs)
        issues.extend(structure_issues)
    return {
        "mode": "vasp",
        "template_dir": str(template_dir),
        "image_root": str(image_root),
        "n_images": len(image_dirs),
        "image_files": image_files if image_dirs else [],
        "atom_counts": atom_counts if image_dirs else [],
        "status": "ok" if not issues else "failed",
        "issues": issues,
    }


def validate_abacus(template_dir: Path, image_root: Path) -> dict:
    issues: list[str] = []
    input_path = template_dir / "INPUT"
    kpt_path = template_dir / "KPT"
    stru_path = template_dir / "STRU"
    for path in [input_path, kpt_path, stru_path]:
        if not path.exists():
            issues.append(f"Missing ABACUS template file: {path}")
        elif path.stat().st_size == 0:
            issues.append(f"Empty ABACUS template file: {path}")

    params = _parse_abacus_input(input_path)
    if input_path.exists():
        calculation = params.get("calculation", "").lower()
        if calculation and calculation != "scf":
            issues.append(f"ABACUS INPUT calculation should be scf for VC-NEB force calls, got {calculation!r}")
        for required in ["basis_type", "cal_force", "cal_stress", "out_stru"]:
            if required not in params:
                issues.append(f"ABACUS INPUT missing recommended key: {required}")
        for required in ["cal_force", "cal_stress", "out_stru"]:
            if required in params and not _abacus_flag_enabled(params[required]):
                issues.append(f"ABACUS INPUT {required} must be enabled for VC-NEB force/stress calls")

    pp_dir = template_dir / params.get("pseudo_dir", ".")
    orbital_dir = template_dir / params.get("orbital_dir", ".")
    pp_files, orbital_files = _abacus_species_files(stru_path)
    for name in pp_files:
        if not (pp_dir / name).exists():
            issues.append(f"Missing pseudopotential referenced by STRU: {pp_dir / name}")
    if params.get("basis_type", "").lower() == "lcao":
        for name in orbital_files:
            if not (orbital_dir / name).exists():
                issues.append(f"Missing orbital referenced by STRU: {orbital_dir / name}")

    image_dirs = _image_dirs(image_root)
    if not image_dirs:
        issues.append(f"No image directories with POSCAR-like files under {image_root}")
    else:
        structure_issues, image_files, atom_counts = _common_structure_report(image_dirs)
        issues.extend(structure_issues)
    return {
        "mode": "abacus",
        "template_dir": str(template_dir),
        "image_root": str(image_root),
        "n_images": len(image_dirs),
        "image_files": image_files if image_dirs else [],
        "atom_counts": atom_counts if image_dirs else [],
        "status": "ok" if not issues else "failed",
        "issues": issues,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["vasp", "abacus"], required=True)
    parser.add_argument("--template-dir", required=True, help="Directory containing VASP or ABACUS template inputs")
    parser.add_argument("--image-root", required=True, help="Directory containing per-image subdirectories")
    parser.add_argument("--json", action="store_true", help="Print machine-readable JSON")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    template_dir = Path(args.template_dir).resolve()
    image_root = Path(args.image_root).resolve()
    if args.mode == "vasp":
        report = validate_vasp(template_dir, image_root)
    else:
        report = validate_abacus(template_dir, image_root)

    if args.json:
        print(json.dumps(report, indent=2, ensure_ascii=False))
    else:
        print(f"mode={report['mode']}")
        print(f"status={report['status']}")
        print(f"n_images={report['n_images']}")
        if report["atom_counts"]:
            print("atom_counts=" + " ".join(str(x) for x in report["atom_counts"]))
        for issue in report["issues"]:
            print(f"issue: {issue}")
    if report["status"] != "ok":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
