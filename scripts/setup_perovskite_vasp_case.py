"""Prepare a five-site ABO3 VASP case, optionally with B-site VCA."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys
import tempfile

from ase import Atoms
from ase.io import read, write


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def potcar_metadata(path: Path) -> dict:
    text = path.read_text(encoding="latin-1")
    title = re.search(r"^\s*TITEL\s*=\s*(.+)$", text, flags=re.MULTILINE)
    zval = re.search(r"ZVAL\s*=\s*([0-9.]+)", text)
    if title is None or zval is None:
        raise ValueError(f"cannot read TITEL/ZVAL from POTCAR dataset: {path}")
    return {"title": title.group(1).strip(), "zval": float(zval.group(1)), "sha256": sha256(path)}


def write_json(path: Path, payload: dict) -> None:
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(path)


def replace_ab_sites(atoms: Atoms, *, template_a: str, template_b: str, a_symbol: str, b_symbol: str) -> Atoms:
    symbols = atoms.get_chemical_symbols()
    a_indices = [index for index, symbol in enumerate(symbols) if symbol == template_a]
    b_indices = [index for index, symbol in enumerate(symbols) if symbol == template_b]
    if len(a_indices) != 1 or len(b_indices) != 1:
        raise ValueError(f"expected one {template_a} A site and one {template_b} B site")
    symbols[a_indices[0]] = a_symbol
    symbols[b_indices[0]] = b_symbol
    result = atoms.copy()
    result.set_chemical_symbols(symbols)
    return result


def expand_b_vca(atoms: Atoms, *, b_symbol: str, components: tuple[str, str]) -> Atoms:
    b_indices = [index for index, symbol in enumerate(atoms.get_chemical_symbols()) if symbol == b_symbol]
    if len(b_indices) != 1:
        raise ValueError(f"expected exactly one physical {b_symbol} B site")
    b_index = b_indices[0]
    symbols, positions = [], []
    for index, (symbol, position) in enumerate(zip(atoms.get_chemical_symbols(), atoms.get_scaled_positions(wrap=False))):
        for expanded_symbol in (components if index == b_index else (symbol,)):
            symbols.append(expanded_symbol)
            positions.append(position)
    return Atoms(symbols, scaled_positions=positions, cell=atoms.cell, pbc=atoms.pbc)


def incar_with_vca(incar: Path, weights: tuple[float, float] | None) -> str:
    lines = [line for line in incar.read_text(encoding="utf-8").splitlines() if not re.match(r"^\s*VCA\s*=", line, re.I)]
    if weights is not None:
        lines.append(f"VCA = 1.00000000 {weights[0]:.8f} {weights[1]:.8f} 1.00000000")
    return "\n".join(lines) + "\n"


def setup(
    *,
    template_initial: Path,
    template_final: Path,
    destination: Path,
    incar: Path,
    kpoints: Path,
    a_potcar: Path,
    b_potcar: Path,
    o_potcar: Path,
    template_a: str = "Ba",
    template_b: str = "Ti",
    a_symbol: str = "Pb",
    b_symbol: str = "Ti",
    vca_b_potcar: Path | None = None,
    vca_b_symbol: str | None = None,
    vca_weights: tuple[float, float] = (0.5, 0.5),
) -> dict:
    required = (template_initial, template_final, incar, kpoints, a_potcar, b_potcar, o_potcar)
    if any(not path.is_file() for path in required):
        raise FileNotFoundError("one or more required templates or POTCAR datasets are missing")
    is_vca = vca_b_potcar is not None or vca_b_symbol is not None
    if is_vca and (vca_b_potcar is None or vca_b_symbol is None):
        raise ValueError("B-site VCA needs both vca_b_potcar and vca_b_symbol")
    if is_vca and (min(vca_weights) < 0 or abs(sum(vca_weights) - 1.0) > 1e-12):
        raise ValueError("VCA weights must be nonnegative and sum to one")
    if destination.exists():
        raise FileExistsError(f"destination exists: {destination}")

    metadata = {"A": potcar_metadata(a_potcar), "B": potcar_metadata(b_potcar), "O": potcar_metadata(o_potcar)}
    if is_vca:
        metadata["B2"] = potcar_metadata(vca_b_potcar)
        if abs(metadata["B"]["zval"] - metadata["B2"]["zval"]) > 1e-12:
            raise ValueError("B-site VCA datasets must have equal ZVAL")

    initial = replace_ab_sites(read(template_initial), template_a=template_a, template_b=template_b, a_symbol=a_symbol, b_symbol=b_symbol)
    final = replace_ab_sites(read(template_final), template_a=template_a, template_b=template_b, a_symbol=a_symbol, b_symbol=b_symbol)
    initial_dir, final_dir = destination / "initial", destination / "final"
    initial_dir.mkdir(parents=True)
    final_dir.mkdir(parents=True)
    for directory, atoms in ((initial_dir, initial), (final_dir, final)):
        write(directory / "CONTCAR", atoms, format="vasp", direct=True, vasp5=True)
        if is_vca:
            write(directory / "POSCAR.vca", expand_b_vca(atoms, b_symbol=b_symbol, components=(b_symbol, vca_b_symbol)), format="vasp", direct=True, vasp5=True)

    (initial_dir / "INCAR").write_text(incar_with_vca(incar, vca_weights if is_vca else None), encoding="utf-8")
    shutil.copy2(kpoints, initial_dir / "KPOINTS")
    with (initial_dir / "POTCAR").open("wb") as handle:
        for path in ((a_potcar, b_potcar, vca_b_potcar, o_potcar) if is_vca else (a_potcar, b_potcar, o_potcar)):
            handle.write(path.read_bytes())
    manifest = {
        "status": "prepared_no_dft", "method": "VASP B-site VCA" if is_vca else "VASP ordered endpoint",
        "composition": f"{a_symbol}({b_symbol}{vca_weights[0]:g}{vca_b_symbol}{vca_weights[1]:g})O3" if is_vca else f"{a_symbol}{b_symbol}O3",
        "physical_atom_count": len(initial), "expanded_vasp_atom_count": len(initial) + int(is_vca),
        "virtual_site": ({"physical_symbol": b_symbol, "components": [b_symbol, vca_b_symbol], "weights": list(vca_weights)} if is_vca else None),
        "potcar_metadata": metadata,
        "inputs": {"incar_sha256": sha256(initial_dir / "INCAR"), "kpoints_sha256": sha256(initial_dir / "KPOINTS"), "potcar_sha256": sha256(initial_dir / "POTCAR")},
    }
    write_json(destination / "setup_manifest.json", manifest)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("template-initial", "template-final", "destination", "incar", "kpoints", "a-potcar", "b-potcar", "o-potcar"):
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument("--template-a", default="Ba")
    parser.add_argument("--template-b", default="Ti")
    parser.add_argument("--a-symbol", default="Pb")
    parser.add_argument("--b-symbol", default="Ti")
    parser.add_argument("--vca-b-potcar")
    parser.add_argument("--vca-b-symbol")
    parser.add_argument("--vca-b-weight", type=float, default=0.5)
    parser.add_argument("--vca-b2-weight", type=float, default=0.5)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        result = setup(
            template_initial=Path(args.template_initial).resolve(), template_final=Path(args.template_final).resolve(),
            destination=Path(args.destination).resolve(), incar=Path(args.incar).resolve(), kpoints=Path(args.kpoints).resolve(),
            a_potcar=Path(args.a_potcar).resolve(), b_potcar=Path(args.b_potcar).resolve(), o_potcar=Path(args.o_potcar).resolve(),
            template_a=args.template_a, template_b=args.template_b, a_symbol=args.a_symbol, b_symbol=args.b_symbol,
            vca_b_potcar=Path(args.vca_b_potcar).resolve() if args.vca_b_potcar else None, vca_b_symbol=args.vca_b_symbol,
            vca_weights=(args.vca_b_weight, args.vca_b2_weight),
        )
    except (OSError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
