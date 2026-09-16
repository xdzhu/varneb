"""Prepare a Ba/Sr VCA BaTiO3 endpoint pair without exposing POTCAR data.

The physical endpoint structures keep the five-site BaTiO3 representation used
by VARNEB.  ``POSCAR.vca`` duplicates the A site for direct VASP smoke tests;
the production calculator adapter performs the same expansion in memory.
"""

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


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_json_atomic(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=path.parent, delete=False) as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.replace(path)


def _potcar_metadata(path: Path) -> dict:
    text = path.read_text(encoding="latin-1")
    title = re.search(r"^\s*TITEL\s*=\s*(.+)$", text, flags=re.MULTILINE)
    zval = re.search(r"ZVAL\s*=\s*([0-9.]+)", text)
    if title is None or zval is None:
        raise ValueError(f"cannot read TITEL/ZVAL from POTCAR dataset: {path}")
    return {"title": title.group(1).strip(), "zval": float(zval.group(1)), "sha256": _sha256(path)}


def expanded_vca_atoms(atoms: Atoms) -> Atoms:
    ba_indices = [index for index, symbol in enumerate(atoms.get_chemical_symbols()) if symbol == "Ba"]
    if len(ba_indices) != 1:
        raise ValueError(f"expected exactly one Ba virtual site, found {len(ba_indices)}")
    ba_index = ba_indices[0]
    order = [ba_index] + [index for index in range(len(atoms)) if index != ba_index]
    scaled = atoms.get_scaled_positions(wrap=False)[order]
    symbols = ["Ba", "Sr"] + [atoms[index].symbol for index in order[1:]]
    positions = [scaled[0], scaled[0], *scaled[1:]]
    return Atoms(symbols, scaled_positions=positions, cell=atoms.cell, pbc=atoms.pbc)


def _incar_with_vca(source: Path, weights: tuple[float, float]) -> str:
    lines = [line for line in source.read_text(encoding="utf-8").splitlines() if not re.match(r"^\s*VCA\s*=", line, re.I)]
    lines.append(f"VCA = {weights[0]:.8f} {weights[1]:.8f} 1.00000000 1.00000000")
    return "\n".join(lines) + "\n"


def setup(
    *,
    initial_path: Path,
    final_path: Path,
    destination: Path,
    incar: Path,
    kpoints: Path,
    ba_potcar: Path,
    sr_potcar: Path,
    ti_potcar: Path,
    o_potcar: Path,
    weights: tuple[float, float] = (0.5, 0.5),
) -> dict:
    sources = (initial_path, final_path, incar, kpoints, ba_potcar, sr_potcar, ti_potcar, o_potcar)
    for path in sources:
        if not path.is_file():
            raise FileNotFoundError(path)
    if destination.exists():
        raise FileExistsError(f"destination already exists: {destination}")
    if min(weights) < 0.0 or abs(sum(weights) - 1.0) > 1e-12:
        raise ValueError("Ba/Sr VCA weights must be nonnegative and sum to one")

    metadata = {
        "Ba": _potcar_metadata(ba_potcar),
        "Sr": _potcar_metadata(sr_potcar),
        "Ti": _potcar_metadata(ti_potcar),
        "O": _potcar_metadata(o_potcar),
    }
    if abs(metadata["Ba"]["zval"] - metadata["Sr"]["zval"]) > 1e-12:
        raise ValueError("Ba and Sr POTCAR datasets must have the same ZVAL for the isovalent baseline")

    initial, final = read(initial_path), read(final_path)
    initial_dir, final_dir = destination / "initial", destination / "final"
    initial_dir.mkdir(parents=True)
    final_dir.mkdir(parents=True)
    for directory, atoms in ((initial_dir, initial), (final_dir, final)):
        write(directory / "CONTCAR", atoms, format="vasp", direct=True, vasp5=True)
        write(directory / "POSCAR.vca", expanded_vca_atoms(atoms), format="vasp", direct=True, vasp5=True)

    (initial_dir / "INCAR").write_text(_incar_with_vca(incar, weights), encoding="utf-8")
    shutil.copy2(kpoints, initial_dir / "KPOINTS")
    with (initial_dir / "POTCAR").open("wb") as output:
        for path in (ba_potcar, sr_potcar, ti_potcar, o_potcar):
            output.write(path.read_bytes())

    manifest = {
        "status": "prepared_no_dft",
        "method": "VASP VCA",
        "composition": f"Ba{weights[0]:g}Sr{weights[1]:g}TiO3",
        "physical_atom_count": len(initial),
        "expanded_vasp_atom_count": len(initial) + 1,
        "virtual_site": {"physical_symbol": "Ba", "components": ["Ba", "Sr"], "weights": list(weights)},
        "potcar_order": ["Ba", "Sr", "Ti", "O"],
        "potcar_metadata": metadata,
        "inputs": {"incar_sha256": _sha256(initial_dir / "INCAR"), "kpoints_sha256": _sha256(initial_dir / "KPOINTS"), "potcar_sha256": _sha256(initial_dir / "POTCAR")},
    }
    _write_json_atomic(destination / "setup_manifest.json", manifest)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--initial", required=True)
    parser.add_argument("--final", required=True)
    parser.add_argument("--destination", required=True)
    parser.add_argument("--incar", required=True)
    parser.add_argument("--kpoints", required=True)
    parser.add_argument("--ba-potcar", required=True)
    parser.add_argument("--sr-potcar", required=True)
    parser.add_argument("--ti-potcar", required=True)
    parser.add_argument("--o-potcar", required=True)
    parser.add_argument("--ba-weight", type=float, default=0.5)
    parser.add_argument("--sr-weight", type=float, default=0.5)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    try:
        payload = setup(
            initial_path=Path(args.initial).resolve(), final_path=Path(args.final).resolve(),
            destination=Path(args.destination).resolve(), incar=Path(args.incar).resolve(),
            kpoints=Path(args.kpoints).resolve(), ba_potcar=Path(args.ba_potcar).resolve(),
            sr_potcar=Path(args.sr_potcar).resolve(), ti_potcar=Path(args.ti_potcar).resolve(),
            o_potcar=Path(args.o_potcar).resolve(), weights=(args.ba_weight, args.sr_weight),
        )
    except (OSError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 2
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
