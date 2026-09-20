"""Command-line entry point for VARNEB.

The import package remains ``vcneb`` for source compatibility; ``varneb`` is
the public distribution/project name.  The legacy ``vcneb`` console alias is
kept so existing scripts continue to work.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import shutil
import sys

from .backends import backend_capability_matrix, get_backend_spec
from .version import __version__


_CONFIG_TEMPLATE = {
    "schema_version": 1,
    "backend": "vasp",
    "initial": "initial/POSCAR",
    "final": "final/POSCAR",
    "workdir": "runs/example",
    "n_images": 7,
    "fmax_ev_per_angstrom": 0.10,
    "calculator": {
        "parameters": {},
        "command": "",
    },
}


def _print_backends(as_json: bool) -> int:
    rows = backend_capability_matrix()
    if as_json:
        print(json.dumps(rows, indent=2, sort_keys=True))
        return 0
    print("backend  executable      variable-cell  status       adapter")
    print("-------  ---------------  -------------  -----------  ----------------")
    for row in rows:
        print(
            f"{row['name']:<8} {row['executable']:<15} "
            f"{'yes' if row['variable_cell'] else 'no':<13} "
            f"{row['status']:<11} {row['python_adapter']}"
        )
    return 0


def _doctor(as_json: bool, selected: str | None) -> int:
    rows = backend_capability_matrix()
    if selected:
        get_backend_spec(selected)
        rows = [row for row in rows if row["name"] == selected.lower()]
    for row in rows:
        module = {
            "abacus": "ase.calculators.abacus",
            "vasp": "ase.calculators.vasp",
            "qe": "ase.calculators.espresso",
            "lammps": "ase.calculators.lammpsrun",
            "cp2k": "ase.calculators.cp2k",
        }[str(row["name"])]
        row["ase_module"] = module
        row["ase_importable"] = importlib.util.find_spec(module) is not None
        row["executable_on_path"] = shutil.which(str(row["executable"])) is not None
    if as_json:
        print(json.dumps(rows, indent=2, sort_keys=True))
        return 0
    for row in rows:
        print(
            f"{row['name']}: ASE={'yes' if row['ase_importable'] else 'no'}, "
            f"PATH executable={'yes' if row['executable_on_path'] else 'no'}; "
            f"{row['notes']}"
        )
    return 0


def _write_template(path: str, *, force: bool) -> int:
    target = Path(path)
    if target.exists() and not force:
        raise FileExistsError(f"refusing to overwrite existing file: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(_CONFIG_TEMPLATE, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {target}")
    return 0


def _validate_config(path: str) -> int:
    source = Path(path)
    data = json.loads(source.read_text(encoding="utf-8"))
    if data.get("schema_version") != 1:
        raise ValueError("config schema_version must be 1")
    backend = data.get("backend")
    get_backend_spec(backend)
    for key in ("initial", "final", "workdir"):
        if not isinstance(data.get(key), str) or not data[key].strip():
            raise ValueError(f"config field {key!r} must be a non-empty path")
    try:
        n_images = int(data.get("n_images"))
        fmax = float(data.get("fmax_ev_per_angstrom"))
    except (TypeError, ValueError) as exc:
        raise ValueError("n_images and fmax_ev_per_angstrom must be numeric") from exc
    if n_images < 3:
        raise ValueError("n_images must include two fixed endpoints and at least one interior image")
    if not 0.0 < fmax:
        raise ValueError("fmax_ev_per_angstrom must be positive")
    calculator = data.get("calculator", {})
    if not isinstance(calculator, dict) or not isinstance(calculator.get("parameters", {}), dict):
        raise ValueError("calculator.parameters must be a mapping")
    print(f"valid VARNEB config: backend={backend}, n_images={n_images}, fmax={fmax:g} eV/A")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="varneb",
        description="VARNEB: calculator-agnostic variable-cell nudged elastic band tools.",
    )
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command")

    backends = subparsers.add_parser(
        "backends", help="list calculator backends and their VARNEB status"
    )
    backends.add_argument("--json", action="store_true", help="emit machine-readable JSON")

    doctor = subparsers.add_parser(
        "doctor", help="check optional ASE adapter imports and local executables"
    )
    doctor.add_argument("--backend", help="inspect one backend only")
    doctor.add_argument("--json", action="store_true", help="emit machine-readable JSON")

    init = subparsers.add_parser(
        "init", help="write a minimal, calculator-agnostic run configuration"
    )
    init.add_argument("path", nargs="?", default="varneb.json")
    init.add_argument("--force", action="store_true", help="replace an existing template")

    validate = subparsers.add_parser(
        "validate-config", help="validate a varneb.json without launching a calculator"
    )
    validate.add_argument("path", nargs="?", default="varneb.json")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.command == "backends":
            return _print_backends(args.json)
        if args.command == "doctor":
            return _doctor(args.json, args.backend)
        if args.command == "init":
            return _write_template(args.path, force=args.force)
        if args.command == "validate-config":
            return _validate_config(args.path)
    except (FileExistsError, OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"varneb: {exc}", file=sys.stderr)
        return 2
    return 0
