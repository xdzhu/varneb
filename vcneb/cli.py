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
from .config import RunConfig, prepare_run
from .optimizer_registry import optimizer_capability_matrix
from .version import __version__


_CONFIG_TEMPLATE = {
    "schema_version": 1,
    "backend": "vasp",
    "initial": "initial/POSCAR",
    "final": "final/POSCAR",
    "workdir": "runs/example",
    "n_images": 7,
    "fmax_ev_per_angstrom": 0.10,
    "k": 0.20,
    "pressure_gpa": 0.0,
    "cell_mode": "full",
    "cell_interpolation": "log_strain",
    "mapping": "auto",
    "mic": True,
    "align_translation": True,
    "minimum_distance": None,
    "maximum_deformation": None,
    "climb": False,
    "optimizer": "FIRE",
    "steps": 300,
    "image_workers": 0,
    "calculator": {
        "kind": None,
        "symbol": None,
        "parameters": {},
        "command": "",
        "factory_kwargs": {},
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


def _print_optimizers(as_json: bool) -> int:
    """List path optimizers without importing or selecting a calculator."""

    rows = optimizer_capability_matrix()
    if as_json:
        print(json.dumps(rows, indent=2, sort_keys=True))
        return 0
    print("optimizer          family       backend-independent  notes")
    print("------------------  -----------  -------------------  ------------------------------")
    for row in rows:
        print(
            f"{row['name']:<18}  {row['family']:<11}  "
            f"{row['backend_independent']:<19}  {row['notes']}"
        )
    return 0


def _doctor(as_json: bool, selected: str | None) -> int:
    rows = backend_capability_matrix()
    if selected:
        get_backend_spec(selected)
        rows = [row for row in rows if row["name"] == selected.lower()]
    executable_candidates = {
        "ase": (),
        "abacus": ("abacus",),
        "vasp": ("vasp_std", "vasp_gam", "vasp_ncl"),
        "qe": ("pw.x",),
        "lammps": ("lammps", "lmp", "lmp_mpi"),
        "cp2k": ("cp2k_shell", "cp2k_shell.psmp"),
        "abinit": ("abinit",),
    }
    for row in rows:
        module = {
            "ase": "ase.calculators.calculator",
            "abacus": "ase.calculators.abacus",
            "vasp": "ase.calculators.vasp",
            "qe": "ase.calculators.espresso",
            "lammps": "ase.calculators.lammpsrun",
            "cp2k": "ase.calculators.cp2k",
            "abinit": "ase.calculators.abinit",
        }[str(row["name"])]
        row["ase_module"] = module
        row["ase_importable"] = importlib.util.find_spec(module) is not None
        candidates = executable_candidates[str(row["name"])]
        resolved = None
        for candidate in candidates:
            resolved = shutil.which(candidate)
            if resolved is not None:
                break
        row["executable_candidates"] = list(candidates)
        row["resolved_executable"] = resolved
        row["executable_on_path"] = None if not candidates else resolved is not None
    if as_json:
        print(json.dumps(rows, indent=2, sort_keys=True))
        return 0
    for row in rows:
        available = row["executable_on_path"]
        executable_state = "n/a" if available is None else "yes" if available else "no"
        print(
            f"{row['name']}: ASE={'yes' if row['ase_importable'] else 'no'}, "
            f"PATH executable={executable_state}; "
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
    config = RunConfig.from_file(path)
    print(
        f"valid VARNEB config: backend={config.backend}, "
        f"n_images={config.n_images}, fmax={config.fmax_ev_per_angstrom:g} eV/A, "
        f"mapping={config.mapping}, cell_interpolation={config.cell_interpolation}"
    )
    return 0


def _prepare_config(path: str) -> int:
    config, report_path = prepare_run(path)
    print(f"prepared calculator-free path: {config.workdir / 'initial-vcneb.traj'}")
    print(f"wrote preflight report: {report_path}")
    return 0


def _run_config(path: str) -> int:
    """Execute the packaged material runner from an explicit reviewed config."""

    config = RunConfig.from_file(path)
    calculator = config.calculator
    if set(calculator) - {"kind", "symbol", "parameters", "command", "factory_kwargs"}:
        raise ValueError("calculator contains unsupported fields")
    kind = calculator.get("kind")
    symbol = calculator.get("symbol")
    if (kind not in {"ase_class", "factory"} or not isinstance(symbol, str)
            or ":" not in symbol or not all(symbol.split(":", 1))):
        raise ValueError("calculator.kind must be ase_class or factory, with an explicit module:attribute symbol")
    command = calculator.get("command")
    if command is not None and not isinstance(command, str):
        raise ValueError("calculator.command must be a string")
    if config.backend != "ase" and not (command and command.strip()):
        raise ValueError("named DFT backends require an explicit calculator.command")
    factory_kwargs = calculator.get("factory_kwargs", {})
    if set(factory_kwargs) & {"parameters", "command"}:
        raise ValueError("factory_kwargs must not override parameters or command")
    if config.workdir.exists() and (
        any(config.workdir.glob("image_[0-9]*"))
        or any((config.workdir / name).exists() for name in (
            "vcneb.opt.log", "vcneb.traj", "vcneb_summary.json", "ase_static_summary.json",
            "image_cache", "snapshots",
        ))
    ):
        raise FileExistsError("workdir contains a previous execution; choose a new workdir or use the explicit resume driver")
    prepare_run(path)
    argv = [
        "--initial", str(config.initial), "--final", str(config.final),
        "--workdir", str(config.workdir), "--n-images", str(config.n_images),
        "--fmax", str(config.fmax_ev_per_angstrom), "--steps", str(config.steps),
        "--k", str(config.k), "--pressure-gpa", str(config.pressure_gpa),
        "--cell-mode", config.cell_mode, "--optimizer", config.optimizer,
        "--image-workers", str(config.image_workers),
        "--image-retries", str(config.image_retries),
        "--candidate-step-retries", str(config.candidate_step_retries),
        "--cell-interpolation", config.cell_interpolation,
        "--mapping", config.mapping,
        "--backend-label", config.backend,
        "--config-source", str(Path(path).resolve()),
        "--parameters-json", json.dumps(calculator.get("parameters", {})),
        "--factory-kwargs-json", json.dumps(factory_kwargs),
    ]
    argv += ["--factory" if kind == "factory" else "--calculator", symbol]
    if command:
        argv += ["--command", command]
    if config.mic:
        argv.append("--mic")
    if config.align_translation:
        argv.append("--align-translation")
    if config.climb:
        argv.append("--climb")
    for name, value in (
        ("--climb-after", config.climb_after),
        ("--minimum-distance", config.minimum_distance),
        ("--maximum-deformation", config.maximum_deformation),
        ("--minimum-endpoint-separation", config.minimum_endpoint_separation),
        ("--maximum-cell-step", config.maximum_cell_step),
        ("--maxstep", config.maxstep),
        ("--endpoint-static-summary", config.endpoint_static_summary),
    ):
        if value is not None:
            argv += [name, str(value)]
    from .material_runner import main as run_material

    run_material(argv)
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

    optimizers = subparsers.add_parser(
        "optimizers", help="list calculator-independent path optimization strategies"
    )
    optimizers.add_argument("--json", action="store_true", help="emit machine-readable JSON")

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

    prepare = subparsers.add_parser(
        "prepare", help="build a calculator-free initial path and geometry preflight report"
    )
    prepare.add_argument("path", nargs="?", default="varneb.json")
    run = subparsers.add_parser(
        "run", help="execute a reviewed varneb.json in the current process (no scheduler submission)"
    )
    run.add_argument("path", nargs="?", default="varneb.json")
    run.add_argument("--execute", action="store_true", help="acknowledge this may launch external calculators")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    try:
        if args.command == "backends":
            return _print_backends(args.json)
        if args.command == "optimizers":
            return _print_optimizers(args.json)
        if args.command == "doctor":
            return _doctor(args.json, args.backend)
        if args.command == "init":
            return _write_template(args.path, force=args.force)
        if args.command == "validate-config":
            return _validate_config(args.path)
        if args.command == "prepare":
            return _prepare_config(args.path)
        if args.command == "run":
            if not args.execute:
                raise ValueError("run requires --execute; use prepare for calculator-free validation")
            return _run_config(args.path)
    except (FileExistsError, OSError, ValueError, json.JSONDecodeError) as exc:
        print(f"varneb: {exc}", file=sys.stderr)
        return 2
    return 0
