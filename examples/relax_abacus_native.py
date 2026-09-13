"""Run an ABACUS-native atomic/cell relaxation and export a VASP endpoint.

Unlike ``relax_abacus_endpoint.py``, this driver does not attach an ASE
optimizer.  ABACUS receives ``calculation cell-relax`` and performs the BFGS
updates itself, writing ``OUT.ABACUS/STRU_ION_D`` and its native relaxation log.
ASE is used only for input/output conversion and lightweight geometry checks.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
from ase.data import atomic_masses, atomic_numbers
from ase.io import read, write


BOHR_PER_ANGSTROM = 1.8897261258369282


def _parse_species_files(values: list[str], option: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"{option} expects SPECIES=FILE, got {value!r}")
        species, filename = value.split("=", 1)
        if not species or not filename:
            raise ValueError(f"{option} expects non-empty SPECIES and FILE, got {value!r}")
        result[species] = filename
    return result


def _args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--structure", required=True, help="ASE-readable endpoint structure")
    p.add_argument("--workdir", required=True)
    p.add_argument("--command", required=True, help="ABACUS command, usually srun ... abacus")
    p.add_argument("--pseudo-dir", required=True)
    p.add_argument("--basis-dir", required=True)
    p.add_argument("--pp", action="append", default=[], metavar="SPECIES=FILE")
    p.add_argument("--basis", action="append", default=[], metavar="SPECIES=FILE")
    p.add_argument("--ecutwfc", type=float, default=100.0)
    p.add_argument("--scf-thr", type=float, default=1e-8)
    p.add_argument("--scf-nmax", type=int, default=150)
    p.add_argument("--kpts", type=int, nargs=3, default=[2, 2, 2], metavar=("NX", "NY", "NZ"))
    p.add_argument("--relax-nmax", type=int, default=200)
    p.add_argument("--force-thr", type=float, default=0.02)
    p.add_argument("--stress-thr", type=float, default=0.1, help="ABACUS stress threshold in kbar")
    p.add_argument("--mixing-beta", type=float, default=0.3)
    p.add_argument(
        "--relax-method",
        choices=["bfgs", "bfgs_trad", "cg", "sd", "cg_bfgs", "fire"],
        default="bfgs",
    )
    return p.parse_args()


def _write_stru(atoms, path: Path, pp: dict[str, str], basis: dict[str, str], pseudo_dir: Path, basis_dir: Path) -> None:
    symbols = atoms.get_chemical_symbols()
    species = list(dict.fromkeys(symbols))
    missing = [s for s in species if s not in pp or s not in basis]
    if missing:
        raise ValueError(f"Missing pseudopotential/orbital mapping for {missing}")
    lines = ["ATOMIC_SPECIES"]
    for s in species:
        source = pseudo_dir / pp[s]
        if not source.is_file():
            raise FileNotFoundError(source)
        lines.append(f"{s} {atomic_masses[atomic_numbers[s]]:.6f} {source.name}")
    lines += ["", "NUMERICAL_ORBITAL"]
    for s in species:
        source = basis_dir / basis[s]
        if not source.is_file():
            raise FileNotFoundError(source)
        lines.append(source.name)
    lines += ["", "LATTICE_CONSTANT", f"{BOHR_PER_ANGSTROM:.15f}", "", "LATTICE_VECTORS"]
    for vector in np.asarray(atoms.cell.array, dtype=float):
        lines.append(" ".join(f"{x:.15f}" for x in vector))
    lines += ["", "ATOMIC_POSITIONS", "Direct", ""]
    scaled = atoms.get_scaled_positions(wrap=False)
    for s in species:
        indices = [i for i, value in enumerate(symbols) if value == s]
        lines += [s, "0.0", str(len(indices))]
        for i in indices:
            q = scaled[i]
            lines.append(" ".join(f"{x:.15f}" for x in q) + " 1 1 1 mag 0.0")
        lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")
    for s in species:
        pseudo_source = (pseudo_dir / pp[s]).resolve()
        pseudo_target = (path.parent / Path(pp[s]).name).resolve()
        if pseudo_source != pseudo_target:
            shutil.copy2(pseudo_source, pseudo_target)
        basis_source = (basis_dir / basis[s]).resolve()
        basis_target = (path.parent / Path(basis[s]).name).resolve()
        if basis_source != basis_target:
            shutil.copy2(basis_source, basis_target)


def _write_input(path: Path, args: argparse.Namespace, ntype: int) -> None:
    text = f"""INPUT_PARAMETERS
calculation       cell-relax
stru_file         STRU
basis_type        lcao
ks_solver         genelpa
dft_functional    pbe
ecutwfc           {args.ecutwfc:g}
scf_thr           {args.scf_thr:g}
scf_nmax          {args.scf_nmax}
mixing_type       pulay
mixing_beta       {args.mixing_beta:g}
smearing_method   gauss
smearing_sigma    0.008
symmetry          0
ntype             {ntype}
cal_force         1
cal_stress        1
out_stru          1
out_chg           0
relax_method      {args.relax_method}
relax_new         1
relax_scale_force 0.5
relax_bfgs_rmax    0.2
relax_bfgs_rmin    1e-05
relax_bfgs_init    0.5
relax_nmax        {args.relax_nmax}
force_thr_ev      {args.force_thr:g}
stress_thr        {args.stress_thr:g}
fixed_axes        None
"""
    path.write_text(text, encoding="utf-8")


def _write_kpt(path: Path, kpts: list[int]) -> None:
    path.write_text("K_POINTS\n0\nGamma\n" + " ".join(str(x) for x in kpts) + " 0 0 0\n", encoding="utf-8")


def _parse_stru(path: Path):
    """Read the small subset of ABACUS STRU written by this driver."""
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    lc_idx = lines.index("LATTICE_CONSTANT")
    lattice_constant_bohr = float(lines[lc_idx + 1].split()[0])
    lv = lines.index("LATTICE_VECTORS")
    # ABACUS stores lattice vectors as dimensionless rows multiplied by the
    # lattice constant (in Bohr).  Return ASE's Angstrom cell.
    cell = np.array([[float(x) for x in lines[lv + i + 1].split()[:3]] for i in range(3)], dtype=float)
    cell *= lattice_constant_bohr / BOHR_PER_ANGSTROM
    ap = lines.index("ATOMIC_POSITIONS")
    mode = lines[ap + 1].lower()
    species = []
    positions = []
    i = ap + 2
    while i < len(lines):
        if not lines[i]:
            i += 1
            continue
        # ABACUS writes labels as ``Hf #label`` and includes a magnetic
        # moment line before the atom count.  Keep only the species token,
        # then locate the first integer count below it.
        name = lines[i].split()[0]
        if name in {"LATTICE_CONSTANT", "LATTICE_VECTORS"}:
            break
        j = i + 1
        count = None
        while j < len(lines):
            tokens = lines[j].split()
            if tokens:
                try:
                    count = int(tokens[0])
                    break
                except ValueError:
                    pass
            j += 1
        if count is None:
            break
        for row in lines[j + 1 : j + 1 + count]:
            values = [float(x) for x in row.split()[:3]]
            positions.append(values)
            species.append(name)
        i = j + 1 + count
    positions = np.asarray(positions, dtype=float)
    if mode.startswith("cart"):
        positions = positions @ np.linalg.inv(cell)
    from ase import Atoms

    return Atoms(symbols=species, scaled_positions=positions, cell=cell, pbc=True)


def _read_structure(path: Path):
    """Read either an ASE structure or an ABACUS STRU/STRU_ION_D restart."""
    text = path.read_text(encoding="utf-8", errors="replace")
    if "ATOMIC_POSITIONS" in text and "LATTICE_VECTORS" in text:
        return _parse_stru(path)
    return read(path)


def _last_force_stress(log_text: str, natoms: int) -> tuple[float | None, float | None]:
    """Extract max |force| and |stress| from the last ABACUS output blocks."""
    force_max = None
    force_mark = "TOTAL-FORCE (eV/Angstrom)"
    force_start = log_text.rfind(force_mark)
    if force_start >= 0:
        block = log_text[force_start + len(force_mark) :]
        values = []
        for line in block.splitlines():
            fields = line.split()
            if len(fields) < 4:
                if values:
                    break
                continue
            try:
                values.append([float(fields[1]), float(fields[2]), float(fields[3])])
            except (TypeError, ValueError):
                if values:
                    break
            if len(values) >= natoms:
                break
        if values:
            force_max = float(np.linalg.norm(np.asarray(values), axis=1).max())

    stress_max = None
    stress_mark = "TOTAL-STRESS (KBAR)"
    stress_start = log_text.rfind(stress_mark)
    if stress_start >= 0:
        block = log_text[stress_start + len(stress_mark) :]
        rows = []
        for line in block.splitlines():
            fields = line.split()
            if len(fields) < 3:
                if rows:
                    break
                continue
            try:
                rows.append([float(fields[0]), float(fields[1]), float(fields[2])])
            except (TypeError, ValueError):
                if rows:
                    break
            if len(rows) >= 3:
                break
        if len(rows) == 3:
            stress_max = float(np.max(np.abs(np.asarray(rows))))
    return force_max, stress_max


def _summary(workdir: Path, args: argparse.Namespace, returncode: int, atoms) -> dict:
    log_candidates = [workdir / "OUT.ABACUS" / "running_cell-relax.log", workdir / "running_cell-relax.log"]
    log = next((p for p in log_candidates if p.is_file()), None)
    log_text = log.read_text(encoding="utf-8", errors="replace") if log else ""
    final_stru = workdir / "OUT.ABACUS" / "STRU_ION_D"
    converged = "Lattice relaxation is converged" in log_text or "Relaxation is converged" in log_text
    force_matches = re.findall(r"Largest gradient in force is\s+([-+0-9.eE]+)\s+eV/A", log_text)
    pressure_matches = re.findall(r"TOTAL-PRESSURE:\s*([-+0-9.eE]+)\s*KBAR", log_text)
    block_force, block_stress = _last_force_stress(log_text, len(atoms))
    result = {
        "driver": "abacus-native-cell-relax",
        "returncode": returncode,
        "converged": bool(returncode == 0 and converged and final_stru.is_file()),
        "workdir": str(workdir),
        "input_structure": str(Path(args.structure).resolve()),
        "final_stru": str(final_stru),
        "optimizer": args.relax_method,
        "calculation": "cell-relax",
        "ecutwfc_Ry": args.ecutwfc,
        "scf_thr": args.scf_thr,
        "kpoints": args.kpts,
        "relax_nmax": args.relax_nmax,
        "force_thr_eV_per_A": args.force_thr,
        "stress_thr_kbar": args.stress_thr,
        "last_largest_gradient_eV_per_A": float(force_matches[-1]) if force_matches else None,
        "last_max_atom_force_eV_per_A": block_force,
        "last_max_stress_kbar": block_stress,
        "last_total_pressure_kbar": float(pressure_matches[-1]) if pressure_matches else None,
        "natoms": len(atoms),
        "composition": {s: int(atoms.get_chemical_symbols().count(s)) for s in sorted(set(atoms.get_chemical_symbols()))},
        "log": str(log) if log else None,
    }
    if final_stru.is_file():
        final = _parse_stru(final_stru)
        write(workdir / "CONTCAR", final, format="vasp", direct=True, vasp5=True)
        result["final_natoms"] = len(final)
        result["final_composition"] = {s: int(final.get_chemical_symbols().count(s)) for s in sorted(set(final.get_chemical_symbols()))}
        result["final_volume_A3"] = float(final.get_volume())
    (workdir / "native_relax_summary.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return result


def main() -> None:
    args = _args()
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    atoms = _read_structure(Path(args.structure).resolve())
    pp = _parse_species_files(args.pp, "--pp")
    basis = _parse_species_files(args.basis, "--basis")
    _write_stru(atoms, workdir / "STRU", pp, basis, Path(args.pseudo_dir), Path(args.basis_dir))
    _write_input(workdir / "INPUT", args, len(set(atoms.get_chemical_symbols())))
    _write_kpt(workdir / "KPT", args.kpts)
    command = args.command
    with (workdir / "native_cell_relax.stdout").open("w", encoding="utf-8") as stdout, (workdir / "native_cell_relax.stderr").open("w", encoding="utf-8") as stderr:
        completed = subprocess.run(command, cwd=workdir, shell=True, stdout=stdout, stderr=stderr, check=False)
    result = _summary(workdir, args, completed.returncode, atoms)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    raise SystemExit(completed.returncode)


if __name__ == "__main__":
    main()
