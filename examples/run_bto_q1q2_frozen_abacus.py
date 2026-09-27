"""Validate or evaluate frozen-C-cell BTO Q1/Q2 points with proven ABACUS settings.

The default is input generation only. ``--run-dft`` requires an hfacnormal01
Slurm allocation and a fresh workdir; it never modifies the archived BTO NEB
calculation. These points are a diagnostic frozen atomic-mode slice, not a
variable-cell T-to-C barrier or a conditionally relaxed PES.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys

import numpy as np
from ase.io import read

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb.abacus import make_ase_abacus_factory


PP = {"Ba": "Ba.upf", "Ti": "Ti.upf", "O": "O.upf"}
BASIS = {
    "Ba": "Ba_gga_10au_100Ry_4s2p1d.orb",
    "Ti": "Ti_gga_10au_100Ry_4s2p2d1f.orb",
    "O": "O_gga_10au_100Ry_2s2p1d.orb",
}
PARAMETERS = {
    "calculation": "scf", "basis_type": "lcao", "dft_functional": "pbe",
    "cal_force": 1, "cal_stress": 1, "out_stru": 1,
    "ecutwfc": 100.0, "scf_thr": 1e-8, "scf_nmax": 150,
    "mixing_beta": 0.3, "kpts": (4, 4, 4),
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_input(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("#", 1)[0].split()
        if len(fields) >= 2:
            key = fields[0].lower()
            if key in values:
                raise ValueError(f"duplicate INPUT key {key} in {path}")
            values[key] = " ".join(fields[1:])
    expected = {
        "calculation": "scf", "basis_type": "lcao", "dft_functional": "pbe",
        "cal_force": "1", "cal_stress": "1", "out_stru": "1",
        "scf_nmax": "150",
    }
    for key, value in expected.items():
        if values.get(key, "").lower() != value:
            raise ValueError(f"INPUT {key} differs from audited BTO protocol: {path}")
    for key, value in (("ecutwfc", 100.0), ("scf_thr", 1e-8), ("mixing_beta", 0.3)):
        if key not in values or not np.isclose(float(values[key]), value, rtol=1e-12, atol=0.0):
            raise ValueError(f"INPUT {key} differs from audited BTO protocol: {path}")
    return values


def _validate_written_case(directory: Path) -> dict[str, str]:
    for name in ("INPUT", "KPT", "STRU"):
        if not (directory / name).is_file() or (directory / name).stat().st_size == 0:
            raise ValueError(f"missing nonempty {name} in {directory}")
    _read_input(directory / "INPUT")
    kpt = (directory / "KPT").read_text(encoding="utf-8")
    if not re.search(r"\bGamma\s+4\s+4\s+4\s+0\s+0\s+0\b", kpt, re.IGNORECASE):
        raise ValueError(f"KPT differs from audited 4x4x4 Gamma mesh: {directory}")
    stru = (directory / "STRU").read_text(encoding="utf-8")
    for filename in [*PP.values(), *BASIS.values()]:
        if filename not in stru:
            raise ValueError(f"STRU is missing the reviewed pseudo/orbital {filename}: {directory}")
    return {name: _sha256(directory / name) for name in ("INPUT", "KPT", "STRU")}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points-dir", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True, help="fresh directory only")
    parser.add_argument("--pseudo-dir", type=Path, required=True)
    parser.add_argument("--basis-dir", type=Path, required=True)
    parser.add_argument("--abacus-bin", type=Path, required=True)
    parser.add_argument("--mpi-ranks", type=int, default=32)
    parser.add_argument("--run-dft", action="store_true", help="otherwise generate and audit inputs only")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.workdir.exists():
        raise FileExistsError(f"refusing to overwrite existing workdir: {args.workdir}")
    if not args.abacus_bin.is_file() or args.mpi_ranks < 1:
        raise ValueError("ABACUS binary must exist and mpi-ranks must be positive")
    if args.run_dft:
        if not os.environ.get("SLURM_JOB_ID") or os.environ.get("SLURM_JOB_PARTITION") != "hfacnormal01":
            raise RuntimeError("--run-dft requires an hfacnormal01 Slurm allocation")
        if int(os.environ.get("SLURM_JOB_NUM_NODES", "0")) != 1 or int(os.environ.get("SLURM_NTASKS", "0")) != 1:
            raise RuntimeError("run requires one Slurm task owning the MPI worker CPUs")
        if int(os.environ.get("SLURM_CPUS_PER_TASK", "0")) < args.mpi_ranks:
            raise RuntimeError("Slurm cpus-per-task is smaller than mpi-ranks")
    manifest_path = args.points_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    supported = {
        "inert_frozen_cubic_cell_gamma_mode_points_not_T_to_C_barrier",
        "inert_frozen_cubic_cell_transverse_soft_increment20_not_T_to_C_barrier",
        "inert_frozen_cubic_cell_transverse_soft_centers16_not_T_to_C_barrier",
        "inert_frozen_cubic_cell_transverse_soft_edges18_not_T_to_C_barrier",
    }
    if manifest.get("kind") not in supported:
        raise ValueError("unrecognized or non-frozen BTO point manifest")
    if manifest.get("n_points") != len(manifest.get("points", [])) or not manifest["points"]:
        raise ValueError("point manifest has inconsistent point count")
    if manifest["kind"] == "inert_frozen_cubic_cell_transverse_soft_increment20_not_T_to_C_barrier":
        if (manifest["n_points"] != 20
                or manifest.get("no_dft_launched") is not True
                or manifest.get("grid") is not None
                or not manifest.get("source_full_manifest_sha256")
                or not manifest.get("shared_axis_reuse_audit_sha256")
                or manifest.get("axis_definition")
                   != "Ti_z-minus-Ba_z_and_Ti_x-minus-Ba_x_projected_to_Gamma_unstable_triplet"
                or len(manifest.get("axis_mode_weights", [])) != 15
                or any(point.get("q2") == 0.0 for point in manifest["points"])):
            raise ValueError("transverse-soft increment lacks its 20-point off-axis provenance")
    if manifest["kind"] == "inert_frozen_cubic_cell_transverse_soft_centers16_not_T_to_C_barrier":
        if (manifest["n_points"] != 16
                or manifest.get("no_dft_launched") is not True
                or manifest.get("grid") is not None
                or not all(manifest.get(key) for key in (
                    "source_coarse_manifest_sha256", "source_dense_manifest_sha256",
                    "source_coarse_assembled_sha256"))
                or manifest.get("axis_definition")
                   != "Ti_z-minus-Ba_z_and_Ti_x-minus-Ba_x_projected_to_Gamma_unstable_triplet"
                or len(manifest.get("axis_mode_weights", [])) != 15
                or any(len(point.get("grid_index_q1_q2", [])) != 2
                       or not all(type(index) is int and index % 2 == 1
                                  for index in point["grid_index_q1_q2"])
                       for point in manifest["points"])):
            raise ValueError("transverse-soft center increment lacks its 16-point provenance")
    if manifest["kind"] == "inert_frozen_cubic_cell_transverse_soft_edges18_not_T_to_C_barrier":
        if (manifest["n_points"] != 18
                or manifest.get("no_dft_launched") is not True
                or manifest.get("grid") is not None
                or not all(manifest.get(key) for key in (
                    "source_dense_manifest_sha256", "source_41_analysis_sha256"))
                or manifest.get("axis_definition")
                   != "Ti_z-minus-Ba_z_and_Ti_x-minus-Ba_x_projected_to_Gamma_unstable_triplet"
                or len(manifest.get("axis_mode_weights", [])) != 15
                or any(len(point.get("grid_index_q1_q2", [])) != 2
                       or not all(type(index) is int for index in point["grid_index_q1_q2"])
                       or (point["grid_index_q1_q2"][0] + point["grid_index_q1_q2"][1]) % 2 != 1
                       for point in manifest["points"])):
            raise ValueError("transverse-soft edge increment lacks its 18-point provenance")
    asset_paths = {**{f"pp_{symbol}": args.pseudo_dir / name for symbol, name in PP.items()},
                   **{f"basis_{symbol}": args.basis_dir / name for symbol, name in BASIS.items()}}
    for label, path in asset_paths.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"missing reviewed ABACUS asset {label}: {path}")
    points = []
    for point in manifest["points"]:
        name = point["name"]
        if not (re.fullmatch(r"image-\d\d-projected", name)
                or re.fullmatch(r"grid-q1-\d\d-q2-\d\d", name)):
            raise ValueError(f"unsafe or unexpected point name: {name}")
        structure = args.points_dir / point["structure"]
        if structure.resolve().parent != (args.points_dir / name).resolve():
            raise ValueError(f"point structure escapes its declared directory: {structure}")
        if _sha256(structure) != point["structure_sha256"]:
            raise ValueError(f"point structure hash changed: {structure}")
        atoms = read(str(structure))
        if atoms.get_chemical_formula() != "BaO3Ti" or len(atoms) != 5:
            raise ValueError(f"point is not one 5-atom BaTiO3 formula unit: {structure}")
        distances = atoms.get_all_distances(mic=True)
        np.fill_diagonal(distances, np.inf)
        if float(np.min(distances)) < 1.6:
            raise ValueError(f"unsafe atom separation in {structure}")
        points.append((point, structure, atoms))
    parameters = dict(PARAMETERS)
    parameters.update({"pseudo_dir": str(args.pseudo_dir), "basis_dir": str(args.basis_dir),
                       "pp": PP, "basis": BASIS})
    command = f"mpirun -np {args.mpi_ranks} {args.abacus_bin}"
    factory = make_ase_abacus_factory(parameters=parameters, command=command)
    args.workdir.mkdir(parents=True)
    result_manifest = {
        "kind": "bto_frozen_cubic_cell_abacus_static_canary_not_T_to_C_barrier",
        "status": "input_preflight_in_progress", "source_manifest_sha256": _sha256(manifest_path),
        "source_kind": manifest["kind"], "dft_requested": args.run_dft,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "abacus_binary": str(args.abacus_bin), "abacus_binary_sha256": _sha256(args.abacus_bin),
        "mpi_ranks": args.mpi_ranks, "calculator_parameters": parameters,
        "asset_sha256": {label: _sha256(path) for label, path in asset_paths.items()},
        "points": [],
    }
    for index, (point, structure, atoms) in enumerate(points):
        directory = args.workdir / point["name"]
        directory.mkdir()
        shutil.copy2(structure, directory / "source_POSCAR")
        calculator = factory(index, atoms, directory)
        calculator.write_inputfiles(atoms, properties=["energy", "forces", "stress"])
        input_hashes = _validate_written_case(directory)
        record = {"name": point["name"], "q1": point["q1"], "q2": point["q2"],
                  "source_structure_sha256": point["structure_sha256"],
                  "input_sha256": input_hashes, "status": "input_validated"}
        result_manifest["points"].append(record)
        if args.run_dft:
            atoms.calc = calculator
            energy = float(atoms.get_potential_energy())
            forces = np.asarray(atoms.get_forces(), dtype=float)
            stress = np.asarray(atoms.get_stress(voigt=False), dtype=float)
            log_path = directory / "OUT.ABACUS" / "running_scf.log"
            log = log_path.read_text(encoding="utf-8", errors="replace")
            if "charge density convergence is achieved" not in log or "!FINAL_ETOT_IS" not in log:
                raise ValueError(f"ABACUS static calculation lacks a converged final marker: {log_path}")
            mpi_sizes = re.findall(r"\bDSIZE\s*=\s*(\d+)", log)
            if len(mpi_sizes) != 1 or int(mpi_sizes[0]) != args.mpi_ranks:
                raise ValueError(f"ABACUS did not use the requested {args.mpi_ranks} MPI ranks: {log_path}")
            if not np.isfinite(energy) or not np.all(np.isfinite(forces)) or not np.all(np.isfinite(stress)):
                raise ValueError(f"ABACUS returned non-finite results: {directory}")
            record.update({"status": "converged", "energy_eV": energy,
                           "max_atomic_force_eV_per_A": float(np.max(np.linalg.norm(forces, axis=1))),
                           "stress_eV_per_A3": stress.tolist(), "mpi_dsize": int(mpi_sizes[0]),
                           "log_sha256": _sha256(log_path)})
        (args.workdir / "result_manifest.json").write_text(
            json.dumps(result_manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
    result_manifest["status"] = "converged" if args.run_dft else "inputs_validated_no_dft"
    if args.run_dft:
        reference = next((item for item in result_manifest["points"]
                          if abs(item["q1"]) < 1e-12 and abs(item["q2"]) < 1e-12), None)
        if reference is not None:
            for item in result_manifest["points"]:
                item["relative_energy_to_C_eV_per_formula_unit"] = item["energy_eV"] - reference["energy_eV"]
    (args.workdir / "result_manifest.json").write_text(
        json.dumps(result_manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps({"workdir": str(args.workdir), "status": result_manifest["status"],
                      "n_points": len(points)}, indent=2))


if __name__ == "__main__":
    main()
