"""Prepare the first PRL 2023 HfO2 Fig. 2(a) VASP/VARNEB cases.

The public transition files contain twenty concatenated POSCAR records whose
comment lines are ``Image 1`` ... ``Image 20``.  ``seeds`` preserves those
records and writes path-specific endpoint inputs.  ``correct`` transfers the
non-linear displacement of the published chain onto endpoints relaxed with
the selected local PAW-PBE potential.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from ase.io import read, write


CASES = {
    "t_to_po": "T-pca21",
    "po_to_m": "pca21-M",
}


INCAR = """SYSTEM = HfO2 PRL2023 Fig2a VARNEB
PREC = Accurate
ENCUT = 600
EDIFF = 1E-6
ALGO = Normal
NELM = 160
ISPIN = 1
ISMEAR = 0
SIGMA = 0.05
GGA = PE
LASPH = .TRUE.
LREAL = .FALSE.
ADDGRID = .TRUE.
LMAXMIX = 4
IBRION = -1
NSW = 0
ISIF = 2
ISYM = -1
SYMPREC = 1E-4
LWAVE = .FALSE.
LCHARG = .FALSE.
KPAR = 4
NCORE = 8
"""


KPOINTS = """HfO2 PRL2023 4x4x4 MP
0
Monkhorst-Pack
4 4 4
0 0 0
"""


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_concatenated_poscars(path: Path) -> list:
    lines = path.read_text(encoding="utf-8").splitlines()
    starts = [index for index, line in enumerate(lines) if line.strip().startswith("Image ")]
    if not starts:
        raise ValueError(f"no Image records in {path}")
    starts.append(len(lines))
    images = []
    for begin, end in zip(starts[:-1], starts[1:]):
        block = "\n".join(lines[begin:end]) + "\n"
        images.append(read(io.StringIO(block), format="vasp"))
    if len(images) != 20:
        raise ValueError(f"expected 20 images in {path}, found {len(images)}")
    symbols = images[0].get_chemical_symbols()
    if any(image.get_chemical_symbols() != symbols for image in images):
        raise ValueError(f"species/order changes within {path}")
    return images


def minimum_distance(atoms) -> float:
    distances = atoms.get_all_distances(mic=True)
    distances[distances == 0.0] = np.inf
    return float(np.min(distances))


def write_seed_cases(source_root: Path, output_root: Path) -> None:
    output_root.mkdir(parents=True, exist_ok=True)
    manifests = {}
    for case_name, source_name in CASES.items():
        source = source_root / source_name
        images = read_concatenated_poscars(source)
        case_root = output_root / case_name
        initial = case_root / "input" / "initial"
        final = case_root / "input" / "final"
        initial.mkdir(parents=True, exist_ok=True)
        final.mkdir(parents=True, exist_ok=True)
        for directory, image in ((initial, images[0]), (final, images[-1])):
            write(directory / "POSCAR", image, format="vasp", direct=True, vasp5=True)
            (directory / "INCAR").write_text(INCAR, encoding="utf-8")
            (directory / "KPOINTS").write_text(KPOINTS, encoding="utf-8")
        write(case_root / "author_path.traj", images)
        manifest = {
            "case": case_name,
            "source": str(source),
            "source_sha256": sha256(source),
            "n_images_total": len(images),
            "n_images_interior": len(images) - 2,
            "formula": images[0].get_chemical_formula(),
            "n_formula_units": len(images[0]) // 3,
            "minimum_image_distance_A": min(minimum_distance(image) for image in images),
            "initial_cell_A": np.asarray(images[0].cell).tolist(),
            "final_cell_A": np.asarray(images[-1].cell).tolist(),
        }
        (case_root / "seed_manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        manifests[case_name] = manifest
    (output_root / "seed_manifest.json").write_text(
        json.dumps(manifests, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def unwrap_scaled(images: list) -> np.ndarray:
    scaled = [np.asarray(image.get_scaled_positions(wrap=False), dtype=float) for image in images]
    unwrapped = [scaled[0]]
    for current in scaled[1:]:
        delta = current - unwrapped[-1]
        delta -= np.rint(delta)
        unwrapped.append(unwrapped[-1] + delta)
    return np.asarray(unwrapped)


def endpoint_near(reference: np.ndarray, endpoint) -> np.ndarray:
    scaled = np.asarray(endpoint.get_scaled_positions(wrap=False), dtype=float)
    return scaled + np.rint(reference - scaled)


def corrected_chain(case_root: Path, initial_path: Path, final_path: Path, output: Path) -> None:
    author = read(case_root / "author_path.traj", index=":")
    initial = read(initial_path)
    final = read(final_path)
    if len(author) != 20:
        raise ValueError("author trajectory must contain 20 images")
    symbols = author[0].get_chemical_symbols()
    if initial.get_chemical_symbols() != symbols or final.get_chemical_symbols() != symbols:
        raise ValueError("relaxed endpoint species/order differs from the published chain")

    published_scaled = unwrap_scaled(author)
    relaxed_initial = endpoint_near(published_scaled[0], initial)
    relaxed_final = endpoint_near(published_scaled[-1], final)
    cell0 = np.asarray(initial.cell, dtype=float)
    cell1 = np.asarray(final.cell, dtype=float)
    old_cell0 = np.asarray(author[0].cell, dtype=float)
    old_cell1 = np.asarray(author[-1].cell, dtype=float)
    corrected = []
    for index, old in enumerate(author):
        t = index / (len(author) - 1)
        old_cell_linear = (1.0 - t) * old_cell0 + t * old_cell1
        cell = (1.0 - t) * cell0 + t * cell1 + np.asarray(old.cell) - old_cell_linear
        old_scaled_linear = (1.0 - t) * published_scaled[0] + t * published_scaled[-1]
        scaled = (
            (1.0 - t) * relaxed_initial
            + t * relaxed_final
            + published_scaled[index]
            - old_scaled_linear
        )
        image = old.copy()
        image.set_cell(cell, scale_atoms=False)
        image.set_scaled_positions(scaled)
        if np.linalg.det(image.cell.array) <= 0.0:
            raise ValueError(f"non-positive cell at image {index}")
        corrected.append(image)
    corrected[0] = initial.copy()
    corrected[-1] = final.copy()
    min_distance = min(minimum_distance(image) for image in corrected)
    if min_distance < 1.40:
        raise ValueError(f"corrected chain has an unsafe minimum distance: {min_distance:.6f} A")
    output.parent.mkdir(parents=True, exist_ok=True)
    write(output, corrected)
    summary = {
        "status": "completed",
        "method": "published nonlinear deviation transferred to locally relaxed endpoints",
        "n_images_total": len(corrected),
        "minimum_distance_A": min_distance,
        "initial": str(initial_path.resolve()),
        "final": str(final_path.resolve()),
        "output": str(output.resolve()),
    }
    output.with_suffix(".json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="stage", required=True)
    seeds = subparsers.add_parser("seeds")
    seeds.add_argument("--source-root", required=True, type=Path)
    seeds.add_argument("--output-root", required=True, type=Path)
    correct = subparsers.add_parser("correct")
    correct.add_argument("--case-root", required=True, type=Path)
    correct.add_argument("--initial", required=True, type=Path)
    correct.add_argument("--final", required=True, type=Path)
    correct.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.stage == "seeds":
        write_seed_cases(args.source_root.resolve(), args.output_root.resolve())
    else:
        corrected_chain(
            args.case_root.resolve(), args.initial.resolve(), args.final.resolve(), args.output.resolve()
        )


if __name__ == "__main__":
    main()
