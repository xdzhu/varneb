"""Build a 12-atom HfO2 tetragonal-to-polar-orthorhombic VC-NEB fixture."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import sys

import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.units import Bohr
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vcneb import interpolate_vcneb


ABACUS_ROOT = Path(r"D:\Work\Code\abacus-agent")
T_VASP = ABACUS_ROOT / (
    r"朱旭东-曙光算例36个\BORN电荷_3+6=9"
    r"\T相HfO2的BORN有效电荷计算-VASP\2-scf\POSCAR"
)
O_ABACUS = ABACUS_ROOT / (
    r"朱旭东-曙光算例36个\BORN电荷_3+6=9"
    r"\O相HfO2的BORN有效电荷计算-ABACUS\1-cell-relax\STRU"
)
OUTDIR = ROOT / "validation" / "hfo2_t_to_po"
SOURCE_DIR = OUTDIR / "source_structures"
T_SOURCE_COPY = SOURCE_DIR / "T_HfO2_primitive.vasp"
O_SOURCE_COPY = SOURCE_DIR / "PO_HfO2_12.STRU"


def read_abacus_stru(path: Path) -> Atoms:
    lines = path.read_text(encoding="utf-8").splitlines()
    lattice_constant = 1.0
    lattice_vectors = []
    species_order = []
    positions = []
    counts = []

    index = 0
    while index < len(lines):
        line = lines[index].strip()
        upper = line.upper()
        if upper == "LATTICE_CONSTANT":
            lattice_constant = float(lines[index + 1].split()[0])
            index += 2
            continue
        if upper == "LATTICE_VECTORS":
            lattice_vectors = [
                [float(x) for x in lines[index + 1 + row].split()[:3]]
                for row in range(3)
            ]
            index += 4
            continue
        if upper == "ATOMIC_POSITIONS":
            coord_type = lines[index + 1].strip().lower()
            if not coord_type.startswith("direct"):
                raise ValueError(f"Only Direct STRU positions are supported: {path}")
            index += 2
            while index < len(lines):
                symbol = lines[index].strip().split()[0] if lines[index].strip() else ""
                if not symbol:
                    index += 1
                    continue
                if symbol.upper() in {"NUMERICAL_ORBITAL", "LATTICE_CONSTANT", "LATTICE_VECTORS"}:
                    break
                species_order.append(symbol)
                _mag = lines[index + 1]
                count = int(lines[index + 2].split()[0])
                counts.append(count)
                index += 3
                for _ in range(count):
                    fields = lines[index].split()
                    positions.append((symbol, [float(x) for x in fields[:3]]))
                    index += 1
            continue
        index += 1

    if not lattice_vectors or not positions:
        raise ValueError(f"Failed to parse STRU: {path}")
    cell = np.asarray(lattice_vectors, dtype=float) * lattice_constant * Bohr
    symbols = [symbol for symbol, _ in positions]
    scaled = np.asarray([pos for _, pos in positions], dtype=float)
    return Atoms(symbols=symbols, scaled_positions=scaled, cell=cell, pbc=True)


def sort_by_species(atoms: Atoms) -> Atoms:
    scaled = atoms.get_scaled_positions(wrap=True)
    order = sorted(
        range(len(atoms)),
        key=lambda i: (
            str(atoms.symbols[i]),
            round(float(scaled[i, 0]), 10),
            round(float(scaled[i, 1]), 10),
            round(float(scaled[i, 2]), 10),
        ),
    )
    return atoms[order]


def tetragonal_12_atom_cell(t_primitive: Atoms) -> Atoms:
    transform = np.array([[1, 1, 0], [-1, 1, 0], [0, 0, 1]], dtype=int)
    inv_transform = np.linalg.inv(transform)
    symbols = []
    scaled = []
    seen = set()
    primitive_scaled = t_primitive.get_scaled_positions(wrap=True)
    for atom_index, (symbol, q_old) in enumerate(zip(t_primitive.get_chemical_symbols(), primitive_scaled)):
        for tx in range(-1, 3):
            for ty in range(-1, 3):
                for tz in range(0, 1):
                    q_new = (q_old + np.array([tx, ty, tz], dtype=float)) @ inv_transform
                    q_wrapped = q_new - np.floor(q_new)
                    if np.all(q_new >= -1e-10) and np.all(q_new < 1.0 - 1e-10):
                        key = (atom_index, tuple(np.round(q_wrapped, 10)))
                        if key in seen:
                            continue
                        seen.add(key)
                        symbols.append(symbol)
                        scaled.append(q_wrapped)
    if len(symbols) != 2 * len(t_primitive):
        raise RuntimeError(f"Expected {2 * len(t_primitive)} atoms in tetragonal supercell, got {len(symbols)}")
    atoms = Atoms(symbols=symbols, scaled_positions=np.asarray(scaled), cell=transform @ t_primitive.cell.array, pbc=True)
    a = t_primitive.cell.lengths()[0] * np.sqrt(2.0)
    c = t_primitive.cell.lengths()[2]
    atoms.set_cell(np.diag([a, a, c]), scale_atoms=True)
    return sort_by_species(atoms)


def mic_delta(frac_a: np.ndarray, frac_b: np.ndarray) -> np.ndarray:
    delta = frac_a - frac_b
    return delta - np.rint(delta)


def assignment_cost(initial: Atoms, final: Atoms, translation: np.ndarray) -> tuple[float, list[int], np.ndarray]:
    final_scaled = (final.get_scaled_positions(wrap=False) + translation) % 1.0
    average_cell = 0.5 * (initial.cell.array + final.cell.array)
    mapping = [-1] * len(initial)
    total_cost = 0.0
    for symbol in sorted(set(initial.get_chemical_symbols())):
        init_idx = [i for i, s in enumerate(initial.symbols) if s == symbol]
        final_idx = [i for i, s in enumerate(final.symbols) if s == symbol]
        cost = np.zeros((len(init_idx), len(final_idx)))
        for row, i in enumerate(init_idx):
            for col, j in enumerate(final_idx):
                dfrac = mic_delta(initial.get_scaled_positions(wrap=False)[i], final_scaled[j])
                dcart = dfrac @ average_cell
                cost[row, col] = np.linalg.norm(dcart)
        rows, cols = linear_sum_assignment(cost)
        total_cost += float(cost[rows, cols].sum())
        for row, col in zip(rows, cols):
            mapping[init_idx[row]] = final_idx[col]
    return total_cost, mapping, final_scaled


def reorder_final_to_initial(initial: Atoms, final: Atoms) -> tuple[Atoms, dict]:
    candidates = []
    init_hf = [i for i, s in enumerate(initial.symbols) if s == "Hf"]
    final_hf = [i for i, s in enumerate(final.symbols) if s == "Hf"]
    qi = initial.get_scaled_positions(wrap=False)
    qf = final.get_scaled_positions(wrap=False)
    for i in init_hf:
        for j in final_hf:
            candidates.append(qi[i] - qf[j])
    candidates.append(np.zeros(3))

    best = None
    for translation in candidates:
        cost, mapping, final_scaled = assignment_cost(initial, final, translation)
        if best is None or cost < best[0]:
            best = (cost, translation, mapping, final_scaled)
    assert best is not None
    cost, translation, mapping, final_scaled = best
    reordered = final[mapping]
    reordered.set_scaled_positions(final_scaled[mapping])
    info = {
        "translation": [float(x) for x in translation],
        "mapping_initial_index_to_original_final_index": [int(x) for x in mapping],
        "total_assignment_distance_A": float(cost),
        "mean_assignment_distance_A": float(cost / len(initial)),
    }
    return reordered, info


def path_distances(images: list[Atoms]) -> list[float]:
    values = []
    for left, right in zip(images[:-1], images[1:]):
        average_cell = 0.5 * (left.cell.array + right.cell.array)
        d = mic_delta(right.get_scaled_positions(wrap=False), left.get_scaled_positions(wrap=False))
        values.append(float(np.sqrt((d @ average_cell * (d @ average_cell)).sum())))
    return values


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DIR.mkdir(exist_ok=True)
    t_source = T_VASP if T_VASP.exists() else T_SOURCE_COPY
    o_source = O_ABACUS if O_ABACUS.exists() else O_SOURCE_COPY
    if not t_source.exists():
        raise FileNotFoundError(f"No tetragonal source found: {T_VASP} or {T_SOURCE_COPY}")
    if not o_source.exists():
        raise FileNotFoundError(f"No orthorhombic source found: {O_ABACUS} or {O_SOURCE_COPY}")
    if t_source != T_SOURCE_COPY:
        shutil.copy2(t_source, T_SOURCE_COPY)
    if o_source != O_SOURCE_COPY:
        shutil.copy2(o_source, O_SOURCE_COPY)

    t_primitive = read(t_source)
    t_12 = tetragonal_12_atom_cell(t_primitive)
    o_12 = sort_by_species(read_abacus_stru(o_source))
    o_mapped, mapping = reorder_final_to_initial(t_12, o_12)

    images = interpolate_vcneb(t_12, o_mapped, n_images=7, align_cells=True, mic=True)
    write(OUTDIR / "T_HfO2_12.vasp", t_12, format="vasp", direct=True, vasp5=True)
    write(OUTDIR / "PO_HfO2_12_mapped.vasp", o_mapped, format="vasp", direct=True, vasp5=True)
    write(OUTDIR / "initial_vcneb_path.traj", images)
    for index, image in enumerate(images):
        image_dir = OUTDIR / f"image_{index:02d}"
        image_dir.mkdir(exist_ok=True)
        write(image_dir / "POSCAR", image, format="vasp", direct=True, vasp5=True)

    report = {
        "sources": {
            "tetragonal_vasp_poscar": str(t_source),
            "orthorhombic_abacus_stru": str(o_source),
            "portable_tetragonal_copy": str(T_SOURCE_COPY),
            "portable_orthorhombic_copy": str(O_SOURCE_COPY),
        },
        "tetragonal_12_cellpar": [float(x) for x in t_12.cell.cellpar()],
        "orthorhombic_12_cellpar": [float(x) for x in o_mapped.cell.cellpar()],
        "symbols": t_12.get_chemical_symbols(),
        "mapping": mapping,
        "path_segment_distances_A": path_distances(images),
    }
    (OUTDIR / "mapping_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    (OUTDIR / "README.md").write_text(
        "# HfO2 T to PO VC-NEB fixture\n\n"
        "This fixture uses a 12-atom sqrt(2) x sqrt(2) x 1 tetragonal HfO2 cell "
        "and a mapped 12-atom polar/orthorhombic HfO2 endpoint from the local "
        "ABACUS Born-charge examples.\n\n"
        "Generated files:\n"
        "- `T_HfO2_12.vasp`\n"
        "- `PO_HfO2_12_mapped.vasp`\n"
        "- `initial_vcneb_path.traj`\n"
        "- `image_00` ... `image_06` POSCAR files\n"
        "- `mapping_report.json`\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
