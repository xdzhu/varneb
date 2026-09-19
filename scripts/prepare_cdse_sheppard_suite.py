"""Audited 8-atom CdSe seeds and mapped paths; no DFT on a login node.

Seed motifs are idealized explicitly from the official TSASE example, then
relaxed with BFGS. A seed is never represented as a converged endpoint.
"""
import argparse
import hashlib
from io import StringIO
import json
from pathlib import Path
import shutil
import sys
import tarfile

import numpy as np
from ase import Atoms
from ase.io import read, write
import spglib

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from vcneb import endpoint_structure_record, validate_path_geometry

ARCHIVE_SHA = "c59d3a6b3599a4704d0b9093d32d249c3e841ae617949bade7bee14e6a0f166f"
ROUTES = ("cell_mapping", "atomic_mapping")
N_IMAGES = 17  # 15 interior images; endpoints are not workers.


def symmetry(atoms):
    records = []
    for tol in (1e-5, 1e-4, 1e-3):
        ds = spglib.get_symmetry_dataset((atoms.cell.array, atoms.get_scaled_positions(), atoms.numbers), symprec=tol)
        if ds is None:
            raise ValueError("symmetry identification failed")
        records.append({"symprec_A": tol, "number": int(ds.number), "international": ds.international})
    return records


def ideal_seeds(archive):
    if hashlib.sha256(archive.read_bytes()).hexdigest() != ARCHIVE_SHA:
        raise ValueError("official reference archive fingerprint mismatch")
    with tarfile.open(archive) as tar:
        structures = {}
        for label, member in (("rs", "CdSe_sq"), ("wz", "CdSe_hex")):
            text = tar.extractfile(f"tsase-799/tsase/examples/{member}").read().decode()
            structures[label] = read(StringIO(text), format="vasp")
    # Exact motifs remove the empirical relaxation's ~1e-4 A noise only in
    # initial seed construction. Production images are never symmetrized.
    a = float(structures["rs"].cell[2, 2])
    aw = float(structures["wz"].cell[0, 0]) / 2
    cw = float(structures["wz"].cell[2, 2])
    z = float(structures["wz"].get_scaled_positions()[4, 2])
    rs_q = np.array([[0, 0, 0], [.5, 0, 0], [.25, 0, .5], [.75, 0, .5],
                     [.25, 0, 0], [.75, 0, 0], [0, 0, .5], [.5, 0, .5]])
    wz_q = np.array([[0, 0, 0], [.5, 0, 0], [1/6, 1/3, .5], [2/3, 1/3, .5],
                     [1/6, 1/3, z], [2/3, 1/3, z], [0, 0, .5+z], [.5, 0, .5+z]])
    seeds = {"rs": Atoms("Cd4Se4", cell=[[2*a, 0, 0], [a/2, a/2, 0], [0, 0, a]], scaled_positions=rs_q, pbc=True),
             "wz": Atoms("Cd4Se4", cell=[[2*aw, 0, 0], [aw/2, np.sqrt(3)*aw/2, 0], [0, 0, cw]], scaled_positions=wz_q, pbc=True)}
    for label, expected in (("rs", 225), ("wz", 186)):
        if any(r["number"] != expected for r in symmetry(seeds[label])):
            raise ValueError(f"{label} ideal seed has incorrect phase")
    return seeds


def atomic_endpoint(rs):
    """Unimodular basal change followed by a rigid orientation, not strain."""
    transform = np.array([[1, -1, 0], [0, 1, 0], [0, 0, 1]])
    q = rs.get_scaled_positions(wrap=False) @ np.linalg.inv(transform)
    h = transform @ rs.cell.array
    rotation, triangular = np.linalg.qr(h.T)
    signs = np.diag(np.where(np.diag(triangular) < 0, -1., 1.))
    rotation = rotation @ signs
    out = Atoms(rs.symbols, cell=h @ rotation, scaled_positions=q, pbc=True)
    if not np.isclose(out.get_volume(), rs.get_volume(), atol=1e-9) or not np.allclose(
        np.sort(out.get_all_distances(mic=True).ravel()), np.sort(rs.get_all_distances(mic=True).ravel()), atol=1e-9
    ):
        raise ValueError("basis/orientation change altered the physical periodic structure")
    return out


def mapped_chain(rs, wz, route, n_images=N_IMAGES):
    if n_images < 3:
        raise ValueError("at least one interior is required")
    if route not in ROUTES:
        raise ValueError("unknown mapped route")
    if rs.get_chemical_symbols() != ["Cd"]*4+["Se"]*4 or wz.get_chemical_symbols() != rs.get_chemical_symbols():
        raise ValueError("ordered equal 8-atom Cd4Se4 cells required")
    start = atomic_endpoint(rs) if route == "atomic_mapping" else rs.copy()
    q0, q1 = start.get_scaled_positions(wrap=False), wz.get_scaled_positions(wrap=False)
    # The published sq/hex motif correspondence defines zero winding in
    # these unwrapped representatives.  Do NOT reselect a nearest image:
    # atomic_mapping has exact half-cell ties, and tiny origin noise changes
    # np.rint's branch, producing a different, colliding initial path.
    # No atom-index reassignment or runtime pathway repair is performed.
    images = [Atoms(start.symbols, cell=(1-t)*start.cell.array+t*wz.cell.array,
                    scaled_positions=(1-t)*q0+t*q1, pbc=True) for t in np.linspace(0, 1, n_images)]
    validate_path_geometry(images, minimum_distance=1.8, maximum_deformation=.8)
    return images


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2)+"\n", encoding="utf-8")


def prepare_seeds(archive, suite, pseudo_root):
    if (suite / "endpoint_seeds").exists():
        raise FileExistsError("seed directory exists; do not overwrite")
    datasets = [pseudo_root / "PBE" / element / "POTCAR" for element in ("Cd", "Se")]
    for path in datasets:
        if not path.is_file():
            raise FileNotFoundError(path)
    potcar = b"".join(p.read_bytes() for p in datasets)
    for label, atoms in ideal_seeds(archive).items():
        target = suite / "endpoint_seeds" / label
        target.mkdir(parents=True)
        write(target / "POSCAR", atoms, format="vasp", direct=True, vasp5=True)
        (target / "POTCAR").write_bytes(potcar)
        (target / "INCAR").write_text("SYSTEM = CdSe Sheppard comparison PBE\nENCUT = 455\nPREC = Accurate\nEDIFF = 1e-8\nNELM = 180\nALGO = Normal\nISMEAR = 0\nSIGMA = 0.02\nLREAL = .FALSE.\nLASPH = .TRUE.\nISYM = -1\nSYMPREC = 1e-4\nIBRION = -1\nNSW = 0\nISIF = 2\nLWAVE = .FALSE.\nLCHARG = .FALSE.\nNCORE = 4\n", encoding="ascii")
        (target / "KPOINTS").write_text("Fixed 10x10x10 MP mesh\n0\nMonkhorst-Pack\n10 10 10\n0 0 0\n", encoding="ascii")
        write_json(target / "seed.json", {"seed_only_not_relaxed": True, "structure": endpoint_structure_record(atoms),
                   "symmetry": symmetry(atoms), "reference_archive_sha256": ARCHIVE_SHA,
                   "reference_url": "https://theory.cm.utexas.edu/code/tsase.tgz",
                   "initial_seed_idealization": "exact 225/186 motif from empirical TSASE endpoints; lattice lengths and WZ z retained",
                   "pressure_gpa": 0., "potcar_sha256": hashlib.sha256(potcar).hexdigest(),
                   "datasets": [{"label": p.parent.name, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()} for p in datasets]})


def prepare_path(suite, route, endpoint_candidates=None):
    root = suite / route
    if root.exists():
        raise FileExistsError(root)
    endpoints = []
    for label, expected in (("rs", 225), ("wz", 186)):
        directory = (endpoint_candidates / label if endpoint_candidates is not None
                     else suite / "endpoint_relax" / label)
        if endpoint_candidates is None:
            result = json.loads((directory / "endpoint_relax_summary.json").read_text())
            if not result["converged"] or result["external_pressure_gpa"] != 0.:
                raise ValueError("endpoint not converged at zero pressure")
        else:
            result = json.loads((directory / "phase_candidate_audit.json").read_text())
            if not result["candidate_only_not_stationarity_verified"] or result["expected_spacegroup"] != expected:
                raise ValueError("explicit phase candidate record required")
        atoms = read(directory / "CONTCAR")
        atoms.set_constraint()
        if any(r["number"] != expected for r in symmetry(atoms)):
            raise ValueError("relaxed endpoint has wrong phase")
        endpoints.append(atoms)
    images = mapped_chain(*endpoints, route)
    source = endpoint_candidates / "rs" if endpoint_candidates is not None else suite / "endpoint_seeds" / "rs"
    for label, atoms in (("initial", images[0]), ("final", images[-1])):
        target = root / "input" / label
        target.mkdir(parents=True)
        for name in ("INCAR", "KPOINTS", "POTCAR"):
            shutil.copy2(source / name, target / name)
        write(target / "CONTCAR", atoms, format="vasp", direct=True, vasp5=True)
    write(root / "initial.traj", images)
    write_json(root / "preparation.json", {"route": route, "pressure_gpa": 0., "n_atoms": 8, "n_formula_units": 4,
               "n_images_total": N_IMAGES, "n_interiors": N_IMAGES-2, "optimization_constraints": "none",
               "endpoint_stationarity_requires_static_gate": endpoint_candidates is not None,
               "coordinate_policy": "explicit species-preserving index mapping; reference-motif zero winding, no nearest-image reselection",
               "periodic_winding": np.zeros((8,3), dtype=int).tolist(),
               "atomic_route_status": "equivalent-basis reconstructed candidate; Fig5 mechanism identity requires result audit",
               "comparison": {"paper_xc": "PW91", "current_xc": "PBE", "paper_encut_eV": 455,
                              "paper_kmesh": [10,10,10], "paper_dft_barrier_meV_per_atom": 2.4,
                              "caution": "not empirical-potential barriers; 0.10 force threshold does not certify meV-scale barrier accuracy"},
               "geometry": validate_path_geometry(images, minimum_distance=1.8, maximum_deformation=.8),
               "endpoint_symmetry": {"initial": symmetry(images[0]), "final": symmetry(images[-1])},
               "trajectory_sha256": hashlib.sha256((root / "initial.traj").read_bytes()).hexdigest()})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("seeds", "path"), required=True)
    parser.add_argument("--suite-root", type=Path, required=True)
    parser.add_argument("--archive", type=Path)
    parser.add_argument("--pseudo-root", type=Path)
    parser.add_argument("--route", choices=ROUTES)
    parser.add_argument("--endpoint-candidates", type=Path,
                        help="explicit phase candidates; production requires a separate raw static force/virial gate")
    args = parser.parse_args()
    if args.stage == "seeds":
        if args.archive is None or args.pseudo_root is None:
            parser.error("seeds requires --archive and --pseudo-root")
        prepare_seeds(args.archive, args.suite_root, args.pseudo_root)
    else:
        if args.route is None:
            parser.error("path requires --route")
        prepare_path(args.suite_root, args.route, args.endpoint_candidates)
    print(f"Prepared {args.stage}; no DFT was executed")


if __name__ == "__main__":
    main()
