"""Register two HfO2 inversion variants using an explicit fluorite scaffold.

The parent scaffold is geometric and not a new relaxed cubic calculation.
Mode coefficients below are projections, not energies or Berry polarization.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import io
import json
from pathlib import Path
import warnings

from ase.io import read, write
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from scripts.audit_hfo2_static_replica import sha256
from vcneb.reference_variants import apply_parent_operation


def fluorite_scaffold(tetragonal):
    if tetragonal.get_chemical_symbols() != ["Hf"] * 4 + ["O"] * 8:
        raise ValueError("requires the historical ordered Hf4O8 T cell")
    sites = tetragonal.get_scaled_positions()
    ideal = np.rint(sites * 4) / 4
    delta = sites - ideal
    delta -= np.rint(delta)
    if np.max(np.abs(delta[:4])) > 1e-8 or np.max(np.abs(delta[4:, :2])) > 1e-8:
        raise ValueError("T sites do not have the documented fluorite xy/Hf embedding")
    if not np.allclose(np.sort(ideal[4:, 2]), [.25] * 4 + [.75] * 4, atol=1e-8):
        raise ValueError("unexpected parent oxygen quarter sites")
    parent = tetragonal.copy()
    parent.calc = None
    parent.set_scaled_positions(ideal)
    return parent


def distortion(atoms, parent):
    delta = atoms.get_scaled_positions() - parent.get_scaled_positions()
    delta -= np.rint(delta)
    values = delta @ parent.cell.array
    return values - values.mean(axis=0)


def translation_sectors(parent, values):
    """Four real translational characters of the fluorite conventional cell.

    Wavevectors use its reciprocal basis; these are q sectors, not X irreps.
    Coordinates/metric remain those of the declared strained scaffold.
    """
    translations = ((0, 0, 0), (0, .5, .5), (.5, 0, .5), (.5, .5, 0))
    characters = {"Gamma": (1, 1, 1, 1), "X_x": (1, 1, -1, -1),
                  "X_y": (1, -1, 1, -1), "X_z": (1, -1, -1, 1)}
    permutations = [apply_parent_operation(parent, parent, np.eye(3), t).source_to_target
                    for t in translations]
    transformed = [values[np.argsort(p)] for p in permutations]
    sectors = {label: sum(sign * vector for sign, vector in zip(c, transformed)) / 4
               for label, c in characters.items()}
    if not np.allclose(sum(sectors.values()), values, atol=1e-10, rtol=0):
        raise ValueError("translational sector reconstruction failed")
    gram = np.asarray([[np.sum(a * b) for b in sectors.values()] for a in sectors.values()])
    if np.max(np.abs(gram - np.diag(np.diag(gram)))) > 1e-10:
        raise ValueError("translational sectors are not orthogonal")
    return sectors


def symmetry_report(atoms):
    structure = AseAtomsAdaptor.get_structure(atoms)
    distances = atoms.get_all_distances(mic=True)
    np.fill_diagonal(distances, np.inf)
    if not np.isfinite(atoms.positions).all() or atoms.get_volume() <= 0 or distances.min() < 1.6:
        raise ValueError("invalid HfO2 geometry")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        records = []
        for tolerance in (1e-4, 1e-3, 1e-2):
            analyzer = SpacegroupAnalyzer(structure, symprec=tolerance, angle_tolerance=1)
            records.append({"symprec_A": tolerance, "angle_tolerance_deg": 1,
                            "number": analyzer.get_space_group_number(),
                            "symbol": analyzer.get_space_group_symbol()})
    return {"volume_A3": atoms.get_volume(), "minimum_distance_A": float(distances.min()),
            "periodicity": atoms.pbc.tolist(), "coordinate_mode": "fractional",
            "disordered": not structure.is_ordered, "oxidation_states_assigned": False,
            "symmetry": records, "warnings": sorted(set(str(x.message) for x in caught))}


def write_clean_poscar(path, atoms):
    buffer = io.StringIO()
    write(buffer, atoms, format="vasp", direct=True, sort=False, vasp5=True)
    path.write_text("\n".join(line.rstrip() for line in buffer.getvalue().splitlines()) + "\n",
                    encoding="utf-8", newline="\n")


def prepare(t_path: Path, po_path: Path, m_path: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError("refusing existing variant namespace")
    t, po, monoclinic = [read(p, format="vasp") for p in (t_path, po_path, m_path)]
    if any(a.get_chemical_symbols() != t.get_chemical_symbols() for a in (po, monoclinic)):
        raise ValueError("T/PO/M ordered species differ")
    parent = fluorite_scaffold(t)
    base = distortion(t, parent)
    base /= np.linalg.norm(base)
    po_d = distortion(po, parent)
    po_sectors = translation_sectors(parent, po_d)
    # Inversion 1 is a symmetry of the actual T parent. Inversion 2 is an
    # operation of ideal fluorite that reverses this T distortion instead.
    operations = {"T_preserving_inversion": [.5, 0, .5], "T_reversing_inversion": [0, 0, 0]}
    variants = {}
    for name, translation in operations.items():
        result = apply_parent_operation(parent, po, -np.eye(3), translation)
        t_variant = apply_parent_operation(parent, t, -np.eye(3), translation).atoms
        d = distortion(result.atoms, parent)
        sectors = translation_sectors(parent, d)
        variants[name] = {"atoms": result.atoms, "operation": {
            "rotation_fractional": result.rotation.tolist(), "translation_fractional": result.translation.tolist(),
            "source_to_parent_target": list(result.source_to_target)},
            "T_projection_A": float(np.sum(d * base)),
            "transformed_T_overlap": float(np.sum(distortion(t_variant, parent) * base) / np.linalg.norm(distortion(t, parent))),
            "oxygen_minus_hafnium_mean_displacement_A": (d[4:].mean(0) - d[:4].mean(0)).tolist(),
            "translation_sectors": {k: {"norm_A": float(np.linalg.norm(v)),
                "overlap_with_PO_sector": float(np.sum(v * po_sectors[k]) / np.linalg.norm(v) / np.linalg.norm(po_sectors[k]))}
                for k, v in sectors.items()},
            "structure_audit": symmetry_report(result.atoms)}
    summary = {"schema_version": 1, "purpose": "ordered_parent_operation_variant_registration",
        "sources": {k: {"path": str(p), "sha256": sha256(p)} for k, p in
                    (("T", t_path), ("PO", po_path), ("M_seed", m_path))},
        "software": {p: importlib.metadata.version(p) for p in ("ase", "pymatgen", "spglib")},
        "parent": {"definition": "ideal fluorite quarter-site scaffold in the unmodified T cell; not relaxed cubic HfO2",
                   "cell_A": parent.cell.array.tolist(), "fractional_sites": parent.get_scaled_positions().tolist(),
                   "distortion_basis": "normalized translation-free T minus ideal scaffold, Cartesian in T metric",
                   "basis_cartesian": base.tolist(), "irrep_label": None},
        "PO": {"T_projection_A": float(np.sum(po_d * base)),
               "translation_sectors": {k: {"norm_A": float(np.linalg.norm(v)), "vector_cartesian_A": v.tolist()}
                                       for k, v in po_sectors.items()},
               "oxygen_minus_hafnium_mean_displacement_A": (po_d[4:].mean(0) - po_d[:4].mean(0)).tolist()},
        "structure_audits": {k: symmetry_report(a) for k, a in (("T", t), ("PO", po), ("M_seed", monoclinic))},
        "variants": {k: {x: y for x, y in v.items() if x != "atoms"} for k, v in variants.items()},
        "limitations": ["geometric inversion-related candidates only; electronic Berry branches pending",
                        "T-pattern projection is not an X2- irrep/eigenvector claim",
                        "q sectors Gamma/(1,0,0)/(0,1,0)/(0,0,1) in scaffold reciprocal axes; not individual irreps",
                        "distinct ordered pathways, not proof of topological inequivalence",
                        "M geometry is a literature seed, not an ABACUS endpoint or energy"]}
    output.mkdir(parents=True)
    for name, atoms in [("T", t), ("PO", po), ("M_seed", monoclinic), ("parent_scaffold", parent)]:
        write_clean_poscar(output / f"{name}.vasp", atoms)
    for name, record in variants.items():
        write_clean_poscar(output / f"PO_minus_{name}.vasp", record["atoms"])
    summary["generated_structure_sha256"] = {p.name: sha256(p) for p in sorted(output.glob("*.vasp"))}
    (output / "variant_manifest.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("t", "po", "m", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.t, args.po, args.m, args.output)
    print(json.dumps({"status": "geometric_candidates_only", "variants": result["variants"]}))


if __name__ == "__main__":
    main()
