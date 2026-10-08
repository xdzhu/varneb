"""Stage 16 T-cell Gamma forces and three HfO2 endpoint candidate statics.

All geometry files are new. Electronic INPUT/KPT/pseudo/orbital files remain
byte-identical to the audited historical HfO2 contract. There is no expansion
and no automatic endpoint relaxation or NEB submission in this script.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
from pathlib import Path
import shutil
import warnings

from ase import Atoms
from ase.io import read, write
import numpy as np
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

from scripts.audit_hfo2_static_replica import INPUT_FILES, audited_results, sha256
from scripts.prepare_hfo2_channel_work_probes import write_structure


def gamma_family(atoms, distance_A):
    """Use Angstrom geometries and eV/Angstrom forces (not Bohr interface)."""
    if atoms.get_chemical_symbols() != ["Hf"] * 4 + ["O"] * 8:
        raise ValueError("ordered twelve-atom Hf4O8 required")
    if distance_A not in (.01, .02):
        raise ValueError("only the two declared finite displacement sizes allowed")
    phonon = Phonopy(PhonopyAtoms(symbols=atoms.get_chemical_symbols(), cell=atoms.cell.array,
                                 scaled_positions=atoms.get_scaled_positions(wrap=False)),
                      supercell_matrix=np.eye(3, dtype=int), primitive_matrix=np.eye(3),
                      symprec=1e-4, is_symmetry=True)
    phonon.generate_displacements(distance=distance_A, is_plusminus=True, is_diagonal=False)
    if len(phonon.supercells_with_displacements) != 8:
        raise ValueError("T symmetry/displacement count differs from reviewed eight-point family")
    for shifted, item in zip(phonon.supercells_with_displacements, phonon.dataset["first_atoms"]):
        if shifted is None or len(shifted) != len(atoms) or shifted.symbols != atoms.get_chemical_symbols():
            raise ValueError("Phonopy altered atom count/order")
        delta = shifted.scaled_positions - atoms.get_scaled_positions(wrap=False)
        delta -= np.rint(delta)
        expected = np.zeros((12, 3))
        expected[item["number"]] = item["displacement"]
        if not np.allclose(delta @ atoms.cell.array, expected, atol=1e-10, rtol=0):
            raise ValueError("Phonopy atom/displacement correspondence changed")
    return phonon


def geometry_roundtrip(path, geometry):
    restored = read(path, format="abacus")
    delta = restored.get_scaled_positions() - geometry.get_scaled_positions()
    delta -= np.rint(delta)
    if (restored.get_chemical_symbols() != geometry.get_chemical_symbols()
            or not np.allclose(restored.cell.array, geometry.cell.array, atol=1e-11, rtol=0)
            or np.max(np.abs(delta @ geometry.cell.array)) > 1e-11):
        raise ValueError("actual ASE-ABACUS STRU roundtrip mismatch")
    distances = restored.get_all_distances(mic=True)
    np.fill_diagonal(distances, np.inf)
    if restored.get_volume() <= 0 or distances.min() < 1.6:
        raise ValueError("invalid staged geometry")
    return float(distances.min())


def prepare(source: Path, variants: Path, output: Path):
    if output.exists():
        raise FileExistsError("refusing existing Gamma namespace")
    if sha256(source / "INPUT") != "dc6684ffa709bf3d953c647272589b05d844c1293477bd144226d495207cb9df":
        raise ValueError("historical 100-Ry INPUT differs")
    if sha256(source / "KPT") != "92b917107d9df11da28a465cc4900f3574956a057ddf9f7ad125e75c956ec508":
        raise ValueError("historical Gamma 2x2x2 electronic KPT differs")
    raw = audited_results(source)
    t = read(source / "STRU", format="abacus")
    recorded_t = read(variants / "T.vasp", format="vasp")
    if (not np.allclose(t.cell.array, recorded_t.cell.array, atol=1e-9, rtol=0)
            or np.max(np.abs((t.get_scaled_positions() - recorded_t.get_scaled_positions())
                             - np.rint(t.get_scaled_positions() - recorded_t.get_scaled_positions()))) > 1e-9):
        raise ValueError("T static source differs from the ordered channel endpoint")
    if np.max(np.linalg.norm(raw["forces"], axis=1)) > .01:
        raise ValueError("T source not stationary enough for the declared local Gamma analysis")
    output.mkdir(parents=True)
    records, software_warnings = [], []

    def stage(geometry, metadata):
        index = len(records)
        directory = output / "points" / f"{index:02d}"
        directory.mkdir(parents=True)
        for name in INPUT_FILES:
            if name != "STRU":
                shutil.copyfile(source / name, directory / name)
        write_structure(directory / "STRU", geometry)
        minimum = geometry_roundtrip(directory / "STRU", geometry)
        records.append({"index": index, **metadata, "minimum_distance_A": minimum,
                        "input_sha256": {name: sha256(directory / name) for name in INPUT_FILES}})

    for distance in (.01, .02):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            phonon = gamma_family(t, distance)
            phonon.save(filename=str(output / f"phonopy_T_d{distance:.2f}.yaml"))
        software_warnings.extend(str(x.message) for x in caught)
        for local_index, (shifted, item) in enumerate(zip(phonon.supercells_with_displacements,
                                                         phonon.dataset["first_atoms"])):
            geometry = Atoms(symbols=shifted.symbols, cell=shifted.cell,
                             scaled_positions=shifted.scaled_positions, pbc=True)
            stage(geometry, {"kind": "T_Gamma", "distance_A": distance, "local_index": local_index,
                             "atom_index": int(item["number"]), "displacement_A": item["displacement"]})
    for name in ("PO_minus_T_preserving_inversion", "PO_minus_T_reversing_inversion", "M_seed"):
        path = variants / f"{name}.vasp"
        stage(read(path, format="vasp"), {"kind": "endpoint_candidate_static", "name": name,
                                         "seed_sha256": sha256(path)})
    write(output / "T_reference.vasp", t, format="vasp", direct=True, sort=False, vasp5=True)
    manifest = {"schema_version": 1, "purpose": "HfO2_T_Gamma_two_step_sizes_and_three_seed_statics",
                "source_directory": str(source), "source_inputs": {n: sha256(source / n) for n in INPUT_FILES},
                "source_log_sha256": sha256(source / "OUT.ABACUS/running_scf.log"),
                "T_reference_sha256": sha256(output / "T_reference.vasp"),
                "source_T_results": {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in raw.items()},
                "variants_manifest_sha256": sha256(variants / "variant_manifest.json"),
                "supercell_matrix": np.eye(3, dtype=int).tolist(), "primitive_matrix": np.eye(3).tolist(),
                "phonopy_geometry_unit": "angstrom", "force_unit": "eV/angstrom", "force_constant_unit": "eV/angstrom^2",
                "Gamma_boundary": "analytic fixed-cell atomic q=0; no NAC, dispersion, or variable-cell TS certificate",
                "software": {k: importlib.metadata.version(k) for k in ("ase", "phonopy", "spglib")},
                "warnings": sorted(set(software_warnings)), "n_static_points": len(records), "points": records}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "variants", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    manifest = prepare(args.source, args.variants, args.output)
    print(json.dumps({"status": "staged_no_DFT", "n_points": manifest["n_static_points"],
                      "software": manifest["software"], "warnings": manifest["warnings"]}))


if __name__ == "__main__":
    main()
