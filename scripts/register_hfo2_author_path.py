"""Align a published PO->M seed once to the actual ordered PO+ endpoint.

The fixed z reflection selects its opposite polar variant; a single parent
endpoint assignment is then applied to every image. No per-image relabelling,
SCF, endpoint relaxation, or claim of a new physical mechanism is made here.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write
from scipy.optimize import linear_sum_assignment

from scripts.audit_hfo2_static_replica import sha256
from scripts.prepare_hfo2_reference_variants import symmetry_report, write_clean_poscar


def register(source: Path, target: Path, output: Path):
    if output.exists():
        raise FileExistsError("refusing existing path registration namespace")
    images = read(source, index=":")
    actual = read(target, format="vasp")
    if len(images) != 20 or any(a.get_chemical_symbols() != actual.get_chemical_symbols() for a in images):
        raise ValueError("reviewed twenty-image ordered Hf4O8 author seed required")
    rotation = np.diag([1., 1., -1.])
    transformed = images[0].get_scaled_positions(wrap=False) @ rotation.T
    target_sites = actual.get_scaled_positions(wrap=False)
    shift = target_sites[0] - transformed[0]
    shift -= np.rint(shift)
    differences = transformed[:, None, :] + shift - target_sites[None, :, :]
    differences -= np.rint(differences)
    distances = np.linalg.norm(differences @ actual.cell.array, axis=2)
    distances[images[0].numbers[:, None] != actual.numbers[None, :]] = 1e6
    source_indices, target_indices = linear_sum_assignment(distances)
    if not np.array_equal(source_indices, np.arange(12)) or max(distances[source_indices, target_indices]) > .03:
        raise ValueError("published endpoint is not close to actual PO+ under the declared z reflection")
    # Reject ambiguous nearest targets instead of inventing a preferred oxygen
    # path identity from a nearly tied assignment.
    margins = np.sort(distances, axis=1)[:, 1] - np.sort(distances, axis=1)[:, 0]
    if np.min(margins) < .1:
        raise ValueError("endpoint site assignment is ambiguous")
    registered = []
    for image in images:
        geometry = image.copy()
        geometry.calc = None
        # Fractional and Cartesian reflections together retain a right-handed
        # lattice, including off-diagonal components; no cell symmetrization.
        geometry.set_cell(rotation @ image.cell.array @ rotation.T, scale_atoms=False)
        sites = image.get_scaled_positions(wrap=False) @ rotation.T + shift
        ordered = np.empty_like(sites)
        ordered[target_indices] = sites
        geometry.set_scaled_positions(ordered % 1.0)
        symmetry_report(geometry)
        registered.append(geometry)
    output.mkdir(parents=True)
    shutil.copyfile(source, output / "published_seed.traj")
    write(output / "registered_author_path.traj", registered)
    write_clean_poscar(output / "M_seed_registered.vasp", registered[-1])
    write_clean_poscar(output / "PO_author_registered.vasp", registered[0])
    report = {"schema_version": 1, "source": str(source), "source_sha256": sha256(source),
              "target_PO": str(target), "target_PO_sha256": sha256(target),
              "rotation_fractional": rotation.tolist(), "rotation_cartesian": rotation.tolist(),
              "translation_fractional": shift.tolist(), "source_to_target_indices": target_indices.tolist(),
              "assignment_max_distance_A": float(np.max(distances[source_indices, target_indices])),
              "assignment_RMS_distance_A": float(np.sqrt(np.mean(distances[source_indices, target_indices]**2))),
              "minimum_assignment_margin_A": float(np.min(margins)),
              "same_permutation_all_twenty_images": True,
              "generated_sha256": {p.name: sha256(p) for p in output.iterdir()},
              "endpoint_audits": {"PO_author_registered": symmetry_report(registered[0]),
                                  "M_seed_registered": symmetry_report(registered[-1])},
              "limitations": "coordinate/species gauge registration only; replace seed endpoints by same-contract relaxed endpoints before production"}
    report["public_source"] = {"paper": "Ma and Liu, Physical Review Letters 130, 096801 (2023)",
                               "doi": "https://doi.org/10.1103/PhysRevLett.130.096801",
                               "repository": "https://github.com/sliutheorygroup/supplementary-material/tree/2804a091ad7d36d3f6d068b9e58da5139ce9f8e8/L23",
                               "raw_record": "figure2/figure2a/transition-poscars/pca21-M",
                               "scope": "atomic structures only; neither VASP POTCAR nor electronic output redistributed"}
    (output / "registration.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "target", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    report = register(args.source, args.target, args.output)
    print(json.dumps({k: report[k] for k in ("assignment_RMS_distance_A", "assignment_max_distance_A", "minimum_assignment_margin_A")}))


if __name__ == "__main__":
    main()
