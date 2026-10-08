"""Prepare the first common-PO+ decay channel using a registered author seed."""

import argparse
import json
from pathlib import Path

from ase.io import read, write
import numpy as np

from examples.hfo2_fixed_input_factory import CONTRACT, same_ordered_geometry
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.prepare_hfo2_reference_variants import write_clean_poscar
from vcneb import endpoint_structure_record, validate_path_geometry, minimum_image_path_lift, validate_periodic_path_lift
from vcneb.analysis import path_reaction_coordinate


def uniform_seed(controls, initial, final, n_total=9):
    """Ordered piecewise-linear geometry seed; no new atom assignment or DFT."""
    controls = [a.copy() for a in controls]
    controls[0], controls[-1] = initial.copy(), final.copy()
    controls, _ = minimum_image_path_lift(controls)
    for a in controls:
        a.calc = None
    arc, _ = path_reaction_coordinate(controls)
    if arc[-1] <= 0 or np.any(np.diff(arc) <= 0):
        raise ValueError("registered control chain has a degenerate segment")
    images = [initial.copy()]
    for s in np.linspace(0, arc[-1], n_total)[1:-1]:
        i = min(int(np.searchsorted(arc, s, side="right") - 1), len(controls) - 2)
        fraction = (s - arc[i]) / (arc[i+1] - arc[i])
        left, right = controls[i], controls[i+1]
        q = left.get_scaled_positions(wrap=False)
        dq = right.get_scaled_positions(wrap=False) - q
        a = left.copy()
        a.set_cell((1-fraction)*left.cell.array + fraction*right.cell.array, scale_atoms=False)
        a.set_scaled_positions(q + fraction*dq)
        images.append(a)
    images.append(controls[-1].copy())
    for a in images:
        a.calc = None
    validate_periodic_path_lift(images)
    return images


def prepare(author_chain, po_path, po_static, m_root, m_audit_path, output):
    if output.exists():
        raise FileExistsError("refusing existing PO-M namespace")
    m_audit = json.loads(m_audit_path.read_text(encoding="utf-8"))
    if (m_audit["status"] != "raw_endpoint_gates_passed"
            or m_audit["symmetry_audit"]["symmetry"][-1]["number"] != 14):
        raise ValueError("M does not pass raw/declared phase-identity gate")
    po, m = read(po_path, format="vasp"), read(m_root / "CONTCAR", format="vasp")
    m_static = Path(m_audit["final_static_source"])
    endpoint_results = []
    for atoms, source in ((po, po_static), (m, m_static)):
        if not same_ordered_geometry(atoms, read(source / "STRU", format="abacus")):
            raise ValueError("PO or M raw endpoint differs from requested geometry")
        if any(sha256(source / n) != h for n, h in CONTRACT.items()):
            raise ValueError("PO/M endpoint electronic contract differs")
        raw = audited_results(source)
        endpoint_results.append({"energy_eV_cell": float(raw["energy"]),
                                 "structure": endpoint_structure_record(atoms),
                                 "source": str(source), "raw_log_sha256": sha256(source / "OUT.ABACUS/running_scf.log")})
    controls = read(author_chain, index=":")
    if len(controls) != 20:
        raise ValueError("expected registered twenty-image public seed")
    images = uniform_seed(controls, po, m)
    geometry = validate_path_geometry(images, minimum_distance=1.6, maximum_deformation=.25)
    output.mkdir(parents=True)
    write(output / "seed.traj", images)
    write_clean_poscar(output / "initial.vasp", images[0])
    write_clean_poscar(output / "final.vasp", images[-1])
    parameters = {"source_directory": str(po_static),
                  "seed_static_directories": [str(po_static)] + [None]*7 + [str(m_static)]}
    (output / "factory_parameters.json").write_text(json.dumps(parameters, indent=2) + "\n", encoding="utf-8")
    arc, segments = path_reaction_coordinate(images)
    manifest = {"purpose": "same_POplus_to_M_unconstrained_decay_channel_pilot",
                "n_total_images": 9, "n_fixed_endpoints": 2, "n_active_images": 7,
                "first_optimization_step_cap": 10, "fmax_eV_A": .10, "climb": False,
                "optimizer": "FIRE", "maxstep_A": .02, "spring_eV_A2": .2,
                "pressure_GPa": 0, "physical_inputs_changed": False,
                "endpoints": endpoint_results, "M_raw_audit_sha256": sha256(m_audit_path),
                "author_chain_sha256": sha256(author_chain),
                "initial_geometry_gate": geometry, "initial_arc_A": arc.tolist(),
                "initial_segment_lengths_A": segments.tolist(),
                "mapping": "registered fixed permutation retained; no per-image relabel; endpoints retain exact ordered periodic geometries",
                "periodic_lift_audit": validate_periodic_path_lift(images),
                "limitations": "published geometries guide only the starting chain; neither source energies nor VASP parameters are reused"}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("author-chain", "po", "po-static", "m-root", "m-audit", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    a = parser.parse_args()
    report = prepare(a.author_chain, a.po, a.po_static, a.m_root, a.m_audit, a.output)
    print(json.dumps({k: report[k] for k in ("purpose", "n_total_images", "n_active_images")}))


if __name__ == "__main__":
    main()
