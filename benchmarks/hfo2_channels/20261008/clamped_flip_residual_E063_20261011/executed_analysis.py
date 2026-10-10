"""Zero-DFT decomposition of an audited ordinary clamped G2 residual.

This is a snapshot diagnosis, not a changed optimizer, a new barrier label,
or a prediction of the improvement obtainable by changing spring constants.
The generalized atomic block is conjugate to q @ H_reference, not the raw
Cartesian force; the cell block uses the registered length scale.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read

from scripts.audit_hfo2_static_replica import sha256
from scripts.export_hfo2_clamped_observation import CELL_SCALE
from vcneb import VCNEB, clamped_plane_vcneb_boundary


def max_vector(vector):
    return float(np.linalg.norm(np.asarray(vector).reshape(-1, 3), axis=1).max())


def decompose(chain):
    """Replay the implemented ordinary force, then split its vector blocks.

    Private coordinate/tangent access is intentionally local to this analysis
    tool. Production VCNEB and its diagnostics remain unchanged; the explicit
    vector sum must reproduce the public force before anything is reported.
    """
    if chain.climb or chain.dynamic_relaxation != 1.0 or chain.n_images < 3:
        raise ValueError("only unweighted ordinary NEB snapshots are supported")
    residuals = chain.get_forces().reshape(chain.n_images-2, chain.image_ndofs)
    mask = chain._active_x_mask()
    coordinates = [chain._image_x(i)*mask for i in range(chain.n_images)]
    split = 3*chain.n_atoms
    records = []
    for i in range(1, chain.n_images-1):
        tangent = chain._tangent(i, chain.enthalpies, coordinates)
        true = chain._project_constraint_force(chain._last_true_forces_x[i]*mask)*mask
        perpendicular = true - np.dot(true, tangent)*tangent
        d_minus = float(np.linalg.norm(coordinates[i]-coordinates[i-1]))
        d_plus = float(np.linalg.norm(coordinates[i+1]-coordinates[i]))
        spring = chain.k[i-1]*(d_plus-d_minus)*tangent
        perpendicular = chain._project_constraint_force(perpendicular)*mask
        spring = chain._project_constraint_force(spring)*mask
        residual = residuals[i-1]
        error = float(np.max(np.abs(perpendicular+spring-residual)))
        if not np.isfinite(residual).all() or error > 1e-10:
            raise ValueError("component sum differs from the implemented ordinary residual")
        vectors = residual.reshape(-1, 3)
        winner = int(np.argmax(np.linalg.norm(vectors, axis=1)))
        norm = float(np.linalg.norm(vectors[winner]))
        direction = vectors[winner]/norm if norm else np.zeros(3)
        record = dict(image_index=i, worst_vector_index=winner,
                      worst_vector_block="atoms" if winner < chain.n_atoms else "open_cell",
                      residual_max_vector_eV_A=norm,
                      component_sum_max_abs_error_eV_A=error,
                      spacing_minus_A=d_minus, spacing_plus_A=d_plus,
                      spring_signed_scalar_eV_A=float(chain.k[i-1]*(d_plus-d_minus)),
                      perpendicular_spring_dot_eV2_A2=float(np.dot(perpendicular, spring)),
                      perpendicular_projection_on_worst_residual_eV_A=float(
                          np.dot(perpendicular.reshape(-1, 3)[winner], direction)),
                      spring_projection_on_worst_residual_eV_A=float(
                          np.dot(spring.reshape(-1, 3)[winner], direction)))
        for name, value in (("residual", residual), ("perpendicular", perpendicular), ("spring", spring)):
            record[name+"_atom_max_vector_eV_A"] = max_vector(value[:split])
            record[name+"_open_cell_max_vector_eV_A"] = max_vector(value[split:])
            record[name+"_euclidean_eV_A"] = float(np.linalg.norm(value))
            record[name+"_generalized_vector_eV_A"] = value.tolist()
        records.append(record)
    winner = max(records, key=lambda r: r["residual_max_vector_eV_A"])
    return dict(images=records, limiting_image_index=winner["image_index"],
                limiting_vector_block=winner["worst_vector_block"],
                fmax_eV_A=winner["residual_max_vector_eV_A"],
                perpendicular_only_fmax_same_geometry_eV_A=max(
                    max(r["perpendicular_atom_max_vector_eV_A"],
                        r["perpendicular_open_cell_max_vector_eV_A"]) for r in records),
                spring_only_fmax_same_geometry_eV_A=max(
                    max(r["spring_atom_max_vector_eV_A"],
                        r["spring_open_cell_max_vector_eV_A"]) for r in records))


def analyze(observation_directory):
    record_path = observation_directory/"observation.json"
    record = json.loads(record_path.read_text(encoding="utf-8"))
    trajectory = observation_directory/"evaluated_chain.traj"
    if (record["new_DFT_calls"] != 0 or record["climb"] or record["pressure_GPa"] != 0
            or record["cell_scale_A"] != CELL_SCALE or record["k_eV_A2"] != .2
            or record["n_total_images"] != 9 or record["n_active_images"] != 7
            or sha256(trajectory) != record["evaluated_chain_sha256"]):
        raise ValueError("audited G2 observation or pinned trajectory contract differs")
    mechanical = record["mechanical_boundary"]
    if (mechanical["kind"] != "clamped_plane" or mechanical["allow_tilt"] is not True
            or mechanical["fixed_cell_rows"] != [0, 1] or mechanical["cell_dofs"] != 3):
        raise ValueError("registered open-cell boundary required")
    images = read(trajectory, index=":")
    if len(images) != 9 or len(record["raw_image_evaluations"]) != 9:
        raise ValueError("nine evaluated images required")
    for i, (image, raw) in enumerate(zip(images, record["raw_image_evaluations"])):
        if not isinstance(image.calc, SinglePointCalculator) or raw["image_index"] != i:
            raise ValueError("single-point-only ordered cache replay required")
        for key, expected in (("energy", raw["energy_eV_cell"]), ("forces", raw["forces_eV_A"]),
                              ("stress", raw["stress_ASE_voigt_eV_A3"])):
            if not np.allclose(image.calc.results[key], expected, atol=1e-12, rtol=0):
                raise ValueError("trajectory differs from audited native E/F/stress")
    boundary = clamped_plane_vcneb_boundary(12, np.array(mechanical["reference_cell_A"]), allow_tilt=True)
    chain = VCNEB(images, cell_scale=CELL_SCALE, k=.2, climb=False, pressure=0,
                  **boundary.vcneb_kwargs(images))
    result = decompose(chain)
    if abs(result["fmax_eV_A"]-record["replayed_fmax_eV_A"]) > 1e-10:
        raise ValueError("snapshot residual no longer matches its audited replay")
    return dict(format_version=1, source_job_id=record["source_job_id"],
                snapshot_step=record["snapshot_step"], new_DFT_calls=0,
                observation_sha256=sha256(record_path), analyzed_trajectory_sha256=sha256(trajectory),
                ordinary_residual_pass=bool(result["fmax_eV_A"] <= .10),
                **result, limitations=[
                    "instantaneous generalized-coordinate residual, not raw Cartesian atomic force",
                    "cell-block norms depend on the unchanged registered metric",
                    "component maxima do not add; signed projections at one vector do add",
                    "removing a spring at fixed geometry is not a reoptimized path benchmark",
                    "not a converged barrier, certified TS, G3 selection or holdout prediction"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("observation_directory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = analyze(args.observation_directory)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2); stream.write("\n")
    print(json.dumps({k:result[k] for k in ("snapshot_step", "fmax_eV_A", "limiting_vector_block")}))


if __name__ == "__main__":
    main()
