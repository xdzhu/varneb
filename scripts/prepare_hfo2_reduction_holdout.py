"""Freeze eight independent tests of the atomic harmonic response prediction."""

import argparse
import json
from pathlib import Path
import shutil

from ase.io import read
import numpy as np

from examples.hfo2_fixed_input_factory import CONTRACT, same_ordered_geometry
from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.prepare_hfo2_channel_work_probes import write_structure
from scripts.prepare_hfo2_gamma_and_seeds import geometry_roundtrip


def holdout_vectors(model):
    q = np.array(model["retained_basis_cartesian"], dtype=float)
    response = np.array(model["families"]["d0.01"]["full_atomic_response"], dtype=float)
    if q.shape != (36, 3) or response.shape != q.shape or not np.isfinite(response).all():
        raise ValueError("invalid response/basis dimensions or values")
    np.testing.assert_allclose(q.T @ q, np.eye(3), atol=1e-10)
    np.testing.assert_allclose(q.T @ response, 0, atol=1e-10)
    vectors = {"frozen": q[:, 0], "responded": q[:, 0] + response[:, 0]}
    for vector in vectors.values():
        np.testing.assert_allclose(q.T @ vector, [1, 0, 0], atol=1e-10)
        np.testing.assert_allclose(vector.reshape(12, 3).sum(axis=0), 0, atol=1e-10)
    return vectors


def prepare(source, variants, model_path, output):
    if output.exists():
        raise FileExistsError("refusing existing holdout namespace")
    model = json.loads(model_path.read_text(encoding="utf-8"))
    if sha256(variants / "T.vasp") != model["T_reference_sha256"]:
        raise ValueError("model T geometry changed")
    if any(sha256(source / n) != h for n, h in CONTRACT.items()):
        raise ValueError("historical electronic contract changed")
    t = read(variants / "T.vasp", format="vasp")
    if not same_ordered_geometry(t, read(source / "STRU", format="abacus")):
        raise ValueError("T baseline/source geometry mismatch")
    raw = audited_results(source)
    vectors = holdout_vectors(model)
    family = model["families"]["d0.01"]
    predicted = {"frozen": family["frozen_curvature_eV_A2"][0][0],
                 "responded": family["relaxed_curvature_eV_A2"][0][0]}
    output.mkdir(parents=True)
    points = []
    for kind, vector in vectors.items():
        for h in (.05, .10):
            for sign in (-1, 1):
                index = len(points)
                shifted = t.copy()
                shifted.positions += sign * h * vector.reshape(12, 3)
                directory = output / "points" / f"{index:02d}"
                directory.mkdir(parents=True)
                for name in CONTRACT:
                    shutil.copyfile(source / name, directory / name)
                write_structure(directory / "STRU", shifted)
                minimum = geometry_roundtrip(directory / "STRU", shifted)
                points.append({"index": index, "kind": kind, "amplitude_A": h, "sign": sign,
                               "minimum_distance_A": minimum,
                               "input_sha256": {n: sha256(directory / n) for n in (*CONTRACT, "STRU")}})
    record = {"purpose": "HfO2_T_atomic_reduction_8_independent_holdouts", "model_sha256": sha256(model_path),
              "prediction_fixed_before_DFT": predicted, "direction_vectors": {k: v.tolist() for k, v in vectors.items()},
              "T_energy_eV_cell": float(raw["energy"]), "T_baseline_log_sha256": sha256(source / "OUT.ABACUS/running_scf.log"),
              "relative_curvature_acceptance": .10, "points": points,
              "retained_basis_cartesian": model["retained_basis_cartesian"],
              "coordinate": "Qx projection in Cartesian Angstrom; response vector not renormalized; Qy/Qz/cell held fixed",
              "limitations": "linear harmonic response line, not fully minimized conditional surface or joint-cell result"}
    (output / "manifest.json").write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "variants", "model", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    report = prepare(args.source, args.variants, args.model, args.output)
    print(json.dumps({"n_independent_points": len(report["points"]), "fixed_prediction": report["prediction_fixed_before_DFT"]}))


if __name__ == "__main__":
    main()
