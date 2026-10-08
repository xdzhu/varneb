"""Audit every HfO2 Gamma/seed static before constructing force constants."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase.io import read

from scripts.audit_hfo2_static_replica import audited_results, sha256
from scripts.prepare_hfo2_gamma_and_seeds import gamma_family
from vcneb.phonons import identify_acoustic_modes, phonopy_gamma_eigenpairs


def point_audit(root, index):
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    if manifest["purpose"] != "HfO2_T_Gamma_two_step_sizes_and_three_seed_statics":
        raise ValueError("unexpected static manifest")
    record = manifest["points"][index]
    if record["index"] != index:
        raise ValueError("manifest index/order changed")
    directory = root / "calculations" / f"{index:02d}"
    if any(sha256(directory / n) != h for n, h in record["input_sha256"].items()):
        raise ValueError("effective static input differs from manifest")
    raw = audited_results(directory)
    result = {"index": index, "status": "passed", "kind": record["kind"],
              "input_sha256": record["input_sha256"], "raw_log_sha256": sha256(directory / "OUT.ABACUS/running_scf.log"),
              "energy_eV_cell": float(raw["energy"]), "forces_eV_A": raw["forces"].tolist(),
              "stress_eV_A3": raw["stress"].tolist(),
              "max_force_eV_A": float(np.max(np.linalg.norm(raw["forces"], axis=1))),
              "max_stress_kbar": float(np.max(np.abs(raw["stress"])) * 1602.176634),
              "net_force_eV_A": raw["forces"].sum(axis=0).tolist()}
    (directory / "point_audit.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    return result


def audit(root, output):
    if output.exists():
        raise FileExistsError("refusing existing Gamma audit")
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    if manifest["n_static_points"] != 19 or len(manifest["points"]) != 19:
        raise ValueError("incomplete declared nineteen-point batch")
    if sha256(root / "T_reference.vasp") != manifest["T_reference_sha256"]:
        raise ValueError("reference geometry changed")
    reference = read(root / "T_reference.vasp", format="vasp")
    points = [point_audit(root, i) for i in range(19)]
    families, arrays = {}, {}
    for distance in (.01, .02):
        phonon = gamma_family(reference, distance)
        records = [r for r in manifest["points"] if r["kind"] == "T_Gamma" and r["distance_A"] == distance]
        if len(records) != 8 or [r["local_index"] for r in records] != list(range(8)):
            raise ValueError("Gamma family count/order changed")
        for r, expected in zip(records, phonon.dataset["first_atoms"]):
            if r["atom_index"] != expected["number"] or not np.allclose(r["displacement_A"], expected["displacement"], atol=1e-10):
                raise ValueError("displacement contract differs from regenerated family")
        forces = np.asarray([points[r["index"]]["forces_eV_A"] for r in records])
        phonon.produce_force_constants(forces=forces, show_drift=False)
        raw = phonon.force_constants.copy()
        raw_pairs = phonopy_gamma_eigenpairs(phonon)
        # Preserve raw diagnostics; ASR/permutation projection is a reported
        # post-processing step, not evidence that the underlying DFT is exact.
        phonon.symmetrize_force_constants(level=1)
        sym = phonon.force_constants.copy()
        pairs = phonopy_gamma_eigenpairs(phonon)
        acoustic, overlap = identify_acoustic_modes(pairs.eigenvectors_mass_weighted, phonon.primitive.masses)
        label = f"T_d{distance:.2f}"
        arrays[label] = {"force_constants_raw_eV_A2": raw,
                         "force_constants_symmetrized_eV_A2": sym,
                         "frequencies_raw_thz": raw_pairs.frequencies_thz,
                         "frequencies_thz": pairs.frequencies_thz,
                         "eigenvectors_mass_weighted": pairs.eigenvectors_mass_weighted,
                         "masses_amu": np.asarray(phonon.primitive.masses)}
        families[label] = {"distance_A": distance, "n_displacements": len(records),
            "frequencies_thz": pairs.frequencies_thz.tolist(), "frequencies_raw_thz": raw_pairs.frequencies_thz.tolist(),
            "acoustic_mode_indices": acoustic.tolist(), "minimum_translation_principal_overlap": overlap,
            "max_abs_acoustic_frequency_thz": float(np.max(np.abs(pairs.frequencies_thz[acoustic]))),
            "negative_optical_indices": [int(i) for i in range(36) if i not in acoustic and pairs.frequencies_thz[i] < -.05],
            "max_raw_ASR_drift_eV_A2": float(np.max(np.abs(raw.sum(axis=1)))),
            "max_symmetrization_change_eV_A2": float(np.max(np.abs(sym - raw))),
            "max_raw_permutation_asymmetry_eV_A2": float(np.max(np.abs(raw - raw.transpose(1, 0, 3, 2)))),
            "max_displaced_net_force_eV_A": float(np.max(np.linalg.norm(forces.sum(axis=1), axis=1))),
            "interpretation": "negative-mode list uses a reporting cutoff only; significance needs two-step/error checks"}
    small, large = arrays["T_d0.01"], arrays["T_d0.02"]
    report = {"schema_version": 1, "status": "nineteen_statics_audited_Gamma_requires_physical_interpretation",
              "manifest_sha256": sha256(root / "manifest.json"), "families": families,
              "step_size_comparison": {
                  "max_raw_fc_difference_eV_A2": float(np.max(np.abs(small["force_constants_raw_eV_A2"] - large["force_constants_raw_eV_A2"]))),
                  "max_sym_fc_difference_eV_A2": float(np.max(np.abs(small["force_constants_symmetrized_eV_A2"] - large["force_constants_symmetrized_eV_A2"]))),
                  "max_sorted_frequency_difference_thz": float(np.max(np.abs(small["frequencies_thz"] - large["frequencies_thz"])))} ,
              "candidate_statics": {r["name"]: points[r["index"]] for r in manifest["points"] if r["kind"] != "T_Gamma"},
              "points": points, "limitations": manifest["Gamma_boundary"],
              "units": {"positions": "angstrom", "forces": "eV/angstrom", "force_constants": "eV/angstrom^2"},
              "postprocessing": "raw force constants retained; level-1 Phonopy ASR/permutation symmetrization reported separately"}
    output.mkdir(parents=True)
    for label, values in arrays.items():
        np.savez_compressed(output / f"{label}.npz", **values)
    (output / "audit.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.index is not None:
        result = point_audit(args.root, args.index)
        print(json.dumps({"index": args.index, "status": result["status"]}))
    else:
        if args.output is None:
            parser.error("--output required for complete Gamma analysis")
        result = audit(args.root, args.output)
        print(json.dumps({"status": result["status"], "step_size_comparison": result["step_size_comparison"],
                          "families": {k: {x: v[x] for x in ("negative_optical_indices", "max_abs_acoustic_frequency_thz")}
                                       for k, v in result["families"].items()}}))


if __name__ == "__main__":
    main()
