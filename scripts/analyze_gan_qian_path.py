"""Offline GaN path audit and per-formula literature comparison.

No DFT, geometry repairs, parameter changes or symmetry projection are performed.
Space groups and the hexagonal-branch hint are diagnostics, not saddle proofs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import sys

import ase
from ase import units
from ase.io.trajectory import Trajectory
import numpy as np
from scipy.signal import find_peaks
import spglib

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from vcneb import endpoint_structure_record, validate_path_geometry

REFERENCES = {"tetragonal": (45.7, .34), "hexagonal": (45.7, .39), "b3_diagonal": (45.0, .57)}
POTCAR_SHA256 = "f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168"


def latest_evaluated_chain(path, n_images):
    # The resume helper intentionally removes calculators. Analysis must keep
    # the saved SinglePoint results, never trigger a new electronic calculation.
    with Trajectory(str(path), "r") as trajectory:
        if len(trajectory) < n_images or len(trajectory) % n_images:
            raise ValueError("completed trajectory has missing or partial chain frames")
        return [trajectory[index] for index in range(len(trajectory) - n_images, len(trajectory))]


def formula_units(atoms):
    symbols = atoms.get_chemical_symbols()
    if set(symbols) != {"Ga", "N"} or symbols.count("Ga") != symbols.count("N"):
        raise ValueError("expected a stoichiometric GaN cell")
    return symbols.count("Ga")


def peak_segments(enthalpies_per_formula, prominence=.02):
    """Return resolved peaks and barriers relative to each local segment start."""
    values = np.asarray(enthalpies_per_formula, dtype=float)
    if values.ndim != 1 or len(values) < 3 or not np.all(np.isfinite(values)):
        raise ValueError("a finite single enthalpy chain is required")
    if not np.isfinite(prominence) or prominence <= 0:
        raise ValueError("peak prominence must be finite and positive")
    peaks, properties = find_peaks(values, prominence=prominence)
    valleys = [int(left + np.argmin(values[left:right + 1]))
               for left, right in zip(peaks[:-1], peaks[1:])]
    boundaries = [0, *valleys, len(values) - 1]
    segments = []
    for start, peak, stop, height in zip(boundaries[:-1], peaks, boundaries[1:], properties["prominences"]):
        segments.append({"start_index": start, "peak_index": int(peak), "stop_index": stop,
                         "forward_barrier_eV_per_GaN": float(values[peak] - values[start]),
                         "reverse_barrier_eV_per_GaN": float(values[peak] - values[stop]),
                         "reaction_enthalpy_eV_per_GaN": float(values[stop] - values[start]),
                         "peak_prominence_eV_per_GaN": float(height)})
    return {"peak_indices": peaks.astype(int).tolist(), "valley_indices": valleys,
            "segments": segments, "prominence_threshold_eV_per_GaN": prominence}


def structure_row(atoms, index):
    volume = float(atoms.get_volume())
    cell = (atoms.cell.array, atoms.get_scaled_positions(wrap=True), atoms.numbers)
    symmetry = {}
    for tolerance in (1e-5, 1e-4, 1e-3, 1e-2):
        dataset = spglib.get_symmetry_dataset(cell, symprec=tolerance, angle_tolerance=-1)
        symmetry[str(tolerance)] = None if dataset is None else {
            "number": int(dataset.number), "international": dataset.international}
    row = {"image_index": index, "volume_A3": volume,
           "cell_lengths_A": atoms.cell.lengths().tolist(),
           "cell_angles_deg": atoms.cell.angles().tolist(), "symmetry": symmetry}
    if atoms.get_chemical_symbols() == ["Ga", "Ga", "N", "N"]:
        q = atoms.get_scaled_positions(wrap=False)
        offsets = np.mod(q[2:, 2] - q[:2, 2], 1)
        row.update(u_fractional_pairs=offsets.tolist(),
                   u_half_deviation=float(np.max(np.abs(offsets - .5))),
                   basal_angle_deg=float(atoms.cell.angles()[2]),
                   paper_complementary_gamma_deg=float(180 - atoms.cell.angles()[2]))
    return row


def analyze(summary, images, route, generic_audit, manifest_records, prominence=.02):
    pressure_gpa, literature = REFERENCES[route]
    if summary.get("status") != "completed" or not summary.get("converged"):
        raise ValueError("path must be completed and converged")
    if generic_audit.get("status") != "ok" or generic_audit.get("issues"):
        raise ValueError("generic VCNEB result audit must pass first")
    if summary["final_max_generalized_force_eV_per_A"] > .10 + 1e-12:
        raise ValueError("path exceeds the requested 0.10 eV/A residual")
    parameters = summary["calculator_parameters"]
    expected = {"encut": 600., "ediff": 1e-7, "isym": -1, "symprec": 1e-4,
                "nsw": 0, "ibrion": -1, "isif": 2, "xc": "PBE", "pp": "PBE", "gamma": True}
    if any(parameters.get(key) != value for key, value in expected.items()):
        raise ValueError("fixed GaN calculation contract changed")
    if parameters.get("kpts") != [8, 8, 6] or parameters.get("setups", {}).get("Ga") != "_d":
        raise ValueError("Ga_d or k-point contract changed")
    if summary["licensed_input_fingerprints"]["POTCAR"]["sha256"] != POTCAR_SHA256:
        raise ValueError("Ga_d+N PAW fingerprint changed")
    if summary.get("runtime_parameter_changes_allowed") is not False:
        raise ValueError("runtime input policy is not locked")
    if summary.get("endpoint_evaluation_policy") != "fixed_cached_once":
        raise ValueError("fixed endpoint cache policy is missing")
    n_images = summary["n_images"]
    if len(images) != n_images or n_images != 29:
        raise ValueError("expected the complete 29-image literature chain")
    if any(not np.all(atoms.pbc) for atoms in images):
        raise ValueError("GaN analysis requires three-dimensional periodicity")
    count = formula_units(images[0])
    expected_count = 4 if route == "b3_diagonal" else 2
    if count != expected_count or any(formula_units(atoms) != count for atoms in images):
        raise ValueError("formula-unit count changed along the chain")
    if any(atoms.get_chemical_symbols() != images[0].get_chemical_symbols() for atoms in images):
        raise ValueError("ordered species changed along the chain")
    pressure = summary["path_diagnostics"]["pressure_eV_per_A3"]
    if not np.isclose(pressure / units.GPa, pressure_gpa, rtol=0, atol=1e-5):
        raise ValueError("path pressure differs from its literature reference")
    geometry = validate_path_geometry(images, minimum_distance=1.4, maximum_deformation=.50)
    enthalpies = np.asarray(summary["image_enthalpies_eV"], dtype=float)
    calculated = np.array([atoms.get_potential_energy() + pressure * atoms.get_volume() for atoms in images])
    if not np.allclose(calculated, enthalpies, rtol=0, atol=1e-8):
        raise ValueError("trajectory snapshot and summary enthalpies differ")
    for label, atoms in (("initial", images[0]), ("final", images[-1])):
        record = endpoint_structure_record(atoms)
        # A wrapped coordinate at a cell boundary can round to 0 or 1 at file
        # precision. Compare periodically, while retaining ordered species.
        expected_record = summary["endpoint_structures"][label]
        difference = np.asarray(record["fractional_positions_wrapped"]) - expected_record["fractional_positions_wrapped"]
        difference -= np.rint(difference)
        if (record["species_order"] != expected_record["species_order"]
                or not np.allclose(record["cell_A"], expected_record["cell_A"], atol=1e-9, rtol=0)
                or not np.allclose(difference, 0, atol=1e-9, rtol=0)):
            raise ValueError("trajectory and recorded endpoint identity differ")
    diagnostics = summary["path_diagnostics"]["images"]
    if any(record.get("is_climbing_image") for record in diagnostics):
        raise ValueError("this goal requires ordinary VCNEB without CI")
    if not manifest_records or any(record.get("status") != "ok"
            or record.get("image_indices") != list(range(1, n_images - 1)) for record in manifest_records):
        raise ValueError("worker manifest does not prove successful interior-only evaluation")
    rows = [structure_row(atoms, index) for index, atoms in enumerate(images)]
    values = (enthalpies - enthalpies[0]) / count
    peaks = peak_segments(values, prominence)
    for segment in peaks["segments"]:
        segment["difference_from_literature_eV_per_GaN"] = segment["forward_barrier_eV_per_GaN"] - literature
    peak = int(np.argmax(values))
    hex_visits = [row["image_index"] for row in rows[:peak + 1]
                  if row.get("u_half_deviation", 1) < .015 and row.get("basal_angle_deg", 0) > 110]
    return {"status": "numerically_audited_structure_diagnostics_require_interpretation",
            "route_requested": route, "external_pressure_gpa": pressure_gpa,
            "n_formula_units": count, "barrier_eV_per_GaN": float(values[peak]),
            "reaction_enthalpy_eV_per_GaN": float(values[-1]), "highest_image_index": peak,
            "literature_barrier_eV_per_GaN": literature,
            "difference_from_literature_eV_per_GaN": float(values[peak] - literature),
            "relative_difference_percent": float(100 * (values[peak] / literature - 1)),
            "resolved_extrema": peaks, "structures": rows, "geometry": geometry,
            "worker_batches": len(manifest_records),
            "hexagonal_branch_hint": {"prepeak_h_MgO_like_indices": hex_visits,
                "criterion": "paired fractional u within 0.015 of 0.5 and basal angle >110 degrees before highest image",
                "not_a_space_group_or_saddle_certificate": True},
            "symmetry_parameters": {"symprec_A": [1e-5, 1e-4, 1e-3, 1e-2], "angle_tolerance_deg": -1},
            "method_difference": "VASP/PBE/Ga_d+N PAW vs literature QE/PW91/ultrasoft; no target-value tuning",
            "relative_enthalpies_eV_per_GaN": values.tolist()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--trajectory", required=True, type=Path)
    parser.add_argument("--audit", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--route", required=True, choices=REFERENCES)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--peak-prominence", type=float, default=.02, help="Diagnostic eV/GaN cutoff, not NEB convergence")
    args = parser.parse_args()
    if args.output.exists() or args.output.with_suffix(".csv").exists():
        raise FileExistsError("use a new output path to preserve previous analyses")
    summary = json.loads(args.summary.read_text())
    images = latest_evaluated_chain(args.trajectory, summary["n_images"])
    records = [json.loads(line) for line in args.manifest.read_text().splitlines() if line.strip()]
    report = analyze(summary, images, args.route, json.loads(args.audit.read_text()), records, args.peak_prominence)
    report["input_sha256"] = {name: hashlib.sha256(path.read_bytes()).hexdigest()
        for name, path in (("summary", args.summary), ("trajectory", args.trajectory),
                           ("audit", args.audit), ("manifest", args.manifest))}
    report["software_versions"] = {"python": platform.python_version(), "ase": ase.__version__,
                                   "numpy": np.__version__, "spglib": spglib.__version__}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    with args.output.with_suffix(".csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["image_index", "relative_enthalpy_eV_per_GaN", "volume_A3", "basal_angle_deg", "u_half_deviation"])
        for row, enthalpy in zip(report["structures"], report["relative_enthalpies_eV_per_GaN"]):
            writer.writerow([row["image_index"], enthalpy, row["volume_A3"], row.get("basal_angle_deg"), row.get("u_half_deviation")])
    print(json.dumps({key: report[key] for key in ("status", "route_requested", "barrier_eV_per_GaN", "resolved_extrema", "hexagonal_branch_hint")}, indent=2))


if __name__ == "__main__":
    main()
