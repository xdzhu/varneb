"""Reconcile the eight frozen G2 work reports with manuscript section 3.8.

Reparse portable native E/F/stress and independently contract segment work.
This does not recreate unpublished mode arrays, audit undistributed runtime
files, certify saddles, select a model, or access the reserved strain labels.
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from scripts.analyze_hfo2_clamped_reference import pinned_json, terminal_images
from scripts.audit_hfo2_static_replica import sha256
from vcneb.strain_work import configuration_work

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "benchmarks/hfo2_channels/20261008"
# Manuscript Table 7 order. Only already admitted training receipts are read.
CASES = (
    ("PO-to-M, zero strain", "clamped_M_terminal_E068_20261011",
     "56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c",
     "clamped_path_work_E092_20261011/zero_M.json",
     "a15e2374d2971fec1b3dcfe88c972408c69341ce5312d051b8417b87b3f8233b"),
    ("PO-to-T, zero strain", "clamped_T_terminal_E079_20261011",
     "c5653b4759b202199cb75739f99038aaa52b3ccc220ef226efbacc980e5d716e",
     "clamped_path_work_E092_20261011/zero_T.json",
     "3f18d539a0d6c0585e3daebf96745c4a26f8cf405d995fa2d893ed8294266567"),
    ("Reversing flip, zero strain", "clamped_reversing_terminal_E083_20261011",
     "632e0e9b2b821b3bfc7d66b361c44c7db734988fb5d018ea2f31441fab14ed41",
     "clamped_path_work_E092_20261011/zero_reversing.json",
     "e6a7a71e362dcfa628ab59c88c1a0ac6f879898577887f64c727ede70e540275"),
    ("PO-to-T, +1% strain", "clamped_plus1_T_terminal_E087_20261011",
     "32a6e49a618e6e62632ae71a5e0ed093095254f6858221be2ade3d83cc194a28",
     "clamped_path_work_E092_20261011/plus1_T.json",
     "0001414317551871351bbf2b241708817076be88ddab3b62f8b27ed3168df90a"),
    ("Preserving flip, zero strain", "clamped_held_flips_E100_20261011/zero_strain_preserving",
     "46145da14f24b2b8bd46aa4656f0982650d6b2dc661282f58d809aab8a6d3202",
     "clamped_preserving_reference_E101_20261011/path_work.json",
     "8177333af0f86146919ef8a39f8fbfcc66835d0bbe9be01dfbf2e7bb6994d56e"),
    ("PO-to-M, +1% strain", "clamped_reversing_release_E104_20261012/plus1_M_terminal",
     "b86cdb72eab9e3d8ff2125e16e6fbb2d0954b34d54c0f05099a581c066108389",
     "clamped_M_strain_reference_E105_20261012/path_work.json",
     "8765c96c7149061465bd5ccba9d496382ad86fb44c55723212e4d3a70167d2a3"),
    ("Preserving flip, +1% strain", "clamped_plus1_preserving_terminal_E145_20261012/completed_HF",
     "0c8d91532bc9051decdfd600eae12ce409fd84a86c677158077948fa352c9d8a",
     "clamped_plus1_preserving_terminal_E145_20261012/analysis/path_work.json",
     "135825034f86a6030da44003e413275cbc28c6c4f968e0687178d05e14efff78"),
    ("Reversing flip, +1% strain", "clamped_plus1_reversing_terminal_E148_20261014/completed_HF",
     "450dca48d9c8fec24dd649df868c47b122b340a88eaf89df941811855ef590e7",
     "clamped_plus1_reversing_terminal_E148_20261014/analysis/path_work.json",
     "bc5af060455c507cea4942944abef36bddb7d9bf176c9588d74c25ac128ef084"),
)
REPRESENTATIONS = ("T_geometric_triplet", "T_two_lowest_complete_optical_doublets",
                   "Cmma_frame_103", "Cmma_frame_106", "Cmma_frame_114", "Cmma_frame_127")


def audit_cohort():
    rows, identities = [], defaultdict(list)
    max_defect = 0.
    for label, case, audit_pin, report_path, report_pin in CASES:
        stored = pinned_json(BASE / report_path, report_pin)
        images, native, residual, receipt = terminal_images(BASE / case, audit_pin)
        if (len(images) != 9 or stored["portable_native_frames_reparsed"] != 9
                or stored["terminal_audit_sha256"] != audit_pin
                or stored["source_job_id"] != receipt["job_id"]
                or stored["new_DFT_calls"] != 0 or stored["scheduler_mutations"] != 0
                or stored["G3_selection_or_holdout_access"] is not False
                or stored["complete_six_physical_bytes_rechecked_locally"] is not False
                or list(stored["representations"]) != list(REPRESENTATIONS)):
            raise ValueError("frozen training cohort or provenance scope differs")
        energy = np.array([image.get_potential_energy() for image in images])
        np.testing.assert_allclose(energy, stored["absolute_energy_eV_cell"], atol=1e-12, rtol=0)
        peak = int(np.argmax(energy))
        if peak != stored["sampled_peak_image"] or stored["formula_units"] != 4:
            raise ValueError("sampled peak or formula-unit conversion differs")
        direct = []
        for left, right in zip(images[:-1], images[1:]):
            dh = right.cell.array - left.cell.array
            ds = right.get_scaled_positions(wrap=False) - left.get_scaled_positions(wrap=False)
            work = [configuration_work(image.cell.array, image.get_stress(voigt=False), dh,
                    forces=image.get_forces(apply_constraint=False), fractional_direction=ds)
                    for image in (left, right)]
            direct.append({key: .5 * (work[0][key] + work[1][key]) for key in
                           ("total_eV_per_control", "atomic_eV_per_control", "cell_eV_per_control")})
        for key, span in stored["representations"].items():
            if (span["retained_rank"] != (3 if key == REPRESENTATIONS[0] else 4)
                    or span["source_modes_or_forces_ASR_corrected"] is not False
                    or span["physical_error_bound_measured"] is not False):
                raise ValueError("registered rank or accuracy scope differs")
            parts = span["segment_partition_work"]
            np.testing.assert_allclose(np.sum(list(parts.values()), axis=0),
                                       span["segment_full_work"], atol=1e-12, rtol=0)
            for values, raw in ((span["segment_full_work"], "total_eV_per_control"),
                                (parts["open_cell"], "cell_eV_per_control")):
                defect = float(np.max(np.abs(np.asarray(values) - [r[raw] for r in direct])))
                max_defect = max(max_defect, defect)
                if defect > 1e-12:
                    raise ValueError("stored work differs from native physical contractions")
            np.testing.assert_allclose(span["cumulative_energy_change"], energy-energy[0], atol=1e-12, rtol=0)
        for index, source in enumerate(native["raw_image_evaluations"]):
            identities[(source["raw_log_sha256"], source["input_sha256"]["STRU"])].append([label, index])
        spans = stored["representations"]
        triplet = spans[REPRESENTATIONS[0]]
        cmma = [spans[key]["cumulative_partition_work"]["retained_atomic"][peak]*250
                for key in REPRESENTATIONS[2:]]
        rows.append(dict(label=label, source_job_id=receipt["job_id"], terminal_audit_sha256=audit_pin,
            work_report_sha256=report_pin, ordinary_fmax_eV_A=residual["fmax_eV_A"], peak_image=peak,
            table_meV_fu=[(energy[peak]-energy[0])*250, triplet["cumulative_full_work"][peak]*250,
                          triplet["cumulative_energy_minus_work"][peak]*250,
                          triplet["endpoint_energy_minus_work"]*250],
            peak_triplet_parts_meV_fu={key: values[peak]*250 for key, values in
                                       triplet["cumulative_partition_work"].items()},
            peak_optical_retained_meV_fu=spans[REPRESENTATIONS[1]]["cumulative_partition_work"]["retained_atomic"][peak]*250,
            peak_Cmma_retained_range_meV_fu=[min(cmma), max(cmma)]))
    return dict(status="eight_completed_training_bands_native_work_and_manuscript_cohort_audit",
        bands=len(rows), frame_records=sum(len(v) for v in identities.values()),
        distinct_raw_log_STRU_identities=len(identities),
        reused_identities=[dict(raw_log_sha256=k[0], STRU_sha256=k[1], occurrences=v)
                           for k, v in identities.items() if len(v) > 1],
        rows=rows, maximum_native_dual_work_defect_eV_cell=max_defect,
        maximum_absolute_peak_translation_meV_fu=max(abs(r["peak_triplet_parts_meV_fu"]["rigid_translation"]) for r in rows),
        new_DFT_calls=0, scheduler_mutations=0, holdout_or_G3_access=False,
        unpublished_reference_arrays_recreated=False, full_native_runtime_or_assets_reaudited_locally=False,
        physical_error_bound_or_TS_certified=False)


def audit_manuscript(text, result):
    section = text.split("### 3.4.", 1)[1].split("### 3.5.", 1)[0]
    required = ("eight completed clamped training bands", "72 frame records", "66 distinct raw-log/geometry")
    if (result["bands"], result["frame_records"], result["distinct_raw_log_STRU_identities"]) != (8, 72, 66):
        raise ValueError("actual cohort no longer matches the declared manuscript snapshot")
    if any(phrase not in section for phrase in required):
        raise ValueError("manuscript cohort coverage is stale")
    for row in result["rows"]:
        line = "| " + row["label"] + " | " + " | ".join(f"{v:.6f}" for v in row["table_meV_fu"]) + " |"
        if line not in section:
            raise ValueError("manuscript Table 3 differs from frozen native evidence")
        values = [row["peak_triplet_parts_meV_fu"][key] for key in
                  ("retained_atomic", "omitted_atomic", "open_cell")]
        values += [row["peak_optical_retained_meV_fu"], *row["peak_Cmma_retained_range_meV_fu"]]
        if any(f"{value:.6f}" not in section for value in values):
            raise ValueError("manuscript signed-work components omit a cohort member")
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("fresh cohort audit output required")
    result = audit_cohort()
    manuscript = ROOT / "paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md"
    result["manuscript_sha256"] = sha256(manuscript)
    result["manuscript_cohort_check"] = audit_manuscript(manuscript.read_text(encoding="utf-8"), result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
