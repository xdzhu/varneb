"""Reparse the three existing training bands; no SCF or scheduler call."""
from pathlib import Path
import json

import numpy as np

from scripts.analyze_hfo2_clamped_reference import terminal_images

ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / "benchmarks/hfo2_channels/20261008"
CASES = [
    ("clamped_M_terminal_E068_20261011", "56b92e8b15195d6f619a8402c7eee83a774b2e364e012478faa3fdc5d2da431c"),
    ("clamped_T_terminal_E079_20261011", "c5653b4759b202199cb75739f99038aaa52b3ccc220ef226efbacc980e5d716e"),
    ("clamped_reversing_terminal_E083_20261011", "632e0e9b2b821b3bfc7d66b361c44c7db734988fb5d018ea2f31441fab14ed41"),
]


def verify():
    first, first_raw = None, None
    rows = []
    for case, digest in CASES:
        images, observation, _, receipt = terminal_images(BASE / case, digest)
        image = images[0]
        raw = observation["raw_image_evaluations"][0]
        if first is None:
            first, first_raw = image, raw
        if (raw["input_sha256"] != first_raw["input_sha256"]
                or raw["raw_log_sha256"] != first_raw["raw_log_sha256"]
                or image.get_potential_energy() != first.get_potential_energy()
                or image.get_chemical_symbols() != first.get_chemical_symbols()
                or not np.array_equal(image.get_forces(), first.get_forces())
                or not np.array_equal(image.get_stress(), first.get_stress())):
            raise ValueError("shared PO native SCF bytes/EFS differ")
        dp = float(np.max(np.abs(image.positions - first.positions)))
        dh = float(np.max(np.abs(image.cell.array - first.cell.array)))
        if max(dp, dh) > 1e-14:
            raise ValueError("ordered initial geometry differs beyond serialization")
        rows.append(dict(case=case, terminal_audit_sha256=digest, job=receipt["job_id"],
                         initial_native_input_and_log_hashes_equal=True,
                         initial_EFS_exact=True, ordered_position_max_difference_A=dp,
                         cell_max_difference_A=dh))
    return dict(scope="same ordered cached PO initial point only; no material prediction",
                reports=rows, native_unique_terminal_frames_reparsed=27,
                shared_input_sha256=first_raw["input_sha256"],
                shared_log_sha256=first_raw["raw_log_sha256"],
                physical_geometry_repaired=False, new_DFT_calls=0, scheduler_mutations=0)


if __name__ == "__main__":
    result = verify()
    output = Path(__file__).with_name("shared_initial.json")
    with output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, allow_nan=False); stream.write("\n")
    print(json.dumps(result))
