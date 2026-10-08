"""Report true energies inside a sparse historical segment, not a new MEP."""

import argparse
import json
from pathlib import Path

import numpy as np

from scripts.audit_hfo2_static_replica import audited_results, sha256


def audit_point(root, index):
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    if manifest["purpose"] != "HfO2_guided_step300_gap01_02_three_statics":
        raise ValueError("unexpected gap manifest")
    record = manifest["points"][index]
    directory = root / "calculations" / f"{index:02d}"
    if record["index"] != index or any(sha256(directory / n) != h for n, h in record["input_sha256"].items()):
        raise ValueError("gap index/input contract changed")
    results = audited_results(directory)
    report = {"index": index, "fraction": record["fraction"], "status": "passed",
              "input_sha256": record["input_sha256"], "raw_log_sha256": sha256(directory / "OUT.ABACUS/running_scf.log"),
              "results": {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in results.items()},
              "relative_energy_meV_fu": (float(results["energy"]) - manifest["T_energy_eV_cell"]) * 1000 / 4}
    (directory / "point_audit.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--index", type=int)
    args = parser.parse_args()
    if args.index is not None:
        report = audit_point(args.root, args.index)
        print(json.dumps({k: report[k] for k in ("index", "status", "relative_energy_meV_fu")}))
    else:
        points = [audit_point(args.root, i) for i in range(3)]
        report = {"status": "three_gap_statics_audited", "points": points,
                  "manifest_sha256": sha256(args.root / "manifest.json"),
                  "highest_new_sample_relative_energy_meV_fu": max(p["relative_energy_meV_fu"] for p in points),
                  "limitations": "sampled straight reconstruction; does not determine the globally optimal MEP barrier"}
        (args.root / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps({k: v for k, v in report.items() if k != "points"}))


if __name__ == "__main__":
    main()
