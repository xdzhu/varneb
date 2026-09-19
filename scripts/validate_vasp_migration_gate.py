"""Compare converged, fixed-input VASP probes across compute environments."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


def compare_probes(reference, candidate, energy_tolerance=1e-4, force_tolerance=1e-4):
    if not all(math.isfinite(x) and x > 0 for x in (energy_tolerance, force_tolerance)):
        raise ValueError("comparison tolerances must be finite and positive")
    for report in (reference, candidate):
        if report.get("mode") != "converged_static_contract_probe" or report.get("all_passed") is not True:
            raise ValueError("migration requires successful converged static probes")
        if not report.get("records"):
            raise ValueError("migration probe has no records")
    if any(reference.get(key) != candidate.get(key) for key in ("isym", "symprec")):
        raise ValueError("fixed symmetry policies differ")

    def by_structure(report):
        records = {}
        for record in report["records"]:
            name = Path(record["structure"]).name
            if name in records:
                raise ValueError("duplicate probe structure")
            if record.get("status") != "converged_static_passed" or record.get("repeated_input_write_passed") is not True:
                raise ValueError("probe did not pass converged SCF and repeated input writing")
            for key in ("energy_eV", "maximum_force_eV_per_A"):
                if not math.isfinite(float(record[key])):
                    raise ValueError("nonfinite probe result")
            records[name] = record
        return records

    old, new = by_structure(reference), by_structure(candidate)
    if old.keys() != new.keys():
        raise ValueError("probe structures differ")
    comparisons = []
    for name in sorted(old):
        de = abs(float(new[name]["energy_eV"]) - float(old[name]["energy_eV"]))
        df = abs(float(new[name]["maximum_force_eV_per_A"]) - float(old[name]["maximum_force_eV_per_A"]))
        comparisons.append({"structure": name, "energy_difference_eV": de,
                            "maximum_atomic_force_difference_eV_per_A": df,
                            "passed": de <= energy_tolerance and df <= force_tolerance})
    return {"status": "ok" if all(x["passed"] for x in comparisons) else "failed",
            "scope": "static_energy_and_maximum_atomic_force_only_not_full_array_or_path_certification",
            "energy_tolerance_eV": energy_tolerance, "force_tolerance_eV_per_A": force_tolerance,
            "comparisons": comparisons}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    result = compare_probes(json.loads(args.reference.read_text()), json.loads(args.candidate.read_text()))
    result["input_sha256"] = {key: hashlib.sha256(path.read_bytes()).hexdigest()
                              for key, path in (("reference", args.reference), ("candidate", args.candidate))}
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
