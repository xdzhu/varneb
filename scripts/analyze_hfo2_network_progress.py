"""Dated G1 progress using the unchanged, source-bound physical analysis.

The old figure/replay entry point intentionally describes its original frozen
stage. This adapter replaces its historical one-source caption, not physical
fields, thresholds, source directions, geometry or the scientific G1 gate.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.analyze_hfo2_network_update import analyze_case
from scripts.audit_hfo2_static_replica import sha256


def analyze_progress(specification: Path) -> dict:
    result = analyze_case(Path(specification))
    historical = "only the ordinary T-PO residual passed; no final competing barrier or strain prediction verdict"
    if result["limitations"].count(historical) != 1:
        raise ValueError("review the source-bound historical caption before adapting its schema")
    passed = [c["name"] for c in result["channels"] if c["ordinary_residual_passed"]]
    result["physical_analysis_source_sha256"] = result["analysis_source_sha256"]
    result["analysis_source_sha256"] = sha256(Path(__file__))
    result["format_version"] = 2
    result["status"] = "G1_dated_progress_not_final_mechanisms"
    result["ordinary_residual_passed_channels"] = passed
    result["ordinary_residual_passed_sources"] = len(passed)
    result["unconverged_channels"] = [c["name"] for c in result["channels"]
                                      if not c["ordinary_residual_passed"]]
    result["historical_caption_replaced"] = historical
    result["limitations"] = [x for x in result["limitations"] if x != historical]
    result["limitations"].append(
        f"{len(passed)} of {len(result['channels'])} sources pass the ordinary residual criterion; "
        "this does not certify G1 phase/polarization/sampling gates, stationary TSs, "
        "a final channel ranking or independent predictions")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--specification", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("refusing an existing progress report")
    report = analyze_progress(args.specification)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"ordinary_passed": report["ordinary_residual_passed_channels"],
                      "unconverged": report["unconverged_channels"],
                      "G1_certified": report["full_G1_passed"], "new_DFT_calls": 0}))


if __name__ == "__main__":
    main()
