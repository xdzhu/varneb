"""Freeze registered B0/B1 controls before reading complete holdout path labels.

This reads only declared training reports and optional endpoint-only features.
The explicit label-access attestation is recorded, not independently proven.
Current free-cell/unconverged G1 observations cannot serve as G2 training.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

from examples.hfo2_fixed_input_factory import CONTRACT
from scripts.audit_hfo2_static_replica import sha256
from vcneb.prediction_controls import simple_prediction_controls


REQUIRED = {"PO_to_T": "decay", "PO_to_M": "decay",
            "PO_flip_T_pattern_preserving": "switching", "PO_flip_T_pattern_reversing": "switching"}


def freeze(lower_path, upper_path, output, *, endpoint_features_path=None,
           holdout_path_labels_unread=False, protocol_root=None):
    output = Path(output)
    if output.exists():
        raise FileExistsError("refusing an existing prediction namespace")
    if holdout_path_labels_unread is not True:
        raise ValueError("explicit caller attestation of unread holdout path labels required")
    paths = [Path(lower_path), Path(upper_path)]
    reports = [json.loads(p.read_text(encoding="utf-8")) for p in paths]
    for report, strain in zip(reports, (0., .01)):
        if (report.get("required_channels") != REQUIRED or report.get("formula_units") != 4
                or report.get("pressure_eV_A3") != 0.
                or report.get("mechanical_family") != "same_substrate_tilt_open"
                or report.get("physical_contract") != CONTRACT
                or report.get("mechanical_parameters", {}).get("biaxial_strain") != strain
                or report.get("mechanical_parameters", {}).get("external_field_V_A") != 0.
                or isinstance(report.get("mechanical_parameters", {}).get("external_field_V_A"), bool)
                or report.get("mechanical_parameters", {}).get("allow_tilt") is not True
                or report.get("fmax_target_eV_A") != .10):
            raise ValueError("registered Hf4O8 clamped 0/+1% training networks required; G1 free cell is not training")
    feature_path = Path(endpoint_features_path) if endpoint_features_path is not None else None
    features = json.loads(feature_path.read_text(encoding="utf-8")) if feature_path is not None else None
    result = simple_prediction_controls(*reports, parameter="biaxial_strain", target=.005, endpoint_features=features)
    root = Path(protocol_root) if protocol_root is not None else Path(__file__).resolve().parents[1]
    protocols = [root / "docs/HFO2_PREDICTION_PROTOCOL.md", root / "docs/HFO2_PREDICTION_PROTOCOL_V2_2026-10-09.md"]
    result.update({"frozen_at_UTC": datetime.now(timezone.utc).isoformat(),
                   "caller_attests_complete_holdout_path_labels_unread": True,
                   "independent_label_blinding_proven": False,
                   "training_source_sha256": [sha256(p) for p in paths],
                   "endpoint_feature_source_sha256": sha256(feature_path) if feature_path is not None else None,
                   "protocol_sha256": {p.name: sha256(p) for p in protocols},
                   "source_sha256": {"freezer": sha256(Path(__file__)),
                                      "controls": sha256(root / "vcneb/prediction_controls.py")},
                   "B2_to_B5_material_predictions_frozen": False,
                   "complete_forecast_batch_ready_for_holdout_paths": False})
    output.mkdir(parents=True, exist_ok=False)
    prediction = output / "prediction_controls.json"
    prediction.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    receipt = {"prediction_file": prediction.name, "prediction_sha256": sha256(prediction),
               "complete_forecast_batch_ready_for_holdout_paths": False,
               "new_DFT_calls": 0, "publish_before_labels": "commit this bundle and all B2-B5 predictions before reading complete holdout paths"}
    (output / "freeze_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--training-zero", type=Path, required=True)
    parser.add_argument("--training-one-percent", type=Path, required=True)
    parser.add_argument("--endpoint-features", type=Path)
    parser.add_argument("--holdout-path-labels-unread", action="store_true", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = freeze(args.training_zero, args.training_one_percent, args.output,
                    endpoint_features_path=args.endpoint_features,
                    holdout_path_labels_unread=args.holdout_path_labels_unread)
    print(json.dumps({"output": str(args.output), "selected": [r["selected_channel"] for r in result["predictions"]],
                      "new_DFT_calls": 0, "full_prediction_batch_ready": False}))


if __name__ == "__main__":
    main()
