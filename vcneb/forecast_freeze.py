"""Byte-bound, calculator-free storage of a complete registered forecast panel.

This closes an artifact boundary, not a material/model/stability gate. It
does not fit models, audit DFT, prove label blinding or authorize holdouts.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Mapping

from .prediction_evaluation import FrozenChannelForecast


_UNIT = "eV/formula_unit"
_BUNDLE_ATTRIBUTES = (b"# Preserve every frozen byte across Git checkouts and archives.\n"
                      b"* -text -eol whitespace=-blank-at-eol,-blank-at-eof,cr-at-eol\n")
_INPUT_KINDS = {"training", "visible_endpoint", "source", "protocol",
                "reference_definition", "preparation_cost"}
_REGISTRATION_KEYS = {"format_version", "energy_unit", "comparison_contract_sha256",
                      "candidate_model", "baseline_models", "case_ids"}
_MODEL_KEYS = {"format_version", "energy_unit", "model_id", "comparison_contract_sha256",
               "model_definition", "input_artifact_sha256", "forecasts"}
_MANIFEST_KEYS = {"format_version", "frozen_at_UTC", "comparison_contract_sha256",
                  "model_files", "input_files", "file_sha256", "n_forecasts",
                  "caller_attests_complete_holdout_path_labels_unread",
                  "independent_label_blinding_proven", "material_gates_certified",
                  "authorizes_holdout_calculation", "new_DFT_calls"}


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _hash(value) -> str:
    if (not isinstance(value, str) or len(value) != 64
            or any(c not in "0123456789abcdef" for c in value)):
        raise ValueError("lowercase SHA256 required")
    return value


def _name(value) -> str:
    if not isinstance(value, str) or not value.strip() or value != value.strip():
        raise ValueError("nonempty trimmed artifact/model/case identity required")
    return value


def _names(values) -> tuple[str, ...]:
    if not isinstance(values, list) or not values:
        raise ValueError("nonempty list of registered identities required")
    result = tuple(_name(value) for value in values)
    if len(set(result)) != len(result):
        raise ValueError("duplicate registered identity")
    return result


def _keys(value, expected) -> None:
    if not isinstance(value, dict) or set(value) != expected:
        raise ValueError("exact registered schema required; no extra label fields")


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate JSON key")
        result[key] = value
    return result


def _nonfinite(value):
    raise ValueError(f"nonfinite JSON constant rejected: {value}")


def _json(data: bytes):
    document = json.loads(data.decode("utf-8"), object_pairs_hook=_pairs,
                          parse_constant=_nonfinite)
    # Exponent overflow (e.g. 1e999) is not passed to parse_constant.
    json.dumps(document, allow_nan=False)
    return document


def _encode(value) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")


def _version(document) -> None:
    if type(document.get("format_version")) is not int or document["format_version"] != 1:
        raise ValueError("format_version=1 required")
    if document.get("energy_unit") != _UNIT:
        raise ValueError("explicit eV/formula_unit required; no implicit unit conversion")


def _assemble(registration_data, contract_data, model_data, input_data):
    registration = _json(registration_data)
    _keys(registration, _REGISTRATION_KEYS)
    _version(registration)
    contract = _json(contract_data)
    if not isinstance(contract, dict) or not contract:
        raise ValueError("nonempty comparison-contract document required")
    contract_sha = _hash(registration["comparison_contract_sha256"])
    if contract_sha != _sha(contract_data):
        raise ValueError("comparison-contract bytes do not match registration")
    candidate = _name(registration["candidate_model"])
    baselines = _names(registration["baseline_models"])
    cases = _names(registration["case_ids"])
    if candidate in baselines:
        raise ValueError("candidate cannot be its own baseline")
    models = (candidate, *baselines)
    if set(model_data) != set(models):
        raise ValueError("complete registered model set required")
    input_sha = {name: _sha(data) for name, data in input_data.items()}
    forecasts = []
    for model in models:
        document = _json(model_data[model])
        _keys(document, _MODEL_KEYS)
        _version(document)
        if document["model_id"] != model or document["comparison_contract_sha256"] != contract_sha:
            raise ValueError("model identity/comparison contract mismatch")
        if not isinstance(document["model_definition"], dict) or not document["model_definition"]:
            raise ValueError("explicit frozen model definition required")
        if document["input_artifact_sha256"] != input_sha:
            raise ValueError("all models must bind the same declared visible input bytes")
        rows = document["forecasts"]
        if not isinstance(rows, list) or len(rows) != len(cases):
            raise ValueError("complete case records required, including abstentions")
        by_case = {}
        for row in rows:
            _keys(row, {"case_id", "prediction_eV_fu", "abstention_reason"})
            case = _name(row["case_id"])
            if case not in cases or case in by_case:
                raise ValueError("extra/duplicate forecast case")
            by_case[case] = FrozenChannelForecast(
                model, case, row["prediction_eV_fu"], row["abstention_reason"],
                contract_sha, _sha(model_data[model]),
            )
        forecasts.extend(by_case[case] for case in cases)
    return registration, tuple(forecasts)


def freeze_forecast_panel(registration_path, comparison_contract_path,
                          model_artifact_paths: Mapping, input_artifacts: Mapping,
                          output, *, holdout_path_labels_unread=False):
    """Snapshot a declared full panel, without reading a reference-label API.

    Each input identity maps to ``{"path": ..., "kind": ...}``; permitted
    kinds describe training/visible features/source/protocol/reference/cost.
    Each model binds every visible input SHA and contains every registered
    case (including abstentions). Source bytes are read once, validated and
    copied from that same snapshot. An external study audit must certify
    physical data, budgets, strong controls and prospective registration.
    """
    destination = Path(output)
    if destination.exists() or destination.is_symlink():
        raise FileExistsError("refusing an existing forecast namespace")
    if holdout_path_labels_unread is not True:
        raise ValueError("explicit caller unread-label attestation required")
    if not isinstance(model_artifact_paths, Mapping) or not model_artifact_paths:
        raise ValueError("explicit model artifact mapping required")
    if not isinstance(input_artifacts, Mapping) or not input_artifacts:
        raise ValueError("explicit shared input artifact mapping required")
    models = {_name(name): Path(path).read_bytes() for name, path in model_artifact_paths.items()}
    inputs, kinds = {}, {}
    for name, record in input_artifacts.items():
        _name(name)
        _keys(record, {"path", "kind"})
        if not isinstance(record["kind"], str) or record["kind"] not in _INPUT_KINDS:
            raise ValueError("registered input kind required; held-out path labels are not input features")
        inputs[name] = Path(record["path"]).read_bytes()
        kinds[name] = record["kind"]
    registration_data = Path(registration_path).read_bytes()
    contract_data = Path(comparison_contract_path).read_bytes()
    registration, forecasts = _assemble(registration_data, contract_data, models, inputs)
    model_files = {name: f"models/{i:04d}.json" for i, name in enumerate(sorted(models))}
    input_files = {name: {"path": f"inputs/{i:04d}", "kind": kinds[name]}
                   for i, name in enumerate(sorted(inputs))}
    files = {".gitattributes": _BUNDLE_ATTRIBUTES,
             "registration.json": registration_data, "comparison_contract.json": contract_data}
    files.update({model_files[name]: data for name, data in models.items()})
    files.update({input_files[name]["path"]: data for name, data in inputs.items()})
    manifest = dict(format_version=1, frozen_at_UTC=datetime.now(timezone.utc).isoformat(),
                    comparison_contract_sha256=registration["comparison_contract_sha256"],
                    model_files=model_files, input_files=input_files,
                    file_sha256={name: _sha(data) for name, data in files.items()},
                    n_forecasts=len(forecasts),
                    caller_attests_complete_holdout_path_labels_unread=True,
                    independent_label_blinding_proven=False, material_gates_certified=False,
                    authorizes_holdout_calculation=False, new_DFT_calls=0)
    # Validation precedes namespace creation. Never clean up/overwrite a partial bundle.
    destination.mkdir(parents=True, exist_ok=False)
    for name, data in files.items():
        path = destination / name
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(data)
    manifest_data = _encode(manifest)
    with (destination / "manifest.json").open("xb") as stream:
        stream.write(manifest_data)
    receipt = dict(manifest_sha256=_sha(manifest_data), n_forecasts=len(forecasts),
                   authorizes_holdout_calculation=False, material_gates_certified=False,
                   independent_label_blinding_proven=False, new_DFT_calls=0)
    with (destination / "freeze_receipt.json").open("xb") as stream:
        stream.write(_encode(receipt))
    return receipt


def _bundle_read(root: Path, relative: str) -> bytes:
    path = root / relative
    if path.is_symlink() or not path.resolve().is_relative_to(root.resolve()):
        raise ValueError("bundle artifact must not escape through a symlink")
    return path.read_bytes()


def load_frozen_forecast_panel(bundle, *, expected_manifest_sha256: str):
    """Verify against an EXTERNALLY pinned manifest, then return score records.

    Do not derive the expected hash from this bundle's mutable receipt.
    Hashes attest bytes, not chronological blindness or physical correctness.
    Only manifest-listed fixed relative paths are read, never executable code.
    """
    _hash(expected_manifest_sha256)
    root = Path(bundle)
    if root.is_symlink() or not root.is_dir():
        raise ValueError("ordinary frozen-bundle directory required")
    manifest_data = _bundle_read(root, "manifest.json")
    if _sha(manifest_data) != expected_manifest_sha256:
        raise ValueError("manifest differs from externally frozen hash")
    manifest = _json(manifest_data)
    _keys(manifest, _MANIFEST_KEYS)
    if (type(manifest["format_version"]) is not int or manifest["format_version"] != 1
            or manifest["caller_attests_complete_holdout_path_labels_unread"] is not True
            or any(manifest[key] is not False for key in ("independent_label_blinding_proven",
                                                         "material_gates_certified",
                                                         "authorizes_holdout_calculation"))
            or type(manifest["new_DFT_calls"]) is not int or manifest["new_DFT_calls"] != 0):
        raise ValueError("invalid freeze status; artifact storage is not a material gate")
    _name(manifest["frozen_at_UTC"])
    model_files, input_files = manifest["model_files"], manifest["input_files"]
    if not isinstance(model_files, dict) or not model_files or not isinstance(input_files, dict) or not input_files:
        raise ValueError("nonempty frozen model/input mappings required")
    expected_models = {_name(name): f"models/{i:04d}.json"
                       for i, name in enumerate(sorted(model_files))}
    if model_files != expected_models:
        raise ValueError("fixed model paths required; no path traversal")
    for i, name in enumerate(sorted(input_files)):
        _name(name)
        _keys(input_files[name], {"path", "kind"})
        if (input_files[name]["path"] != f"inputs/{i:04d}"
                or not isinstance(input_files[name]["kind"], str)
                or input_files[name]["kind"] not in _INPUT_KINDS):
            raise ValueError("fixed paths/registered input kinds required")
    required_files = {".gitattributes", "registration.json", "comparison_contract.json", *model_files.values(),
                      *(record["path"] for record in input_files.values())}
    if not isinstance(manifest["file_sha256"], dict) or set(manifest["file_sha256"]) != required_files:
        raise ValueError("complete exact file-hash graph required")
    data = {}
    for name in sorted(required_files):
        data[name] = _bundle_read(root, name)
        if _sha(data[name]) != _hash(manifest["file_sha256"][name]):
            raise ValueError(f"frozen artifact bytes changed: {name}")
    if data[".gitattributes"] != _BUNDLE_ATTRIBUTES:
        raise ValueError("fixed byte-preservation contract required")
    registration, forecasts = _assemble(
        data["registration.json"], data["comparison_contract.json"],
        {name: data[path] for name, path in model_files.items()},
        {name: data[record["path"]] for name, record in input_files.items()},
    )
    if (manifest["comparison_contract_sha256"] != registration["comparison_contract_sha256"]
            or type(manifest["n_forecasts"]) is not int or manifest["n_forecasts"] != len(forecasts)):
        raise ValueError("manifest contradicts registered forecast panel")
    return forecasts
