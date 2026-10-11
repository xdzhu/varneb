"""Synthetic storage/provenance checks, never HfO2 prediction evidence."""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile

import pytest

from vcneb.forecast_freeze import freeze_forecast_panel, load_frozen_forecast_panel
from vcneb.prediction_evaluation import AuditedBarrierReference, evaluate_frozen_panel


MODELS = ["B5", "B0", "B1_fixed", "B1_follow", "B2_T",
          *(f"B2_Cmma_gauge{i}" for i in range(1, 5)), "B3", "B4"]
CASES = ["switch", "escape"]


def put(path, value):
    path.write_text(json.dumps(value, sort_keys=True, allow_nan=False), encoding="utf-8")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def inputs(tmp_path):
    source = tmp_path / "source"
    source.mkdir()
    contract = source / "contract.json"
    put(contract, {"synthetic_only": True, "formula_units": 4})
    registration = source / "registration.json"
    put(registration, dict(format_version=1, energy_unit="eV/formula_unit",
        comparison_contract_sha256=sha(contract), candidate_model="B5",
        baseline_models=MODELS[1:], case_ids=CASES))
    artifacts = {}
    for kind in ("training", "source", "protocol", "preparation_cost"):
        path = source / (kind + ".json")
        put(path, {"kind": kind, "synthetic_only": True, "DFT_calls": 0})
        artifacts[kind] = {"path": path, "kind": kind}
    hashes = {name: sha(record["path"]) for name, record in artifacts.items()}
    models = {}
    for i, model in enumerate(MODELS):
        path = source / f"model{i}.json"
        put(path, dict(format_version=1, energy_unit="eV/formula_unit", model_id=model,
            comparison_contract_sha256=sha(contract),
            model_definition={"kind": "synthetic constant; no material model fitted"},
            input_artifact_sha256=hashes,
            forecasts=[dict(case_id=case, prediction_eV_fu=value, abstention_reason=None)
                       for case, value in zip(CASES, (.2, .4))]))
        models[model] = path
    return dict(registration=registration, contract=contract, models=models,
                artifacts=artifacts, output=tmp_path / "frozen")


def freeze(inputs, **kwargs):
    return freeze_forecast_panel(inputs["registration"], inputs["contract"],
        inputs["models"], inputs["artifacts"], inputs["output"],
        holdout_path_labels_unread=kwargs.get("attest", True))


def edit(path, change):
    document = json.loads(path.read_text(encoding="utf-8"))
    change(document)
    put(path, document)


def test_complete_panel_roundtrip_preserves_four_gauges_and_harmonic_null(inputs):
    receipt = freeze(inputs)
    forecasts = load_frozen_forecast_panel(inputs["output"],
        expected_manifest_sha256=receipt["manifest_sha256"])
    assert len(forecasts) == 22
    assert {f.model_id for f in forecasts} == set(MODELS)
    assert [(f.model_id, f.case_id) for f in forecasts] == [(m, c) for m in MODELS for c in CASES]
    assert all(len({f.frozen_model_sha256 for f in forecasts if f.model_id == model}) == 1 for model in MODELS)
    references = [AuditedBarrierReference(case, value, .005, sha(inputs["contract"]), "a"*64)
                  for case, value in zip(CASES, (.2, .4))]
    scored = evaluate_frozen_panel(forecasts, references, candidate_model="B5",
                                  baseline_models=MODELS[1:])
    assert not scored["full_panel_strict_advantage"]
    assert all(row["best_available_baseline_oracle_gain_interval_eV_fu"] == (0, 0)
               for row in scored["comparisons"])
    assert receipt["new_DFT_calls"] == 0
    assert not receipt["independent_label_blinding_proven"]
    assert not receipt["material_gates_certified"]
    assert not receipt["authorizes_holdout_calculation"]


def test_signed_abstention_survives_freeze_load_and_scoring(inputs):
    def change(model):
        model["forecasts"][1].update(prediction_eV_fu=-.02, abstention_reason="synthetic branch failure")
    edit(inputs["models"]["B5"], change)
    receipt = freeze(inputs)
    forecasts = load_frozen_forecast_panel(inputs["output"], expected_manifest_sha256=receipt["manifest_sha256"])
    failed = next(f for f in forecasts if f.model_id == "B5" and f.case_id == "escape")
    assert failed.prediction_eV_fu == -.02 and failed.abstention_reason == "synthetic branch failure"


@pytest.mark.parametrize("fault", ["missing_model", "extra_model", "missing_case", "duplicate_case",
    "extra_case", "unit", "contract", "wrong_model", "missing_definition", "hidden_label_field",
    "changed_input", "missing_input_hash", "negative_without_reason", "missing_without_reason"])
def test_invalid_panel_rejected_before_namespace_creation(inputs, fault):
    path = inputs["models"]["B5"]
    def change(model):
        if fault == "missing_case": model["forecasts"].pop()
        elif fault == "duplicate_case": model["forecasts"][1] = model["forecasts"][0]
        elif fault == "extra_case": model["forecasts"][1]["case_id"] = "posthoc"
        elif fault == "unit": model["energy_unit"] = "eV/cell"
        elif fault == "contract": model["comparison_contract_sha256"] = "a"*64
        elif fault == "wrong_model": model["model_id"] = "B0"
        elif fault == "missing_definition": model["model_definition"] = {}
        elif fault == "hidden_label_field": model["heldout_barrier"] = .2
        elif fault == "missing_input_hash": model["input_artifact_sha256"].pop("protocol")
        elif fault == "negative_without_reason": model["forecasts"][0]["prediction_eV_fu"] = -.1
        elif fault == "missing_without_reason": model["forecasts"][0]["prediction_eV_fu"] = None
    if fault == "missing_model": inputs["models"].pop("B2_Cmma_gauge4")
    elif fault == "extra_model": inputs["models"]["posthoc"] = path
    elif fault == "changed_input": put(inputs["artifacts"]["training"]["path"], {"changed": True})
    else: edit(path, change)
    with pytest.raises(ValueError): freeze(inputs)
    assert not inputs["output"].exists()


@pytest.mark.parametrize("fault", ["attestation", "duplicate_model", "duplicate_case", "candidate_baseline", "version_bool", "extra_label"])
def test_invalid_registration_or_attestation_rejected(inputs, fault):
    def change(registration):
        if fault == "duplicate_model": registration["baseline_models"].append("B0")
        elif fault == "duplicate_case": registration["case_ids"].append("escape")
        elif fault == "candidate_baseline": registration["baseline_models"].append("B5")
        elif fault == "version_bool": registration["format_version"] = True
        elif fault == "extra_label": registration["holdout_reference"] = .2
    edit(inputs["registration"], change)
    with pytest.raises(ValueError): freeze(inputs, attest=fault != "attestation")
    assert not inputs["output"].exists()


def test_declared_path_label_input_kind_is_rejected_without_reading_it(inputs, monkeypatch):
    label = inputs["artifacts"]["training"]["path"]
    inputs["artifacts"]["training"]["kind"] = "holdout_path_labels"
    original = Path.read_bytes
    def checked_read(path):
        assert path != label, "label input must be rejected before reading"
        return original(path)
    monkeypatch.setattr(Path, "read_bytes", checked_read)
    with pytest.raises(ValueError, match="held-out path labels"):
        freeze(inputs)
    assert not inputs["output"].exists()


@pytest.mark.parametrize("file", [".gitattributes", "registration.json", "comparison_contract.json", "models/0000.json", "inputs/0000"])
def test_raw_artifact_tampering_fails_before_loading_score_records(inputs, file):
    receipt = freeze(inputs)
    path = inputs["output"] / file
    path.write_bytes(path.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="bytes changed"):
        load_frozen_forecast_panel(inputs["output"], expected_manifest_sha256=receipt["manifest_sha256"])


def test_rewriting_manifest_and_receipt_does_not_defeat_external_pin(inputs):
    receipt = freeze(inputs)
    manifest_path = inputs["output"] / "manifest.json"
    edit(manifest_path, lambda manifest: manifest.update(frozen_at_UTC="rewritten"))
    edit(inputs["output"] / "freeze_receipt.json", lambda record: record.update(manifest_sha256=sha(manifest_path)))
    with pytest.raises(ValueError, match="externally frozen hash"):
        load_frozen_forecast_panel(inputs["output"], expected_manifest_sha256=receipt["manifest_sha256"])


@pytest.mark.parametrize("fault", ["traversal", "extra_file_hash", "claim_material_gate", "forecast_count"])
def test_malformed_manifest_is_not_a_material_gate_even_if_caller_pins_its_hash(inputs, fault):
    freeze(inputs)
    path = inputs["output"] / "manifest.json"
    def change(manifest):
        if fault == "traversal": manifest["model_files"]["B5"] = "../../outside.json"
        elif fault == "extra_file_hash": manifest["file_sha256"]["../../outside.json"] = "a"*64
        elif fault == "claim_material_gate": manifest["material_gates_certified"] = True
        elif fault == "forecast_count": manifest["n_forecasts"] = 1
    edit(path, change)
    with pytest.raises(ValueError): load_frozen_forecast_panel(inputs["output"], expected_manifest_sha256=sha(path))


def test_original_input_edit_after_freeze_does_not_mutate_snapshot(inputs):
    receipt = freeze(inputs)
    original_model_sha = sha(inputs["models"]["B5"])
    put(inputs["models"]["B5"], {"later": "changed"})
    forecasts = load_frozen_forecast_panel(inputs["output"], expected_manifest_sha256=receipt["manifest_sha256"])
    assert next(f for f in forecasts if f.model_id == "B5").frozen_model_sha256 == original_model_sha
    assert next(f for f in forecasts if f.model_id == "B5").prediction_eV_fu == .2


def test_existing_namespace_not_overwritten(inputs):
    receipt = freeze(inputs)
    before = (inputs["output"] / "manifest.json").read_bytes()
    with pytest.raises(FileExistsError): freeze(inputs)
    assert (inputs["output"] / "manifest.json").read_bytes() == before
    assert sha(inputs["output"] / "manifest.json") == receipt["manifest_sha256"]


@pytest.mark.parametrize("payload", [b'{"format_version":1,"format_version":1}', b'{"x":NaN}'])
def test_duplicate_or_nonfinite_JSON_cannot_hide_in_artifacts(inputs, payload):
    inputs["registration"].write_bytes(payload)
    with pytest.raises(ValueError): freeze(inputs)
    assert not inputs["output"].exists()


def test_exponent_overflow_in_model_metadata_is_rejected(inputs):
    path = inputs["models"]["B5"]
    edit(path, lambda model: model["model_definition"].update(coefficient="OVERFLOW"))
    path.write_bytes(path.read_bytes().replace(b'"OVERFLOW"', b'1e999'))
    with pytest.raises(ValueError, match="Out of range float"):
        freeze(inputs)
    assert not inputs["output"].exists()


@pytest.mark.parametrize("kind", [None, True, ["training"], {"kind": "training"}])
def test_malformed_input_kind_is_an_explicit_validation_error(inputs, kind):
    inputs["artifacts"]["training"]["kind"] = kind
    with pytest.raises(ValueError, match="input kind"):
        freeze(inputs)
    assert not inputs["output"].exists()


def test_snapshot_is_read_once_and_not_reread_while_writing(inputs, monkeypatch):
    source_paths = {inputs["registration"], inputs["contract"], *inputs["models"].values(),
                    *(record["path"] for record in inputs["artifacts"].values())}
    original_bytes = {path: path.read_bytes() for path in source_paths}
    reads = {path: 0 for path in source_paths}
    original_read = Path.read_bytes
    def changing_read(path):
        data = original_read(path)
        if path in reads:
            reads[path] += 1
            # A subsequent reread would see different bytes; the captured
            # original snapshot, not those changed bytes, must be stored.
            path.write_bytes(b"changed after snapshot")
        return data
    monkeypatch.setattr(Path, "read_bytes", changing_read)
    receipt = freeze(inputs)
    forecasts = load_frozen_forecast_panel(inputs["output"], expected_manifest_sha256=receipt["manifest_sha256"])
    assert all(count == 1 for count in reads.values())
    assert len(forecasts) == 22
    manifest = json.loads((inputs["output"] / "manifest.json").read_text(encoding="utf-8"))
    for name, path in inputs["models"].items():
        assert (inputs["output"] / manifest["model_files"][name]).read_bytes() == original_bytes[path]


@pytest.mark.parametrize("relative", [None, "models/0000.json"])
def test_symlink_guard_branches_reject_before_following_target(inputs, monkeypatch, relative):
    receipt = freeze(inputs)
    unsafe = inputs["output"] if relative is None else inputs["output"] / relative
    original = Path.is_symlink
    monkeypatch.setattr(Path, "is_symlink", lambda path: path == unsafe or original(path))
    # This is a synthetic branch check, not creation of a native symlink.
    with pytest.raises(ValueError):
        load_frozen_forecast_panel(inputs["output"], expected_manifest_sha256=receipt["manifest_sha256"])


def test_identifiers_are_metadata_not_filesystem_paths(inputs):
    alias = "../B5-metadata-only"
    model_path = inputs["models"].pop("B5")
    inputs["models"][alias] = model_path
    edit(model_path, lambda model: model.update(model_id=alias))
    edit(inputs["registration"], lambda registration: registration.update(candidate_model=alias))
    receipt = freeze(inputs)
    forecasts = load_frozen_forecast_panel(inputs["output"], expected_manifest_sha256=receipt["manifest_sha256"])
    assert forecasts[0].model_id == alias
    assert not (inputs["output"].parent / "B5-metadata-only").exists()


def test_git_cross_platform_archive_preserves_frozen_raw_bytes(inputs, tmp_path):
    git = shutil.which("git")
    if git is None:
        pytest.skip("Git is required only for the raw-byte transport regression")
    training = inputs["artifacts"]["training"]["path"]
    training.write_bytes(b"synthetic training row1\r\nsynthetic row2\r\n")
    for path in inputs["models"].values():
        edit(path, lambda model: model["input_artifact_sha256"].update(training=sha(training)))
    receipt = freeze(inputs)
    manifest = json.loads((inputs["output"] / "manifest.json").read_text(encoding="utf-8"))
    empty_config = tmp_path / "empty-git-config"
    empty_config.write_bytes(b"")
    env = {key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
    env.update(GIT_CONFIG_NOSYSTEM="1", GIT_CONFIG_GLOBAL=str(empty_config))
    command = [git, "-C", str(inputs["output"]), "-c", "core.autocrlf=true",
               "-c", "core.safecrlf=false", "-c", "core.attributesFile=" + str(empty_config)]
    def run(*args):
        return subprocess.run([*command, *args], env=env, check=True,
                              capture_output=True, text=True).stdout.strip()
    run("init", "--quiet")
    run("add", ".")
    tree = run("write-tree")
    archive = tmp_path / "raw-byte-transport.tar"
    # Stage under Windows-like CRLF conversion, export with Linux-like LF
    # settings. A same-platform smudge can hide normalization in the index.
    run("-c", "core.autocrlf=false", "archive", "--format=tar", "--output=" + str(archive), tree)
    with tarfile.open(archive) as stream:
        with stream.extractfile("manifest.json") as manifest_member:
            assert hashlib.sha256(manifest_member.read()).hexdigest() == receipt["manifest_sha256"]
        for relative, expected in manifest["file_sha256"].items():
            member = stream.extractfile(relative)
            assert member is not None
            with member:
                assert hashlib.sha256(member.read()).hexdigest() == expected, relative


def test_arbitrary_git_attribute_rules_are_not_admitted_by_a_new_manifest_pin(inputs):
    freeze(inputs)
    attribute = inputs["output"] / ".gitattributes"
    attribute.write_bytes(b"* filter=not-a-frozen-byte-contract\n")
    manifest = inputs["output"] / "manifest.json"
    edit(manifest, lambda data: data["file_sha256"].update({".gitattributes": sha(attribute)}))
    with pytest.raises(ValueError, match="byte-preservation contract"):
        load_frozen_forecast_panel(inputs["output"], expected_manifest_sha256=sha(manifest))
