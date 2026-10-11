# Freeze and verify a complete declared forecast panel

`vcneb.forecast_freeze` stores predictions already computed by a caller's
registered models. It is calculator-independent artifact plumbing, not a
prediction theory, material audit or label-blinding proof. It connects the
existing [local response controls](NESTED_RESPONSE_CONTROLS.md) and
[complete-panel scorer](PREDICTION_EVALUATION.md) without requiring test labels
to construct or validate the prediction panel.

## Inputs and exact byte identities

The registration JSON has exactly these fields:

```json
{
  "format_version": 1,
  "energy_unit": "eV/formula_unit",
  "comparison_contract_sha256": "<actual lowercase SHA256 of contract file>",
  "candidate_model": "<registered candidate>",
  "baseline_models": ["<every registered baseline>"],
  "case_ids": ["<every registered case>"]
}
```

Placeholders are explanatory, not executable hashes. The separate comparison
contract is a nonempty JSON object with the caller-audited physical/mechanical
definitions, normalization, information/cost policy and domain. Its bytes are
bound exactly; this module does **not** validate those scientific declarations.
Only `eV/formula_unit` is accepted: no implicit conversion from eV/cell occurs.

Each model artifact has exactly `format_version`, `energy_unit`, `model_id`,
`comparison_contract_sha256`, `model_definition`, `input_artifact_sha256` and
`forecasts`. The definition is a nonempty JSON object retaining the actual
chosen basis, coefficients, domain, calibration and failure rules. Each
forecast row has exactly `case_id`, `prediction_eV_fu`, `abstention_reason`.
A missing value requires an explicit reason; a negative signed value must
retain its abstention reason rather than be clipped to zero. Complete records
are required even if every model abstains. These definitions and outputs
are preserved, not recomputed or certified by the freezer.

`input_artifact_sha256` maps every declared visible-input identity to the
actual captured byte hash. All models bind the **same full visible input set**,
including unused preparations. Model-specific use/subsets remain in their
definitions. Input entries use `{"path": ..., "kind": ...}`; accepted kinds
are `training`, `visible_endpoint`, `source`, `protocol`,
`reference_definition`, `preparation_cost`. A path-label kind is refused.
The cost file remains an externally audited record, not a cost calculation
by this interface. Preparation is not made free by the freezer's zero calls.

Inputs are read once, validated before output creation, and copied from that
same captured snapshot. Source edits after capture cannot change the stored
bytes. An existing output namespace is refused, including a symlink; partial
output after an I/O failure is retained for investigation, not overwritten or
automatically deleted. Model/input identities are metadata, never filenames.
Only fixed ordinal paths are emitted inside the new bundle.

The bundle emits its own `.gitattributes`, also covered by the hash graph,
to preserve **all frozen bytes** through Git newline conversion and archives.
This rule applies only inside that new namespace, not to global Git settings
or other project files. The loader requires the fixed byte-preservation rule;
arbitrary Git attribute/filter rules are not admitted. A real local regression
stages CRLF inputs under `core.autocrlf=true` and exports with it false, then
checks the external manifest pin and every raw file hash. No remote Git or
DFT operation occurs in that temporary-repository test.

## Freeze first; pin externally; score only after authorized labels

```python
from vcneb.forecast_freeze import freeze_forecast_panel, load_frozen_forecast_panel

receipt = freeze_forecast_panel(
    "registration.json", "comparison_contract.json",
    model_artifact_paths={...},
    input_artifacts={...},
    output="new_forecast_namespace",
    holdout_path_labels_unread=True,
)
# Retain receipt["manifest_sha256"] outside this mutable directory, e.g.
# in the committed preregistration. Do not regenerate it from a later receipt.
forecasts = load_frozen_forecast_panel(
    "new_forecast_namespace",
    expected_manifest_sha256=externally_pinned_manifest_sha256,
)
```

The loader validates the external manifest pin, exact file-hash graph,
registration, complete rectangular case/model panel, source links, units,
signed abstentions and unchanged model artifact per model. It rejects changed
bytes, extra/missing records, unexpected paths, escaping symlinks and false
material-gate flags before returning `FrozenChannelForecast` records. It
does not execute model/source files. Rewriting both the manifest and its
local receipt cannot defeat the **original external pin**.

Use returned records with `evaluate_frozen_panel` only after independently
audited test references are authorized and available. Identical predictions
retain exactly zero shared-reference advantage. Snapshotting successful or
abstaining records does not establish TS identity, model validity, numerical
error, equal preparation cost or an improvement over strong controls.

## HfO₂ production gates remain separate

This is not a new HfO₂ study adapter or a complete material forecast. Its
receipt explicitly leaves `material_gates_certified`,
`independent_label_blinding_proven` and `authorizes_holdout_calculation` false.
It does not run G3, select a bottleneck/reference or launch/read a +0.5% chain.

For the actual study, first finish and audit the complete 0/+1% G2 matrix,
then pass the registered training-only G3, stationary/domain, reference,
sampling/error and budget gates. Freeze **all** v1/v2 controls: B0, both B1
arms, valid T controls, all four Cmma registrations and the registered
atomic/joint/branch arms, retaining failures and any missing required data.
The existing [B0/B1 material adapter](PREDICTION_CONTROLS.md) alone still does
not satisfy that full batch. A smaller declared panel accepted by the generic
API is not evidence that the HfO₂ protocol is fulfilled.

Commit the complete actual models, registration, inputs/cost audits and
external manifest pin before authorized complete holdout labels are accessed.
A timestamp/hash and caller attestation cannot prove that a person never
saw labels, detect labels misdeclared as training or concealed in arbitrary
metadata, or prove that output predictions were computed by the recorded
model. Those provenance/chronology checks remain external requirements.
No physical parameter, DFT budget, optimizer or ordinary threshold changes.

## Software checks, not material prediction accuracy

```sh
python -m pytest -q tests/test_forecast_freeze.py tests/test_prediction_evaluation.py
```

Fixtures use synthetic sources, eleven model IDs (including four distinct
Cmma-like labels), two artificial cases and zero DFT. They check roundtrip,
complete coverage, abstentions, input snapshots, malformed schemas, tampering,
external pinning, fixed paths and no manufactured harmonic advantage. The
symlink branch tests simulate the guard; they do not claim creation or testing
of native Windows symlinks. No HfO₂ training/holdout artifact enters these tests.
