# Executable prospective B0/B1 controls

The calculator-free `vcneb.prediction_controls.simple_prediction_controls`
implements the conventional controls already registered in
[HfO2 protocol v1](HFO2_PREDICTION_PROTOCOL.md). These controls are **not** a
new theory, an acceleration method or material prediction evidence by themselves.

## Training-only selection and one mechanical family

Both training inputs are audited outputs of `summarize_competing_paths`.
They require complete declared switching/decay coverage, ordinary residual
pass, sampling audit, measured directional-barrier error and provenance.
The API also checks the numeric residual, relative profile, reaction energy
and interval rather than trusting a top-level ready flag alone. Physical
input hashes, formula normalization, force target, pressure and boundary
family must match. Only the declared condition parameter may change.

For each role at the lower training condition, let U_min be the smallest
channel upper bound. Every channel with L_j<=U_min remains a possible lowest
candidate under the supplied bounds. Select the lexical-first name from this
set, as preregistered, and report all possible candidates. This is not proof
of a unique minimum. Keep this representative even if another edge is lowest
at the upper condition; record the difference rather than selecting on the
held-out result. The prediction covers two selected edges, not the minimum
over the entire network at the unseen condition.

Record the invariant T substrate **reference** in the mechanical definition
and the strain as its changing parameter. Actual strain-dependent cell vectors
remain in the raw geometry audit; they are not a second independently changed
parameter hidden under a fixed family name.

## Values, units and limits

For x0<x*<x1 and w=(x*-x0)/(x1-x0), all B below are in eV/f.u.:

```text
B0_j(x*) = (1-w) B_j(x0) + w B_j(x1)

delta_IS = [E_IS(x*) - E_IS(x0)] / formula_units
delta_FS = [E_FS(x*) - E_FS(x0)] / formula_units
B1_fixed_j  = B_j(x0) - delta_IS
B1_follow_j = B_j(x0) + delta_FS - delta_IS
```

B1_fixed assumes an unchanged absolute bottleneck energy; B1_follow assumes
that its change follows the final well. Both are zero-fit null models, not
assumptions that these behaviours physically occur. Endpoint energies are raw
eV per simulation cell. Their anchor values must match the chosen training
path and common initial well. Missing features remain explicitly unavailable,
not zero energy changes or failed performance.

This raw-E-only endpoint feature schema is restricted to **P=0**. At nonzero
pressure, shifts in E alone omit P times the endpoint volume change; audited
volumes/enthalpies would be required. Such feature use is rejected explicitly.
B0 without these features can interpolate the audited training H=E+PV barriers
at a matched finite pressure, but does not obtain an unseen error guarantee.

The B0 training interval is propagated through the same weighted expression.
It bounds the interpolation of supplied training-reference values **only**:
it does not bound curvature/nonlinearity or model error at x*. Accordingly,
`prediction_error_bound_eV_fu` remains null. No force threshold becomes an
energy uncertainty, confidence interval or prospective error guarantee.

When a raw prediction is negative, or below a visible positive target
reaction energy, it violates the endpoint-inclusive barrier definition.
Preserve the signed prediction and record an abstention; do not clip it to
zero or move it to an endpoint after seeing the test result. Abstention is
coverage/assumption failure, not a successful accurate prediction.

## Optional visible endpoint features

The feature document contains **only**:

- `physical_contract`, `mechanical_family`, `formula_units`, `pressure_eV_A3`;
- `anchor_parameters` and `target_parameters` matching the registered conditions;
- `initial` and `final_by_channel` endpoint records.

Each record contains `anchor_energy_eV_cell`, `target_energy_eV_cell`,
`anchor_audit_sha256`, `target_audit_sha256` and positive integer
`target_endpoint_DFT_calls`. Audit the actual endpoint geometry, phase,
atomic force and open traction before supplying these features. Hash strings
do not perform that external audit. Complete target path labels are not
accepted as fields or used as input to this API.

Sum the preparation calls for **all** provided endpoint records, including
unused final features; do not report only the selected records' lower cost.
Supply the same visible features and actual preparation cost to B0--B5 in
the eventual fair comparison. The API itself calls no calculator. Feature
availability is not a claim of zero-DFT prediction or independent label blinding.

## HfO2 freezing entry point

`scripts.freeze_hfo2_prediction_controls` additionally fixes HfO2's original
ABACUS input/orbital hashes,4f.u.,P=0,E=0, tilt-open common substrate,0/+1%
training conditions,0.10 ordinary target and +0.5% unseen condition.
It refuses free-cell G1 reports, incomplete training and existing output
namespaces. The numeric module remains backend-independent; this adapter
applies the current study's preregistered material contract.

After **actual G2 training and audit**, not now, run from the repository root:

```powershell
python -m scripts.freeze_hfo2_prediction_controls --training-zero /actual/audited/strain0_network.json --training-one-percent /actual/audited/strain1_network.json --endpoint-features /actual/audited/visible_endpoint_features.json --holdout-path-labels-unread --output /new/prediction_namespace
```

Omit `--endpoint-features` if they are not yet available, and retain B1 as
unavailable. The fresh bundle hashes training inputs, optional features,
source code, both protocol versions and the prediction file. Commit it
before complete holdout path labels are read. The unread-label flag records
the caller's attestation; a checksum/timestamp cannot prove that no person
or other process has previously seen labels.

This bundle freezes B0/B1 **only**. It does not freeze B2--B5 or authorize
holdout production/reading. The complete registered forecast batch remains
required. No target error, model advantage, H1 window, TS or JCTC readiness
is evaluated by this entry point.

The dated E025 [G1 evidence](../benchmarks/hfo2_channels/20261008/network_update_20261009/README.md)
is deliberately rejected: it is free-cell and lacks three ordinary passes
and measured sampling/barrier-error evidence. The
[E026 delivery](../benchmarks/hfo2_channels/20261008/prediction_controls/README.md)
contains the actual negative-gate evidence and synthetic algorithm tests,
not invented G2 training or material forecasts.

The later [fixed-control stationary-branch prerequisite](CONTROLLED_STATIONARY_BRANCH.md)
retains affine offsets and distinguishes a minimum from a first-order saddle
inside the actual mechanical internal space. It does not freeze B2--B5,
choose material branches, authorize holdout or turn the unconverged G1
observations into predictions. Actual same-contract training and all
registered uncertainty/stability/coverage gates remain necessary.

The [paired stationary-gap step](STATIONARY_GAP_RESPONSE.md) computes the
initial-versus-bottleneck response with explicit reference corrections and
per-centre scale factors. It preserves frozen-space residuals and does not
replace actual G2 training, B2--B5 selection/freezing or independent labels.
