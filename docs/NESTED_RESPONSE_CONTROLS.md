# Nested, measured-column stationary-response controls

This implementation makes the local B2--B5 **algebra and assumptions**
comparable. It does not yet freeze or validate the registered HfO2 forecasts.
Schur reduction and harmonic continuation are standard mathematics, not our
novelty claim. The prospective scientific question remains whether measured
mode--strain release improves *unseen channel-response predictions* over
strong, fairly resourced references.

## One dataset; explicit missing information

`vcneb.MeasuredStationaryCentre` records an orthonormal measured basis M,
directional physical-gradient derivatives D=H M, the full physical gradient,
raw potential, retained coordinates Q, admissible internal space B, atomic
space A, prescribed control direction c, and source/measurement identities.
All coordinates have one explicit Angstrom metric; H is in eV/Angstrom².
At nonzero pressure use derivatives of H=E+PV consistently, not E alone.
Never differentiate NEB/spring forces or substitute a parent-phase Gamma
matrix for curvature at a different well or bottleneck.

Q and A must lie within B, c must be orthogonal to B, and M must lie within
the actual internal-plus-control chart. Complex modes need a declared real
representation before this interface; imaginary components are not dropped.
Measured reciprocity is checked against the caller's declared tolerance.
These declarations and hashes are not independent audits of physical data.

The symmetric action extension is

```text
K = sym(M.T D)
D_corrected = D + M (K - M.T D)
H_extension = D_corrected M.T + M D_corrected.T - M K M.T
```

The unmeasured M-orthogonal block is **unknown**, not physically zero. A
control is unavailable unless every retained, released and prescribed-control
direction lies inside M. Rows of D outside M remain available to diagnose
unrepresented movable gradients. Merely having a nonzero projection of B
onto M is not complete coverage of B. This avoids manufacturing stable
complements or apparent stationarity from unsampled data.

## Four local controls and their limits

| Control | Retained coordinates | Released coordinates | Additional requirement |
| --- | --- | --- | --- |
| B2 frozen | Q | none | Q and c measured; omitted movable residual is reported |
| B3 atomic release | Q | orthogonal complement of Q within A | Q must be atomic-only |
| B4 joint release | Q | orthogonal complement of Q within B | released block resolved and positive above the registered floor |
| B5 training instability promotion | Q plus complete nonpositive/unresolved released eigensubspace | remaining positive released subspace | combined internal curvature resolved, correct index |

All four see the same input columns and gradient, including preparation
data not ultimately selected. B3 abstains if Q already contains cell motion;
it does not project Q and quietly change the coordinate definition. B5 only
implements the registered **local increment-dimension arm**. It does not
search anharmonic branches, choose a Cmma gauge or demonstrate improved
material accuracy. True zero curvature, excess negative directions, wrong
internal index, missing data or unresolved release produce explicit failure,
not a pseudoinverse or a zero barrier. External control curvature is not part
of the internal index. Where B4 and B5 are both valid on the same quadratic,
their stationary response must agree; equivalence is not an advantage claim.

## Pair the well and bottleneck; keep offsets and domains

`nested_stationary_gap_responses` uses the existing
[paired response contract](STATIONARY_GAP_RESPONSE.md). The initial centre
has expected represented index0, the bottleneck index1. Each has a separately
declared control scale, parameter interval and displacement radius; both use
the same physical parameter, Hamiltonian, mechanical family and energy zero.

The output retains the raw gap, nonstationary-anchor correction, signed gap,
initial and bottleneck responses, first/second parameter derivatives and
omission/reaction diagnostics. It never clips a signed gap into a barrier.
Out-of-domain responses are unavailable. Radius and interval checks do not
prove probe-hull coverage, anharmonic accuracy, branch connection or TS
identity. Ordinary-NEB convergence alone is not a Hessian certificate.

The v2 strong T/Cmma/gauge-matched reference construction, training-only
basis selection, complete B0--B5 forecast freeze, independent errors and
unseen labels remain external requirements. No held-out structure is
generated or read by this module.

## Record cost without giving missing or failed calls a free pass

`ResponseEvaluationCost` and `recorded_response_dataset_cost` accept actual
source evaluation IDs, audit hashes, completed/failed outcomes, recorded
elapsed times and allocations. Exact shared calls are counted once;
conflicting duplicates or missing required identities are rejected. Include
all visible preparation calls, even unused columns and failures. Cache reads
refer to the original call ID rather than create fictitious new evaluations.
Unknown timing/resources give unknown totals and explicit partial sums,
never zero-cost calls. Wall time times allocated cores is a transport-cost
proxy, not CPU utilization, queue time or total study cost.

The bounded `scripts.export_hfo2_scf_costs` independently reads only the two
completed E053 pilot allocations, all154interior SCF audits, original six
physical files and native32-rank E/F/stress. It excludes earlier endpoint
preparation and ongoing continuations explicitly; therefore it is **not**
the eventual full prediction/preparation cost. It cannot launch DFT.

## Reproduce implementation checks, not a material forecast

```sh
python -m scripts.check_nested_response --output /new/nested-analytic-check.json
python -m pytest -q tests/test_nested_response.py tests/test_response_cost.py tests/test_nested_response_benchmark.py
```

The bounded synthetic benchmark checks48paired model points against
independent full-block solves and finite derivatives, including rotated
charts and distinct local scales. It also checks missing-column abstention
and an unstable cell-direction promotion. These residuals are implementation
errors, not DFT uncertainty, Hessian validation or HfO2 prediction errors.
The actual E055 delivery is recorded separately in
[the case evidence](../benchmarks/hfo2_channels/20261008/nested_response_E055_20261010/README.md).
