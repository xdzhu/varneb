# E152: native raw-output replay, not a complete Hessian or TS

The E150 training batches remain the only two material allocations. This
milestone adds a separate read-only auditor, not a new runtime or SCF. The
deployed E150 source/wrapper, original six physical assets, cached centers,
probe order and sealed physical+0.005condition are unchanged.

`scripts/audit_hfo2_G3_curvature.py` verifies the original center/binary/assets
and declared evaluator contract; snapshots atomic complete `result.json`
files in the original ordered prefix; independently rereads native SCF
energy, force and stress; checks one DSIZE32 SCF, six actual byte hashes,
geometry, call audit and the full physical gradient against its cache. Each
probe log is parsed twice by the existing raw reader and native validator;
these are duplicate reads, not new or double-counted material calculations.
The original center is also reparsed, never recalculated.

## Diagnostics without filling missing directions

At each amplitude, report only complete signed pairs. The existing joint-
stencil assembler produces the measured40-by-k action and k-by-k projection,
raw reciprocity/symmetrization change, energy-derived diagonal and physical-
gradient diagonal disagreement. An incomplete block has no full internal
index or invented missing columns. A complete37-column block separates the
36-dimensional fixed-strain internal space from its prescribed control.

Only directions completed at both amplitudes receive the existing two-step
energy/gradient check and five-point integrated-work residual. Only all148
completed records enable the two-step internal operator/direction/tangent
comparison. Negative eigenvalue counts remain raw observations, not measured
resolution/error-certified signs. No soft direction, plane, TS, frequency,
branch, reaction mechanism or forecast is selected by this audit.

No acceptance tolerance, physical error bound, ASR repair, pseudoinverse or
cutoff/electronic-parameter adjustment is introduced. Incomplete or failed
reservations remain charged; a failure artifact is explicitly reported and
cannot be hidden behind a complete-status label. Native call-time proxy and
allocation cost remain different quantities.

## Reproduction

Run the audited helper from a separate directory, with the immutable E150
source on PYTHONPATH and PYTHONDONTWRITEBYTECODE=1. Example arguments:

```text
--manifest E150_SOURCE/benchmarks/hfo2_channels/20261008/G3_curvature_E150_20261014/manifest.json
--manifest-sha256 5f047f54ca83cb2dceb8cb71687eb9206a85c04205cdf26431565bc8a9b98a5a
--channel T --workdir E150_ROOT/T --output FRESH_EXTERNAL_REPORT.json
```

The output must not already exist or lie inside the live study/source or
physical-assets tree. The auditor neither constructs a CalculatorModeEvaluator
nor submits/releases/cancels jobs. It can be rerun after completion with a
fresh output, without changing or restarting either original allocation.

## Actual implementation checks

95local related checks passed with0failure/error/skip in29.475JUnit seconds.
Independent fresh checkout:95passed/0failure/error/skip in32.278JUnit seconds,
with imports explicitly proved inside that checkout. All6455tracked files
were present before execution; tested tree
09f4aa45e3c823f999ecfc31cea4446c402a0acc, archive SHA256
800aef31eed8801869efa1c440ad7f8e23ac28d0fdd8df2b12b34ed0618626e0.
The two actual source/test texts are normalized-equal to that archive; only
code LF/CRLF normalization is allowed, not physical-evidence normalization.

The first53-test attempt had1test failure: the assertion used a nonexistent
`residuals` field rather than the existing `energy_minus_gradient_work` API.
The test was corrected; the production API/DFT were not changed. Its original
JUnit is retained. The first snapshot shell also mistakenly launched tests
from the original cwd: that95-pass output is retained as a repeat local run,
not independent evidence. A separate actual checkout execution with verified
import roots provides the independent report cited above. No whole E152
regression or new material prediction is claimed.

## Actual native replay

The copied/executed helper has raw SHA256
09ebabb3f74dfbbce00604f32d8bef1a9fc48e65d08812f8cfecffcde524a49b,
equal to the tested working script. It executes from separate native audit
namespace20261014-varneb-G3-curvature-audit-E152-r1 with the unchanged
E150 runtime imports. Each channel's actual snapshotted prefix has34native
results,17complete signed pairs at.01Angstrom and no failure artifacts.
All68distinct probe SCFs pass the original six actual byte hashes, unique
DSIZE32, electronic convergence, full ordered geometry and raw E/F/stress/
gradient replay. This does not mean all148points in either batch finished.

| Existing batch | Reparsed prefix | Measured pairs | Partial reciprocity defect | Maximum energy/gradient diagonal difference(eV/Angstrom²) |
|---|---:|---:|---:|---:|
| T | 34 | 17 | 0.000145486 | 0.0584274 |
| Preserving | 34 | 17 | 0.000119830 | 0.1890415 |

No physical threshold is assigned to those differences. The second amplitude
is not yet available in this snapshot, so finite-step/energy-force effects
are unresolved; no missing internal/cell/control direction is assumed stable.
Full Hessians, branch tests, selection/error resolution and independent
forecasts remain incomplete. The summed call times3509.200/3603.502s are
subsets of the already charged E150 calls, not fresh work or allocation cost.

Reports are native byte copies with matching remote/local SHA256, not
regenerated normalized replacements. At19:42:32CST the same two Slurm handles
28855614/28855615 are RUNNING,32CPU/node10,node27. No new scheduler action,
DFT, live-source change or reserved physical geometry. Goal remains active.
