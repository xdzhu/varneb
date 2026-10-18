# Physical +0.5% holdout: sampled-path response (E268)

## Completed result, 2026-10-18

Both original handles completed0:0 on hfacnormal01:28938772 at08:10:04CST
and28938773 at08:40:05CST. Both meet the registered ordinary criterion.

| Channel | Sampled peak (meV/HfO2) | Change from zero | NEB residual (eV/A) | Updates | Interior SCFs |
| --- | ---: | ---: | ---: | ---: | ---: |
| Preserving flip | 98.530877 | 1.702694 | 0.079773 | 0 | 7 |
| PO-to-T | 90.418178 | 6.001554 | 0.080825 | 2 | 21 |

The selected-pair T-minus-preserving contrast is -8.112699meV/HfO2,
4.298861 above its zero-strain value. It remains negative. B0's signed
errors are+12.595430/+0.171959meV/HfO2; B1-fixed gives-2.933378/-7.232239;
B1-follow gives-1.667843/-5.004639 in preserving/T order. All138frozen
records, including132unavailable transfers, are scored without refitting.
The test establishes no predictive gain for the proposed local reduction.

Native replay checked every28interior call, the six original physical input
files, real32MPI and complete E/F/stress. All11endpoint calls were already
checked before prediction freeze. Total physical-test cost is39SCFs,
47.289888native-call proxy core-hours and48.808889allocation core-hours;
all previous training/curvature/model costs remain additional.

The transferred archive has SHA256
`72276fb0fa799e2c4ce07dc31ad270ca246e60cc88f81bec59a375d1307ee68a`.
Its242payload files and their hashes are in
[completed_paths_20261018/files.json](completed_paths_20261018/files.json).
It contains both terminal observations and all28native interior call logs,
inputs, geometries and receipts; pseudopotential/orbital bytes are excluded.
The existing terminal reader was copied to a fresh postprocess directory
after both jobs ended. The first collection attempt found that reader had
not been deployed at the root and created no output; deploying the unchanged
reader resolved the postprocess invocation without any new DFT.

The bounded warm-start test has finished. A zero-update preserving pass is
reported explicitly, not presented as a fresh path search or a speedup.
No additional response models or calculations are triggered by this result.
The [paper analysis and plots](../../../../paper/VARNEB_JCTC/figures/hfo2_independent_20261018/README.md)
replay the frozen forecasts, selected-pair responses and original-frame modes.

The remaining sections preserve the registration and deployment chronology.

The author delegated protocol decisions on 2026-10-18. The researcher selected
the finite sampled-peak test in [protocol.md](protocol.md). Original certified
selector refusals and the previous T8/P8 evidence remain unchanged.

Before any target geometry, `registration.json` was generated from the pinned
E148 eight-band report. External registration SHA256:
`faebe7588579f5dca9253e672755d4db5c85f8a33165685dc5cdce3cef2899ff`.
It fixes the two representatives, B0, B1 formulas/visible endpoint policy,
all strong references, unavailable transfer models and execution budgets.
It is not the final numeric forecast panel: B1 awaits the three predeclared
endpoint features. No target interior may be evaluated until that panel is
frozen separately.

B0 predicts 111.126307 and 90.590137 meV/HfO2 for preserving and PO-to-T,
respectively. These are sampled-peak predictions without physical error bounds.

`study.py register` is a zero-DFT replay; `study.py endpoints` requires the
external registration pin and creates only three endpoint seeds in a fresh
directory. Neither invokes a calculator or submits a job. Source checkout
with the project's normal scientific dependencies is required.

Related numerical and refusal tests: 40 passed (`test_hfo2_sampled_holdout.py`
and unchanged strict `test_prediction_controls.py`). The first direct-file
CLI attempt lacked the repository import path and wrote no artifact; the
entry point was repaired and all40 rerun successfully. This was a local
invocation issue, not a DFT failure or changed scientific protocol.

At registration: target endpoint/path DFT calls0. Deployment and actual
results will be recorded below with immutable handles and hashes.

## Endpoint deployment

After registration commit `77d1a7ca`, three +.005 endpoint-only seeds were
prepared from the audited zero-strain wells. No old SCF result was transferred
to strained geometry. Actual HF preflight checked all six original input
hashes, the ABACUS binary, all375 exported runtime files, and all three native
STRU roundtrips. No target interior was generated/evaluated.

Bundle pin: `689a26fed0f47dede6d7c2f2230d3181fe5a16bbc5eee0356da6f033d89987cf`.
Source archive: `8e1457e90827696c6d976a1dbca46fb43a5421acda1fdb41d7bd763bf99a00fc`.
Remote root: `/public/home/iai806/abacus/agent-runs/20261018-varneb-physical-E268-r1`.
Archive contains the375 existing runtime files from `77d1a7ca`; the frozen
endpoint worker is in the separate bundle. No running source is overwritten.

57 related local checks passed, including real seed-boundary/hash checks,
numeric controls and scheduler-resource refusal checks. A preceding test
command named a nonexistent module and ran zero tests; the corrected command
is the recorded57, not an additional test success. The native preflight is
actual HF evidence, not a local synthetic DFT validation.

Original submitted endpoint handles are `28938150` (PO+ then preserving PO-)
and `28938151` (T). Both were independently verified held before release,
then started at07:26:26CST on node6, each one node/32CPUs; sharing a node
does not change each native32MPI group. Exact submission receipts are retained.

The follow-on `forecast_and_chain.py` uses the existing complete-panel freezer
and audited endpoint/path readers. It preserves69 registered model/variant
records across two cases (three numerical baselines plus unavailable physical
transfers); these are not69 newly fitted models or successful predictions.
The two endpoint-feature B1 values are substituted only after native endpoint
and identity checks. Synthetic freezer plus existing chain tests:80 passed.
No +.005 interior has been evaluated at this stage.

## Native endpoint results and prospective numeric freeze

Both original jobs completed0:0:28938150 at07:43:17CST and28938151 at07:37:45CST.
All11 actual SCFs replay from native energy/forces/stress logs. PO+, preserving
PO- and T consumed3/3/5SCFs; final atomic forces were0.016430/0.017196/0.0000275
eV/A and open tractions1.1806/1.5335/1.7154kbar, within the unchanged0.03eV/A
and2kbar endpoint criteria. All three symmetry-tolerance screens retained
Pca21/Pca21/P42/nmc, respectively; the two PO polar displacement signs and
the registered preserving pattern remained distinct. These identity screens
are not Hessian-stability or electronic-polarization certificates.

The89-file portable native export excludes pseudopotential/orbital bytes;
their original HF full-byte audit remains preserved. Export archive SHA256:
`3299b20e11355fa468ea7d495a42a41f9140f8a1793db70bb39a2692f6851e18`.
The first local extraction attempt used a newer Python tar API and created
no extraction; explicit member/type/path/hash checks permitted extraction
with the installed Python. No DFT or scientific settings changed.

Audited feature SHA256:
`0441a5e7e63ec86a68807068ea21db6901f4a762d740592473d91bf15bbd8f23`.
Final numeric forecast manifest SHA256, frozen before any target interior
was generated or evaluated:
`1dd8dc4cacb22b073a6402acb38ffb956b3e0dd3738a85db18fab2afb279a88b`.

| Forecast (meV/HfO2) | Preserving flip | PO-to-T |
| --- | ---: | ---: |
| B0 midpoint | 111.126307 | 90.590137 |
| B1 fixed absolute peak | 95.597499 | 83.185940 |
| B1 final-well-following peak | 96.863034 | 85.413539 |

These are six numeric records, with132 explicit unavailable-transfer records
retained, not discarded competitors or zero errors. No physical error bounds
are inferred from the ordinary NEB force tolerance. Preparation cost:
13.703824 native-SCF core-hours,15.022222 allocated core-hours. All previous
training, curvature and P8 costs remain additional.

Independent scoped regression:162 passed,0failed/0errors/0skipped in39.811s,
covering the new sampled protocol, unchanged strict prediction controls,
forecast freezer, clamped endpoints/chains and prediction evaluation. Raw
JUnit evidence is `independent-tests.xml`; this is not a full-suite claim.
The tested staged tree was`e0af9051951b5d3cd404314d6f0f3f3b658ebebe`;
the scoped supplement archive SHA256 was
`79f5500d7f3ebbb2b4ff5d786756908e1c01fafde8d55d7ab65814b1a4680edd`.

## Actual two-path deployment

The complete prediction panel was pushed in `a8dcca7d` before constructing
the target interiors. The separate path bundle retains the original runtime
and physical inputs; no source used by the completed endpoint jobs was edited.
Its external SHA256 is
`ce9bfed78e7b3e8d9c1708f4cc72879e977ef1c54653f7b45e04675611edfe7b`.
Actual native preflight validated all 18 input structures and four endpoint
cache uses, with fresh DFT disabled. All 96 deployed source/forecast/seed
files also reproduce byte-for-byte from the `165fd143` Git archive, as do
all 89 portable endpoint evidence files. The remaining generated config is
the two-channel list plus the already recorded forecast manifest pin.

Original path handles, submitted held and resource-checked before release:

| Job | Channel | Resources | Fixed limit |
| --- | --- | --- | --- |
| 28938772 | PO+ to preserving PO- | hfacnormal01, one node, 32 MPI | 80 updates / 24 h |
| 28938773 | PO+ to T | hfacnormal01, one node, 32 MPI | 80 updates / 24 h |

Both started at 07:53:24 CST on 2026-10-18, on node6 in separate allocations.
At 07:54:33 both were RUNNING and each had begun one fresh interior SCF;
no complete chain or convergence claim follows. The seven interior images
are evaluated sequentially within each allocation; the two paths run in
parallel. Each image uses native `mpirun -np 32`, with one thread per rank.
The new converged endpoints are cached, so they are not recalculated each
iteration. Ordinary NEB uses 0.10 eV/A and no climbing image.

`inspect_paths.py` makes one read-only snapshot of these exact handles.
Its optional prefix audit reparses only the first completed interior call
per channel; this is a launch-health check, not the full terminal audit.
`audit_paths.py` will retain converged versus bounded-step outcomes, all
native calls, point errors, unavailable predictions and actual costs.
No restart, follow-up submission, automatic parameter change, or model
refit is performed by either inspector.

At 07:57:02 CST, both jobs remained RUNNING and each had completed its first
interior SCF and begun the second. Native replay confirmed electronic
convergence, complete energy/forces/stress, the six unchanged physical input
hashes and genuine MPI32 (no singleton fallback) in both first calls.
They took 131.078 s and 181.611 s, respectively. Evidence is retained in
`path_first_call_snapshot.json`. This is a two-call startup audit, not an
all-call audit or a converged path. If subsequent calls have similar cost,
one seven-interior chain evaluation is roughly 15--22 minutes; total NEB
completion cannot be inferred before observing its residual trajectory.
