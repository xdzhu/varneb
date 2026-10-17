# E266: execute the frozen remaining P8 candidate tests

**Latest E267 result:** both original jobs completed normally;8/8candidates,
24native SCFs audited and all23frozen arms scored. See the complete-result
section below. Earlier preparation/startup paragraphs are historical, not
instructions to submit again. Physical+.005remains sealed.

2026-10-18. This implements the local-test stage of the author-approved
[blueprint Section5C](../../../../docs/HFO2_PAPER_BLUEPRINT_2026-10-17.md):
freeze the revised local models before P8; **afterwards** complete the
separate physical+0.005 forecast and at most two paths. The latter remains
sealed here. No new training family or calibration can use the P8 labels.

The exact E265 [protocol](protocol.md) and [forecasts](forecasts.json) retain
four surrogate families/steps (31values/one refusal), all19 original controls,
the old2Drefusal and the unverified-stability warning. This is a separately
registered gradient-candidate energy target, not certification of a stable
conditional PES, a TS or an advantage of the nonlinear fit.

## Immutable inputs and scope

- Exactly8 original half-grid Q points, IDs1--8, original E264 training-only
  starts/parents2/4/7/8/2/2/4/7. No independent-to-independent continuation.
- At most12 new attempted SCFs per point,96total, no rollover/retry/new start.
- Same E260 BFGS/positive numerical seed/trust0.05A/max200 updates; both
  released norm0.03eV/A and zeta-gradient1e-4eV/A; no post-label loosening.
- Original0.5A total internal-coordinate radius about zero and1.6A distance.
- ABACUS/PBE/100Ry/full10auDZP/original Gamma2x2x2/binary/six physical hashes/
  P0/E0 unchanged. Old380-file runtime and original E260 eight-file bundle
  are read-only; the external wrapper never overwrites them.
- Two disjoint four-point shards[1,2,3,4]/[5,6,7,8], each hf/hfacnormal01/
  1node/32genuine MPI/mpirun/threads1/4h. Count pending/held in the two-allocation
  limit. Endpoints/full paths are not recomputed by this local test.

Config SHA256`2ca8fd47f8a53208501c046a0afe38c3d42f8792db0ac989c1139839b6ab22ae`;
six-file manifest SHA256`100f1a55d4c18f54a82412b346d26085aa33f1c57e3f43bbb956ac95fff548dc`.
The copied forecast/start/protocol bytes retain their original external pins.
`freeze_local.py` refuses existing artifact names; do not rerun it for status.

## Zero-DFT native preflight

The new, previously absent HF namespace is
`/public/home/iai806/abacus/agent-runs/20261018-varneb-G3-P8-E266-r1`.
Only this external bundle was copied. `bash -n` passes; the original runtime,
native origin and physical assets revalidate. All8 static input decks use
the native embedding of frozen internal starts, preserving the existing
2e-12 arithmetic-reconstruction check, not changing physical tolerances.
All six physical hashes match at each deck; no calculator was invoked.
[native_preflight.json](native_preflight.json) SHA256 is
`96a45eef626f0c0b4a2f0d05580be1cb2a6145761dbeeea499d50bc4a19ae8be`.
The queue was empty at that snapshot, not proof of a later scheduler state.

Source checks cover immutable target/config replay, budget/gradient/MPI
mutation refusal, unallocated execution, actual Slurm one-node range syntax,
wrong held resources and a synthetic full four-target controller using an
analytic oracle. The latter verifies reader routing/point accounting without
claiming real MPI, physical geometry or material convergence. Actual native
input preflight is a separate check; completed SCF/independent errors require
later native audit. No full regression or independent material success is
claimed from these software tests.

`dispatch.py` is a one-shot controller, not a status command. It submits held,
records intents/handles, verifies exact resources, records acceptance, then
releases the **same IDs once**. An interrupted/unknown outcome must be
inspected, never blindly resubmitted. Wrong resources remain held. The worker
checks its accepted Slurm identity before constructing any calculator.

Predictions/model/domain are fixed before any label. Keep all failed attempts,
missing labels, refusals and costs; the ordinary0.10NEB target is unchanged.
The seven E260 gradient candidates are not stable branches; the same limitation
applies to any new local label. All prior training/curvature/path preparation
costs remain additional, and T8 remain already seen.

Preparation verification passes30local/30independent related checks, three
modules/13new checks, no failures/errors/skips; not fullsuite. The full-source
export and all XML/config/native-input hashes are recorded in
[validation_preparation.json](validation_preparation.json). These checks
precede submission and do not claim a completed independent test.

## Actual deployment and native startup (E266, Oct18)

Prepared source was pushed as `1ef27c1c` before submission. Original handles
`28935677`/`28935679` were submitted held, verified at one node/32CPUs/four
hours, accepted and then released exactly once. They started at01:56:11/12
CST on node13/node40 respectively. The held/accepted/released JSON receipts
are historical scheduler snapshots, not current status commands. Do not rerun
`dispatch.py` or submit a replacement because a status observation is missing.

The read-only [startup audit](native_startup.json) reparses the first actual
SCF of each shard, confirms full E/F/stress, genuine MPI32 and all original
physical hashes. Durations139.7664/123.5693s are two completed-call timings,
not a forecast of optimizer convergence. Its SHA256 is
`a215c09e8f0ae40591a75512f8ea66106e333a0d9795ebd4c0728b62c1c131bf`.
This is a two-call prefix, **not** a complete eight-target audit or stability
certificate. P8 is now being evaluated; physical+.005 remains sealed.

[postprocess.py](postprocess.py) lives outside the immutable running bundle.
It never invokes DFT or changes the scheduler. `audit` requires both original
jobs COMPLETED0:0 and reparses every charged attempt, native source, starting
geometry, budget and dual-gradient result. `score` requires an externally
pinned complete audit and the exact pre-label forecast; it reports all23
arms/all184 point records, refusals, missing labels, raw and origin-aligned
errors and comparable subsets without refitting or changing the denominator.
The revised candidate target does not retroactively validate the original
stable conditional-plane target.

The combined software checks pass37local/37independent, four modules/seven
new postprocess checks, no failures/errors/skips; not a full regression.
The independent source is the previously tested full Git export plus an
exact two-file Git-archive supplement, not a newly claimed full export.
Hashes and verification scope are in [validation_execution.json](validation_execution.json).
Synthetic scoring/routing tests are not material accuracy or real-MPI evidence.

The manuscript Methods2.4/Results3.3--3.4 now integrate the completed ten-start
training result and the fixed fit comparison, without reporting pending P8
errors. All original tables/figures and four Results subsections remain.
Related checks pass138locally/138in the independent source, twelve modules,
not fullsuite; two intermediate failures were stale status-text guards, not
physical/data failures. The existing LuaLaTeX build succeeds,40pages; every
page rendered, five contact sheets and new pages6/17/18/19 actually inspected.
No new clipping/notation/layout error found; older sparse figure layout is
retained, not final editorial approval. Raw XML and
[manuscript_integration.json](manuscript_integration.json) preserve scope/hashes.
No production source, parameter, model or threshold changed in this integration.

## E267: complete native cohort and prospective scores

Both original handles are COMPLETED0:0:28935677 in1069s and28935679 in2075s,
ending02:14:00/02:30:47CST. Every one of24new SCFs was reparsed for ordered
geometry, exact inputs/native E/F/stress, genuine MPI32, source identity and
budget. All8candidates meet the original released-norm0.03eV/A and directional
1e-4eV/A tests; SCFcounts are2/2/2/2/6/4/3/3. No new start, retry or changed
threshold. [Complete native audit](native_complete_E267.json) SHA256:
`dd245ac8f7335442ef5355b484b5c01dcba03ae3201ec1ad02131a3ae16ac654`.

The unchanged read-only scorer used the exact pre-label forecasts and this
external audit pin. [All23-arm score](score_E267.json) SHA256:
`c51dec44d4e92e05ad824bd7952a9e1cee16d662476684675c9ae0c7489c2607`.
Four new arms share targets1--7; harmonic0.02still abstains at8although its
DFT candidate converged. On that common set, harmonic/nonlinear maximum
errors are0.033195/0.035625meV/fu at0.01and0.033928/0.036390at0.02;
nonlinear RMSE is slightly smaller. There is no metric-independent or
physically resolved nonlinear accuracy advantage. All new raw/aligned
predictions coincide because the same native origin was fixed, not refitted.

All19old controls and their new-target warnings remain. Tensor's raw maximum
is0.029103meV/fu but it has a different registered representation/information
set; a numerical comparison does not revalidate its old stable-plane claim.
The two restricted T-triplet arms have no supported target; each Cmma arm
supports2/8. The all23-arm common set is empty, not a zero-error success.
Full tables and all184signed records are retained, not only favorable arms.

Actual new cost is26.469276completed-call proxy and27.946667allocation
core-hours. Training81SCFs/85.6allocation core-hours, original296+272
derivatives, seen T8 and all other preparations stay additional. Source-point
R33positivity and these candidate energies do not certify a stable conditional
sheet/TS. The original2Drefusal remains; neither a new Hessian nor another
training family is scheduled. Physical+.005/private arrays remain unopened.

The [public standard-library replay](../../../../examples/cases/hfo2_competing_channels/README.md)
now exports all8P candidates,23arms/184predictions alongside unchanged T/G2
tables. It verifies the exact pre-label predictions and recomputes signed-error
units/statistics; it performs no new native audit or DFT. Paper Results3.4/
Table5 and AppendixE.9/TableA19 report these data, without claiming JCTC
readiness. The next physical-path stage awaits the author's clarification on
sampled-peak versus strict error-certified selection; the old selector is not
silently relaxed and the goal is not declared complete.

Final related checks pass149local/149independent, thirteen modules/seven new,
zero failures/errors/skips; not fullsuite. The independent source is the
recorded E266 full export/supplements plus an exact thirteen-file E267 Git
supplement, not a new full-source export. Actual `python -S` CLI replays
produce seven byte-identical files. Initial146pass/2failure XML retains the
two stale182-versus213table-count guards; prior row-value guards remain intact.
Both new tables are source-derived checked, including coverage and missing
values. [validation_E267.json](validation_E267.json) records all source/data/
XML/export hashes and verification scope.

The existing LuaLaTeX build succeeds to42pages without compiler/layout
warnings. All final pages are rendered. Seven contact sheets were inspected;
after a caption-only repair, first40page pixels are unchanged and final
pages19/41/42 individually reviewed. New TableA19's caption now stays with
its table; no new clipping/overlap/broken notation found. Pre-existing sparse
figure pages are retained, not final editorial sign-off or submission readiness.
