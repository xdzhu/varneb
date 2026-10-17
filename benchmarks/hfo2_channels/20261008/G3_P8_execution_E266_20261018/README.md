# E266: execute the frozen remaining P8 candidate tests

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
