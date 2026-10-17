# E260: bounded preserving-branch experiment

2026-10-17. The author approved the [paper blueprint](../../../../docs/HFO2_PAPER_BLUEPRINT_2026-10-17.md),
the new third-coordinate protocol, at most10 training candidates and the
subsequent independent-prediction route. The exact pre-label experiment is
in [protocol.md](protocol.md). This first stage computes training only.

Five Q pairs, each from independently fixed positive/negative0.10A starts
along the original E255 omitted direction, give10 candidates. Q is fixed;
the explicit internal zeta coordinate can relax with the other allowed
directions and can cross sign. No other point's result supplies a start.
The norm target remains0.03eV/A; a separately declared1e-4eV/A directional
test prevents an unresolved shallow direction from passing immediately.
Neither criterion changes the ABACUS SCF or ordinaryNEB0.10 settings.

The numerical seed is positive on the free space: measured first-step R33
plus the outer stiffness0.04986124402846963eV/A2 from all four seen gradient
probes. Its zero cross blocks are a declared preconditioner choice, not
modified physical curvature. Original two-step negative signs and the2D
refusal stay in the source. No stable branch or stationary TS is presumed.

## Exact scope and accounting

- Original ABACUS/PBE/100Ry/full10auDZP/Gamma2x2x2/P0/E0/binary/six assets.
- Training commitments T34/P39 -> T34/P49, each cap81; unused slots are caps.
- Up to12 new attempted SCFs/candidate,120 total, no retry/rollover.
- Two allocations maximum including held/pending; each one node/32 genuine
  mpirun ranks, OMP/MKL/OpenBLAS1,hfacnormal01,4h maximum.
- Original T8 are seen; preserving8 remain unseen during training. A revised
  target/domain/model/control forecast freezes before the remaining8 labels.
- Physical+0.005 geometry/labels and private author arrays remain sealed.
- All296 original derivatives,272 post-release derivatives, training and
  other preparation remain additional cost. No acceleration claim here.

The old380-file HF runtime is read-only. The new portable
[branch optimizer](../../../../vcneb/branch_relaxation.py) is separately
pinned in this small external bundle, so no running source is overwritten.
The two controller shards contain five starts each, not endpoints or full NEB.

## Validation and reproducibility

The first analytic/source run had41passes/1failure: a hand-entered expected
preconditioner number was inaccurate. It was replaced by the value from
the unchanged four-probe formula; the original XML is retained. Final
optimizer/source checks passed43locally and43in an independent Git export.
Eight additional dispatcher checks use the actual archived Slurm schema,
including `CPUs/Task` and exact `NumNodes=1-1`, with no scheduler mutation.
These are implementation checks, not material validation or a full-suite run.

The original inert r1 preflight refused exact equality of arithmetic arrays.
Its read-only [diagnosis](arithmetic_diagnosis_r1.json) found1382 last-bit
differences, maximum1.0658141036401503e-14. Physical inputs, sources, target
Q/sign/amplitude and policy matched; no inert decks, new DFT or jobs existed.
The original bundle, config, manifest and first failure are preserved.

The r2 check keeps every source, policy and target identity exact, verifies
matrix-product reconstruction with the existing2e-12 arithmetic tolerance,
and uses frozen optimizer seeds/internal starts. Native preflight writes
the same internal-to-full embedding used by the actual evaluator. It passed
all ten static decks and the original native origin/runtime/asset checks.
This does not alter energies, SCF thresholds or scientific error criteria.

Current config SHA256:
`63b1aa0eea88e7cd7c2f1a78ce9f4f3c7f39e16a9caf04d3f6d8c6ef5f84ab80`.
Bundle manifest SHA256:
`4422902d6474afe1056181b978f6098c9fc1284ff6cc66d11ff865d499d87eff`.
The independent preparation CLI reproduces the config byte-for-byte.

```sh
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
python -m pytest -q tests/test_branch_relaxation.py \
 tests/test_hfo2_preserving_branch_trial.py tests/test_hfo2_branch_dispatch.py \
 tests/test_conditional_hessian_seed.py tests/test_hfo2_T_original_independent.py
```

`worker.py prepare` accepts a fresh output config and performs no DFT.
`preflight` requires the pinned HF sources. `dispatch.py` has one-shot held
submission and exact-handle release modes; persistent intents prevent blind
resubmission. Do not rerun it to check progress. Native results and final
curvature/branch identities need a separate read-only audit after calculation.

## Execution location and first handles

HF namespace:
`/public/home/iai806/abacus/agent-runs/20261017-varneb-G3-P-branches-E260-r2`.
Original rejected inert r1 remains at the corresponding `-r1` path.
28934048 owns IDs0/2/4/6/8;28934049 owns1/3/5/7/9.
Both were submitted held at23:23:49CST. Accepted/release/status receipts
record later live states; submission alone is not completion or convergence.
No remaining independent labels or physical-strain holdouts are computed
by this training submission. The full research goal remains incomplete.

At23:32:30CST both jobs were RUNNING,436s elapsed. Six complete evaluator
records and no failure records existed, with no completed point row yet.
The first SCF of each shard was independently reparsed against native
E/F/stress, six input hashes and32 genuine MPI ranks:119.285/109.394s.
This [first native prefix](first_native_prefix_2331.json) is a bounded startup
proof, not a complete audit of the later training dataset.

The final optimizer/source/dispatch/audit-admission tests passed57local and
57independent checks across6modules,32new checks, zero failures/errors/skips;
not a full repository run. Source treee33b523efd584f6a360a773408693ec0fac8a0f0,
archive SHA256e63f51d86086c9fad1b12ca7b570c4c0850b5dab7dc6cbcca80969e47b211cd2.
The [complete read-only auditor](audit.py) refuses live, failed, partial or
wrong-resource jobs and later rechecks every charged native attempt and
both stopping criteria. Its admission unit fixtures are synthetic, not
completed DFT labels. validation.json records the final verification scope.

At the observed109--119s/SCF, four to twelve attempts per point imply
roughly40min--2h for the two simultaneous five-point shards. A planning
window is Oct18 00:15--01:30CST, conditional on healthy calls and scheduler
continuity. The actual4h allocation caps are03:25:14. This is processing-time
estimation, not guaranteed convergence; unsuccessful candidates retain their
original outcomes and costs. No periodic monitor is created.

## E262: first-four native replay, not a completed cohort (Oct18)

At00:10:10CST both original allocations were still RUNNING. The separate
[prefix auditor](prefix_review_E262.py) reparsed all36 native SCFs of the
first four completed candidates, including both unsuccessful starts. It
checked original source/input identities, ordered attempts, actual E/F/stress
and MPI provenance. An ordered replay of the unchanged optimizer consumed
every record with zero coordinate difference on HF and no new DFT. It
recovered accepted states, rather than treating arbitrary trial records as
an optimization trajectory. The complete ten-candidate auditor is unchanged
and still requires both original jobs to finish normally.

| Q / starts | Native SCFs | Outcome at the preset cap |
|---|---:|---|
| (0,0), positive / negative |12 /12|Both capped, not converged|
| (-0.1,0), positive / negative |6 /6|Both pass the declared gradient tests; stability unchecked|

The [read-only diagnosis](prefix_diagnosis_E262.json) identifies a concrete
numerical limitation at Q=0. All12 native records in each start were accepted
energy-decreasing states, not repeated energy-rejected evaluations. Energy
decreases6.710/6.714meV per cell, but the final full released norms remain
0.29791/0.29780eV/A. Total radii reach0.49999059/0.49999942A against the
original0.5A guard. Constant-radius tangent gradient norms remain
0.29559/0.29547eV/A; these are neither stationary candidates nor boundary
KKT solutions. Successively tiny accepted steps are consistent with the
scalar search being limited by the original radius, not an SCF failure.
The radius is centered at internal-coordinate zero, **not** the relaxed
origin. No radius extension, projected-gradient success criterion or retry
is authorized by this diagnosis.

At Q=(-0.1,0), final released norms are0.001411/0.001446eV/A and signed
zeta gradients -7.6003e-5/+7.3270e-5eV/A. Final zeta values are
-0.001461/+0.001292A: opposite initial signs do not prove distinct branches.
The final structures are0.002771A apart in this internal metric, with an
energy difference0.000675meV/cell. Neither positive curvature, uniqueness,
a whole sheet nor independent predictive gain follows from these results.

All36 native calls retain their costs:37.6168 completed-call proxy core-hours,
with median119.58s/SCF (96.63--155.43s). This is not full allocation cost or
the entire preparation cost. The native report SHA256 is
`c122d71eb2dd46f3a3e892e49389b5d38a24a38eb86794ee0dd35ea7bc100772`.

Seven new replay/decomposition checks plus15 related checks pass locally
and in a targeted independent Git export (22each; not fullsuite). The exact
export scope, tree/archive and XML hashes are in [validation_E262.json](validation_E262.json).
One initial diagnostic assertion failed because
local BLAS reconstruction differed by at most2.776e-17A. The raw failed XML
is retained. Exact accepted energy/order identifies the record; the existing
2e-12 matrix-reconstruction check handles arithmetic only. Native record
identities and the original HF replay remain exact; no physical tolerance
changed. The remaining six candidates keep their original runtime, starts,
caps and inputs. P8 independent labels and physical+0.005 remain unopened.
