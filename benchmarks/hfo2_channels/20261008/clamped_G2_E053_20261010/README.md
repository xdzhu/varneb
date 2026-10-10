# G2: matched clamped switching versus escape

These are **new unevaluated starting chains**, not clamped barrier results.
The first pair is zero-strain `PO_flip_T_pattern_preserving` and `PO_to_M`.
The remaining two channels and +1% training condition retain their registered
place in the <=8chain budget; the +0.5% holdout is not generated or read.

The common T-substrate plane, P=0/E=0, released third-vector tilt/length,
100Ry/10au basis, INPUT/KPT/pseudo/orbital bytes and ordered Hf4O8 are fixed.
Nine total images means seven fresh internal SCFs and two fixed wells,
whose actual **clamped** terminal raw SCFs are reused after exact geometry
and input/log hash checks. Free-cell G1 geometry informs the starting polygon,
not any new energy, force or stress. Source index interpolation is only a
starting-path construction, not an arclength/MEP assertion.

The narrow case preparation is `scripts.prepare_hfo2_clamped_chains`.
It records wrapping/lifting of endpoint BFGS raw STRUs, the original ordered
G1 endpoint lift, endpoint relaxation corrections, source polygon resampling,
and hashes. It refuses changed source endpoints, remapping, unsafe distances,
unregistered conditions, stale caches or overwrites. The dedicated factory
does not depend on optional ASE ABACUS I/O and never caches the interiors.

Pilot recipe: FIRE10steps/4h, .10eV/A ordinary threshold/no CI, .02maxstep,
.2spring, candidate geometry/continuous-lift checks, one32MPI HF allocation
per chain/direct mpirun/thread1; at most two concurrent study allocations.
Slurm timeout or healthy step-cap is reconciled with raw observations, not
automatically resubmitted. No running source or endpoint runtime is overwritten.

This is a bounded scientific comparison, not a certified saddle or new
optimizer benchmark. Finishing the matched barriers, sampling/error checks,
strong controls and unseen prediction is still necessary for the research claim.

## Actual launch, 2026-10-10

The immutable executable/observable tree is c3e4bb13, clean full regression
1373passed/2skipped/0failed;140new endpoint observables byte-matched.
The binary-inclusive seed bundle is transported separately and each artifact
hash is checked on HF. Two exact cache/geometry preflights passed with external
executable launch prohibited. Source or physical inputs were not modified to
resolve the login-only environment/NumPy probe checks; [receipt](submission.json)
records these two zero-DFT harness failures and the reviewed correction.

Jobs28574708 (preserving flip, node28) and28574709 (PO→M, node142) started
11:32:38CST, both confirmedRUNNING11:34. Each has one32CPU allocation; only
seven internal images evaluate DFT serially inside its own allocation. This
two-allocation choice follows the finite study concurrency cap, not a claim
that the universal manager requires serial images. Ordinary .10/noCI applies.
The first10step segments have4h caps ending15:32CST if still running. Based
on endpoint SCFs near2minutes, a seven-image round is roughly14minutes;
initial11rounds are about2.5h before variation in intermediate SCF difficulty.
This estimates the pilot segment, **not guaranteed path convergence**.

No automatic continuation, scheduler monitor, release or holdout is created.
These handles must be reconciled before any next submission. Next independent
work is prepare the already registered strong-control analysis and error/cost
tables without reading the +0.5% held-out condition.
