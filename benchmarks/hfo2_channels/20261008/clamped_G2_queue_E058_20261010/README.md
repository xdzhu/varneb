# E058: second zero-strain G2 wave, queued without extra concurrency

This milestone submits two previously prepared, registered starting chains.
It does not add materials, strain conditions, Hessians or new independent
trials beyond the eight-chain G2 training matrix. Submission is not convergence.

## Actual state, 2026-10-10 CST

| Channel | Slurm handle | Observed state | Latest complete ordinary-NEB row |
|---|---|---|---|
| PO flip, T pattern preserving | 28661019 | RUNNING | step7, 0.467214 eV/Angstrom |
| PO to M | 28661020 | RUNNING | step10, 0.355683 eV/Angstrom |
| PO to T | 28692775 | PENDING / Dependency | not started |
| PO flip, T pattern reversing | 28692776 | PENDING / Dependency | not started |

The parent log observation is at22:37:35; the actual child allocation and
dependency check is at22:41:08. Both parents have no actual
`band/vcneb_failure.json`; neither has reached0.10. They are not cancelled
because of a rebound or replaced by the new channels.

Each child requires **both** parents to finish successfully:
`afterok:28661019:28661020`. Thus this wave cannot add running allocations
while either parent is still active. Each requests one node, one controller
task with32CPUs, four hours and at most10optimizer steps; electronic calls
remain serial over seven moving images, using direct `mpirun -np32` and
one thread per rank. Two fixed endpoints are exact caches, not workers.
The accepted resources and dependencies were independently read from
`scontrol`; two allocations run and two starts wait. No unrelated job is changed.

Both dependencies are execution-success dependencies, **not convergence
dependencies**. A healthy step-capped segment may exit0without reaching0.10.
A timeout or real failure leaves children waiting for reconciliation; it is
not permission to duplicate a chain, blindly rerun this script or edit an
active source tree. No requeue loop, daemon or new automation is created.

## Fixed contract and evidence

The original ABACUS100Ry/full10-auDZP INPUT/KPT/PP/orbital bytes, P=0,
common clamped T-plane, third-vector tilt freedom, ordinary0.10 and noCI
remain unchanged. Only the remaining two **zero-strain** channels are queued.
The four+1% starts remain prepared but unsubmitted. The+0.5% held-out
condition is neither generated nor read. Later G3 retains its existing
two-bottleneck budget; no extra full reference/phonon matrix is authorized.

The runtime is the existing immutable source commit
`86d2637ebcd6c502df43b91ca03977fd50120716`, never replaced. All338Python
files are again byte-compared to its original tested archive, whose SHA256
is`cbe6552a3fa4237fec2098c08c8056735a9a95df4e18691d09fed6b4f24a613c`.
The pilot script is unchanged, SHA256
`7bf764238c4cc857ca2b7e3e94642f0517e1f3b84f4e4c5b17833379ac1101ed`.
The original endpoint preparation receipt and ABACUS binary are also pinned.
The earlier full source regression remains1459passed/2skipped; it is not
claimed as a new full-project or HFpytest run in E058.

Fresh HF preflight completed22:37:11: all18geometries and four exact
endpoint caches passed ordered-atom, continuous-lift, clamped-plane,
original-six-file, native32MPI and full raw E/F/stress checks. Fourteen
moving geometries remain genuinely uncached. Preparation invokes no DFT;
the concurrently running parent SCFs are not counted as zero-cost or complete.
The six E056 seed records are preserved as historical preparation evidence.

- `queue_next_pair.py`: explicit preflight/submission actions; no default
  action, automatic retry or future monitoring. Its exclusive, flushed
  submission journal records intent before each submission and a handle
  immediately after acceptance, including an ambiguous-response stop.
- `executed_queue_next_pair.py`, `preflight.json`, `cache_preflight/`:
  exact executed helper, immutable source identities and actual raw cache audits.
- `submission_journal.jsonl`, `submission_receipt.json`:
  actual22:38:42/43accepted handles, not guessed from queue ordering.
- `verify_queued_jobs.sh`, `queued_jobs_verification.json`:
  one real Slurm resource/dependency observation, not a recurrent monitor.
- `snapshot_live_jobs.sh`, `live_jobs_snapshot.json`:
  parent log/SCF timing observation, explicitly not a full raw-chain TS audit.
- JUnit files: actual offline source and material-delivery checks; no HFpytest
  installation or licensed pseudopotential/orbital redistribution.

The initial Windows tests exposed three path-rendering/mock issues before
upload or submission (11passed/3failed); that report is retained. POSIX
argument rendering was fixed and14guards passed both locally and in an
independent extracted two-file source subset. The first remote checksum
command was misquoted in PowerShell and executed no preflight or job;
the absent namespace and actual helper checksum were reconciled before
the correctly quoted command. No physical parameter was changed for either.

Final delivery checks add two actual-HF receipt fixtures:16cases pass locally,
and all45selected material/guard/delivery cases pass in the explicit independent
archive cwd (45.82s, zero failure/error/skip). All19then-present tracked case
files and the test source match that archive. An initial comparison incorrectly
included ignored Python bytecode and stopped before tests; the corrected
tracked-file comparison and actual clean replay are recorded in`validation.json`.
These selected checks are not a new full-project or HFpytest claim. Subsequent
additions are reports and documentation, not executable or scientific data changes.

At22:37the mean fresh SCFs are177.17s for flip and117.95s for M. Remaining
13/10steps imply about4.48/2.29hours to the20-step segment cap, approximately
Oct11 03:06/00:55CST. These are **cap estimates, not convergence promises**.
The next pair can become eligible only after the later successful parent
completion, then remains subject to actual scheduling availability.

## Next research gate

Finish and audit the matched clamped paths, fill the four registered+1%
training paths, select at most two bottlenecks using training data only,
then execute the finite G3/prediction freeze and independent holdout.
No barrier ranking, response superiority, TS certificate or JCTC readiness
is inferred from queue acceptance, endpoint caches or analytic checks.

Offline replay from the repository root:

```sh
python -m pytest -q tests/test_hfo2_G2_dependency_wave.py tests/test_hfo2_clamped_chains.py tests/test_hfo2_G2_training_matrix_delivery.py
```

Do not replay the HF submission actions as an example quick-start: the
remote namespace is already used and the four handles above already exist.
