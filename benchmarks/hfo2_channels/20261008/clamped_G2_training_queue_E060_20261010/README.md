# E060: finish queueing the registered +1% training matrix

The case identifier preserves this turn's 2026-10-10 start date; actual
checks and accepted handles carry their own CST timestamps. This is the
remaining four starts of the **existing eight-chain G2 matrix**, not four
extra experiments. Queue acceptance is not ordinary-NEB convergence.

The proposed finite graph is:

```text
28661019 + 28661020  (existing zero-strain continuations)
          |
          | existing afterok dependency; no mutation
          v
28692775 + 28692776  (existing zero-strain T/reversing starts)
          |
          | new afterok dependency
          v
+1% PO->T + PO->M
          |
          | afterok uses the two actual newly accepted handles
          v
+1% pattern-preserving flip + pattern-reversing flip
```

Each new start retains the existing tested pilot: one hfacnormal01 node,
one controller task with32CPUs,4hwalltime,10optimizer steps,seven moving
images computed serially with direct32MPI and one thread per rank. Two
fixed endpoints use exact terminal caches. At most77fresh interior SCFs
are needed for a healthy complete10-step starter; convergence can stop
earlier. Continuations remain the same independent chain, not new trials.

The P=0/common T-plane/third-vector-tilt/ordinary0.10/noCI and original
ABACUS100Ry/full10auDZP/sixphysical-byte contract are unchanged. The runtime
is the existing immutable86d2637source and tested archive, not an overwrite
of any active run. This batch reads only the four prepared strain_p0100
starts from E056. No+0.5%holdout geometry/label, additional material,
Hessian, G3surface or electronic-parameter variant is created.

Before any submission the helper must read all four known handles from
sacct, inspect the user's queue without touching unrelated pc-* jobs,
verify the actual zero-strain successor dependency/resource/source graph,
recheck338source bytes against the original archive, original binary,
registered seeds,36geometries and8native32MPI/full-EFS endpoint caches.
The preflight cannot launch a DFT/external executable when reading caches.

The exclusive, flushed submission journal records intent before each
sbatch and accepted handles immediately. The second pair depends on the
**returned** first-pair handles. An ambiguous reply, timeout, duplicate,
failed predecessor or changed source stops; it is never blindly retried.
afterok is execution success, not a0.10 certificate: healthy step-capped
paths may exit0without convergence. A failed/timed-out parent requires
authoritative reconciliation, not automatic cancellation or restart.

Actual preflight/submission/resource receipts are added only after those
actions succeed. No daemon, recurring monitor, scheduler/source rewrite,
release or independent forecast is created. The overall research goal
remains incomplete until full matched barriers, finite G3 and the frozen
independent prediction/strong-control gates have actual evidence.

## Actual accepted graph and validation

HF preflight completed2026-10-11 00:04:13CST, checking all338immutable
runtime Python files, the original pilot/binary,36geometries and8exact
native32MPI/rawEFS/sixphysical-input caches. All28moving geometries are
fresh. Preflight launches zeroDFT; concurrently running parent SCFs are
not described as zero-cost or completed. The original six E056 starts and
all previous production sources remain unchanged.

| Registered +1% channel | Accepted handle | Actual dependency | State at00:08:31CST |
|---|---|---|---|
| PO to T | 28709788 | afterok28692775+28692776 | PENDING/Dependency |
| PO to M | 28709789 | afterok28692775+28692776 | PENDING/Dependency |
| Pattern-preserving flip | 28709790 | afterok28709788+28709789 | PENDING/Dependency |
| Pattern-reversing flip | 28709791 | afterok28709788+28709789 | PENDING/Dependency |

All four were actually accepted00:05:58CST. Independent `scontrol` and
all-eight-handle sacct/squeue reads confirm1node/32CPUs/1task/4hfor each
new job,2RUNNING+6PENDING study handles and the exact two-wave dependency
graph. No extra registered G2 start remains without a handle; this does
**not** imply eight completed chains, converged barriers or material
predictions. The helper is single-use and must not be replayed on HF.

`preflight.json`, eight cache audits, exact executed helper,
`submission_journal.jsonl`, `submission_receipt.json` and
`queued_jobs_verification.json` retain actual timestamps, accepted IDs,
raw scheduler fields and hashes. Only13declared files were exported in a
restricted tar (SHA256f858f1dccca6493c0e6aec8970d82dd2bc5bfbca3895bff7c4581c5b2c1e6be9);
no PP/orbital, charge/wavefunction or other licensed file is included.

Initial70guard/material cases pass before adding predecessor-control checks;
77then pass locally42.31s and in actual clean Git-archive cwd47.40s. Tree
1409207c7bdcac5a044acd66828b13a18645da03, archiveSHA256
7e4e65ba7832ae2b6708c3245206e5d95722134d772be0153af20e1879a5ee67;
helper/test/README bytes match. Actual-delivery fixtures are verified in
separate subsequent reports. These are not a fresh full-project or HFpytest
claim. The existing production source's earlier full1459pass/2skip regression
and last E059 full1515pass/2skip regression are separate historical evidence.

Final actual-delivery tests pass79cases locally43.247s and in the explicit
independent full-archive cwd49.464s, with zero failure/error/skip. All21
then-present tracked case/test files match archive bytes. Tested delivery
tree5410b33ead6a3be1722e107f2fd5cb4abb94fa26, archiveSHA256
4adca7418798ebe27960aeb681a6ea3abb015b2e316846c9f8b16a0e77fb0cee.
Subsequent changes only add records/docs; the executed helper and scientific
runtime are unchanged. A final00:20CST eight-handle sacct read still shows
2RUNNING/6PENDING with no terminal failure. Existing user dirty files remain
un-staged and unchanged. This turn isPROGRESS, not overall goal completion.
