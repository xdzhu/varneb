# E075: fresh-starter terminal audit and one-slot admission preparation

This stage is **preparation**, not a completed PO-to-T result or an actual
release. At the 2026-10-11 12:06:46 CST native read, PO-to-T 28692775 and
preserving continuation 28723655 are both RUNNING. T step 6 rebounded from
0.103118 to 0.143524 eV/Angstrom; it is not stopped for that single rebound.
The existing reversing starter 28692776 and four +1% starters remain held.

## Correct fresh-SCF accounting before inspecting the new starter

The prior terminal auditor was exercised on geometry continuations, whose
nine step-zero images are cached. A fresh starter has only its two endpoint
caches: seven step-zero interior SCFs must also be counted. The new check
requires the byte-registered production recipe, exact runtime factory and
matching cache policy, then counts `7*last_step + initial_fresh_interiors`.
`initial_fresh_interiors` is 7 for a new starter and 0 for a continuation.
Unexpected extra or missing calls still reject the audit for investigation.
The native copied starter preflight proves the former policy; no nearby,
permuted or unevaluated geometry is treated as a cached result.

Only the read-only analysis source is changed/deployed. Original production
source, all six physical files, ordinary0.10/P0/clamping/spring/FIRE settings,
running jobs and raw outputs are not edited. Existing M receipts are not
rewritten; their 91 calls remain correct for 13 continuation updates.

## Later admission, deliberately not executed yet

`release_existing_reversing_once.py` will only release the already registered
28692776 after a positively completed, genuinely ordinary-converged T audit,
all fresh native calls/source/inputs/EFS checks, 338 pending runtime hashes,
the original pending seed and exact endpoint preflight. Native accounting
must corroborate one live preserving continuation; the queue must prove one
used slot and all five protected waiters. All four downstream holds must
remain intact. No submission, dependency rewrite, running-job mutation or
automatic retry is present. An existing journal prohibits replay; ambiguous
writes require native reconciliation. Admission is not startup evidence.

Native audit support was deployed to the fresh namespace
`hf:/public/home/iai806/abacus/agent-runs/20261011-varneb-clamped-T-terminal-E075-source-r1`.
Four copied source hashes and the starter runtime-preflight hash match local
bytes exactly. This does not claim the audit/admission was run. Current
`live_optimizer_1206.log` is a read-only nonterminal snapshot, not a full
SCF audit or a converged material label. No holdout was generated or read,
no G3 centre selected, no new independent chain or Hessian budget added.

144 related tests pass locally (12.10s) and in a complete independent Git
archive (14.86s), direct exit0. They cover old material replay, correct initial
call counting, invalid cache policies, native accounting, exact capacity,
synthetic one-release/ambiguous-write safety and the actual copied starter
metadata. Mock admission is not native release evidence; no HF pytest or
new material DFT is claimed. The full scientific goal remains incomplete.
