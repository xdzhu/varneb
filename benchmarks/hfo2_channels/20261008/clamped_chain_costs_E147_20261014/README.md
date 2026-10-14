# E147: cumulative completed-chain cost, not a speed-up benchmark

Seven ordinary-converged G2 chains are reconstructed from their complete
registered production histories, including all capped continuation segments.
The original sixteen Slurm allocations were queried through ssh hf on
2026-10-14. The command/stdout and UTC are retained in native_accounting.json;
all are COMPLETED0:0, one node/32CPU/hfacnormal01/iai806. Start/end, ElapsedRaw
and CPUTimeRAW agree. CPUTimeRAW is allocated core-seconds, not a CPU profiler.

| Chain | Segments | Updates | Fresh interior SCFs | SCF core-hour proxy | Allocation core-hours |
|---|---:|---:|---:|---:|---:|
| PO-to-M, zero strain | 3 | 43 | 308 | 310.297723 | 312.168889 |
| PO-to-T, zero strain | 1 | 10 | 77 | 92.211809 | 92.720000 |
| Reversing flip, zero strain | 1 | 5 | 42 | 47.667216 | 47.884444 |
| PO-to-T, +1% strain | 1 | 5 | 42 | 51.689291 | 51.751111 |
| Preserving flip, zero strain | 5 | 105 | 742 | 791.163330 | 794.782222 |
| PO-to-M, +1% strain | 2 | 32 | 231 | 258.431562 | 259.502222 |
| Preserving flip, +1% strain | 3 | 67 | 476 | 415.495917 | 417.555556 |

Total1918fresh interior SCFs,1966.956847629call-time-proxy core-hours and
1976.364444444allocation core-hours. These include49initial interior
evaluations plus7times267updates. Initial/final endpoint caches and nine-frame
restart cache accesses are not additional SCFs. Different new source calls
are counted separately even if their geometry/log hashes happen to agree.
The +1% preserving terminal's119calls/99.555915core-hours alone would omit
357earlier calls/315.940002core-hours from that same production path.

The existing backend-independent ResponseEvaluationCost and
recorded_response_dataset_cost API is reused, not replaced. The new small
material adapter verifies byte-pinned original per-call receipts, original
six-input declarations, exact7-image call grids, full native segment identity
and serial costs within allocation time. It does not re-audit undistributed
raw SCFs/assets, rerun a calculator, read running jobs or select a G3 model.
The all-call native E/F/stress and genuine DSIZE32 audits are prior evidence,
not silently relabeled as fresh local raw-data audits. Original mixed E054/E055
runtime receipt provenance is retained, not equated by ignoring code hashes.

Run from repository root with an unused output path:

```sh
python -m scripts.audit_hfo2_clamped_chain_costs --output /path/to/new_report.json
```

costs_final.json is the current exact-grid export; costs.json retains the
initial same-number export before adding the strict image/call grid check.
Neither modifies the old receipts. Local related checks:42pass/0fail/error/skip
(24new material/fault checks and18existing cost API checks). Independent
archive replay/testing is recorded only after execution.

Actual independent42pass/0fail/error/skip. Tested tree506390df64fb82f51d6094c3acdd9062224833a0,
tar SHA76404675f3317561b4b97ad658eb43c462ef52961da8adfd42dec76d6b92c27c.
All6323tracked files existed before tests;16receipt/accounting files were
byte-identical and both source texts normalized-identical before execution.
The CLI actually regenerated the entire cost export in that extracted archive.
An initial byte-comparison correctly refused equivalence: only the recorded
executed-script SHA differed. The two original script byte hashes are retained,
their LF/CRLF-normalized content agrees, and every other output field agrees
exactly. No physical input, native receipt or result was normalized or rewritten.
The independent export is retained verbatim; no new full-suite/HFpytest,
paper/PDF, physical parameter, production source or scheduler change occurred.

See validation.json and the original local/independent JUnits for all source
hashes. Only explanatory manual/README, copied independent reports, validation
and ledger follow the tested snapshot; adapter, tests and source receipts stay
the tested versions. At16:48:52CST the separately running28852667 remains
RUNNING/node75/32MPI/hfacnormal01, step2/.106872, stderr0. It is intentionally
not included in a completed-chain denominator. The paper/PDF and five
protected user edits remain unchanged in this cost-accounting milestone.

This is a complete account of the declared completed-production subset,
not the whole study. Endpoint relaxation/cache preparation, G1/initial-seed
preparation, reference phonons/polarization, earlier diagnostics, the unfinished
+1% reversing chain, and future conditional/failed-start/prediction costs are
explicitly outside this table. A cost-complete flag inside this subset does
not turn excluded work into zero. Different channels/seeds and restarted FIRE
histories are not a controlled optimizer comparison. No acceleration,
independent prediction, TS accuracy or G3 gate pass is claimed; G2 remains7/8.
