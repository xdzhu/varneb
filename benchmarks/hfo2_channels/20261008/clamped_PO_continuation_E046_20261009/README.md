# E046: one audited clamped PO+ geometry continuation

Prerequisite E045 job28446324 completed0:0 at22:07:34CST, five actual SCFs
with original six physical-file bytes and genuine32MPI. Its four BFGS steps
retained Pca2_1 at all three registered tolerances but did not converge:
atomic max0.037333755eV/A, open traction6.948494kbar. No phase stability,
spontaneous polarization or G2 path conclusion follows from that canary.

This experiment continues **only this zero-strain common-T-substrate PO+**.
The terminal scf_000004 geometry/E/F/stress is reused after exact ordered
geometry, six physical bytes, STRU and raw-log checksum checks. No free-cell
evaluation is substituted and no atomic reorder/rotation/symmetry restoration
occurs. A geometry change invalidates the normal ASE cache, then the original
fixed-input calculator performs a new SCF in a fresh directory. The BFGS
Hessian is freshly initialized, not resumed or counted as optimizer-history
reuse/an acceleration benchmark.

- At most20new BFGS steps/20fresh SCFs; one cached initial evaluation.
- Same .03eV/A atomic,2kbar open-traction and .02maxstep gates.
- P0/E0, PBE/100Ry/full10auDZP, Gamma2x2x2; all six physical input bytes fixed.
- One HF/hfacnormal01 node32MPI/directmpirun, one thread per rank,2hour cap.
- At most two concurrent study allocations; no matrix, holdout, CI or restart
  submitted automatically. A complete step-limit segment is not a converged
  endpoint. An actual phase change is reported, not repaired into PO.

`seed_local/` is a zero-DFT preparation using the exported observable files,
not a full-six-byte local proof because basis/pseudopotentials are omitted.
The runtime must regenerate a fresh `seed_HF/` after full physical-byte replay
of the original canary. Those parent files and its production source are
read-only. Receipt hashes and actual job states are added only after observed
execution. No result is preregistered as positive.

The case's narrow Direct/Hf4O8 STRU reader supports portable raw replay with
stock ASE; it rejects other units/masses/species/mobility/header formats.
It does not advertise a general ABACUS I/O implementation. The generic public
clamped endpoint entry remains calculator-independent.

## Observed result, 2026-10-09

Job28453655 completed0:0 at22:47:31CST. Only two new SCFs/two fresh BFGS
steps were needed: terminal atomic force0.022284634eV/A, open
traction1.928434990kbar, Pca2_1 at all three registered symmetry tolerances.
The endpoint passes the declared physical screen. Its energy is
-9782.973830941437eV/Hf4O8; this is not a Hessian-minimum certificate,
electronic-polarization measurement or a switching/escape barrier.

The two SCFs used236.309s/2.100523core-hours; the allocation used244s/
2.168889core-hours. The cached initial SCF belongs to E045 and is not charged
again or treated as evidence of a faster optimizer. The original runtime
source and parent canary remain unchanged.

`seed_HF/`, `cache_preflight.json`, `completed_HF/endpoint/continuation_audit.json`
and [the completion receipt](submission_completion.json) preserve the actual
HF full-six-byte evidence. The26 exported observable files permit local raw
E/F/stress replay without redistributing pseudopotentials/orbitals/charge.
Local replay explicitly reports the omitted-basis scope rather than claiming
those absent physical files were locally checked. Registered T, M and PO-
endpoints are tracked in the [finite endpoint matrix](../clamped_endpoint_matrix_20261009/README.md).
