# E066: terminal-audit guard for the existing G2 queue

Slurm `afterok` proves successful execution, not ordinary NEB convergence:
the registered runners also exit successfully at their update limit. Before
the flip segment's expected terminal time, two still-pending successors were
held so that their parent results can be audited and an existing-chain
continuation safely inserted if needed. This guard is an operational step,
not a new scientific acceptance rule or a new material result.

At **2026-10-11 02:54:40 CST**, the tested single-use helper held exactly
`28692775` and `28692776`. A separate 02:55:27 read-only snapshot confirms
both are PENDING/JobHeldUser/Priority=0, with their original AND dependencies
`28661019` and `28722320` unchanged. Both parents remain RUNNING/32 CPU in
hfacnormal01. The four +1% jobs stay pending downstream. No running job,
source, physical input, independent chain, DFT call or holdout was changed.
Other projects are not mutated. This is not a cancellation or an additional
allocation.

Release is deliberately separate: first audit the actual terminal segment,
then either verify a required same-chain continuation/dependency handoff or
prove no continuation is needed. Do not release based solely on exit 0 or
rerun `hold_waiters_HF.py`. The journal records each intent before its one
write; ambiguous dispatch and stale acknowledgement fail closed and require
read-only reconciliation, not repeated mutations. Identity, original command,
output paths, resources, live parents and exact dependencies are checked.

The helper SHA is
`f231b82c3826af4c0b457d5913541df2c65947a21e4938114b0507753d7e7832`;
`executed_hold_waiters_HF.py` retains those actual execution bytes. The full
raw scheduler snapshot, journal and receipt are preserved. The first 22
offline guard tests passed locally and in a separate complete Git archive
before execution (0.16/0.18 s). These are mocked scheduler/refusal tests,
not DFT or fresh HF pytest. Subsequent receipt replay and related tests are
reported only after running them in `validation.json`.

This fixes a real queue handoff gap but does not substitute for the still
missing complete G2 barrier matrix, finite G3 conditional branches,
strong-control material forecasts or independent tests. The original five
user modifications remain untouched. No new Hessian/phonon budget, CI,
release, PyPI publishing or new monitoring automation is introduced.

Final related local and independent-archive tests each pass **82**
(10.94/11.48 s), with zero failed/skipped. Eleven new guard/source/receipt
files byte-match tree `05f35938b93e600d8de29ddf5ca1b6273cad7186`, archive
`47c983782bdd0505ad807258778ab5bd8a34255fdfaec5139d541477465eed45`.
Later edits append records/docs and copy test output only; the tested helper,
library and production drivers are not changed. At 02:57, the latest complete
M step is 7/0.120235 eV/Å and flip step18/0.272906. Neither is yet converged.
