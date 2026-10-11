# E074: four native pending-task holds preserve the two-allocation cap

The +1% E060 dependency DAG named the original zero-strain segment IDs,
not the later preserving-flip continuation 28723655. Once the two old parents
exited, two +1% starters could otherwise start beside that live continuation.
This is a scheduler-capacity risk, not an electronic or material failure.
At discovery there were only two live study allocations; the reversing
zero-strain job was already user-held, so no third study allocation started.

At 2026-10-11 11:54:54–11:55:23 CST the exclusive single-use transaction
verified owner, exact task names, partition, 32-CPU resources, immutable
production script, output namespaces and all original dependencies. It then
held exactly 28709788, 28709789, 28709790 and 28709791, each once. Native
readbacks prove PENDING/JobHeldUser/Priority=0 and unchanged dependencies.
The fsynced journal retains the full before/after native scontrol output.
No running job, dependency, input, production source or unrelated job changed;
no sbatch/release/DFT/holdout operation occurred.

Next admission is explicitly slot-based: audit a finished material segment,
prove a free slot with native queue/accounting, then release at most one exact
registered eligible waiter per spare slot. Do not automatically release all
four after a parent finishes. Existing afterok dependencies remain additional
eligibility requirements, not proof that continuation capacity is free.

Offline tests cover exact writes, wrong identities/states, exclusive replay
rejection and an ambiguous write with no repetition. A separate actual-data
test verifies all four native readbacks and the unchanged dependencies. Passing
mock tests is not evidence that a real hold occurred; `native_transaction/`
is that evidence. This script must **not** be rerun with another receipt path:
on an interrupted operation, inspect the existing journal and native state.

Actual local and independent complete Git-archive runs each pass 9 tests
(0.26/0.22 seconds, directly captured exit 0). Tree, archive and native timing
are recorded in `validation.json`; no HF pytest or additional DFT is claimed.
The 11:57:11 native queue still shows exactly two running study allocations,
and all four protected +1% starters plus the old reversing job are held.

This operational repair does not complete G2/G3, certify a saddle, freeze a
prediction or replace the full HfO2 scientific objective with scheduler work.
