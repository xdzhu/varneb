# E062: audited continuation of the same clamped PO→M chain

This case retains the common T-plane substrate at training strain zero,
P=0, 12 ordered atoms, seven moving images and two fixed cached endpoints.
ABACUS uses the original six physical files, 100 Ry/full 10 au DZP,
32 genuine MPI ranks and one thread. No CI or electronic-parameter tuning.

## Actual terminal segment, not a converged MEP

Parent job `28661020` completed with exit `0:0` at 2026-10-11 00:49:09 CST.
Its 20-update geometry continuation reached the step cap, not the ordinary
0.10 eV/Å target. Same-boundary native E/F/stress replay gives:

| Local segment step | Maximum NEB force (eV/Å) | Discrete forward peak (meV/HfO₂) | Reverse peak (meV/HfO₂) |
|---|---:|---:|---:|
| 0, previous chain step 10 | 0.615276004 | 114.427812 | 206.968306 |
| 20, cumulative chain step 30 | 0.161868645 | 104.109701 | 196.650196 |

These peaks are **unfinished observations**, not converged barrier labels,
certified TS energies or a channel ranking. The complete 140 new interior
SCFs were checked on HF against the six original physical-file digests,
native `DSIZE=32`, electronic convergence and finite full E/F/stress.
Their summed transport wall time is 16801.571431 s, or 149.347302 core-hours
at 32 ranks; this is not CPU profiling. Export and cache preflight used no
new DFT calls and prohibited external execution.

`observations/step_0000` and `step_0020` include exact POSCAR lifts,
evaluated trajectories and portable native INPUT/KPT/STRU/log/audit data.
Proprietary pseudo/orbital/binary files remain on HF. Full six-file physical
verification happened on HF; local replay checks only the exported bytes
and native E/F/stress, not absent proprietary files.

## Same-chain restart and scheduler handoff

`seed/` contains the audited identical step-20 geometry and nine exact
hash-pinned E/F/stress caches. After an interior moves, its cache invalidates
normally. Fixed endpoints remain cached. FIRE velocities/time step are not
restored: this is explicitly a geometry continuation, not an optimizer-state
checkpoint. The segment allows at most 20 updates/6 h and keeps all numerical
and mechanical parameters unchanged.

Exactly one job was accepted: `28722320`, initially held at 01:08:52 CST.
It was released at 01:17:31 and accounting independently confirmed RUNNING
32 CPU at 01:18:27, start 01:17:32. Production still uses E054's immutable
`source-fixed`, previously fully regression-tested; no running source was
overwritten. This adds zero independent chains to the eight-chain G2 matrix.

Independent 01:20:30 startup check finds the step-zero residual exactly
reproduces the audited seed to log precision (0.161869), empty stderr tail,
and the first fresh image-1 native SCF has `DSIZE=32`. It does not claim that
this new SCF has finished. The other running flip is at step14/0.220938.

Both waiting jobs `28692775/28692776` now depend on the existing flip
`28661019` **and** this continuation. The four +1% starters retain their
original downstream dependencies. At most two study allocations can run;
other projects are untouched. A successful step-cap exit is not physical
convergence. The +0.5% holdout remains unopened.

The first immediate post-submission query could not see the accepted held
handle. An immediate post-update query subsequently returned stale dependency
metadata; later independent snapshots showed both writes had succeeded.
This is evidence of read-after-write visibility delay, not an ABACUS failure
or proof of the cluster service's internal root cause. No allocation or
dependency write was repeated. Slurm's canonical comma-separated `afterok`
AND form is now accepted; OR/afterany/unrelated dependencies remain rejected.
Bounded read acknowledgements fail closed; separate verified release only
acts after both actual prerequisites are visible. All attempts, partial
journals and executed source versions are retained. **Do not rerun these
single-use dispatch scripts.**

## Delivery validation

The initial guard test incorrectly searched comments for `while`; it was
replaced by an AST check (one failed/38 passed retained). Initial clean
packaging omitted Git-ignored trajectory/log files (two failed/75 passed
retained). Scoped `.gitignore` exceptions and raw-byte attributes fix actual
delivery, without weakening geometry, physical or numeric gates.

Earlier related local and fixed independent-archive runs each passed 77 tests.
Final delivery results and archive identity are recorded in `validation.json`
only after execution. Current scientific runtime is unchanged; this is not
a fresh full-project or HF pytest claim. No new G3/Hessian budget or paper
claim is inferred from this incomplete chain.

Final local and independent full-archive cwd runs each pass **82 tests**
(54.80/53.96 s), including actual native material replay and dispatch receipts.
149 source/material/receipt files byte-match the tested delivery tree
`bdf47db73ac2d3e094747403efc32d1319eae83b`, archive SHA256
`8c2d64dc4ec2d0271d0581b0f745ab8f11150fabc3087d6476ffb3814b040624`.
Evolving JUnit files are explicitly excluded from the source-byte comparison.
Later changes append records/docs only. Original five user changes stay intact.
