# E182 — integrate the fair complete-reference comparison into the paper

Previous goal turn E181 was PROGRESS, pushed1f50ddc8. This turn integrates
already executed E180 evidence, not a new DFT result or another reference
search. Reading E159 first confirmed that native retained/complement sign
tests were already present; those measurements were not repeated.

## Manuscript argument and evidence

The same complete H/g/Q problem is invariant to its orthonormal reference.
It cannot acquire a prediction advantage by renaming coordinates.

| Paragraph job | Actual change | Evidence / boundary |
|---|---|---|
| Method | Same-observable affine KKT equation and null-space minimum criterion | Standard constrained quadratic; zero raw floor is not physical resolution |
| Result | New3.8/Table8, complete T33+cell3 vs coupled2+released34 | E180's original two centers/two steps/all25Q, max2.64e-16eVcell and1.56e-14A numerical equality, not DFT error |
| Result | Same finite-domain shape discrepancy,24T/16preserving nonzero candidates | Twelve-atom-cell energies; seen origin alignment, not independent accuracy |
| Comparison | AppendixD.7/TableA9, all seven references and four Cmma frames | Empty all-reference nonzero eligible cohort; no ranking by unequal subsets or extrapolated values |
| Interpretation | Discussion separates basis identity from nonlinear/branch evidence | Does not complete native-chart comparisons, stable branches or prospective material validation |

The final methods/results/discussion use the existing terminology ledger,
with three explicitly defined new control terms. No abstract or conclusion
asserting unfinished results is added. No new external source/citation is
invented; all13external citations,9figure paths/bytes, all119previous table
lines in original order and the exact TeX wrapper are preserved.
There are12new quantitative rows plus4header/rule lines:135table lines total.
Paper rawSHA92f152399040a3b70a6596625a042a05e5f50b19f96d385a46591c1a3b043de0.
Original E180 sourceSHA39e31d790592e7a2f640a8428f00efecc28fd14775b7fb4b2b699b7da54c6402
is pinned, not regenerated or refitted.

## Verification and actual rendering

113local/113independent related checks across9modules,20new,0fail/errors/skips;
not full regression or new material validation. Three existing editorial
checks were updated only for the legitimate eighth Results subsection and
two new tables: old119rows remain byte-identical/in-order, with all prior
source/figure/citation checks. Mutations of every new row and false source
coverage/zero-error/certification claims are refused. Independent CLI output
is byte-identicalb3fb3b652d61b6eb854ad6200d3f9a630ad69c1a2eb855968bc2aa52a0a7469d.
Final testedtreeb745699e8768771abbac7370e726406606fe37d1,
archive5dbe91fcbb274603e0b13335fc4eb6cec0786995a54482ba73ea293e3c9e7b7d.

Existing TeXLive/LuaLaTeX compiles the actual Markdown and nine actual figure
files,34pages,0overfull/underfull/missing-character diagnostics. All first
34page renders were inspected; observed Table6/8caption page splits
were fixed without altering numbers or the wrapper. Against the already
reviewed second build, final rendering preserves22pages pixel-byte-identically;
all12changed pages were re-inspected in
contacts, and method8/results20–21/appendix34received detail checks.
The first PDF, render metadata, CLI outputs and successful113-check reports
remain in Cscratch, not overwritten or presented as the final build.

Final internal reading PDF, not a submission-ready article:
`C:/Users/zhu/.codex/tmp/varneb-E182-manuscript-20261015-r3/varneb_JCTC_draft.pdf`,
SHA9d52dba9b9e1988c5a0ee1777395b2816290072eead97382f3a0fc62e95741f1.

## Research state and remaining input

07:59:54CST sacct/squeue confirm **28867177and28868827bothPENDING**,
AssocGrpCpuLimit, no assigned node/elapsed computation/new completed labels.
Two committed allocations remain, no resubmission/cancellation/source overwrite.
No DFT, raw-SCF reparse, physical inputs, radius, threshold or reservations
changed this turn;+.005and8independent face points each remain sealed.
No private author arrays were accessed for this paper integration.

Missing material evidence remains actual B3 treatment labels, released-base
transverse curvature, stable nonlinear branch/domain, numerical resolution,
appropriate strong controls and prospectively frozen independent predictions.
The full JCTC research goal remains active/incomplete. No new material,
whole path, CI, installation, release/PyPI, monitor or submission.

Reproduce from repository root into a **new** output file:

```text
python benchmarks/hfo2_channels/20261008/manuscript_shared_controls_E182_20261015/check_tables.py --output fresh_check.json
python -m pytest -q tests/test_hfo2_manuscript_shared_controls.py
```
