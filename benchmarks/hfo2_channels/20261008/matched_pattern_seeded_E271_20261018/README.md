# E271: bounded fixed-pattern numerical repair

This is a T-only numerical continuation of unavailable E270 treatments,
not another material, path search, physical parameter choice or model test.
See [the frozen protocol](protocol.md). The new worker imports the original
immutable E270 worker and E268 runtime; neither running bundle is edited.

The source of the initial matrix is the already measured zero-strain T
curvature at0.01A. `prepared.json` contains its transported block-Jacobi
numerical matrices for the same three charts. Positive release restrictions
were checked at both old probe amplitudes, but only0.01A is admitted here.
The actual DFT forces/stress retain their coupling. No clipping, damping,
new Hessian or different electronic input is introduced.

Current stage: preparation only,0 new DFT calls and0 submissions. Original
E27028947338 must be COMPLETED0:0 before preflight/held admission. At22:13
it was still RUNNING: zero-strain A and midpoint J had exhausted their12
attempts, midpoint A had passed, and the+1% treatments were still pending.
Successful arms will be reused; failed/previously skipped arms receive at
most12 new attempts once. All old outcomes and preparation costs remain.

`repair.py prepare --output FRESH.json` creates inert matrices locally.
`repair.py preflight --bundle-pin PIN` on HF verifies the terminal original
job, all three frozen sources, positive numerical restrictions and reused
native A candidates without DFT. `run` requires its immutable bundle and
the exact accepted single-node32CPU Slurm allocation. It never submits.

29 scoped tests passed (including existing E270 tests). A synthetic quadratic
check verifies constraint/seed plumbing, not DFT acceleration or stability.
There is no release/CI and no change to the paper's original holdout scores.
