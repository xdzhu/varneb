# E271: bounded fixed-pattern numerical repair

## Complete native result, 22:53:48 CST

Job28947723 is COMPLETED0:0 after27m30s. All14 new native calls passed
six-asset, real32MPI and complete E/F/stress replay. All six F/A/J release
candidates now meet the unchanged0.03eV/A released-gradient target;
the midpoint A is the exact original E270 result, not a new call.

| Strain | E_A-E_F | E_J-E_A | E_J-E_F | New A/J calls |
| --- | ---: | ---: | ---: | --- |
| 0% | -0.057642 | -0.119733 | -0.177375 | 2/4 |
| +0.5% | -0.008858 | -0.077234 | -0.086093 | 0/3 (A reused) |
| +1% | -0.166696 | -0.073329 | -0.240024 | 4/1 |

All energy differences are meV/HfO2. These small same-pattern local
accommodation energies concern native sampled peaks with different held
pattern values across strain. They are not fixed-Q strain derivatives,
corrected barriers or stability/TS certificates. Full allowed gradients
at J are0.02744/0.10813/0.16338eV/A; only the32-dimensional released
subspace was the optimization target.

Additional actual cost is14.266900 native-call proxy and14.666667 allocated
core-hours. Including the preserved E270 failures gives55 new SCFs,
51.536867 proxy and52.373333 allocated core-hours; all older curvature/path
preparation costs remain additional. The14-call repair alone is not an
end-to-end acceleration benchmark. The123-file native export is retained
as `native_terminal.tar.gz`, with its receipt, without scattering raw files.

## Actual admission, 2026-10-18 22:26 CST

Original E27028947338 completed0:0 with41 new SCFs and explicit incomplete
arms. The zero-DFT preflight then passed all three charts, numerical positive
restrictions and the exact reused midpoint A's native E/F/stress. Five arms
can receive fresh calls:0% A/J,+0.5% J,+1% A/J (up to60 total; J requires A).
Midpoint A is reused, not recomputed. See `preflight.json`.

Slurm returned28947723. It was submitted held, inspected as1node/32CPUs,
1task/32CPUs per task,hfacnormal01/6h, then released. At22:27:27 it was
RUNNING onnode3; the preserving E269 original handle remained the only
other committed study allocation. `held.json`, `accepted.json` and
`released.json` preserve the exact handle/resources. This is launch
evidence, not proof of improved convergence or a completed F/A/J treatment.

At22:31:25 two fresh SCFs had completed with no failure file. Zero-strain A
then passed the unchanged target in2 attempts: released gradient0.02101510
eV/A and E_A-E_F=-0.05764191meV/HfO2. Its final native candidate was replayed
locally, four retained-variable defect<1.5e-18A and cell coordinates exactly
zero; full allowed gradient remains0.10574105eV/A, so this is not a full TS
or minimum certificate. `startup_receipt.json`, `zero_A_prefix/` and
`zero_A_review.json` retain the partial result. Only its final one of2 calls
is reparsed locally here; the full batch/cost export awaits completion.
This repairs one actual identity-seed budget failure, but does not establish
a general speedup or erase the prior Hessian-preparation/failure costs.
Forty scoped tests now pass, including the complete339-file original batch
fixture and the new real zero-strain A candidate. They do not certify the
unfinished treatments or the full six-neighborhood physical mechanism.

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

Historical preparation stage: inert bundle deployed,0 new DFT calls and0 submissions. Original
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

Execution source977b9b6f was pushed before deployment. The fresh HF namespace
is `/public/home/iai806/abacus/agent-runs/20261018-varneb-matched-pattern-E271-r1`.
Its four-file bundle and shell syntax passed; see `package_receipt.json`,
`executed_bundle_files.json` and `deployment_receipt.json`. Bundle pin is
`bcdb328562237d84522ddc5aa597b291ef5ba622716907e6692cd530cbda32f6`.
This is not the terminal-baseline/native preflight or admission to submit:
E270 was still RUNNING in the deployment receipt. Do not overwrite this
bundle, create a speculative dependent job or bypass the original-job gate.
