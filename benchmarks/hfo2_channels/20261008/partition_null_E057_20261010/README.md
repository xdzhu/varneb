# E057: full-joint harmonic partition null and real-chart design counts

This check prevents a misleading method claim: when the same complete local
quadratic and mechanical space are used, changing the Q/R partition cannot
improve its exact stationary response. It does not show that all frozen,
conditional, anharmonic or material models are equivalent.

Five retained subspaces of dimensions1/2/3, rotated/unrotated charts, index0
and index1 centres give114stationary and54paired-response checks against
independent original-full-space solves. Two invalid B4 release blocks are
correctly refused; B5 promotes their omitted negative direction and returns
the same full stationary result. Frozen versus joint energies differ. These
are implementation residuals, not DFT uncertainty or prediction accuracy.

The10existing registered Hf4O8 clamped starter geometries independently
give33translation-free atomic plus3open-cell internal directions and one
controlled-strain direction. Complete central-paired actions therefore need
37columns:75uncached evaluations at one step,149at two. This is an explicit
stencil design count, not a new allowance, TS/Hessian certification or the
minimum possible cost of every conceivable method. A frozen2Dplane does
not by itself cover those37directions. Conditional optimization evaluations
and geometric overlap must be recorded separately rather than guessed free.

## Decision

Use one actual bottleneck dataset for all valid T/Cmma/reference partitions,
not duplicate Hessian jobs per gauge. Do not expect an exact harmonic
partition to win on accuracy. Test the originally proposed independent
material response against strong references, and seek nonharmonic/branch or
equal-accuracy cost evidence where an actual incremental claim is possible.
The direct full-space solve is an equal-information harmonic reference even
when a particular stable-R partition fails. A refused Schur partition is not
evidence that this stronger, still computable reference predicts poorly.
The finite G2 matrix must still finish; G3 probes require explicit budget,
sampling, error and centre registration. The +0.5%holdout stays untouched.

Run from the repository root, with a new output namespace:

```sh
python -m scripts.check_response_partition_invariance --output /new/check.json
python -m pytest -q tests/test_response_partition_invariance.py tests/test_nested_response.py tests/test_stationary_gap.py
```

`analytic_local.json` records the original working-checkout numerical result.
The selected80tests pass locally; independent clean-source and HF checks are
separate receipts, not implied by this local run. No production module,
DFT INPUT/KPT, optimizer, running job, published historical source or material
geometry is changed. No new DFT, calculator invocation, submission or monitor.

The independent archived source tree23c7b3a passed80selected tests in3.01s;
the working checkout passed80in2.29s. HF then replayed the114stationary,
54paired checks and all10geometry-chart counts with exactly the same eight
audited runtime/helper source bytes. Its existing Python3.10.9/NumPy1.26.4/
ASE3.23.1b1environment was not changed; HFpytest is neither installed nor
claimed. `executed_source/`, `verification_receipt.json`, `analytic_HF.json`
and the local/clean receipts preserve those scopes. Two fixture checks were
subsequently added, so final delivery checks are separate from the80case run.

The first HF checksum guard correctly rejected an archive still uploading:
64,716,800of89,917,440bytes, wrongSHA, and output namespace absent. The SCP
process was still live; it was allowed to finish, not restarted. Only after
the complete archive matched its exactSHA was the same unused namespace
created. This transport-order failure made zeroDFT/calls/mutations and is
not an electronic failure or a reason to duplicate production.

At22:10:29CST actual jobs28661019/28661020 remain RUNNING/32CPUs, no true
`band/vcneb_failure.json`. Latest complete optimizer rows are flip step5/
0.725003andMstep8/0.403503eV/Angstrom. Mean new SCF timings176.643/116.739s
give about5.15/2.72h to the20step segment cap, aroundOct11 03:20/00:55,
not guaranteed convergence. These are live log observations, not complete
raw-chain/TS certificates. No additional allocation or monitor is created.

Final independent archive-cwd delivery passes82selected tests in3.26s,
zero failures/errors/skips. All20then-present case files and the scientific
script/test/two manuscript sources match archive bytes. The final related
scope is not a claim of a new entire-project/HFpytest run. `validation.json`
pins the tested trees, archives, JUnit and HF observations; later additions
only record these results and do not change executable code or data.
