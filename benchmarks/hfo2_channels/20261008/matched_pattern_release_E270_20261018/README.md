# E270: matched geometric-pattern controls at T's three strain conditions

This finite batch addresses the mechanism, not a new reduced-model contest.
The existing T continuation control is complete; the preserving E269 job is
still live. This batch admits only T, using its original0,+0.5,+1% native
sampled peaks. All old path labels, forecast scores and earlier F/A/J results
are retained without replacement.

The [protocol](protocol.md) fixes three geometric oxygen-pattern amplitudes
and the structural polar z proxy during each release. Those definitions
are the same across strain, but their native values need not be equal.
Consequently the experiment is a within-peak F/A/J intervention, **not** a
fixed-Q strain derivative or a unique partition of the whole barrier.

F is the exact archived static result; A releases29 atomic directions at
the exact frozen cell; J starts from A and adds the three allowed open-cell
directions. Each arm has at most12 new attempts, total at most72, and uses
the existing safeguarded BFGS. The distinct local released-gradient target
is.03eV/Å; ordinary NEB remains.10eV/Å, with no CI. A cap/failure is a
reported outcome, not a request for more attempts or parameter changes.

## Preparation and validation

```text
python benchmarks/hfo2_channels/20261008/matched_pattern_release_E270_20261018/study.py prepare --output NEW_PREPARATION
python -m pytest -q tests/test_hfo2_matched_pattern_release.py
```

`prepared/` contains the frozen source logs/input/geometry and four physical
covectors per peak. Their action was checked against the common original-T
projection for arbitrary small atomic **and cell** changes at all three
strains. They are not Hessian eigenvectors or assigned phonon irreps.

The public two-axis optimizer is reused with two additional held directions
and three translation gauges. A uses an atomic-only36D embedding, so no
tiny SVD cell component can masquerade as fixed-cell relaxation. J reuses
the exact A coordinate/cache result; nearby distinct geometries are never
substituted. Completed A/J outcomes are written immediately and survive a
later calculator failure. Tests also cover the ranks29/32, exact-byte cache,
per-arm attempt cap, no-J-after-failed-A and duplicate submission rejection.

The original six ABACUS/PBE100Ry/full10auDZP/Gamma2³ input assets and E268
375-file runtime are reused read-only. The new bundle is separate; it never
overwrites E268/E269 sources or calculation output. No pseudopotential or
orbital binary is redistributed. Only hf/hfacnormal01 is authorized, with
one node32realMPI, at most6h, and at most two committed study allocations.

At initial preparation no new SCF, job submission or preserving intervention
has occurred. Actual preflight and Slurm handles are recorded only after
they exist; preparation is not a claim of completed mechanism validation.

## Actual admission, 2026-10-18

The prepared source/protocol was pushed in868da9e0 before execution. The
19-file native bundle passed its archive hash, original375-file runtime,
six-asset/binary checks, Bash syntax and all three cached E/F/stress/native
structure roundtrips without new SCF. Its external pin and archive identity
are in `package_receipt.json`; native evidence is in `native_preflight.json`.
Thirty-five scoped local checks passed, not the full repository suite.

Slurm returned28947338. It was submitted held, inspected as one node32CPUs,
one manager task32CPUs/task,6h onhfacnormal01, then released. At21:15:42CST
it was RUNNING onnode51. The original preserving28946180 was still RUNNING
onnode60, step7/.140321. No third study allocation or preserving intervention
was submitted. Both old running sources remain unchanged.

At that snapshot no new E270 native call had completed yet; RUNNING is not
an assertion that all calls, releases or stability checks succeeded. The
initial atomic released-gradient norms, computed from cached native data,
are.050124/.037348/.091682eV/Å for0/+.5/+1%, all above the fixed.03 target.
This is a zero-new-DFT starting diagnostic, not a forecast or result label.

At21:17:53, two new native SCFs had completed without failure files. The
first record, preserved verbatim in `first_native_result.json`, passed the
fixed-runtime validator for genuine32MPI, original six physical inputs and
complete E/F/stress. It took102.523s. Its retained-pattern defect is
4.12e-18Å and all three cell coordinates are exactly zero, confirming the
fixed-cell/physical-pattern contract on an actual DFT trial, not only a toy
function. No treatment was complete at that snapshot. `startup_receipt.json`
is explicitly a prefix observation, not an audit of the unfinished batch.

## Observed bounded failure and root-cause diagnostic

At21:39:14, the zero-strain A arm had consumed its12attempts without a
gradient-converged candidate; its J arm was correctly not attempted. Both
raw outcome files are retained in `partial_outcomes_2139/0/`, with missing
candidates kept null, never zero energy. The original job continued the
other predeclared strains;14native calls had completed and there were no
DFT failure files. This is an optimizer-budget outcome, not SCF failure,
structural divergence or termination of the entire study.

The actual first0.05Å trial raised the energy by8.595753meV/HfO2. Its
energy/force directional secants are29.5114/29.4176eV/Å². This supports a
step-scale/conditioning explanation for repeated backtracking with a unit
initial Hessian. It is not a claim that every SCF was an accepted optimizer
step, nor that an alternative has already accelerated the calculation.

The independent read-only `diagnose.py` reuses the original zero-strain
T two-amplitude matrix at this same native input. The references differ
only by a fixed, explicitly recorded integer lattice offset on a few Hf
coordinates: the native STRU hash and complete zero-coordinate gradient
agree, and the offset remains constant under atomic and cell increments.
No atom remapping, altered structure or looser electronic threshold occurs.

At both original probe steps, the29D atomic block and3D cell block are
positive in the measured symmetric matrix, whereas their coupled32D block
has one weak raw negative eigenvalue(-.003673/-.006216eV/Å²). Removing the
mixed block algebraically removes that negative sign. This is evidence
about the *measured local quadratic*, not a new physical experiment that
turns off coupling. The point is nonstationary; no physical uncertainty
bound or saddle-index certificate follows. The scaled-cell fraction of an
eigenvector is not the applied substrate strain or an energy contribution.

`numerical_seed.py` checks a possible **standard numerical block-Jacobi
initializer**, using those existing positive diagonal blocks and affine
coordinate transport to the three strain charts. The actual DFT mixed
coupling is not removed. Both original probe steps pass positive-definite
restriction checks for A and J at all three target chart definitions,
without eigenvalue clipping. This is not a strained physical-Hessian
prediction, new optimizer theory, production repair or demonstrated speedup.
The numerical arrays can be regenerated; no new calculation is admitted.
If a later targeted numerical repair is justified by the completed batch,
freeze the first probe step(.01Å) before its DFT outcomes, keep the second
as a sensitivity control, preserve all baseline failures/costs, and use a
new namespace. Do not alter the running72-attempt protocol or automatically
submit more jobs from this diagnostic.

```text
python benchmarks/hfo2_channels/20261008/matched_pattern_release_E270_20261018/diagnose.py --output NEW_DIAGNOSTIC.json
python benchmarks/hfo2_channels/20261008/matched_pattern_release_E270_20261018/numerical_seed.py --output NEW_SEED_FEASIBILITY.json
```

The diagnosis and seed feasibility are saved in `initialization_diagnostic.json`
and `numerical_seed_feasibility.json`. They add zero DFT/curvature probes;
no current input, production source, path label or forecast score changes.

## First completed midpoint arm, observed21:50:55

The+0.5% A arm reached its fixed released-gradient target in5new native
SCFs. The final raw INPUT/KPT/STRU/log/call record and accepted outcome are
preserved in `midpoint_A_prefix/` (archive SHA256
5dbbae65bd01e93a9af524c209863f6a2f09479930fe9c298df1b39b9dcb98d7).
`review_arm.py` independently reparses that **one final** portable SCF and
replays its physical coordinate gradient. It does not claim to recheck all
five attempts or locally read undistributed orbital/pseudopotential bytes.

- E_A-E_F=-0.008858301meV/HfO2; this is a candidate-energy contrast, not a
  path-barrier correction or an upper bound on further relaxation.
- Released29D gradient norm=.023923316eV/Å; complete allowed39D gradient
  norm=.129115502eV/Å. The arm is not a full stationary point or TS.
- All three cell coordinates are exactly zero; the four fixed physical
  coordinate defects are at most1.22e-19Å.
- Zero-strain A remains budget-exhausted and its J unavailable. At the
  midpoint snapshot J was still being calculated; no value is imputed.

```text
python benchmarks/hfo2_channels/20261008/matched_pattern_release_E270_20261018/review_arm.py --exported benchmarks/hfo2_channels/20261008/matched_pattern_release_E270_20261018/midpoint_A_prefix --index 1 --arm A --output NEW_ARM_REVIEW.json
```

The result is in `midpoint_A_review.json`. These partial outcomes do not
form a completed three-strain mechanism curve or a revised path table.

At22:00:47 both original handles were still RUNNING; E270 had27completed
native SCFs and no failure files. Midpoint J was unfinished and+1% had not
yet supplied a completed arm. The preserving band had reachedstep10 with
residual.133896eV/Å, still above its ordinary.10target. No job was cancelled,
restarted or given a larger budget. Fifty-four scoped checks passed across
the diagnostics, coordinate transport, native arm replay and prior mechanism/
physical-holdout analyses; this is not the full repository suite or completed
physical research. The live worker and bundle are unchanged.
