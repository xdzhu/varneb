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
