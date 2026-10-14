# E168 — diagnose the bounded conditional-relaxation hole without new SCF

Previous goal turn E167 was PROGRESS, pushed as b39db373: completed local
atomic/cell and fixed-quadratic tests entered the connected paper. This case
addresses the next observed training-domain failure, not a new material,
Hessian center, path, electronic setting or independent-prediction result.

## Question

Preserving point13, Q=(-0.2,-0.1) Angstrom, exhausted twelve attempts without
a final gradient candidate. Its original record contains48 geometry checks,
35 rejections and12 complete SCFs, rather than a backend crash. The counts
alone do not identify which guard fired or whether the actual residual could
be reduced within the original chart domain.

This retrospective diagnostic asks which original guard short-circuit caused
the rejected steps and whether the completed trials approached the registered
0.5-Angstrom total-internal-displacement boundary. The separate minimum
distance is1.6 Angstrom; the released-gradient norm target remains0.03
eV/Angstrom. That is not the ordinary NEB0.10 residual. No policy is changed.

## Method and provenance

The already tested E163 standalone auditor is imported with the immutable
E161 runtime on PYTHONPATH. It rechecks all380 Python hashes, original config,
accepted job handle, Slurm resources, cached native origin, every complete
prefix SCF's E/F/stress/assets, fixed Q and substrate, attempts and final
candidate identity. The unfinished current node is not accepted as complete.

Each complete nonconverged prefix node receives the original exact-order
cache-only optimizer replay. An observer calls the original guard, records
its outcome, then independently classifies the same short-circuit order.
A mismatch fails the diagnostic; the original decision is always returned.
The temporary hook is restored even on exception. No calculator is available
to satisfy a missing cache point. Replay must reproduce original coordinates,
callback and guard counts, final candidate or precise stop reason.

Reasons distinguish nonfinite coordinates, the original total-radius bound,
short-distance geometry, invalid chart/cell, and accepted proposals. Distance
is not evaluated after an earlier radius rejection, matching the original
guard: those proposals cannot be declared physically safe from this report.

Every complete cached trial additionally reports its internal radius and
margin, released-gradient norm, signed energy change from the original center
and the algebraic radial indicator -2r dot g_R. A positive indicator means
the released steepest-descent direction initially increases squared radius.
It does not prove a constrained KKT point, an unreachable unconstrained
minimum, a unique branch or a stable conditional minimum. Trials are not
silently relabeled accepted optimizer iterates or final candidates.

## Decision boundary

A domain-guarded hole must remain a hole. It is not corrected by projecting
away the radial residual, changing the gradient target, inserting the last
trial into a contour, enlarging the live job's domain or retrying its budget.
The original two allocations continue. The prepared E166 frozen controls
remain subject to the normal parent-completion and two-allocation admission
gate. Any later training-domain revision must be declared prospectively,
retain these failed observations and obey the existing finite budget.

Original ABACUS/PBE/100 Ry/full10-au DZP/Gamma2x2x2/P0/E0, all six input
hashes, ordinary NEB0.10, cached origins and physical +0.005 labels remain
unchanged. This audit reads no private author arrays. No new DFT, scheduler
write or mutation of live source occurs. The full goal still requires native
complete-face audits, stable branches and strong-control independent tests.

## Reproduction

```text
python -m pytest -q tests/test_hfo2_G3_bound_diagnostic.py tests/test_hfo2_G3_face_audit.py
```

`diagnostics.py` contains only observation/arithmetic, not a new optimizer.
`native_driver.py` requires the dedicated fresh external HF report path and
exact helper source hashes. Native execution must use the original E161
runtime, not the newly edited checkout. Saved native report and receipt bind
the helper, existing auditor and original runtime. Analytic mock tests are
not material SCFs; native evidence is recorded only after actual execution.

## Native result and next discriminating action

The real HF run audits14 completed preserving nodes, including the cached
origin,13 candidates and one hole. All60 existing nonzero-Q prefix SCFs are reparsed under
the original six assets and genuine32MPI; all13 nonzero-node traces reproduce
their original order and outcome from cached E/g. The original380 Python
hashes and immutable policy also pass. Native report SHA:
7c9581eb0c102f6110415ca590b7f045b4875e3d3498a7af5019d0cb8a14f57b.
The tested helper snapshot tree314cc28e was deployed outside the live study.

For point13, all35 refusals occur at the original total-radius short-circuit;
13 guard calls pass, and12 SCFs are complete. No rejected proposal reaches a
DFT call. The final cached trial has radius0.4998853735 Angstrom and remaining
margin0.0001146265, with released gradient0.0735963713 eV/Angstrom. Its radial
gradient component is-0.0136123401, while its constant-radius tangent gradient
is0.0723265516 eV/Angstrom. The radial decomposition is a fixed-metric scalar
identity computed from the native replay, not an additional force calculation.

| Trial index | Internal radius, Angstrom | Released gradient | Tangent gradient | Energy change, meV/cell |
|---|---:|---:|---:|---:|
| 0 | 0.354451 | 0.307094 | 0.307093 | -52.890984 |
| 5 | 0.397032 | 0.033791 | 0.031118 | -58.986418 |
| 11 | 0.499885 | 0.073596 | 0.072327 | -60.475669 |

The table shows the first, descriptively minimum-gradient and last completed
trials; all12 remain in `analysis.json`. Gradient units are eV/Angstrom,
energy is relative to the original nonstationary center per twelve-atom cell,
not a barrier or per-formula-unit energy. No trial is promoted to an accepted
iterate or final candidate. The minimum observed norm0.0337907153 also exceeds
the original0.03 target. The original optimizer continues through the norm
rebound and stops at its registered twelve-attempt cap.

Decision: the literal cause of the rejected steps is established, but failure
of every feasible relaxation direction is not. The substantial remaining
tangent gradient makes feasible-direction step handling a specific next
optimizer hypothesis within the unchanged domain, rather than justification
to enlarge the radius or weaken convergence. A subsequent method experiment
must use matched cached starts, common targets/domain and explicit attempt
costs; neither a tangent projection alone nor an analytic toy is a material
speedup. No such continuation job is prepared or submitted in this milestone;
the already planned E166 controls retain priority and their resource gate.

`analyze.py --output a-new-output.json` reproduces the arithmetic from the
pinned native report without raw SCF reparse. Analysis SHA:
a4adb82578102dcf987e648d9c8002cfb847a01132de4e3fc5f43e13b7f34935.
These observations do not certify a stable branch, boundary KKT point,
unreachable unconstrained minimum, TS, full surface or independent forecast.

## Verification outcome

Before native execution,41 related checks passed locally and41 in the first
independent archive. Final50 local/50 independent checks include21 new case
checks, native observations, scalar radial/tangent identities and all three
table rows. The two-module scope is not a full-repository regression. Final
tested tree:3bef044dfe823286b4ab2417667bbe703c89a292; archive and exact source
hashes are retained in `validation.json`. Deployed observer/driver bytes are
unchanged from the first tested archive. The independent arithmetic CLI
reproduces the entire saved analysis byte-for-byte. No new material run,
source repair, physical holdout label or extra training reservation is hidden
in these checks.
