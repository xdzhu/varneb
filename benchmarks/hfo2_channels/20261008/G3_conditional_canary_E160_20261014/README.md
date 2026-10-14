# E160: actual same-input HfO2 conditional optimization benchmark

Both E159canaries finished normally (COMPLETED/0:0). They use the same
original ABACUS/PBE/100Ry/full10auDZP/Gamma2³/P0/E0 contract, cached centers,
fixed Q=(0,0), .03eV/Angstrom released-gradient norm target, safeguarded BFGS,
.05Angstrom initial trust radius, .5Angstrom internal-domain bound and
1.6Angstrom minimum distance. Ordinary NEB remains.10eV/Angstrom; no CI,
new Hessian center, cutoff change or physical+.005holdout access.
The retained axes are the frozen E159canary-only coupled internal eigenspace,
not named parent phonons or a final forecast domain. Cell energies refer to
the12atom Hf4O8cell (four HfO2units), not one formula unit.

| Existing center | Identity initialization | Curvature initialization | Final seeded released norm/eV Angstrom−1 |
|---|---|---|---:|
| T,28859983/node40 | 12new SCFs, budget exhausted | 2new SCFs, target reached | .02629995 |
| preserving,28859804/node52 | 12new SCFs, budget exhausted | 5new SCFs, target reached | .01959108 |

Thus both curvature-initialized arms reached the common target within the
same12attempt cap, whereas neither identity arm did. The completed results
are17+14=31new SCFs, not the48attempt upper bound. No SCF failed; a conservative
optimizer-budget stop is not a backend failure. The identities' final cached
trials are NOT accepted optimizer iterates or converged candidates.

Do not report12/5or12/2as measured speedup factors: the identity convergence
cost is censored. Native SCF-call-time core-hour proxies are T13.23573/2.17610
and preserving9.94616/4.02135(identity/seed). Full allocations, which contain
both arms plus overhead, cost16.24/14.75556core hours for T/preserving. The
original296curvature SCFs and earlier G2/center preparation are ADDITIONAL,
not free or included in those per-arm proxies. No net end-to-end advantage,
new optimizer invention, equal final branch or universal acceleration is claimed.

## What is genuinely checked

The separate native reader uses the unchanged378Python-file E159runtime,
exact configs and original native curvature report pins. It reparses every
new SCF, matches geometry/atom order, original six asset hashes, DSIZE32,
electronic convergence, raw E/F/stress, physical gradient, reservation digest,
conservative counters, fixed control and Q. Accounting is queried by the
original accepted Slurm handle, not inferred from an absent queue row.

It then reexecutes the pinned optimizer strictly on the original ordered
cached E/g. Duplicate-coordinate cache hits and exact-center reuse are
replayed; a missing/differently ordered point raises WITHOUT a calculator.
Actual full native replay matches both arms' calls, center hits, geometry
screens, stop reasons and final coordinates at both centers. This replay is
required in the pinned native numerical environment; it does not authorize
approximate-geometry cache reuse on another BLAS/runtime.

The new auditor is outside both production runtimes, in a fresh E160audit
directory; no running-source edits or repeated DFT. Its native source SHA is
51445a26007b225220807d84258fb9220bdc7672f38413fea8fea4b28b2cd6d3.
`*_native_audit.json` are exact native reports; `*_worker_outputs.zip` preserve
the original42/36worker files, including all17/14result JSONs, reservation
markers, contracts, summaries and outcomes. Decompressed bytes are checked.
No pseudo/orbital payload or raw/transformed author reference arrays are included.

## A physical diagnostic, not a new barrier or TS certificate

At the preserving center, conditional release lowers energy by
.0494849300eV/cell and moves.31058Angstrom in the declared internal metric.
The retained two-coordinate gradient is about.00012134eV/Angstrom; the
released gradient is.01959108. Removing only the rigid-translation gauge for
ANALYSIS, its common-chart distance to the original T center decreases from
.36375to.15906Angstrom. Relative to the separately relaxed T candidate the
distance is.19178Angstrom, not zero; their sampled energies differ by
.0008134484eV/cell, with no resolved physical error bound.

This shows sensitivity of the sampled bottleneck geometry/energy to
transverse equilibration and motivates the two local conditional faces.
It does NOT prove channel merging, phase identity, a new activation barrier,
release-block stability or a first-order TS. Do not replace the G2full-path
barrier table with one relaxed local point. The remaining nonlinear faces,
strong controls, independent tests and frozen predictions are still required.

Reproduce the ordered-chart diagnostic from the repository root:

```sh
python benchmarks/hfo2_channels/20261008/G3_conditional_canary_E160_20261014/analyze_preserving_release.py
```

The script compares original atom order without rotation/rematching; reports
the translation-free common-chart metric separately from affine-referenced
ordered minimum-image atomic distances and cell Frobenius norms. It does not
alter structures/forces, impose ASR, run DFT or infer polarization.

Native completion audit entry point: scripts/audit_hfo2_G3_conditional_canary.py.
Use the pinned E159source via PYTHONPATH, exact E159native curvature report,
the selected channel and a fresh output OUTSIDE production/source/assets.
An unfinished/missing handle is not restarted or labeled a failed calculation.

First16replay tests had15passes/1fixture failure: the zero-gradient synthetic
origin expected one callback but the solver deliberately performs a final
verification, hence two. Only the fixture count changed; native data, solver,
physics and final tolerances did not. The original failure report is retained.
Final local/independent related results and source identities are in validation.json.
