# E269: does the hafnia midpoint depend on the seed path?

**20:20 latest:** PO-to-T28946181 completed0:0 at20:12:08, ordinary fmax0.082374;
preserving28946180 remains RUNNING at step3/fmax0.140626. Only T has a terminal
comparison. No replacement paired-gap or final preserving sensitivity is inferred.

Actual HF deployment, 2026-10-18. This is the first finite batch of the
[mechanism-closure plan](../../../../docs/HFO2_MECHANISM_CLOSURE_2026-10-18.md).
The two opposite-strain initializations test the same already-seen +0.5%
condition; this is not a second independent prediction test.

## Actual status

| Job | Channel | Actual start / node | Resources / bound |
| --- | --- | --- | --- |
| 28946180 | preserving reversal | 19:27:47 CST / node60 | 1 node, 32 CPUs, 80 updates / 24 h |
| 28946181 | PO-to-T | 19:27:47 CST / node60 | 1 node, 32 CPUs, 80 updates / 24 h |

At 19:33:20 CST both were RUNNING. Preserving had two completed interior
SCFs and T one, with the next call in progress. Their first completed calls
passed native electronic convergence, genuine DSIZE=32, unchanged six physical
asset hashes and complete energy/forces/stress replay. First-call times were
114.589 s and 148.130 s, respectively. This is a startup/prefix audit, not
an all-call or converged-band claim. See [the original snapshot](first_calls_snapshot.json).

With these few samples, a seven-interior evaluation is roughly 12--18 min
plus overhead; whole-path convergence cannot be estimated before its first
residual sequence. Each path evaluates its interiors serially with native
mpirun -np32; the two paths run concurrently. Sharing node60 does not merge
their independent32CPU allocations. No endpoint is recalculated.

The new remote root is
`/public/home/iai806/abacus/agent-runs/20261018-varneb-midpoint-reverse-E269-r1`.
It reuses the original E268 immutable runtime read-only. Local Git confirms
no vcneb, calculator or native-exporter change since the deployed source;
the original checker also verifies all375 runtime files and physical assets
on HF. The executed bundle and all eighteen input decks passed native preflight
with fresh DFT disabled, including four exact endpoint-cache uses.

Bundle SHA256: `381ae9787cca91b1435b121dafd74ee072d1aa3bd22335f05a7e7ba42bfc7ae1`.
Archive SHA256: `26439ba0d9a0e271f6564b753ed14122b9d2f92c05a706020110ddcc556ec254`.
Protocol, geometry seeds and controls were pushed in `2591d29a` before submission.
Actual held, accepted, release, preflight and prefix receipts are retained here.

## Physics and the next decision

All electronic parameters are unchanged: ABACUS/PBE, 100 Ry, full10auDZP,
original Hf/O datasets, Gamma2x2x2, Hf4O8, P=0 and zero field. The ordinary
NEB target remains0.10 eV/A, no climbing image, same FIRE and step/geometry
limits. Only the seed polygon comes from +1% rather than zero strain.
The three target endpoints and their exact static caches are identical to E268.

The [six-neighborhood table](six_neighborhoods.json) is a zero-DFT replay of
the already completed 0/+0.5/+1% preserving and T bands. It does not provide
new Hessians or change a barrier. The preserving peaks have very small
allowed tangential force, while their remaining physical gradient is mainly
transverse. At the midpoint these are -0.001176 and0.248150 eV/A for the
tangential component and whole allowed-gradient Euclidean norm. The latter
is not the max-vector NEB stopping metric and must not be compared directly
with the0.10 target. Excluded substrate reactions are reported separately
by the core and are not counted as a failure of the released space.

This focuses subsequent intervention on transverse oxygen accommodation,
not merely finer sampling along the original tangent. It does not certify
a stable branch, identify a unique mode cause or establish a TS. Both seed
outcomes, not just the lower, will be compared. A start sensitivity comparable
with the original4.30meV/fu gap improvement prevents attributing that signal
uniquely to intrinsic barrier curvature. Agreement does not supply an error bound.

## Reproduce the zero-DFT preparation

From the repository root, with the normal scientific dependencies:

```sh
python -m scripts.prepare_hfo2_midpoint_reverse --output NEW_SEED_DIRECTORY
python benchmarks/hfo2_channels/20261008/midpoint_reverse_E269_20261018/neighborhoods.py --output NEW_NEIGHBORHOOD_JSON
python -m pytest -q tests/test_hfo2_midpoint_reverse.py tests/test_hfo2_sampled_holdout.py tests/test_hfo2_physical_holdout_results.py
```

The latest scoped run passed26 tests, including the direct read-only CLI,
source/target geometry, exact resources and the six-neighborhood replay.
This is not a full repository test run or evidence of physical convergence.
Do not rerun submission when original handles or intent files exist.

The first postprocessing helper was named `inspect.py`, which shadowed the
Python standard library and failed before reading jobs; no DFT job was affected.
It was preserved as `failed_inspect_shadowing.py.txt` on HF and replaced by
`inspect_paths.py`, outside the frozen executed bundle. A subsequent local
direct-CLI path issue was corrected and regression-tested. The production
controls, source, input files and exact original handles were not changed.

## Complete initial evaluations, not terminal results

Both original jobs completed the first full seven-interior evaluation and
continued optimizing. Initial force residuals were0.135999/0.151237eV/A,
above the unchanged0.10 target. The native first snapshots are exported
read-only under `initial_complete_chains`; all18frame geometries/EFS are
replayed, using14already-completed new interiorSCFs and four cached endpoints.
The export launches no new DFT and does not modify a job or replace a label.

The initial archive SHA256 is
`5aa870888e5bc49da5388ac430ca47e89a439d213b23027d419b64f1cae0b2d0`;
all98payload files have separate hashes. The initial collector and terminal
collector reside in the separate HF `postprocess_v1` directory, not the running
bundle or old375-file runtime. The latter collector refuses live handles.

Preserving's initial peak is138.889103meV/fu,40.358226 above E268's98.530877;
T's is90.503363, only0.085185 above90.418178. The preserving peak also changes
Qx from-0.798917 to-1.039388A and open-y stretch from-0.416145 to+0.654708%.
This is a simultaneous atomic/cell seed contrast, not a unique energy partition
or a final sensitivity after equal relaxation. Do not assign the old midpoint
nonlinearity uniquely to intrinsic barrier curvature before E269 finishes.

[The stage figure and source tables](../../../../paper/VARNEB_JCTC/figures/hfo2_midpoint_initial_20261018/README.md)
retain both optimization stages explicitly. The updated38scoped checks passed,
including real initial-data replay and rejection of promotion to terminal paths.
The figure was visually reviewed, with editablePDF/SVG, larger normal-weight
panel labels, aligned axes and the author's framed legend/top-right axes style.

At19:52:09CST both exact jobs remained RUNNING, each with the complete step-zero
record above; step1SCFs were in progress. No failure, restart or parameter change.

At19:55:26CST the same handles remain RUNNING. Preserving completed step1
with fmax0.136367 after0.135999 at step0; this small rebound does not stop
the path. It has15completed interiorSCFs and T13, with more in progress.
The terminal collector and comparison tool are ready but have not produced
or claimed terminal results. No new job or calculator setting was added.

## Minimal mechanistic follow-up prepared, not submitted

The [frozen atomic–cell intervention](FACTORIAL_INTERVENTION.md) now has a
tested preparer. It preserves internal fractional patterns while exchanging
the two candidate cells' open row. Two existing diagonal statics per channel
are reused; only two crossover SCFs per channel would be new. Missing energies
remain null, so no interaction or causal attribution is reported yet.
Four diagnostic native decks were actually generated in ignored build/ with
zero DFT/submissions and the original INPUT/KPT. They are not run-ready on HF
and do not bypass review of the terminal E269 comparison. The default preparer
requires terminal results; initial snapshots need an explicit diagnostic flag.
No new material, backend, Hessian, NEB chain or model search has been opened.

At20:03:20CST both original jobs were RUNNING. T completed step1 at19:58:35,
with residual falling from0.151237 to0.123360eV/A. Preserving still had its
completed step1 at0.136367 and continued step2SCFs. Neither was converged.

## Completed T response, 20:12:08 CST

The original PO-to-T handle28946181 completed0:0 after44m21s. The two FIRE
updates used21new interiorSCFs, retaining the same two endpoint caches.
The full native audit rechecked all21calls, every original physical asset,
genuine32MPI, complete E/F/stress, exact accounting and the unchanged runtime.
The terminal snapshot is retained alongside all21call logs/inputs/receipts.

| Initialization at the same +0.5% target | Sampled peak (meV/HfO2) | Ordinary fmax (eV/A) | Updates |
| --- | ---: | ---: | ---: |
| from0%, E268 | 90.418178 | 0.080825 | 2 |
| from1%, E269 | 90.286243 | 0.082374 | 2 |

The difference is-0.131935meV/HfO2. Both sampled peaks are close at the
declared ordinary resolution; this is not a continuous-barrier error bound,
a unique-MEP proof or an independently trained predictor. The peak Qx values
remain-0.725241/-0.741484A, with open-y stretches-0.605986/-0.502644%.
Similar peak energies do not imply identical configurations. The original
E268 prediction scores remain unchanged; this is a follow-up initialization test.

New actual cost is23.242526native-call proxy core-hours and23.653333allocated
core-hours, additional to all earlier E268/training costs. Identical21call
counts are not a demonstration of optimizer acceleration.

[Native completed T data](completed_channels/T/audit.json),
[the separate T comparison](T_terminal_comparison/comparison.json), and
[the export receipt](T_terminal_export.json) preserve this result. The archive
hash was checked before safe local extraction. The46scoped tests pass,
including replay of every21native call and refusal to invent a missing P
outcome or draw the two-channel comparison before both paths complete.

`compare_paths.py --channel T` exports this single-channel result without
computing a new T-minus-preserving value (`paired=null`). Its default remains
the complete two-channel comparison. The running preserving path is not
restarted or altered. At20:20 it had28completedSCFs and step3 residual0.140626.
