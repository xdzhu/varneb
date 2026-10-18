# E269: does the hafnia midpoint depend on the seed path?

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
