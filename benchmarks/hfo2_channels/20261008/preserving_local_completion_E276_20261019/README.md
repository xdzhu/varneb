# +1% preserving local-control completion

## Completed outcome

Both original handles completed with exit 0 on 2026-10-19:
minus 28954771 at 10:59:50 CST, plus 28954772 at 11:01:15 CST.
The 27 new SCFs (13/14) all have independently replayed native energy,
forces, stress and geometry. Allocation cost: 27.813333 coreh;
completed-call proxy: 26.824676 coreh. Prior work/cached starts are additional.

| Original sign | A lowering | Additional J lowering | A gradient | J gradient |
|---|---:|---:|---:|---:|
| minus |17.172494|0.484821|0.028933|0.022660|
| plus |17.174898|0.484106|0.018045|0.010588|

Energy units: meV/HfO2; gradient units: eV/A. All four held physical observables
remain fixed; A retains the entire F cell. These are qualified local
gradient candidates, not stable conditional minima or stationary TS.
The old 0.5 A searches remain failures; new chart norms 0.753--0.755 A lie
inside the explicitly revised 1.0 A window. Two probes at each A have
positive energy/gradient secants along the followed direction only.

In the common original-T metric, F-to-A atomic norms are 0.747--0.749 A;
the transverse structural proxy changes by +0.206/-0.205 A while the
geometric Qy remains fixed. Those are different observables (rank 4 to 5).
Oxygen contributes 77.9% of squared displacement, not energy. J contracts
open y by 0.310--0.321 percentage points relative to F, in contrast to the
zero-strain expansion. Most observed local lowering is already atomic.

The +1% J points sit 107.765--107.767 meV/fu above IS, still about 16.1 above
the existing free-band peak 91.699. That free peak also changes Qx by +0.229 A
and has transverse proxy 0.123 A rather than about 0.206. Thus the local contrast
supports internal accommodation but cannot alone reproduce the whole-path
strain response with the original descriptors fixed. No additive
mode-energy partition or new global barrier ordering is inferred.

`completed/` holds both original native archives/receipts;
`complete_analysis/` independently reparses every new call and candidate;
`mechanism_analysis/` supplies all four 0/+1% signed rows and free-path context.
`figures/matched_controls/` contains main Figure 10 (PDF/SVG/PNG) and QA.
The original numerical windows differ across strain and are disclosed.

Reproduction from repository root; outputs must be fresh directories:

```sh
python benchmarks/hfo2_channels/20261008/preserving_local_completion_E276_20261019/analyze_terminal.py --output /fresh/native-replay
python benchmarks/hfo2_channels/20261008/preserving_local_completion_E276_20261019/mechanism_data.py --output /fresh/mechanism-data
python benchmarks/hfo2_channels/20261008/preserving_local_completion_E276_20261019/plot_mechanism.py --source /fresh/mechanism-data/analysis.json --output /fresh/figure
```

The mechanism builder uses the retained complete analysis in this case
and validates its actual archive/source hashes; alternate output paths
do not redirect or replace that source. No new DFT or scheduler write.

## Actual admission and historical startup prefix

Original handles28954771(minus) and28954772(plus) were submitted held,
checked for one-node32CPU/one-task/4h/hfacnormal01 resources, then released
on2026-10-19 around10:34CST. Both run onnode117 with separate32CPU
allocations. Their immutable execution source is ffa272da;25-file bundle
pin`16bd55d10c00d579026f6aff11e4b472db07d4526c177e0f0e8a6a0ae43422d1`.
Root:`/public/home/iai806/abacus/agent-runs/20261019-varneb-preserving-local-E276-r2`.

The native zero-SCF preflight replays F and both exact old starts without
new electronic calls. `admission/` retains preflight, accepted held-resource
rows and release records. The first new native call of each sign has
passed six-asset/MPI/geometry/EFS checks:32MPI,115.843/113.790seconds.
Norms.541657/.541656A are outside the old artificial guard but inside the
new declared window; released gradients drop to.234430/.234468eV/A,
and energy changes from F are-13.943955/-13.943086meV/fu. This is a
validated startup prefix, not qualified A/J, a batch outcome or stability.
The first namespace's configuration-check error was stopped before any
allocation/SCF; see`preflight_r1_failure.json`. No deployed source was overwritten.

Actual terminal exports, not forecasts, will be processed by
`collect_terminal.py` and independently replayed by`analyze_terminal.py`.
They reject live/missing/timeout handles and preserve the known old failure.
The new collector runs as a separate postprocess sidecar, never edits
runtime sources or issues scheduler writes. Collection v1 passed native
replay but failed a source-hash Path/string conversion after copying files;
its partial exports remain untouched. Corrected sidecar v2 used a fresh
export namespace, regression-tested source hashing, and no fresh SCF.
Both v2 receipts bind their collector hash. The immutable DFT runtime is unchanged.
The completed result is integrated into main Figure 10 and the article.

## Preparation and diagnosis

The [new finite protocol](protocol.md) follows explicit author authorization
to fill the previous missing A/J intervention with identical electronic
settings. The preparation records below remain distinct from the actual
admission and native startup records above.

`original_failure_diagnosis/` independently replays all24failed-arm native
SCFs. Near the0.5A guard, BFGS backtracking collapses the step while both
radial and tangential physical gradients remain. The budget failure was
not electronic nonconvergence and is not proof that no candidate exists.

`registered_connection_scale.json` uses one endpoint-verified integer
offset for the whole existing band, then projects its peak to the old
four held observables and fixed F cell. Its0.810A scale supports extending
the numerical exploration window to1.0A. No energy is assigned to that
projection. `connection_scale.json` retains the explicitly unregistered,
inadmissible diagnostic; it is neither used for the bound nor as a seed.

`prepared/` contains the exact original F, original four covectors and
numerical matrix, plus the two original last native evaluated trials.
Those are continuation starts, not qualified old A candidates. Each has
at most32A+2checks+16J fresh calls in a new isolated runtime. The old
protocol/results remain untouched; the new physical Hamiltonian is identical.

This addresses the same scientific question as the completed free bands:
which internal accommodation remains under the clamp, and how much of the
local lowering is already atomic rather than requiring an open-cell change?
The free-band result cannot replace this fixed-observable comparison.
