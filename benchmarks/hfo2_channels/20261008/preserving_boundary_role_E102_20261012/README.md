# E102: the same group label has a different sampled energetic role

Bounded question: does the ordinary zero-strain clamped preserving band retain
the low-energy central region of the already audited free-cell preserving band?
Use only those two completed nine-frame bands; no DFT, scheduler mutation,
G3 bottleneck/reference selection, new phase or held-out label is required.

## Figure contract before rendering

Conclusion: both central samples have the same tolerance-stable Pbcn label,
but one is a sampled energy trough and the other a sampled high point; the
free-cell path's substantial T-chart x stretch is excluded by the clamp.
This describes these paths, not certified minima/saddles or a causal barrier
derivative. Archetype: quantitative two-panel grid, existing Python/matplotlib.
Panel(a) shows raw energy relative to each path's own PO+ well; panel(b) shows
F_xx in the same original T chart. Both use each path's own normalized arc;
equal abscissae are not matched geometries. Original electronic contract
unchanged; free-cell versus clamped is NOT one mechanical boundary.
Final size183x88mm, aligned axes and labels, nonbold13pt(a,b), no panel titles,
all spines/top-right ticks, larger10pt labels and outside white0.85framed legend.
Eighteen existing records are neither new SCFs nor independent replicates.
Native samples/straight connectors only, no fitted/smoothed MEP or error band.
Export PDF/SVG/PNG/600dpiTIFF plus CSV and traceable JSON; final visual QA
recorded only after inspection. Symmetry uses the existing structure audit's
fixed .001/.01/.05Angstrom and1degree sweep, never tolerance selection,
structure standardization, atom reordering or oxidation-state guessing.
Reviewer risks: space-group identity is not a basin identity; a sampled
trough is not a stable phase; a sampled peak is not an index-one TS; nine-point
maxima are not error-bounded barriers. The existing free-path six-point
sampling supplement's higher38.508448meV/HfO2 union maximum is retained
separately, not substituted for the original nine-point profile.

Actual results, tests and delivery evidence follow execution.

## Actual material result

The existing exact observation readers and periodic structure audit are
reused, not duplicated or altered. E045's pinned HF receipt attests the
historical free-band raw audit; locally its ordered snapshots/cache/EFS and
production residual are replayed, not its undistributed original raw logs.
The clamped receipt is rechecked and nine portable native SCFs are reparsed.
Both use the original electronic six-file contract, with different mechanical
boundaries and PO+ initial wells. Original generalized cell scales remain
5.131328250206249/5.12968067458423Angstrom respectively; they are not retuned
for this figure. T geometry uses the prescribed inverse cyclic rotation only.

Both central samples are Pbcn at all three registered tolerances. Free/clamped
central energies are-14.8382007842/+96.8281830028meV/HfO2 relative to their own
well. Their centre-minus-adjacent-mean signs differ. The free centre has
F_xx1.1220219827 and volume147.2026537Angstrom^3; the clamped centre has
F_xx1,volume134.9025096. These are different structures, not the same frozen
point's response. The new band lacks the sampled central trough; no claim
is made that every possible clamped Pbcn branch is absent. Full branch
stability and TS certification remain unmeasured. The first rendering attempt
correctly rejected an incorrectly assumed common numerical cell metric;
the fix retains each production metric instead of changing either run.

Section3.10/Figure8 are added to the working manuscript. The actual PNG has
been inspected: aligned axes/labels, readable unbold(a,b), all spines,
outside framed semiwhite legend, no point/label obstruction. No image copied
from literature, no geometry/basis export, no smoothing or covariance band.
Reproduce the analysis/figure in a fresh directory from repository root:

```sh
python -m scripts.plot_hfo2_preserving_boundary_role --output /fresh/figure
python -m pytest -q tests/test_hfo2_preserving_boundary_role.py
```

The Python runtime retains installed pymatgen2025.10.7/spglib2.7.0/ASE3.28.0.
No upgrade is performed. Dataset counts are observations, not independent
replicates. Original five dirty user files and historical source/figures remain.
Actual test, PDF build/visual checks, native status and clean-archive receipts
are added after verification; the full study goal remains active/incomplete.

## Verified delivery

The117 focused tests pass locally and in a fresh clean Git-tree archive,
zero failures/skips and one existing fontTools warning per run. Eight science
files (code/tests/main manuscript/analysis/CSV/vector/raster) compare byte-equal
to that actually tested archive. This is not a new full-suite or HFpytest run.
The25-page reading PDF compiles with existing LuaLaTeX; final TeX diagnostic
count is zero, while external Perl locale warnings remain. All25-page contacts
and final24/25 detail pages are inspected. Figure8's caption is now complete
on the same page and new math subscripts render correctly. This draft is
not submission-ready; older page22 blank space remains an explicit reading
layout limitation, not final journal typesetting.

D drive filled during r2 preview creation. Only this experiment's generated
r1/r2 review directories were moved recoverably to E; scientific/source files
were untouched. The final r3 review is generated on E. D remains critically
low in free space; later large exports/previews must use E rather than
silently deleting unrelated user data. Exact PDF/archive/data/test hashes and
actual results are in validation.json. Two original32CPU HF allocations stay
live, latest completed steps18/.129055 and8/1.680009. The existing reversing
handle stays held. No scheduler, physics, G3 or holdout change is made here.
