# E101: same-reference analysis of the audited preserving flip

Apply the existing displacement, dual atomic-gradient, external-work and
signed-path-work analyzers to the new admitted ordinary zero-strain preserving
band28727382/step15. Terminal audit SHA256
46145da14f24b2b8bd46aa4656f0982650d6b2dc661282f58d809aab8a6d3202
pins E100's original-input/native32MPI/full-EFS/334-file terminal audit.
Local analysis reparses only nine portable native frames, not undistributed
orbitals or every prior SCF. No new DFT, Hessian, G3/reference selection,
strain condition, holdout label or source-mode array is generated/exported.

The original T chart, rank3 geometric triplet, two complete rank4 optical
doublets and all four preregistered rank4 Cmma frames are kept. Atomic modes
and permitted cell strain remain separate. Geometry and energy covectors use
the prescribed proper cyclic frame and fixed ordered integer lift; there is
no fitted frame, atom remapping, ASR repair or input/production modification.

## Figure contract before rendering

Core conclusion: the two registered ordered zero-strain flip paths have
different sampled energy and reference-coverage profiles under one mechanical
and electronic contract; coverage alone is not energy or forecast accuracy.
Archetype: quantitative comparison grid. Existing Python/matplotlib workflow.
Target: JCTC working manuscript, not a claim of submission readiness.
Final size:183x195mm; six aligned panels, nonbold13pt(a)--(f), all four spines
and top/right ticks, larger readable labels, outside white0.85framed legends.
Panels(a,b): raw energy profiles, preserving/reversing, one native PO+ zero
and common energy scale. Panels(c,d): common-T squared atomic displacement
coverage. Panels(e,f): dual atomic-gradient coverage in that same chart.
Each chain keeps its own normalized arc; equal abscissae are not matched
configurations. All four Cmma frame traces are retained, not an error band.
Hero evidence: nine native samples per path and unchanged terminal audits.
Control evidence: complete preregistered spans and explicit omitted/absolute
gradient norms in source CSV. Statistics:18frame records, not18independent
SCFs/replicates; shared initial cache counted separately. No confidence test.
Integrity: raw samples with straight connectors; no smoothing, crop, copied
published panels, basis-array redistribution or certified TS annotation.
Export: PDF/SVG/PNG/600dpi TIFF, full scalar CSV, hashes and actual visual QA.
Reviewer risk: different ordered endpoints/path candidates do not prove the
global minimum switching mechanism; percentages do not partition energy;
coarse signed-work discrepancy is not a barrier uncertainty estimate.

Fresh reports and actual tests/build/visual checks are recorded after running.
Older analyses/figures and the user's dirty working files remain unchanged.

## Actual result and reproduction

The existing four analyzers produced displacement.json, gradient.json,
external_work.json and path_work.json from the admitted nine native frames.
The unchanged source registration, seed, Gamma and author-source hashes are
recorded in those reports. Original author mode arrays remain outside Git;
their pinned source is needed to repeat projections, not to redraw the figure.
All four analyses add zero DFT, scheduler writes, model selections or labels.

The preserving sampled forward peak is96.828183meV/HfO2, ordinary residual
0.098816809eV/Angstrom. Peak T-triplet squared displacement/gradient fractions
are82.6301%/12.9583%; common atomic-gradient norm0.070196452 and omitted
norm0.065490597eV/(sqrt(amu)Angstrom). All four Cmma frames remain, with
14.8336--14.8337% displacement and52.7203% gradient coverage. None of these
fractions is energy accuracy. The peak/endpoint signed-work residuals are
14.523948/0.032398meV/HfO2, not a physical error bound.

Figure7 in the working manuscript compares this band with the reversing
candidate's338.968965meV/HfO2 sampled peak. Eighteen displayed frame records
contain17 distinct native-log/STRU identities; across five work reports,
45 records contain42 identities. These are provenance counts, not replicates.
The preserving/M sampled difference1.629939meV/HfO2 does not establish robust
channel ordering. Three +1% training candidates, conditional-branch/stability
and prospectively frozen strong-control comparisons remain necessary.

From repository root, redraw into a fresh directory:

```sh
python -m scripts.plot_hfo2_clamped_flips --output /fresh/figure-directory
```

The dated figure folder in paper/VARNEB_JCTC/figures contains PDF/SVG/PNG,
complete source_data.csv and automated qa.json. The command also makes a
regenerable600dpiTIFF; do not overwrite historical outputs. Actual visual
inspection is recorded separately in figure_review.json, not retroactively
claimed by the automatic exporter. The24-page reading build is rendered
with review_pdf.py; local PDF rasters are regenerable/ignored, not publication
data. Its whole-page contacts and changed pages15/16/21--24 were inspected.
Page22 retains a large blank region because the unbreakable figure follows
on page23: acceptable for this reading draft, not final journal typesetting.

The actual verified local and independent r2 archive runs both pass109 tests,
zero failures/skips, with one existing fontTools deprecation warning each.
Only verified_local_tests.xml and verified_independent_tests.xml identify
the final delivered code; other JUnit files preserve intermediate attempts.
The first archive's sole mismatch is source-hash metadata for historical
vcneb/biaxial_curvature.py: localLF and archivedCRLF have identical normalized
code. Exact known byte variants and normalized code hashes are now required;
every scientific field still compares exactly, and raw physical-data hashes
are never normalized. The initial reversed newline-pin diagnosis is also
recorded, not silently erased. Original reports and both archives are kept.

Validation, source hashes, actual PDF/visual checks, native status and
independent delivery checks are recorded in validation.json. The actual
24-page final PDF has no TeX overfull/underfull/glyph/LaTeX/package warning;
external Perl locale warnings remain. No new full suite or HF pytest result
is claimed.27initial science/figure/test/evidence files compare byte-equal
to the tested archive before later evidence/README additions.

At23:40:25CST28727368/28709790 are RUNNING on hfacnormal01, each32CPU,
at new-segment steps16/.168904 and4/2.467425eV/Angstrom; stderr bytes0.
28709791 remains JobHeldUser,0allocatedCPUs. No submit/release/replay,
dependency/input/source change or other scheduler write occurs in E101.
The .10 ordinary target,100Ry/full10auDZP inputs and original user edits
remain. The goal is active and scientifically incomplete, not blocked by
expected live calculations and not completed by these descriptive results.
