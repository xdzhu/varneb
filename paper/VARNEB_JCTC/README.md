# Hafnia JCTC study: evidence-bound working manuscript

This is a research manuscript in development, **not a submission-ready paper**.
The software paper in `../VARNEB_CPC` remains a separate artifact.

The scientific question is whether the same mechanical intervention can lower
homogeneous switching barriers without lowering the easiest examined escape
barrier from the same PO+ well, and whether a branch-aware mode-strain model
can predict that response better than strong declared controls. Neither a
compact projection nor a visually smooth landscape establishes this claim.

## Read the study

- [Main-text draft](MANUSCRIPT_DRAFT.md): connected introduction, theory,
  computational definitions, dated pilot results and their limits.
- [Detailed methods](METHODS_DRAFT.md) and [pilot evidence](PILOT_RESULTS_2026-10-08.md).
- [Research design](RESEARCH_DESIGN.md) and
  [finite execution plan](../../docs/VARNEB_JCTC_HFO2_RESEARCH_PLAN.md).
- [Prospective protocol](../../docs/HFO2_PREDICTION_PROTOCOL.md) and
  [dated v2 addendum](../../docs/HFO2_PREDICTION_PROTOCOL_V2_2026-10-09.md).

The manuscript's numerical tables are dated snapshots, not a live job dashboard.
All four final G1 candidate bands have passed the ordinary residual criterion
and the whole-candidate evidence gate has passed. The gate permits the
registered G2 experiment; it does not certify stationary saddles, exhaustive
escape coverage or distinct winding sectors. All five common-substrate ordered
endpoint representations now pass the declared atomic-force/open-traction
screens. The two minus endpoints are related by a parent-declared translation
and same-species permutation, but are not identical at ordered path indices.
Their phase equivalence does not merge their paths or establish electronic
polarity/whole-path topology. The +1% T-derived/PO+/M training endpoints
28469832/28469833/28472493 have all completed and passed their original
raw E/F/stress and input-byte reviews. The T-derived structure is Ccce already
in the strained input, not forcibly assigned the free-T group. The two +1%
minus endpoints28568571/28568573 also completed and passed all raw screens,
each6BFGSsteps/7SCFs. All ten registered endpoint representations are now
screened. Jobs28574708/28574709 were the historical first segments of the
matched zero-strain preserving-flip/PO-to-M chains, not current live handles.
The later PO-to-M segment28722320 has passed an actual ordinary terminal
audit; Tables5–6 and Figure5 document that first clamped candidate and its
mode/strain trace. The other necessary G2 candidates remain incomplete.
The reserved+.005holdout remains ungenerated/unread.
The complete G2 matched-boundary comparison, G3
conditional-branch validation and independent material predictions remain
missing. An ordinary residual pass is not a certified saddle or a barrier
uncertainty bound. No abstract or conclusion asserting these missing results
is supplied.

The [public Ma-Liu v3 supplement audit](../../outputs/HFO2_MA_V3_MAPPING_AUDIT_2026-10-11.md)
adds mapping-coverage evidence to the introduction, without equating our
candidate to a literature path, changing the finite matrix or claiming a
publisher-SI/full-preprint audit.

Its [dated reading-build receipt](MANUSCRIPT_DRAFT.Ma-mapping.validation.json)
identifies the17-page PDF, the selected actually inspected pages and the
unchanged M numbers. It is not an all-page final submission audit or a new
material result; older compile receipts retain their historical meanings.

The [native external-work trace](../../benchmarks/hfo2_channels/20261008/clamped_external_work_E077_20261011/README.md)
adds a raw-stress controlled-coordinate diagnostic to the clamped M discussion.
Its two frozen charts are not relaxed-barrier predictions or new G2 labels.
The [dated validation](../../benchmarks/hfo2_channels/20261008/clamped_external_work_E077_20261011/validation.json)
identifies the actual17-page build and selected reading-layout checks.

The [later complete-frame residual audit](../../benchmarks/hfo2_channels/20261008/residual_update_20261009_0450/README.md)
tracks the preserving flip's atomic/cell residual crossover without claiming
final channel energies. It is additional
source-verified pilot evidence, not a new stationary-TS or prediction result.

## Historical milestones

The following dated records retain the knowledge and figures available at
their own timestamps. Statements of nonconvergence there are historical,
not the current job state.

The [2026-10-09 morning update](../../benchmarks/hfo2_channels/20261008/morning_update_20261009_0850/README.md)
adds the second ordinary-converged source, PO to M, and the preserving flip's
developing split-peak profile. Both switching chains remain unconverged.
Those morning frames and the earlier step6/20/14/10 figure remain historical
evidence and are not overwritten by later values.

The [12:55 update](../../benchmarks/hfo2_channels/20261008/switching_converged_update_20261009_1255/README.md)
adds preserving step69 ordinary convergence. Table1 and the new Figure2 now
use complete step6/39/69/32 records. The central Pbcn image is lower than
PO+, but not certified stable; physical tangential forces of about
0.24 eV/Angstrom at the sampled side peaks prevent a stationary-TS label.
This is new material evidence, not a completed G1 gate or a prediction result.

## Compile the reading PDF

The scientific text has **one source**, `MANUSCRIPT_DRAFT.md`. The small
[LaTeX wrapper](varneb_JCTC_draft.tex) reads it directly with the TeX `markdown`
package; there is no second, manually synchronized scientific text.

With an existing TeX Live installation providing LuaLaTeX, latexmk, the
`markdown` package and the packages named in the wrapper, run from this folder:

```sh
latexmk -norc -lualatex -interaction=nonstopmode -halt-on-error -synctex=1 -outdir=build varneb_JCTC_draft.tex
```

The output is `build/varneb_JCTC_draft.pdf`. The `-norc` flag avoids loading
latexmk configuration files; shell escape is not required. This is a
multi-file project (Markdown and figure assets), not a standalone TeX document.
No software installation is performed by this recipe. Keep the generated PDF
in `build/`: its relative evidence links resolve against that location and
require a reader that permits local links. DOI links are ordinary web links.

Build products are ignored by Git. The source, reproduction command and
verification receipt are committed instead. PDF bytes may change with engine
versions, timestamps or fonts; the receipt identifies the actual local build,
not a promise of byte-identical PDFs on all machines.

The historical morning-update build was compiled with existing TeX Live
2023/LuaHBTeX1.16.0. All eight pages were rendered and visually reviewed with
existing Poppler; the current table's jobs, steps, numbers and pass labels
were checked against the frozen E029 data. Both unchanged figure assets and
local evidence-link targets were verified. See
[its compile/scientific-check receipt](MANUSCRIPT_DRAFT.morning-update.validation.json).
That is a reading-layout and source-data consistency check, not scientific
completion or journal approval. The earlier seven-page
[compile receipt](MANUSCRIPT_DRAFT.compile.validation.json) is retained as
historical evidence, not the hash or page count of this updated draft.

The earlier [main-text validation receipt](MANUSCRIPT_DRAFT.validation.json)
is preserved as a **historical milestone** for commit `456d376`. Its manuscript
hash predates the current math/table formatting changes and must not be used
as the hash of the current draft.

The current preserving-pass build and its actual rendered-page/numeric checks
are identified by [the dated receipt](MANUSCRIPT_DRAFT.preserving-pass.validation.json).
It includes the new Figure2 without claiming the missing material predictions.

The subsequent six-point sampling audit adds Figure3 and Section3.3:
the preserving reconstruction's highest observed sample is38.508448meV/f.u.,
not the original nine-image32.806023 value. It retains the ordinary0.10
residual target but does not certify a stationary saddle or a relaxed
continuous MEP. The [new build receipt](MANUSCRIPT_DRAFT.sampling-audit.validation.json)
identifies the corresponding PDF;
earlier receipts/figures keep their historical meanings.

The subsequent five-point audit adds Table2: T-to-PO sampling changes by only
0.011360meV/HfO2, while PO-to-M changes by11.164864. The inserted M residual
0.123408 does **not** pass ordinary0.10, so it requires bounded refinement;
these static peaks are not final continuous-MEP or TS barriers. The
[current reading-build receipt](MANUSCRIPT_DRAFT.peak-sampling.validation.json)
identifies the ten-page PDF and all-page/numeric/layout checks. Earlier build
receipts remain historical. Missing G2 and independent predictions are not
filled with placeholder success claims.

The latest source adds the [qualified G1 review and canary](../../benchmarks/hfo2_channels/20261008/G1_review_clamped_canary_20261009/README.md),
[two-SCF PO+ continuation](../../benchmarks/hfo2_channels/20261008/clamped_PO_continuation_E046_20261009/README.md)
and [single-SCF T endpoint](../../benchmarks/hfo2_channels/20261008/clamped_endpoint_matrix_20261009/README.md).
Section3.6 reports a12.359975meV/HfO2 common-substrate well separation,
not a switching or escape barrier. Source checks and the corresponding
new reading build are recorded in `MANUSCRIPT_DRAFT.clamped-endpoints.validation.json`.
That receipt refers to the first two-endpoint build. The subsequent four-well
table and its actual reading-build checks are identified by
`MANUSCRIPT_DRAFT.zero-strain-wells.validation.json`; the earlier receipt is
retained as historical evidence, not the hash of the latest PDF.

The eight-endpoint stage now includes the declared endpoint correspondence and
the three completed +1% training endpoints. Its current13page reading-build,
raw numerical/link and rendered-page checks are recorded in
[the dated receipt](MANUSCRIPT_DRAFT.eight-endpoints.validation.json).
This is a reviewed working manuscript, not a submission-ready JCTC result.

## Computation boundaries

Hafnia uses the registered ABACUS/PBE contract: 100 Ry, full 10-au DZP orbitals,
original input/KPT/UPF/orbital hashes, Gamma-centred 2x2x2, P=0 and E=0.
Ordinary NEB is 0.10 eV/Angstrom, without climbing images. Mechanical active
spaces must match across every image and endpoint of a comparison.

Finish and audit G1 candidate channels before G2; freeze predictions before
reading held-out labels. An unstable or unsupported branch is an abstention,
not an implicit successful prediction. No new functional, cutoff, phase,
domain-wall or MPB matrix is opened to fill waiting time. The execution plan
contains the finite budgets and stop/acceptance gates.

The [E078 native-gradient audit](../../benchmarks/hfo2_channels/20261008/clamped_mode_gradient_E078_20261011/README.md)
complements displacement coverage without modifying references or production
calculations. At the sampled clamped M peak, Cmma's93.51% squared atomic-gradient
coverage contrasts with36.48% squared displacement coverage. Neither is an
energy-accuracy measure. This preserves Cmma as a strong comparison and does
not preselect G3 or claim completed prospective predictions.

The [E079 PO-to-T terminal record](../../benchmarks/hfo2_channels/20261008/clamped_T_terminal_E079_20261011/README.md)
adds the second same-substrate ordinary-converged training band to Table5.
Its77 fresh SCFs and nine portable terminal frames pass their distinct native
and local checks. Its84.416624meV/HfO2 discrete peak is not a continuous TS
certificate or proof of escape ordering relative to the M candidate. Six G2
training bands, conditional branches and prospective independent tests remain.

The [E082 source-reading receipt](../../benchmarks/hfo2_channels/20261008/lee_SI_text_E082_20261011/validation.json)
records the18-page build after a partial Lee supplementary-text reading.
The introduction excludes restricted-release profile lowering as novelty;
Section2.3 distinguishes internal-response coordinates and elastic stiffness.
The existing34 table rows and5 figure references are unchanged. All18 pages
were inspected at contact-sheet scale, with changed pages2/3 reviewed separately.
This is neither a full external-SI graphic/version audit nor a new material
surface, Hessian, prospective prediction or submission-readiness claim.

The [E083 reversing terminal case](../../benchmarks/hfo2_channels/20261008/clamped_reversing_terminal_E083_20261011/README.md)
adds the third ordinary-converged G2 band to the main text and Table5.
Its42 fresh native SCFs and nine portable terminal frames pass distinct
checks; the338.968965meV/HfO2 forward discrete peak is not a continuous-TS
or sampling-error certificate. That dated build is identified by this
case's validation receipt, not the historical E082 PDF hash. One existing
+1% PO-to-T job is admitted into a spare slot without new submissions or
input changes; the other three training starters remain held.

The [E085 fixed-reference analysis](../../benchmarks/hfo2_channels/20261008/clamped_reversing_reference_E085_20261011/README.md)
adds the reversing band's three descriptive Table6 rows using the unchanged
T/Cmma references. Its T-triplet displacement/gradient fractions are26.2915%/
0.0320%; all four Cmma frames retain65.1046% of the squared atomic gradient
at displayed precision. These are different denominators, not energetic
accuracy or stationary-saddle claims. The three completed bands share the
same native cached PO inputs/logs and EFS; geometric serialization differences
are retained, not repaired. The [current reading-build receipt](../../benchmarks/hfo2_channels/20261008/clamped_reversing_reference_E085_20261011/validation.json)
identifies the19-page PDF, actual61-test local/archive runs and inspected
layout. G2 remains3/8; five training bands, conditional branches, strong
controls and frozen independent predictions remain necessary.
