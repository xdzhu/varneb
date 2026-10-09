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
escape coverage or distinct winding sectors. Two common-substrate endpoints
(T and PO+) now pass the declared atomic-force/open-traction screens.
G2 matched-boundary channel barriers, G3
conditional-branch validation and independent material predictions remain
missing. An ordinary residual pass is not a certified saddle or a barrier
uncertainty bound. No abstract or conclusion asserting these missing results
is supplied.

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
