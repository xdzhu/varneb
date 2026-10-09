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

The manuscript's numerical table is a dated snapshot, not a live job dashboard.
One T-PO band has passed the ordinary residual criterion; three other pilot
observations in that table have not. G2 matched-boundary comparisons, G3
conditional-branch validation and independent material predictions remain
missing. An ordinary residual pass is not a certified saddle or a barrier
uncertainty bound. No abstract or conclusion asserting these missing results
is supplied.

The [later complete-frame residual audit](../../benchmarks/hfo2_channels/20261008/residual_update_20261009_0450/README.md)
tracks the preserving flip's atomic/cell residual crossover without replacing
the draft's dated table or claiming final channel energies. It is additional
source-verified pilot evidence, not a new stationary-TS or prediction result.

The [2026-10-09 morning update](../../benchmarks/hfo2_channels/20261008/morning_update_20261009_0850/README.md)
adds the second ordinary-converged source, PO to M, and the preserving flip's
developing split-peak profile. Both switching chains remain unconverged.
These later results are recorded in the pilot text; the compiled main text's
earlier dated numerical table and figures are not silently overwritten.

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

The 2026-10-09 build was compiled with existing TeX Live 2023/LuaHBTeX 1.16.0,
then all seven pages were rendered and visually reviewed with existing Poppler.
Three display equations, four table rows and two embedded figure assets were
checked. See [the compile/visual receipt](MANUSCRIPT_DRAFT.compile.validation.json).
That is a reading-layout check, not scientific completion or journal approval.

The earlier [main-text validation receipt](MANUSCRIPT_DRAFT.validation.json)
is preserved as a **historical milestone** for commit `456d376`. Its manuscript
hash predates the current math/table formatting changes and must not be used
as the hash of the current draft.

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
