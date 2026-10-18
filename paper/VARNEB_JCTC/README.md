# Hafnia pathway study: complete reading manuscript

This is the complete scientific reading manuscript for author review,
including its abstract, conclusions, references and in-article appendices.
The [CPC software paper](../VARNEB_CPC/README.md) is a separate artifact.
One editable source, [MANUSCRIPT_DRAFT.md](MANUSCRIPT_DRAFT.md), supplies the
reading PDF through [the thin TeX wrapper](varneb_JCTC_draft.tex).

[Read the compiled manuscript](reading/VARNEB_HfO2_manuscript.pdf).
The single-column review layout retains the main article, references and
in-article appendices. Historical page counts belong to their dated builds;
the current PDF is regenerated from the single Markdown source. The
zero-strain branch revision has55pages. Figure4 is onpage11, the new
three-strain T release Figure5 onpage12 and the independent strain test
Figure6 onpage14. The separate T initialization control remains in
Section3.4, without changing the original labels or forecast scores.

The October19 update adds the completed zero-minus preserving local
control to Results3.3 and Discussion4. Holding the original three T-pattern
amplitudes and z structural polar proxy still allows0.122A transverse
mean-O-minus-Hf motion and13.599meV/HfO2 fixed-cell lowering. This is a
structural proxy, not electronic polarization. Its reconnected polygon
fails the ordinary NEB target, so no corrected barrier or channel
reordering is claimed. The separately bounded free-interior E273
connection test, not another reduced-model search, addresses that gap.

The [Chinese reading guide](AUTHOR_READING_GUIDE.md) explains the argument,
the supported contributions and the role of each main figure. The main
methods now contain the material/comparison definitions and calculation
recipe; the extended response theory and diagnostic procedures are retained
in Appendix G. All218existing numerical table lines and13prior figure
targets remain; the matched-observable control adds one new main figure.
The abstract has208words and the conclusion about
300words; the bibliography supplies17verified article records.

## Current reading route: mechanism first

The October 18 revision connects the results to one material question:
**how mechanical constraints change competition between homogeneous
switching and escape from the polar well, and which atomic rearrangements
accompany that response.** The main Results retain four sections:
path competition, mechanical boundary, microscopic interpretation, and
independent prediction checks.

- Figure 1 compares all four candidate peaks and their endpoint-only controls.
- Figure 2 retains the free-cell/clamped Pbcn comparison.
- Figure 3 follows signed oxygen shuffles and open-axis stretch on both flips
  and T escape; the +0.5% independent paths are now included.
- Figure 4 separates observed atomic-only and additional joint relaxation
  at all five local coordinates per center. Its added lower row connects
  those same candidates to the T-pattern triplet and its geometric
  complement; most preserving F-to-A motion lies in the latter.
- Figure 5 connects the completed T three-strain same-observable local
  energy release to actual atomic readjustment. It does not turn local
  lowering into a corrected barrier or compare nonnested controls by subtraction.
- Figure 6 completes the independent strain test with actual paths, the
  nonlinear sampled competition and errors of all three numeric forecasts.

The [three new main figures and source data](figures/hfo2_mechanism_20261018/README.md)
have a single replay command:

~~~sh
python -m scripts.plot_hfo2_mechanism --output paper/VARNEB_JCTC/figures/hfo2_mechanism_20261018
~~~

No new DFT, force audit, parameter change, atom rematching or fitted smoothing
is introduced by that command. Previous representation figures and every
numerical table row are retained in the same article: the former main
Figures 2--5 are now Figures A5--A8; former Tables 3--5 are Tables A20--A22.
The Results plus Discussion are shorter even with the new mechanism text.

The bounded ten-start preserving trial and revised P8 test are complete.
Their details and negative model-comparison findings remain in Appendices
E and F; they do not define the visual center of the paper.
The [separately frozen +0.5% physical test](../../benchmarks/hfo2_channels/20261008/physical_holdout_E268_20261018/README.md)
is also complete and included in Results3.4/Table3/Figure6. Peaks are
98.530877/90.418178meV/HfO2 for preserving/T; residuals0.079773/0.080825eV/A.
The T-minus-preserving contrast improves at +0.5% relative to zero but
remains negative, then worsens at +1%. B0 predicts T closely and misses
preserving; endpoint-only B1 models trade those errors differently. This
does not establish a new reduced-model prediction advantage.
The [updated figures and complete scoring](figures/hfo2_independent_20261018/README.md)
replay with:

~~~sh
python -m scripts.analyze_hfo2_physical_holdout --output paper/VARNEB_JCTC/figures/hfo2_independent_20261018
~~~

All six numeric predictions and132unavailable transfer records are retained.
The registered test closes without new fitting or a larger calculation matrix. The
[research blueprint](../../docs/HFO2_PAPER_BLUEPRINT_2026-10-17.md)
remains the bounded completion plan, not a new experiment matrix.

## Scientific question and current evidence

Can one prescribed mechanical intervention reduce homogeneous switching
without reducing the easiest examined escape barrier from the same PO+ well?
Can a branch-aware mode-strain model predict that response better than the
registered strong references and simpler controls?

All eight ordinary clamped training bands and all ten endpoint representations
are audited. At zero/+1% substrate strain the sampled preserving switching peak
increases by 28.596248 meV/HfO2, while the higher reversing peak decreases by
44.024331. The lower sampled escape peak increases by 12.347027. The raw
candidate-set selectivity proxy therefore worsens by 16.249222, not the desired
favorable signature. This is not an error-certified physical ordering: the
strained T/M escape peaks differ by only 0.464068 meV/HfO2 and their continuous
bottleneck errors are unmeasured.

The current paper is organized by the material argument, not job chronology:

1. Complete channel-response matrix and endpoint-only controls.
2. Mechanical boundary and the changed energetic role of the Pbcn snapshot.
3. Signed oxygen-shuffle routes, open-cell response and matched atomic/joint release.
4. Completed independent strain test, full frozen controls and local candidate-energy checks.

Detailed signed-work, cross-strain, curvature and representation checks are
retained in Appendix E, not treated as six additional main discoveries.
The complete matched atomic treatment is in Appendix D.8/Table A10.
All five local preserving coordinates showed 44.06--48.35 meV per twelve-atom
cell of atomic-release lowering; joint candidates were a further 0.80--1.14
meV/cell lower. This local treatment contrast is not a unique mode-energy
partition, whole-path mechanism or stable-surface certificate.

The complete post-release check now supplies an actual preserving-origin
sign-screen refusal, not only a pending warning. Its omitted oxygen-dominated
direction stays raw-negative at both steps and in the fixed-cell restriction.
Four additional same-contract training statics measured only0.0126--0.0172
meV/HfO2 lowering. Appendices E.7 and F.1 distinguish that
shallow response from a corrected barrier, full stationary index or independent
accuracy claim. Both secant definitions and adjacent work discrepancies remain.

Free-cell candidates, finite sampling checks, endpoint screens and complete
reference tables remain **inside the paper's appendices**. No separate
supplementary manuscript is introduced. All nine original figure assets and
every published numerical table row are retained; figure/table numbers have
changed. Dated validation files retain their historical hashes and page counts,
not the current article numbering.

## Terminology

| Term | Meaning and boundary |
|---|---|
| Sampled peak, hat B | Maximum among the retained images, not a certified stationary barrier |
| Ordinary NEB residual | Default 0.10 eV/Angstrom max-vector target, not a barrier-error bound |
| Physical gradient g | Raw energy derivative, not spring or perpendicular NEB force |
| Conditional surface | One identified branch relaxed in its stable permitted complement, not a frozen slice |
| Signed work W | Coordinate/connection-dependent contraction, not a unique phonon energy |
| Prescribed strain epsilon | External control, excluded from the internal instability count |
| Strong T/Cmma references | Complete declared ranks; all four tied Cmma frames remain, not a selected winner |

The two lengths in the response chart are not mass-weighted phonon amplitudes.
Geometry coverage, gradient coverage, integrated work and prediction accuracy
are different measurements.

## Scope closure and remaining manuscript work

The later [finite mechanism-closure plan](../../docs/HFO2_MECHANISM_CLOSURE_2026-10-18.md)
is active and is distinct from the completed forecast protocol. The +0.5%
T continuation from the +1% band has completed, with a sampled peak of
90.286243 meV/HfO2 versus the original 90.418178. Section 3.4 retains both
and explains why agreement of this scalar does not establish identical
structures or a continuous-barrier error bound. The preserving counterpart
is still being evaluated; it is not imputed from the T outcome. None of the
original target labels or frozen prediction scores is replaced.

The bounded ten-start trial, revised P8 test and two-path physical +0.5%
test are complete. Their failures and unavailable forecasts are results,
rather than requests to start another fitting or DFT campaign. The physical
test establishes no prediction advantage for the local branch-adapted model;
that stronger claim is withdrawn from the current paper.

Stable connected conditional surfaces and certified stationary indices have
not been established. They would be required for stronger conditional-PES or
transition-state claims, but are not scheduled as additional experiments
under this completed sampled-path protocol. The original preserving-plane
refusal remains visible.

The abstract, introduction and conclusion now state the measured material
contribution consistently. The introduction synthesizes the relevant
primary literature; the complete references credit both the material
studies and the numerical tools. Journal positioning follows the material
contribution; a new-method JCTC claim cannot be inferred from local model
fits. The author-review PDF uses the existing data and figures throughout.

## Sources and reproduction

- [Detailed methods](METHODS_DRAFT.md), [research design](RESEARCH_DESIGN.md),
  [finite plan](../../docs/VARNEB_JCTC_HFO2_RESEARCH_PLAN.md).
- [Prospective protocol](../../docs/HFO2_PREDICTION_PROTOCOL.md) and its
  [v2 addendum](../../docs/HFO2_PREDICTION_PROTOCOL_V2_2026-10-09.md).
- [Complete eight-band source audit](../../benchmarks/hfo2_channels/20261008/clamped_plus1_reversing_terminal_E148_20261014/README.md).
- [Cross-strain numerical source](../../benchmarks/hfo2_channels/20261008/cross_strain_work_E151_20261014/README.md).
- [Historical pilot record](PILOT_RESULTS_2026-10-08.md) and dated validation
  files in this folder. Earlier nonconvergence statements refer to their own
  timestamps; they are not current job status.

From this directory, with the existing TeX Live/LuaLaTeX/markdown packages:

~~~sh
latexmk -norc -lualatex -interaction=nonstopmode -halt-on-error -synctex=1 -outdir=build varneb_JCTC_draft.tex
~~~

The reading output is build/varneb_JCTC_draft.pdf. No installation, shell
escape or latexmk rc loading is required. This is a multi-file project,
not a standalone TeX file. Source and figure dependencies are read directly;
no manually synchronized scientific TeX body exists. Keep a portable local
reading copy in build so its relative evidence links resolve with a reader
that permits local links. External review-build paths recorded in receipts
are layout snapshots, not portable source-data bundles.

Build products are ignored by Git; the reviewed PDF is retained separately
in reading/ for repository access. The October18 complete reading build
passed119related scientific/manuscript checks. All51pages were inspected
in contact sheets and the abstract, conclusions and bibliography at reading
resolution; no overfull/underfull boxes or missing characters remain.
All218numerical table lines and13figure targets are retained.
Receipts identify the actual compiler,
PDF bytes and inspected pages, not byte-identical PDFs across environments.
Existing published assets and raw outputs are not overwritten by a new build.

The Oct17 E257 build has38 pages including appendices. Its actual native
four-point result,157 local/157 independent related checks (not a full-suite
rerun), source preservation and rendered-page review are recorded in the
[E257 receipt](../../benchmarks/hfo2_channels/20261008/G3_soft_probe_E256_20261016/results_validation_E257.json).
This supersedes the older36-page build status, not the scientific completion
gates above. The previous generated reading copy was retained on C:.

The Oct17 E259 build updates Section3.4/Table4 and AppendixE.8/TableA18 with
the first original joint-space T independent cohort:8 labels,18 nativeSCFs,
all19 frozen arms and17.833999/19.173333 call-proxy/allocation coreh. The
40-page working PDF has been rendered and checked; it remains unfinished
research, not a submission-ready paper. See the [complete T test](../../benchmarks/hfo2_channels/20261008/G3_T_independent_E258_20261017/README.md)
and its separate results_validation_E259.json receipt. E195/E257 receipts
are retained as historical structures/builds rather than rewritten.

The E261 reading build remains40pages. Section2.4 now records the approved
fixed-Q signed-start branch experiment, its explicit directional stopping
test and numerical-preconditioner boundary. The original182 table lines,
nine figure assets and13 external citations are unchanged; training outcomes
are still pending. The [E261 methods/build receipt](../../benchmarks/hfo2_channels/20261008/G3_preserving_branches_E260_20261017/manuscript_methods_E261.json)
records67related manuscript checks and visual review of all40pages, including
the new paragraph on page6. This is a reproducibility update, not a new
material result or a submission-readiness claim.

## Fixed computation scope

ABACUS/PBE, Hf4O8/four formula units, 100 Ry/full 10 au DZP/original six file
hashes/Gamma2x2x2/P0/E0. No cutoff, orbital, pseudopotential, SCF, CI or
unregistered material/functional/orientation matrix change. The same
mechanical active space applies at every image and endpoint of a comparison.
The physical +0.005 path predictions were frozen before the two bands were
launched; both completed terminal paths now appear in the figures and scoring.
