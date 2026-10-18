# Hafnia JCTC study: evidence-bound working manuscript

This is a research manuscript in development, **not submission-ready**.
The [CPC software paper](../VARNEB_CPC/README.md) is a separate artifact.
One editable source, [MANUSCRIPT_DRAFT.md](MANUSCRIPT_DRAFT.md), supplies the
reading PDF through [the thin TeX wrapper](varneb_JCTC_draft.tex).

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
  and T escape at both strain conditions.
- Figure 4 separates observed atomic-only and additional joint relaxation
  at all five local coordinates per center.

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
has been launched for preserving switching and T escape. Its terminal
results are not included in this manuscript revision. The
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
4. Prediction controls, the completed T local test and the still-unfinished mechanical-response prediction.

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

## What remains before the intended JCTC claim

- The original two-center joint-curvature batches are complete at both probe
  steps. A measured stability-resolution floor and the appropriate stationary
  fixed-strain internal index remain separate; highest images are not automatically TSs.
- Stable, connected conditional branches and covered sampling domains.
- The preserving2D origin's admission remains refused; a changed retained
  coordinate, branch or test definition needs a separately recorded pre-label
  protocol, not an overwrite of the original16-point forecast record.
- The original T8 local test is now natively audited and scored, not pending.
  All five standard maximum aligned errors are below0.257meV/HfO2, but
  unaligned maxima change the ranking and the full T reference equals fullH.
  This does not prove a new method's superiority or complete the preserving8
  and physical+0.005 mechanical-response tests.
- Training-only model/representative selection and frozen predictions before
  any physical +0.5% geometry or path labels are read.
- Independent comparisons against B0--B5, including strong fixed-reference,
  endpoint-information and interpolation controls, with actual preparation,
  SCF and allocation costs. An unavailable model abstains.
- A demonstrated nonlinear/finite-representation/branch or equal-accuracy-cost
  gain: changing a complete harmonic basis alone cannot improve its exact
  stationary response.
- Final numerical/figure/layout/citation audit and permitted reproducibility
  package. This structural revision is not scientific completion.

The completed [matched-ablation analysis](../../benchmarks/hfo2_channels/20261008/G3_frozen_release_E173_20261015/README.md)
and [connected manuscript validation](../../benchmarks/hfo2_channels/20261008/G3_release_manuscript_E174_20261015/README.md)
supersede the earlier preparation-only status. They do not close the stable-
branch or independent strong-control prediction gates.
The first bounded native local-curvature batches were
[registered in E150](../../benchmarks/hfo2_channels/20261008/G3_curvature_E150_20261014/README.md).
The [E152 read-only replay](../../benchmarks/hfo2_channels/20261008/G3_curvature_audit_E152_20261014/README.md)
checks actual completed raw probe prefixes without filling missing columns,
new DFT or a live-source change. Neither is a completed conditional-face or
prediction result. Job progress belongs to those dated records, not the paper.

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

Build products are ignored by Git. Receipts identify the actual compiler,
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
launched; no terminal labels from those bands are used in the present figures.
