# Hafnia pathways: scientific reading manuscript

The current paper is **Atomic accommodation reshapes strain responses of
competing hafnia pathways**. The October 21 scientific review separates
completion of the calculation protocol from journal readiness.

Read the [Chinese guide](AUTHOR_READING_GUIDE.md), the
[compiled manuscript](reading/VARNEB_HfO2_manuscript.pdf), and the
[explicit pre-submission recommendation](../../docs/HFO2_PRESUBMISSION_VERDICT_2026-10-21.md).
The recommendation is that this is a complete, reviewable material-study
draft; JCTC readiness has not passed the scientific acceptance gate.
The outstanding risk is the increment over known route concepts, not
missing code tests or a need for an ever-larger calculation matrix.

## Scientific result and reading route

The material question is how atomic accommodation permitted by one partial
clamp changes reversal relative to escape from the same polar well.
The lower preserving candidates have sampled peaks of 83.166, 87.719 and
91.699 meV/HfO2 at 0, +0.5% and +1% strain. Their strain penalty is 8.533
rather than 28.596. The matched T-escape peaks are 84.417, 90.418 and 96.764;
the small gaps do not certify global barrier ordering.

Start with **main Figure 4**: it shows mixed oxygen/Hf rearrangement, all
twelve displacement vectors, the candidate strain response and the final
local states. Figure 9 gives the complete three-strain comparison;
Figures 10--11 give matched local and crossed controls.
The former two-coordinate release overview remains intact as Figure A9
in in-article Appendix F.3. No old numerical table or figure source is deleted.

The final cross changes the pattern contrast from 31.342 to 15.369 meV/HfO2.
At +1%, retaining the zero-strain vector gives a local candidate only
0.699 above the lower connected peak. These interventions connect patterns,
transverse motion and the open row; they are not a unique mode-energy
partition or a new theory. The [closest-primary-source comparison](../../benchmarks/hfo2_channels/20261008/mechanism_synthesis_20261021/literature_comparison.md)
states what prior work already establishes and where our atom-resolved
correspondence remains partial.

## Source and build

[MANUSCRIPT_DRAFT.md](MANUSCRIPT_DRAFT.md) is the **only editable scientific
reading source**. The [TeX wrapper](varneb_JCTC_draft.tex) reads it directly.
Other planning/draft Markdown files in this directory are historical working
documents, not competing manuscript sources. The [CPC software manuscript](../VARNEB_CPC/README.md)
is a separate artifact.

From this directory, with the existing LuaLaTeX/markdown toolchain:

~~~sh
latexmk -norc -lualatex -interaction=nonstopmode -halt-on-error -outdir=build varneb_JCTC_draft.tex
~~~

The build output is copied to reading/VARNEB_HfO2_manuscript.pdf after visual
inspection. The [scientific-review build receipt](reading/scientific_acceptance_20261021.json)
records the actual source/data/PDF hashes, page anchors and relevant tests.
Earlier dated receipts record their own historical builds, not the current PDF.

## Reproduce the figures and tables offline

~~~sh
python -m scripts.plot_hfo2_mechanism_synthesis --output NEW_FIGURE_DIRECTORY
python -m examples.hfo2_mechanism_results --output NEW_TABLE_DIRECTORY
~~~

Run these from the source checkout root. Neither command starts DFT or SSH.
The plotting command needs the established Python/ASE/NumPy/matplotlib
environment; the table export needs only the standard library.
All output directories must be fresh.

Key evidence cases:

- [Atom mechanism and literature correspondence](../../benchmarks/hfo2_channels/20261008/mechanism_synthesis_20261021/README.md).
- [Final two crosses: fourteen new native SCFs](../../benchmarks/hfo2_channels/20261008/strain_pattern_cross_E277_20261021/README.md).
- [Same-source three-strain complete bands](../../benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019/README.md).
- [Matched preserving local completion](../../benchmarks/hfo2_channels/20261008/preserving_local_completion_E276_20261019/README.md).
- [Lower connected zero-strain candidate](../../benchmarks/hfo2_channels/20261008/preserving_connection_E273_20261019/README.md).
- [Midpoint initialization control](../../benchmarks/hfo2_channels/20261008/midpoint_reverse_E269_20261018/README.md).

All calculations keep the declared ABACUS/PBE 100 Ry, full10auDZP,
pseudopotential/orbital identities, k points and electronic settings.
The current synthesis uses existing data only. Original independent
predictions, failures and source labels remain in the manuscript appendices
and archived cases. No further calculation is launched by the scientific
review, and no submission or journal acceptance is represented as completed.
