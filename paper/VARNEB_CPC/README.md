# VARNEB CPC manuscript package

This is the working Computer Physics Communications manuscript package for
VARNEB. It follows the supplied ZStar CPC package only for its `elsarticle`
layout, Program Summary convention, and self-contained submission structure.
The wording, scientific claims, source data, and figures are specific to this
repository.

## Build

From this directory, run:

```text
latexmk -pdf -interaction=nonstopmode -halt-on-error varneb_CPC.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error varneb_CPC_supplement.tex
```

The outputs are `varneb_CPC.pdf` and `varneb_CPC_supplement.pdf`. The main
source uses the two-column `elsarticle` 5p layout; the supplement is a separate
two-page article. The package uses `elsarticle-num`, supplied by standard
Elsevier/TeX Live installations. `figures/` contains the PDF source used by
the manuscript and its CSV source data.

## Evidence policy

The manuscript distinguishes three evidence levels:

1. Analytic and fixed-cell regression tests establish the generalized
   coordinate, force, CI, and restart mechanics.
2. ABACUS/PBE material calculations establish the BTO and HfO2 examples.
   BTO's cubic Gamma projection uses the five-atom `1x1x1` phonon cell;
   its 59-point contour is a *frozen-cell* slice, not a relaxed conditional
   surface. The main text now shows the separately audited `Q_y=0`
   symmetry-restricted variable-cell sheet, with nine measured training nodes
   and two prospective DFT holdouts below the predeclared 2-meV/BTO gate. It
   does not claim an unrestricted conditional PES or a finite-temperature
   free-energy surface. The exploratory frozen cut and four-point conditional
   analysis remain staged evidence. The five
   hash-matched input JSON files are under `evidence/`; the plotting command
   and claim limits are in `FIGURE_LOGIC_AND_STYLE.md`.
   A separately sourced 79-point path-adapted frozen plane is retained as a
   supplementary candidate; it combines an unstable and a stable $\Gamma$
   triplet and must not be substituted for the main symmetry-restricted sheet.
   Regenerate it from the repository root with
   `python scripts/plot_bto_paper_path_adapted_79.py`; the committed QA
   distinguishes 79 DFT nodes from interpolated pixels and visual guides.
3. ABACUS, VASP, QE, ABINIT, and CP2K have each completed a 45.7-GPa GaN
   B4-to-B1 material path. The five paths support a common dominant barrier
   topology, not identical calculator protocols or barriers. Additional VASP
   GaN and CdSe paths remain mapping and functional-sensitivity controls.

The HfO2 literature bar is external comparison data with a different
functional, code, image count, and force criterion. The BTO literature number
is a locally restrained PBEsol DFT distortion-cost comparison, not a BTO
NEB benchmark. Neither is represented as a
statistical replica.

## Scope and length

VARNEB is a compact software paper, not a feature-for-feature counterpart to
the ZStar reference package. The working target is approximately ten typeset
pages including the Program Summary, figures, tables, and references. Add
material only when it closes a reproducibility, method, or verification gap;
do not pad the manuscript with duplicate workflow descriptions or unsupported
benchmark claims. The 2026-09-28 draft compiles to **10 main-text pages plus
three supplement pages/figures**. Author metadata and journal-specific final checks
remain submission gates. The main evidence set contains the architecture
figure, BTO restricted-sheet/path projection, GaN five-backend comparison,
the same-600-eV local GaN joint-coordinate cut, and the validated **central**
GaN atomic-transverse enthalpy cut. The BTO/HfO2 material controls and the
all-image GaN endpoint-Gamma mode/strain projection, and the supporting
GaN Hessian/endpoint-mode chart are Supplementary Figs. S1--S3. The central
GaN cut is not a full-path or
conditionally relaxed two-mode landscape, and the local candidate is not a
certified variable-cell transition state.
Per-path literature figures remain archived under
`outputs/neb_literature_benchmarks/` rather than being duplicated in the main
text.

The panel-by-panel argument, data normalization, and plotting conventions are
recorded in `FIGURE_LOGIC_AND_STYLE.md`. Plot generators remain in `scripts/`;
the manuscript reads the files under `figures/` explicitly, so stale root-level
exports cannot silently replace a revised panel.

The CP2K GaN path is now part of the converged five-backend figure. Its
45.7-GPa, 29-image chain reached `fmax=0.09409 eV/Å` and a
`0.29276 eV/GaN` barrier at image 15. The source-data CSV, exact evaluated
chain, and endpoint/geometry evidence are kept with the manuscript package.

The acceleration table uses Slurm-accounted, consecutively numbered ABACUS
process launches through the common `0.10 eV/Å` force threshold, not estimated
manager cache misses. Six first-crossing chains and common-arc profile/geometry
diagnostics are archived in `benchmarks/convergence/`; the source-bound claim
limits are in `ACCELERATION_EVIDENCE_AUDIT_2026-09-27.md`. These results do
not establish identical saddle basins or a universal speedup.

## Author metadata

The author and affiliation block is deliberately marked as a draft. Confirm
author order, affiliations, funding, CRediT roles, and the CPC Library link
before any submission.
