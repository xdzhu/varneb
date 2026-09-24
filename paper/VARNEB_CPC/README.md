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
```

The main output is `varneb_CPC.pdf`. The package uses `elsarticle-num`, which
is supplied by standard Elsevier/TeX Live installations. `figures/` contains
the PDF source used by the manuscript and its CSV source data.

## Evidence policy

The manuscript distinguishes three evidence levels:

1. Analytic and fixed-cell regression tests establish the generalized
   coordinate, force, CI, and restart mechanics.
2. ABACUS/PBE material calculations establish the BTO and HfO2 examples.
3. The VASP adapter has completed material paths on the Hefei Slurm system.
   The GaN tetragonal and hexagonal routes reproduce the published barrier
   scale; the B3 route and CdSe route are retained as explicit mapping and
   functional-sensitivity controls.

The HfO2 literature bar is external comparison data with a different
functional, code, image count, and force criterion. The BTO literature number
is a restrained-NEB energy-scale comparison. Neither is represented as a
statistical replica.

## Scope and length

VARNEB is a compact software paper, not a feature-for-feature counterpart to
the ZStar reference package. The working target is approximately ten typeset
pages including the Program Summary, figures, tables, and references. Add
material only when it closes a reproducibility, method, or verification gap;
do not pad the manuscript with duplicate workflow descriptions or unsupported
benchmark claims. The main evidence set now contains the architecture figure,
the BTO mode decomposition, the ABACUS BTO/HfO2 validation figure, and the GaN
multi-backend comparison. Per-path literature figures remain archived under
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

## Author metadata

The author and affiliation block is deliberately marked as a draft. Confirm
author order, affiliations, funding, CRediT roles, and the CPC Library link
before any submission.
