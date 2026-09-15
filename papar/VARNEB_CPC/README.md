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
3. The VASP adapter is a single-image calculator smoke test only, not a
   production cross-engine VC-NEB benchmark.

The HfO2 literature bar is external comparison data with a different
functional, code, image count, and force criterion. The BTO literature number
is a restrained-NEB energy-scale comparison. Neither is represented as a
statistical replica.

## Scope and length

VARNEB is a compact software paper, not a feature-for-feature counterpart to
the ZStar reference package. Keep the submitted article to approximately eight
typeset pages including the Program Summary, figures, tables, and references.
The current draft is intentionally shorter. Add material only when it closes a
specific reproducibility or verification gap; do not pad the manuscript with
duplicate workflow descriptions, redundant plots, or unsupported benchmark
claims. The single four-panel material-validation figure is the primary
results figure unless new evidence materially changes the scientific case.

## Author metadata

The author and affiliation block is deliberately marked as a draft. Confirm
author order, affiliations, funding, CRediT roles, and the CPC Library link
before any submission.
