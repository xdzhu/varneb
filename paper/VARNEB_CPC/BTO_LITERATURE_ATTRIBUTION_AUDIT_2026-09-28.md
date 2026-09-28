# BTO literature-attribution audit — 2026-09-28

## Citation resolution

One targeted reference was checked: `BTOReaxFF2019`, DOI
`10.1039/C9CP02955A`. The BibTeX title, first author, journal, and year
match the [publisher record](https://pubs.rsc.org/en/content/articlelanding/2019/cp/c9cp02955a):
metadata status **verified**. The method attribution in the prior manuscript
and source-data CSV was a **mismatch** and has been corrected.

The [author-version article](https://pubs.rsc.org/en/content/getauthorversionpdf/C9CP02955A),
pages 7--8, distinguishes two calculations. Its NEB method and residual-force
criterion belong to *oxygen-vacancy migration*. For a BaTiO3 unit-cell
tetragonal-to-cubic distortion, it reports a 2.1-kcal/mol PBEsol DFT cost
under **local restraints**; its Figure 2b places the cubic state at the
high-energy end of that conversion. The latter calculation used VASP/PBEsol,
a 6-by-6-by-6 electronic mesh, and a 600-eV plane-wave cutoff. The cited
paper does not claim that its T-to-C conversion is an NEB calculation.

Our separate ABACUS/PBE 100-Ry/10-au-DZP, 4-by-4-by-4 electronic-mesh,
five-atom variable-cell path has a monotonic 0.0871289017-eV/BTO endpoint
rise, equal to 2.00924 kcal/mol per mole of BTO formula units. The numerical
scale is similar, but different functionals, restraints, geometries, and
path algorithms preclude a like-for-like NEB benchmark or an activation-
barrier claim. The corrected manuscript, figure source-data generator,
committed CSV, package note, and regression test preserve that distinction.
