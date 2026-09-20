# VARNEB path-by-path literature benchmarks

This directory contains one standalone comparison figure for every distinct
physical path that had a completed, auditable VARNEB summary on 2026-09-20.
PNG files are previews; SVG and PDF are the manuscript-quality versions.  Each
figure has a same-stem `*_source_data.csv` containing every plotted value and
its provenance.

## Barrier comparison

| Path | VARNEB | Literature | Interpretation |
|---|---:|---:|---|
| BaTiO3 T -> C | 0.08713 eV/f.u. | 0.09106 eV/f.u. | Close energy scale; literature evidence archived here contains only the 2.1 kcal/mol endpoint rise, not a full path table. |
| HfO2 T -> PO | 0.002973 eV/f.u. | 0.002500 eV/f.u. | Full 20-image curves agree closely in barrier and overall downhill profile. |
| HfO2 PO -> M | 0.08362 eV/f.u. | 0.10508 eV/f.u. | Full curves reproduce the same single-barrier topology; the VARNEB barrier is about 20% lower. |
| GaN B4 -> B1, tetragonal | 0.33849 eV/GaN | 0.34 eV/GaN | Excellent barrier agreement at 45.7 GPa; the path shapes are also consistent. |
| GaN B4 -> B1, hexagonal | 0.38523 eV/GaN | 0.39 eV/GaN | Excellent barrier agreement at 45.7 GPa; the saddle position shifts because each curve uses its own normalized coordinate. |
| GaN B3 -> B1 | 0.95529 eV/GaN | about 0.57 eV/GaN | Negative benchmark: the present diagonal mapping collapses to one high barrier and does not reproduce the three-event literature mechanism. |
| CdSe rock-salt -> wurtzite, cell mapping | 7.216 meV/atom | 2.4 meV/atom | Same small-barrier/downhill energy scale, but not quantitative agreement; current PBE and literature PW91 paths are not like-for-like. |

All VARNEB entries satisfy the project default maximum generalized-force
criterion of `0.10 eV/A`.  Repeated optimizer runs of the same physical path
are not counted as additional benchmark paths.  The unfinished CdSe atomic
mapping branch is intentionally excluded.

## Provenance and normalization

- BaTiO3 is the completed ABACUS/PBE 100 Ry, 10 au DZP, 7-total-image
  static-audited chain.  The literature comparison is the PBEsol restrained
  NEB value reported as `2.1 kcal/mol`; it is shown as a horizontal reference
  rather than a fabricated path curve.
- HfO2 VARNEB curves are the completed VASP/PBE 20-total-image chains copied
  from the Hefei `hfacnormal01` calculations.  Literature points are exact
  author-workbook values from the source data behind PRL 130, 226801 (2023),
  sheets `T-pca21` and `pca21-M`.  Cell energies are divided by four HfO2
  formula units.
- GaN VARNEB curves are the completed VASP/PBE Hefei results at 45.7 GPa
  (B4 -> B1) or 45.0 GPa (B3 -> B1).  Cell energies are divided by two or four
  GaN formula units as appropriate.  Literature curves were manually
  digitized from Figs. 4, 6, and 8(a) of Qian et al.; the paper-stated barrier
  values, rather than digitization maxima, are used in the comparison text.
- CdSe VARNEB energies are divided by all eight atoms and converted to
  meV/atom.  The literature composite `a -> c -> e` route was manually
  digitized from Fig. 11 of Sheppard et al. and therefore supports a
  qualitative path-shape comparison, not sub-meV numerical claims.
- The blue VARNEB abscissa is normalized geometric arc length.  Literature
  abscissae are normalized published image index or digitized path distance.
  Overlaying them compares topology and energy scale; it does not assert an
  identical reaction-coordinate metric.

The exact input-summary hashes, plotted barriers, and generated filenames are
recorded in `manifest.json`.  Recreate the figures with:

```bash
python scripts/plot_literature_path_benchmarks.py
```

The default paths assume the completed Hefei summaries have been staged under
`tmp/plot_data`.  The command-line options allow explicit archived summary
locations to be supplied instead.

## Literature

1. G.-R. Qian et al., *Computer Physics Communications* **184**, 2111-2118
   (2013), DOI: 10.1016/j.cpc.2013.04.004.
2. G. Sheppard et al., *Journal of Chemical Physics* **136**, 074103 (2012),
   DOI: 10.1063/1.3684549.
3. HfO2 Fig. 2(a) author-source workbook accompanying *Physical Review
   Letters* **130**, 226801 (2023).
4. BaTiO3 restrained-NEB comparison: *Physical Chemistry Chemical Physics*
   (2020), DOI: 10.1039/C9CP02955A.
