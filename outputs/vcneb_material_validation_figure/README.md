# VC-NEB material validation figure

## Figure conclusion

The BaTiO3 T-to-C path is monotonic at 5, 7, and 9 total images, whereas the
HfO2 T-to-PO path has a reproducible interior barrier.  For HfO2, the
7-image CI refinement gives `32.29 meV/f.u.`, close to the `32 meV/f.u.`
reported by Liu *et al.* for their 40-image LDA/QE/USPEX VC-NEB calculation.

## Figure contents

- **a**: BaTiO3 ordinary 5/7/9-image paths, normalized by each archived
  geometric path length.  No line has an interior peak, so CI is intentionally
  absent.
- **b**: HfO2 ordinary log-strain 7/9-image paths, the 7-image linear-cell
  control, and the accepted 7-image CI refinement.  The open circle marks the
  climbing image.
- **c**: HfO2 barrier values in meV per HfO2 formula unit.  The literature bar
  is a reference value, not a fourth replica or an uncertainty estimate.
- **d**: Relative volume evolution for the HfO2 paths, showing the cell degree
  of freedom rather than treating the calculation as a fixed-cell NEB.

`vcneb_path_comparison_source_data.csv` is the long-form source-data table;
`vcneb_barrier_literature_comparison.csv` gives all plotted barrier values and
their provenance.  SVG and PDF retain editable text; PNG is only a preview.

## Comparison limits

The BTO literature result is a PBEsol calculation with local restraints, so
its `2.1 kcal/mol` endpoint rise is only an energy-scale comparison to the
unconstrained ABACUS/PBE VC-NEB result (`2.01 kcal/mol`).  The HfO2 literature
path instead uses LDA/QE/USPEX, 40 images, and an RMS force criterion of
`0.025 eV/A`; the present calculations use ABACUS/PBE, 7 or 9 total images,
and a maximum generalized-force criterion.  The close CI barrier is therefore
strong validation evidence, but it is not a claim of exact method identity.

No statistical error bars are shown: these are deterministic convergence and
control calculations, not independent stochastic replicas.  The loose-threshold
endpoint-displacement mode-guided HfO2 branch is intentionally excluded because
it is not comparable to the strict ordinary/CI evidence table.

## Literature

1. Liu *et al.*, *Phys. Rev. Materials* **3**, 054404 (2019),
   DOI: [10.1103/PhysRevMaterials.3.054404](https://doi.org/10.1103/PhysRevMaterials.3.054404).
2. The BTO restrained-NEB comparison is reported in *Phys. Chem. Chem. Phys.*
   (2020), DOI: [10.1039/C9CP02955A](https://doi.org/10.1039/C9CP02955A).
