# Atomic mechanism and pre-submission synthesis

This zero-DFT case connects existing ordered atom motion to complete paths
and the final strain--pattern cross comparison. It introduces no material,
strain condition, force-constant calculation, phonon assignment or path search.

Read [the primary-source comparison](literature_comparison.md) for the
distinction between a supported material result and a new-path claim.
The four-panel main figure shows the original/lower center projection,
all twelve displacement vectors, three-strain path peaks and local
transverse/open-axis states. Atom arrows use the original T metric and
their actual magnitude, with common translation removed. They exclude
affine cell deformation; panel (d) separately shows the open row.

```sh
python -m scripts.plot_hfo2_mechanism_synthesis --output FRESH_DIRECTORY
```

The source reports and their hashes are recorded in analysis.json. They
are reused, not recomputed. Points are finite sampled candidates; local
energies are not new complete-path barriers. No uncertainty interval is
inferred from a force tolerance. Existing labels and forecast scores remain.
