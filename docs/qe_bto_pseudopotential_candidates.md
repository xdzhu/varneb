# BTO QE pseudopotential candidate register

This register records a candidate only. It is not a source of UPF files, does
not authorize a calculation, and must not be passed to a production template
until the user accepts the selected pseudopotential family.

`examples/bto_qe_sssp_1_3_pbe_precision_candidate.json` transcribes the Ba,
Ti and O entries from the official SSSP 1.3.0 PBE precision collection. The
collection is intentionally mixed-family (Ba PseudoDojo NC, Ti GBRV USPP, and
O PSLibrary PAW), as selected by SSSP for the respective elements. Its largest
reported recommendation is 75 Ry for wavefunctions and 600 Ry for charge
density. The project candidate retains 100/600 Ry: 100 Ry is the requested
BTO wavefunction cutoff, while 600 Ry preserves the oxygen density-cutoff
recommendation.

Before activation, all four requirements in the JSON must be met. In
particular, the 100 Ry choice is a conservative starting point, not a BTO
cutoff-convergence result. A single-image convergence check must precede the
five-worker VCNEB production submission.
