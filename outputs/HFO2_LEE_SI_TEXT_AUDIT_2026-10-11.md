# Lee supplementary-text audit: scope, units and the method comparison

Reading date: 2026-10-11. Source: the publisher-linked supplement to
[Lee, Lee and Yu, npj Quantum Materials 11, 34 (2026)](https://doi.org/10.1038/s41535-025-00841-9).
[Supplementary PDF](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41535-025-00841-9/MediaObjects/41535_2025_841_MOESM1_ESM.pdf).
This is a new text reading, not a completed PDF/figure/version audit.
Exact access limitations are in the [provenance record](HFO2_LEE_SI_TEXT_AUDIT_2026-10-11.provenance.md).

Later evidence, 2026-10-12: [E108](HFO2_LEE_SI_VISUAL_AUDIT_2026-10-12.md)
obtains and visually checks the four-page publisher PDF. The historical
limitations below remain the actual E082 record, not current access claims.

## What the newly accessible text establishes

- Part 1 relates cubic/tetragonal mode descriptions in different primitive
  cells. Its notation does not assign irreps to our registered T patterns.
- Part 2 describes a coupled three-mode interpolation and restricted
  relaxation. Release-induced lowering of a profile is therefore prior work,
  not an innovation claim for VARNEB.
- Part 3 develops a long-wavelength elastic argument. Its printed elastic
  notation must not be silently imported as our internal response matrix.

Figures S1-S3 and the numerical entries of Table S1 have **not** been visually
verified in this turn. No plotted curves, thresholds, Landau coefficients or
mode vectors are copied, digitized or used as new material inputs.

## Our dimensional and scientific distinction

For the declared length-valued response chart, q and t are in Angstrom,
t=L(epsilon-epsilon0), and the same-centre quadratic A/b/d are in
eV/Angstrom^2 (cell energy). Its stationary internal response is

```text
dq*/d epsilon = -L solve(A,b)       [Angstrom per unit dimensionless strain]
```

An elastic stiffness is energy per volume. A fitted mode-amplitude slope
has displacement per strain units (or the explicitly converted mass-weighted
equivalent); neither is automatically an elastic modulus or a mixed Hessian.
Percentage strain and dimensionless strain also differ by a factor of 100.
This is our dimensional inference, not a numerical correction of the source.

A continuum acoustic construction alone does not determine the actual
multi-atom internal-coordinate response. The latter needs the registered
mixed derivatives, reference transport and mechanical freedoms. Our existing
[stationary-response implementation](../docs/CONTROLLED_STATIONARY_BRANCH.md)
already supplies the relevant standard algebra; it is not new theory.
No source coefficient is substituted for an ABACUS/PBE measurement.

## Consequence for the study

The main-text introduction now explicitly disclaims relaxation-induced
profile lowering as novelty; Section 2.3 distinguishes length-valued response,
mass-weighted projection and elastic stiffness. The meaningful proposed test
is still matched-boundary switching versus escape and a resolved prospective
advantage over the strongest valid registered controls.

This does not select G3 modes, revise B0-B5, change thresholds, add a training
condition, or authorize a full Hessian/extra matrix. G2 remains incomplete and
the +0.005 control is still unread. The original electronic inputs and two-
allocation limit remain unchanged. Behara's subscription-only supplement and
the historical Lee reference-PDF version discrepancy remain unresolved.
