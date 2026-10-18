# E272: preserving-mode branch preparation, not a submitted job

The main question is whether oxygen rearrangements outside the displayed
patterns account for the preserving bottleneck's residual accommodation.
This differs from T's positive-atomic-block conditioning problem. Keep the
same physical Qx/Qz/Qy and polar structural proxy as the T interventions,
not merely similarly named orthonormal axes.

Reuse of the original zero-strain P native input and both completed raw
curvature probes shows:

| Original probe | Minimum released atomic curvature (eV/A²) | Gradient along that direction (eV/A) |
| --- | ---: | ---: |
| 0.01A | -0.447449 | 9.47e-7 |
| 0.02A | -0.455586 | -8.01e-7 |

The sign-fixed directions overlap0.99999975. Oxygen supplies78.6% of its
squared physical displacement norm, Hf the remainder; this is **not** an
elemental energy fraction, phonon occupation or finite-temperature result.
The original frozen geometry is not stationary (released atomic norm0.2881
eV/A), so these matrices do not certify a physical saddle index.

T's whole-atomic-block Jacobi initializer is not admissible here: that block
itself has a negative released direction. The prepared numerical on-site
3x3-block initializer is positive without clipping, while the true DFT
couplings and measured negative curvature remain present.

`prepared/` contains only the original0/+1% P peaks, cached F E/F/stress,
matched covectors, signed0.10A atomic starts and transported numerical seeds.
Midpoint admission remains false pending the E269 initialization control.
These are geometric proposals, not evaluated candidates. No new Hessian,
DFT, deployment or submission occurred. See [the finite protocol](protocol.md):
at most96 optimization attempts plus8 single-direction checks if executed
later; no whole Hessian, model fitting, altered parameters or unlimited retry.

```text
python benchmarks/hfo2_channels/20261008/preserving_pattern_branches_E272_20261018/prepare.py --output FRESH_DIRECTORY
python -m pytest -q tests/test_hfo2_preserving_pattern_branch_preparation.py
```

Six focused preparation tests and related numerical/export tests passed
(30 total). They test the actual native inputs, direction signs/ranks,
exact atomic-only cell embedding, positive numerical seed and nonexecution;
they do not prove successful branch relaxation or a completed mechanism.
