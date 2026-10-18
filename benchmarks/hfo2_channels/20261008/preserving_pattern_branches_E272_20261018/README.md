# E272: preserving-mode branch controls

## Actual native start, 2026-10-18 23:13 CST

Execution source31c31178 was pushed before deployment. The fresh HF root is
`/public/home/iai806/abacus/agent-runs/20261018-varneb-preserving-branches-E272-r1`;
its17-file immutable bundle pin is
`d12af56eac84edf28e536c69b524d12681f876db1918c889fbf4d89eb59e06c2`.
Zero-DFT native preflight passed both original sources and all four signed
geometry/input roundtrips after T28947723 completed. Slurm returned28948425;
held1node/32CPU/1task/6h resources were checked before release.

Job28948425 started23:10:00 onnode1/hfacnormal01. At23:13:10 its first
SCF passed real32MPI, the same six physical inputs, full E/F/stress and
physical-gradient replay,92.861s. The zero-strain minus0.10A seed lowers
energy by0.540401meV/HfO2, retaining four-pattern defect<9e-18A and exactly
zero cell coordinates. This is one evaluated seed, not a completed A/J
treatment, two-side equivalence, stable minimum or a corrected barrier.
See the preflight/held/accepted/released/startup receipts. E269 preserving
remains the only other study allocation; no midpoint branch is admitted.

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
itself has a negative released direction. The prepared numerical initializer
uses default unit BFGS scale along that branch and on-site3x3-block scale
on its complement. It is positive without clipping, while the true DFT
couplings and measured negative curvature remain present. A bounded quartic
regression rejected a pure on-site seed before deployment; it is not a
DFT acceleration benchmark. The initial proposal remains in Git4cd02905.

`prepared/` contains only the original0/+1% P peaks, cached F E/F/stress,
matched covectors, signed0.10A atomic starts and transported numerical seeds.
Midpoint admission remains false pending the E269 initialization control.
These are geometric proposals, not evaluated candidates. At preparation,
no new Hessian, DFT, deployment or submission occurred. See [the finite protocol](protocol.md):
at most96 optimization attempts plus8 single-direction checks if executed
later; no whole Hessian, model fitting, altered parameters or unlimited retry.

```text
python benchmarks/hfo2_channels/20261008/preserving_pattern_branches_E272_20261018/prepare.py --output FRESH_DIRECTORY
python -m pytest -q tests/test_hfo2_preserving_pattern_branch_preparation.py
```

Eight focused preparation tests and five execution-plumbing tests now pass,
alongside numerical/export checks (34 total in the admission run). They test the actual native inputs, direction signs/ranks,
exact atomic-only cell embedding, positive numerical seed and nonexecution;
they do not prove successful branch relaxation or a completed mechanism.

`worker.py` and `run.slurm` implement the finite prospective batch. They
import only pinned original E270/E271 helpers and the immutable E268 runtime,
and never overwrite those sources or submit jobs. Each signed arm has its
own evaluator/cache and at most26 fresh calls for A/checks/J; four signed
arms give104 total. A missing A blocks J and direction checks. A known
direction check can be positive while whole stability remains uncertified.
The native preflight and held32CPU/6h admission must occur after T28947723
has completed. Do not submit a speculative third dependent allocation.
