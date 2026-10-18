# E272: preserving-mode branch controls

## What is preserved, and what is not

The registered name describes a signed reference distortion, **not** the
coordination identity of each oxygen. A zero-DFT geometric check of the
original0/+0.5/+1% endpoint pairs finds the same ordered role exchange at
all three strains. With full-cell zero-based indices4--11, the initial
oxygen Hf-neighbor counts are[4,4,3,3,4,4,3,3] and the final counts are
[3,3,4,4,3,3,4,4]. Every assignment is unchanged across geometric distance
cutoffs2.4/2.5/2.6/2.7A. This is postprocessing, not a changed electronic
cutoff or input. Atom order and periodic images are retained.

`oxygen_roles.json` gives each endpoint's third/fourth/fifth periodic Hf
neighbor distances and count sensitivity. `inspect_oxygen_roles.py` and
`tests/test_hfo2_preserving_oxygen_roles.py` reproduce this source-pinned
three-condition check. No E/F/stress or new DFT is claimed by this
geometry-only calculation. Coordination exchange alone does not identify
plane crossing, an irrep, electronic polarization or a complete prior
path-family label.

### Literature boundary checked on2026-10-19

| Primary source | What is already established | Consequence for our interpretation |
| --- | --- | --- |
| Qi, Singh and Rabe,2025, [lattice-mode analysis](https://doi.org/10.1103/PhysRevB.111.134106) | Polar and antipolar changes distinguish switching variants. | Coupled modes and multiple switching paths are not our new claim; our geometric T descriptors are not automatically their irreps. |
| Mu et al.,2025, [bulk and electrode/ferroelectric/electrode switching](https://www.nature.com/articles/s41524-025-01633-2) | Their pathway definitions include transverse oxygen displacement and three/four-fold coordination exchange. | Our transverse structural proxy or endpoint role exchange alone is not discovery of a new family. Ordered trajectory and boundary matching are still required. |
| Fan, Zhu and Liu,2025, [functional dependence of Pbcn](https://www.nature.com/articles/s41524-025-01647-w) | Mechanical constraints and strain reference affect phase/path comparisons. | Their hydrostatic, PO-referenced strain is not our two-vector T-referenced clamp; do not equate thresholds or change our PBE settings to match them. |

The two open-access npj publisher HTML texts were read directly. For Qi,
the2025 publisher abstract was checked; the accessible
[2021v1 manuscript](https://arxiv.org/pdf/2108.12538) is an earlier version,
not a substitute for claiming to have read the complete2025 published text.
OpenAlex discovery returned429 and was not repeatedly queried. Only primary
texts/metadata support the table, not third-party paper summaries.

The candidate paper-level question remains whether admitting this mixed
atomic accommodation changes reversal-versus-escape competition under
the **same specified mechanical ensemble**, and whether its strain response
survives initialization and the tested finite resolution. The present
transverse geometry is an observation, not yet a new switching mechanism,
a polarization-rotation measurement or a JCTC-readiness certificate.

## Complete original batch, 2026-10-19 00:36:18 CST

Job28948425 finished COMPLETED0:0 after1h26m18s. All54 native SCFs,
including both budget-failed branches, were exported with450 files in
`native_terminal.tar.gz` (SHA256
`fee229e74ed9132edf7ec2a38288893a4e66bdb6b113d3ace8ce6b0fb97681e9`).
All six assets/MPI/EFS/chart gradients were rechecked on HF; all54 native
EFS/physical-gradient records were independently replayed locally. Local
exports omit pseudo/orbital binaries, so their full bytes are not claimed
to have been rehashed locally. See `native_terminal_receipt.json` and
`complete_analysis/analysis.json`.

| Strain | Signed start | A-F (meV/HfO2) | J-F (meV/HfO2) | Fresh A/check/J |
| --- | --- | ---: | ---: | --- |
| 0 | minus | -13.599258 | -13.617082 | 12/2/1 |
| 0 | plus | -13.599564 | -13.617889 | 12/2/1 |
| +1% | minus | unavailable | unavailable | 12/0/0 |
| +1% | plus | unavailable | unavailable | 12/0/0 |

Both zero-strain signs pass the local released-gradient target. In the
common original-T metric, negative A has0.476962A atomic displacement,
80.364% oxygen squared norm and a+0.121857A change in the transverse y
mean-O-minus-Hf structural proxy, while the three displayed T patterns
and z proxy stay fixed. The other A sign gives-0.121879A transverse
change. These are geometric observations, not electronic polarization,
energy fractions or a proof of global branch equivalence. The separate
known-direction energy/gradient secants at A are also independently
reproduced; no full stability or J Hessian is certified.

At+1%, the two atomic arms exhausted12 attempts each. Their final trial
norms0.499979/0.499972A approached the predeclared0.5A neighborhood bound;
released atomic gradients0.377865/0.377769eV/A remained above0.03.
No qualified A, J or direction-check value is imputed. This is a bounded
optimization outcome, not an electronic-SCF failure or proof that no
other minimum exists. No bound/parameter change or repeated budget follows.

Actual new cost is45.648128native-proxy/46.026667allocated core-hours;
all earlier curvature/path/preparation costs remain additional. Original
labels and forecast scores are unchanged. The next scientific test is
the separately frozen E273 zero-minus free-interior connection, not an
automatic extension of these local failures.

The complete54-call analysis is reproducible locally without DFT or SSH:

```text
python benchmarks/hfo2_channels/20261008/preserving_pattern_branches_E272_20261018/analyze_terminal.py --output FRESH_DIRECTORY
python -m pytest -q tests/test_hfo2_preserving_terminal_analysis.py
```

## Zero-new-DFT connection result, 2026-10-19

The four selected native SCFs in the zero-minus prefix have now been
independently reparsed locally: energy, forces, stress and physical chart
gradients agree, the four held observables are preserved, and the separate
energy/gradient directional secants at A are reproduced. The archive is
still a partial batch snapshot, not all attempts or whole-job completion.

Replacing only image4 with the native J candidate while preserving the
original fixed per-atom integer lift gives the following full-polygon
replay. Every tangent, spring and allowed atomic/cell force is recomputed:

| Polygon | Sampled peak (meV/HfO2) | Peak image(s) | Ordinary NEB fmax (eV/A) |
| --- | ---: | --- | ---: |
| Original zero-strain preserving | 96.828183 | 4 | 0.098817 |
| Original chain with free J image4 | 85.499889 | local peaks3/5, maximum3 | 0.125611 |

The new center itself has residual0.018500eV/A; its neighbors3/5 have
0.125611/0.124071. The central turn has cosine-0.4560 and the adjacent
segments are longer. Thus the local13.617082meV/fu lowering does not
automatically produce a smooth converged path, and subtracting it from
the old barrier would give the wrong sampled maximum. This directly
motivates a free-interior connection test, not another local Hessian.

See `zero_minus_connection/analysis.json` and `analyze_prefix.py`.
[The finite E273 follow-up](../preserving_connection_E273_20261019/protocol.md)
admits only the first declared zero-minus branch, no midpoint or new strain
matrix. Its inert exact-cache seed has all seven interiors free and the
original endpoints fixed; it has **not** been submitted or admitted to run.
Wait for E272 terminal data and a free study slot; no third dependent job.

At00:15CST, the whole E272 job is still RUNNING with42 completed SCFs and
no failure file. Both zero-strain signed A/J controls passed (12/1 calls
each). Their J energies differ by0.000808meV/fu; this is a two-start
observation, not a proof of symmetry or global equivalence. At+1%, minus A
exhausted its12-attempt budget, so J/checks were not attempted and no
candidate energy is reported. The original plus arm continues. This is
an optimization-budget outcome, not an electronic-SCF failure or a
permission to alter the physical settings.

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

## First completed signed branch, partial batch snapshot

At23:42:46 the whole job was still RUNNING,21 fresh calls completed and no
failure file. Zero-strain minus A passed in12 attempts and J in1; two
additional probes at **A** completed the one-direction check. The A/F and
J/F changes are-13.599258/-13.617082meV/HfO2; J-A is-0.017824meV/HfO2.
A/J released-gradient norms are0.008776/0.011617eV/A. At A only, the known
direction has positive finite-step gradient/energy secants1.247123/0.965819
eV/A². These are separate raw estimates, not a full Hessian or a TS certificate;
J's full stability has not been tested.

`prefix_0_minus.tar.gz` retains the four selected native records for A,J
and its two direction probes, plus outcomes. The36-file byte manifest was
verified locally. It is **not** all15 attempts or a terminal whole-job
export; remaining signs/strain are still running. See the explicit scope
in `prefix_0_minus_receipt.json`. Independent local physical replay follows.
Do not overwrite or relaunch the original batch.

This local lowering is comparable with the original12.4meV/HfO2 channel
contrast, motivating a connection check, not subtraction from a barrier.
The next bounded step is a zero-new-DFT replay of a complete nine-image
polygon using original native images and this J candidate. All tangents,
spring and allowed atomic/cell forces must be recomputed; copying the old
NEB residual is invalid. Preserve endpoints/order/periodic gauge and all
physical settings. Keep J as a free interior candidate, never a certified
TS endpoint. Missing geometry/EFS evidence stops replay without imputation.
Any later real continuation requires a free allocation and a separately
frozen finite protocol. Original labels and prediction scores stay unchanged.

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
