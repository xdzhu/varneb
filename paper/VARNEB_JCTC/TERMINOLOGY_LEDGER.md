# Manuscript terminology and claim boundaries

Updated 2026-10-18. Scope: HfO2 manuscript and its source-to-table checks.

| Term | Meaning used here | Does not establish |
|---|---|---|
| PO, T, M | Polar orthorhombic Pca21, tetragonal P42/nmc and monoclinic P21/c structural endpoints | Equivalent atom mappings or one unique pathway for each phase pair |
| Signed geometric Qx/Qy/Qz | Unit-Cartesian-norm rotated T oxygen-shuffle projections in the original T lattice metric, with one fixed periodic gauge per ordered chain | Phonon occupations, assigned irreps, local joint Q1/Q2, or a complete atomic representation |
| Polar displacement proxy | Mean O minus mean Hf displacement along original z in the same reference metric | Electronic polarization or switching charge |
| Open-axis stretch | 100(Fyy-1) in the original T chart; substrate fixes original x and z | Green strain, elastic-energy partition or an isolated causal mechanism |
| Sampled peak | Highest evaluated image of a declared ordinary NEB band | Continuous maximum, stationary TS or barrier error |
| Ordinary-converged band | Meets the declared 0.10 eV/Angstrom NEB residual | Zero physical gradient at each image |
| Internal local curvature | Finite differences in 36 translation-free joint atomic-cell length coordinates at fixed substrate | Mass-weighted phonon frequencies or an uncertainty-resolved stationary index |
| Prescribed control | Substrate strain, excluded from internal instability counting | A free unstable mode |
| Coupled local direction | A direction of the joint atomic-cell matrix in the declared chart | An endpoint phonon eigenmode |
| Gradient-converged candidate | Fixed retained coordinates, released-gradient norm at most 0.03 eV/Angstrom | Stable conditional minimum or unique branch |
| Conditional candidate sheet | Continued candidates at prescribed coordinates, with holes retained | Globally minimized PES or a certified smooth interpolation domain |
| Curvature seed | Standard measured-Hessian quasi-Newton initialization | A new optimizer theorem, guaranteed convergence or whole-path speedup |
| Raw negative count | Negative eigenvalues of a finite-step symmetrized matrix | Physical error bound or certified TS index |
| Censored cost | Calls spent before a common attempt cap without reaching the target | Completed time-to-target or a denominator for a speedup ratio |
| Open-cell trace | Tr(P_a P_-) for the entire raw negative subspace at the declared cell length scale | Barrier-energy fraction, phonon amplitude or arbitrary-metric invariance |
| Origin-aligned shape discrepancy | Model and candidate energy changes after separately subtracting their seen Q=0 values | Calibration-free forecast, independent error or whole-path barrier uncertainty |
| Matched frozen-R control | Same retained Q with R=0 about the original sampled center | Fully fixed cell when Q contains cell motion, or a relaxed conditional minimum |
| Released endpoint work | Signed endpoint-gradient trapezoid contracted with the released displacement between candidates | Missing energy contribution assigned uniquely to a mode, or an energy-error bound |
| Closed-loop endpoint circulation | Sum of endpoint-gradient work around a declared sampled grid loop | Nonconservative DFT forces, interpolation error or a branch-jump certificate |
| Missed release-shape change | Model-minus-DFT Q-dependent release change after each subtracts its own seen origin | A unique microscopic cause, whole-path barrier error or independent predictive advantage |
| Shared-observable control | Same physical covectors q=C.T x, H, affine gradient and mechanical boundary for all arms | Equal native parent-phase coordinates or equal preparation cost |
| Complete-space equivalence | Constrained reference and Schur solves represent the same complete quadratic | Independent accuracy, a new theorem or a basis-induced gain |
| Geometry-eligible comparison set | Coordinates with an observed candidate and a model geometry passing the original guard | Stable released branch, a global domain or the error of a missing point |
| Post-hoc coarsening | Geometric coarse/fine split made after the training labels and relaxation history were seen | A blinded preparation workflow, independent prediction or DFT-call saving |
| Tensor-quadratic interpolation | Energy interpolation of degree at most two in each retained coordinate, including mixed total-degree-four terms | A harmonic Hessian model, a new algorithm or a validated conditional branch |
| Atomic-only treatment | Release of 31 atomic fixed-Q directions at each frozen-Q seed cell | One common cell across Q, joint stationarity or transverse stability |
| Additional joint-release difference | Signed E_J-E_A between the joint and atomic-only candidates at the same Q | A unique strain-energy contribution, global nesting or a whole-path mechanism |
| Directional amplitude zeta | Internal length along the frozen omitted direction; case field t_A | Prescribed substrate-strain coordinate t or a phonon wave vector |
| Finite-amplitude energy/gradient secants | Separately defined kappa_E and kappa_g at the declared amplitude | Equal values outside a quadratic profile, total error or stationary-index certification |
| Formal third-coordinate retention | Removing the measured omitted direction from the released block at the same point | A computed stable three-dimensional nonlinear surface or new independent prediction |
| Original T local independent errors | Frozen E185/E188 predictions compared with eight later single-start gradient-candidate energies | An unseen-strain path test, full transverse stability or a new-method advantage |
| Two-sided preserving training trial | Ten registered starts at five Q positions, with both released-norm and directional-gradient tests | Ten minima, distinct stable branches or a complete three-dimensional PES |
| Retained-direction nonlinear fit | Registered 109-coefficient polynomial/harmonic-bath ansatz fitted to 81 native physical point observations | A new theorem, independently improved accuracy or a validated force field |
| Revised preserving P8 target | Eight original Q positions with frozen training-parent starts and dual-gradient candidate energies | Retroactive success of the refused two-coordinate sheet or unseen-strain path barriers |
| Independent physical-strain test | Two preselected nine-image ordinary bands at +0.5%, with all forecasts frozen before target-interior evaluation | Four-channel coverage at +0.5%, a fresh search over starting paths or a stationary-barrier error bound |
| Reverse-continuation control | A separate +0.5% ordinary band initialized from the same channel's +1% training band, using the original target endpoint cache | A new independent prediction target, replacement prospective label, unique MEP or continuous-barrier error bound |
| Nonmonotonic sampled competition | T-minus-preserving peak contrast of -12.411559/-8.112699/-28.660781 meV/HfO2 at 0/+0.5/+1% | Positive selectivity, a physical bifurcation or resolved continuous-path curvature |

Notation is not interchangeable: Figure 3 uses the original T Cartesian
metric, whereas the appendix representation-coverage figures use the declared
mass-weighted atomic T chart; local joint-curvature coordinates have length
units with L=5.12968067458423 Angstrom. Reported E160 total energies and energy
changes use the twelve-atom cell. Existing G2 barrier tables use their stated
per-formula-unit normalization. No ASR, atom rematching, new reference selection
or held-out labels are introduced by the manuscript integration.

The local tests ask whether omitted coordinates are equilibrated well enough
for a reduced description. The completed +0.5% physical test answers a
different question about the sampled response under a new substrate control;
it did not establish a new reduced-model prediction advantage.
Full T optical and all four Cmma registrations remain
strong controls; unequal ranks are disclosed rather than called equivalent.
