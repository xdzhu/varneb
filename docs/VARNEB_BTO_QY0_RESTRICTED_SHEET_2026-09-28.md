# BTO: fixed-\(Q_y=0\) sheet for an endpoint-connected two-coordinate figure

This is a deliberately narrower physical construction than the fully relaxed
\(\min_{Q_y,\mathrm{other}} E(Q_z,Q_x,\ldots)\) envelope. The latter cannot
contain the cubic C endpoint at its reference energy: cubic BTO has a third
independent unstable Γ mode, \(Q_y\). The code here fixes the two plotted
cubic-soft-mode coordinates \((Q_z,Q_x)\), the rigid-translation gauge, **and
\(Q_y=0\)**. Other atomic coordinates and all six symmetric cell strains may
relax. This is a *mode-restricted sheet*, not a proof of the full mirror-group
symmetry, the global lower envelope, a VCNEB barrier, or a finite-temperature
T/C free-energy surface. The omitted transverse instability must be disclosed
beside any contour.

The source remains the archived five-atom **1×1×1 Γ** force-constant/eigenvector
calculation. Every new static uses ABACUS 3.10.0 LTS, PBE, `ecutwfc=100 Ry`,
Ba/Ti/O **10 au DZP**, the original pseudopotentials, the **4×4×4 electronic**
k mesh, `scf_thr=1e-8`, and the same cubic energy zero. No cutoff or orbital
parameter was raised. Conditional-point optimization targets an open-subspace
gradient of `0.003 eV/(sqrt(amu) Å)` and raw maximum stress of `2 kbar`;
these are not the ordinary NEB `0.10 eV/Å` stopping criterion.

## First endpoint canary: completed and raw-audited

- Preflight: `outputs/batio3_t_to_c_pbe100_dzp10au/bto_qy0_q000_q000_preflight_2026-09-28.json`
  (SHA-256 `bca070dea043423966cae79bbecbfe97aeb2474fec1a6d6b45f6320f41be5a73`).
  Its cubic starting structure has minimum periodic atomic distance `2.01411 Å`.
- hf `hfacnormal01` Slurm `27794614`: `COMPLETED 0:0`, elapsed `00:01:25`,
  one 32-MPI ABACUS static and no duplicate endpoint optimization.
- Independent raw-output audit:
  `/public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot/audit-q000-27794614.json`
  (SHA-256 `8916f75827a40061cb82a291443e48b5634fb68298c7532b701d1b20cfd1ba50`).
  It matched the cached point to raw SCF energy, all atomic forces, all stress
  components, `INPUT/KPT/STRU` hashes and `DSIZE=32`.
- At `(Q_z,Q_x)=(0,0)`: `Q_y=0`, cell strain zero, maximum atomic force zero,
  open-subspace gradient `0.00010423 eV/(sqrt(amu) Å)`, maximum raw stress
  `0.09074 kbar`, and `E−E_C=4.2141e-8 eV/BTO` (reference rounding).
  Thus the C endpoint is representable **under this explicit restriction**.
  Orthogonal curvature and branch/global-minimum status remain unchecked.

The earlier open-\(Q_y\) pilot already contains `Q_y≈0` stationary branches
at `(0.6,0)`, `(0.6,0.3)`, `(0.9,0)`, and `(0.9,0.3)`. A new independent
read-only cross-audit, `scripts/audit_bto_qy0_reuse.py`, matched each branch
replay to its original cached point and raw SCF energy/force/stress, checked
the common 100-Ry INPUT/KPT hashes and Γ basis, verified `Q_y≈0`, and rebuilt
the 15-dimensional restricted gradients. All four pass the `0.003`
gradient/`2 kbar` stress gate **without new DFT**. Their energies are
`−0.0447322552`, `−0.0549838599`, `−0.0740973712`, and
`−0.0816209937 eV/BTO`; the restricted gradient norms are `0.00120035`,
`0.00179149`, `0.00266678`, and `0.00173965`. The hf source-path
manifest including the independent center holdout is
`benchmarks/numerical_integrity/bto_qy0_reuse_sources_with_holdout_2026-09-28.json`;
the compact five-point report includes all 21-component terminal geometries
and original log hashes.
They are now **eligible measured nodes** of the restricted sheet, not
curvature-certified minima. The previously selected open-\(Q_y\) lower
branches lie `28.308–57.467 meV/BTO` below them and must never be silently
substituted into this sheet.

The **independently pre-existing center** `(0.75,0.15)` also has a `Q_y≈0`
branch. The same read-only raw-log/metric audit passes it without new DFT:
`E−E_C=−0.0630474180 eV/BTO`, restricted gradient `0.00219180`, maximum
stress `0.27434 kbar`. The five-branch report is
`benchmarks/numerical_integrity/bto_qy0_reuse_with_holdout_audit_20260928.json`
(SHA-256 `20a25940ee4abc860102b88472d9eda4b9f3d5f2038218ff02758febc880919a`),
with its own hf source-path manifest. Four-corner bilinear interpolation on
`(0.6,0)–(0.9,0.3)` predicts `−0.0638586200 eV/BTO` at that center;
the independent measured energy differs by `+0.81120 meV/BTO`. The
full atom-plus-strain metric-coordinate discrepancy is
`0.06591 sqrt(amu) Å` (`0.00960` atomic, `0.06520` strain). This validates
**one local cell at one center**, not the entire restricted sheet; the
structure discrepancy and branch continuity still require an explicit
figure-level interpretation.

## Finite next sampling block, not a dense blind grid

Slurm array `27794649` submitted four **new** fixed-\(Q_y=0\) points:
`(0,0.3)`, `(0.3,0)`, `(0.3,0.3)`, `(1.2,0.3)` in `sqrt(amu) Å`, at most two
one-node/32-MPI tasks concurrently on `hfacnormal01`. Their local geometry
preflights all pass the `1.6 Å` minimum-distance guard (minimum `1.706 Å`).
Each result is a candidate only until the original SCF/forces/stresses,
fixed-\(Q_y\) coordinate, gradient, stress, branch identity and input hashes
have been audited. A Slurm `COMPLETED` state alone does not pass the gate.

The first array element `(0,0.3)` finished seven valid ABACUS evaluations and
reached open-subspace gradient `0.0023853`, but maximum raw stress was
`2.3176 kbar`, above the predeclared `2 kbar` gate. Its Slurm state is
`FAILED 1:0` because the runner deliberately rejected this **stress-only**
failure; it is not an SCF or geometry crash. The seven raw points passed the
independent audit at
`/public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot/audit-q000_q030-27794649_0-v2.json`.
Warm continuation `27794672` lowered only the optimizer-coordinate gradient
target from `0.003` to `0.001` while keeping the same Hamiltonian, Q
constraints, `2 kbar` physical stress gate and existing completed cache.
It `COMPLETED 0:0` after ten new evaluations: gradient `0.00084591`,
maximum raw stress `0.53564 kbar`, and `E−E_C=−0.0137253015 eV/BTO`.
All 17 original-plus-new DFT points passed the independent raw audit at
`/public/home/iai806/abacus/agent-runs/20260928-varneb-bto-qy0-pilot/audit-q000_q030-27794672.json`.
This is a targeted response to the measured residual, not a retrospective
loosening of the stress threshold. It still lacks a restricted-subspace
curvature and grid-interpolation certificate.

The `(0.3,0)` and `(0.3,0.3)` elements had the same stress-only outcome after
seven raw-audited points each: maximum stress `2.32295` and `3.40490 kbar`,
respectively. Their cache-preserving, `0.001`-gradient continuations were Slurm
array `27794708`; both `COMPLETED 0:0`. Independent raw audits cover all 17 and
19 completed DFT evaluations, respectively, and verify the same fixed $Q_y=0$
and INPUT/KPT contract. At `(0.3,0)`, the final gradient is `0.00084038`,
stress `0.54363 kbar`, and `E−E_C=−0.0137356970 eV/BTO`; its audit is
`audit-q030_q000-27794708_0.json` (SHA-256
`24cc7a026dca36700b356a90e1577da6d1eb8f3c5de6a43aa2775fb647cb0ac1`).
At `(0.3,0.3)`, the corresponding values are `0.00086182`, `0.31064 kbar`,
and `−0.0263718031 eV/BTO`; its audit is
`audit-q030_q030-27794708_1.json` (SHA-256
`6ce52295156a45c2e9f31a3e87f02849d07a7b218aacad772366b65496a9521b`).
These are gradient/stress-eligible local candidates, not curvature-certified
conditional minima.

The `(1.2,0.3)` initial element `27794649_3` terminated after 22 **valid**
ABACUS SCFs. The safeguarded optimizer reached the predeclared orthogonal
coordinate amplitude bound `3.0 sqrt(amu) Å`; the final open-subspace gradient
was still `0.0591044 eV/(sqrt(amu) Å)` and maximum raw stress `18.1299 kbar`.
This is an optimizer-coordinate boundary hit, **not** DFT convergence, an
acceptable sheet node, or an electronic failure. A dedicated independent audit
of all 22 raw energies, forces, stresses, geometry, fixed $Q_y$, input hashes,
and MPI records is `audit-q120_q030-27794649_3.json` (SHA-256
`d92f7c9de77f9c78f691092b40f94c66d6d3b9a6ccddcfdd718d8d644b76e722`).
The active orthogonal coordinate was `2.99999993`; the nearby actual T endpoint
has a larger $\eta_{zz}$ than the failed point, so the `3.0` optimizer bound
was not a physically justified cap on this T-side boundary. The cache-preserving
continuation `27794797` used an independently audited cached starting geometry
and a `4.0` optimizer-coordinate bound. It changed **no ABACUS Hamiltonian,
100 Ry cutoff, 10 au DZP orbital, k mesh, fixed order parameter, or 2 kbar
physical stress target**. It `COMPLETED 0:0` after 17 new DFT points. The
independent audit of **all 39** old-plus-new raw SCFs, energies, forces,
stresses, fixed Q values, and input hashes passes:
`E−E_C=−0.0921388405 eV/BTO`, open-subspace gradient `0.00279650`,
maximum stress `0.42050 kbar`, and $Q_y≈0$. Audit:
`audit-q120_q030-27794797.json` (SHA-256
`0aab778d111913404f28cd485cc467684b737e96b7c27510af660af5fb811c03`).
The frozen runner version mislabeled this warm start as `Q_y=0_frozen` in its
human-facing branch label; the independent auditor required that exact runner
hash and the failed-cache-audit hash before recording the label correction.
The calculation was seeded from the audited terminal cache geometry, with no
repeat of the first 22 DFT evaluations. This is now an eligible T-side
measured node, still not a curvature-certified minimum.

The five new pilot nodes, four older grid nodes, and the older independent
center holdout have been assembled **without interpolation** into
`benchmarks/numerical_integrity/bto_qy0_ten_measured_nodes_20260928.json`
(SHA-256 `c04a2c0927f5d57279f1ae21405eddd3b97ec36d40541f6f26a03a912df73a07`)
and a compact CSV sibling. The assembler verifies each new summary against
its independent raw audit and its SHA-256, rejects duplicate Q coordinates,
and preserves all ten 21-component terminal atom-plus-strain geometries.
There are **nine grid nodes and one independent center holdout**; this is a
measured patch, not yet an accepted smooth conditional-PES contour. The
archived source manifest also records the remote raw-audit paths and hashes.

Before calculating any further center energies, we froze a three-point
holdout plan at `(0.15,0.15)`, `(0.45,0.15)`, and `(1.05,0.20)` in
`benchmarks/numerical_integrity/bto_qy0_three_holdout_plan_20260928.json`
(SHA-256 `b3bef1ea6360c662d4913c1d8c4d8dcf4c09d5ec9efd687c67869b3fdf9bfced`).
The first two are centers of previously measured rectangular cells; the last
is inside the T-side measured triangle. It records the already-computed
bilinear/barycentric energy and full-structure predictions, plus prospective
`2 meV/BTO` energy and `0.10 sqrt(amu) Å` full atom-plus-strain coordinate
error screens. These numerical screens are local figure gates, not rigorous
global error bounds. All three fresh Qy=0 preflights passed the same source
and minimum-distance checks; Slurm array `27794991` runs at most two 32-MPI
tasks concurrently on `hfacnormal01`. No holdout result is accepted merely
because its Slurm task finishes; each will be independently raw-audited.

Before plotting an interpolated restricted surface, assemble a measured
path-covering grid from these points plus **audited, nonduplicated** existing
stationary branches; fill any missing T-side boundary point, measure interior
holdouts chosen before fitting, and compare energy **and full atom-plus-strain
coordinates**. Show every actual DFT node. If branch continuity or holdout
tests fail, keep a measured-point/triangulated stage figure rather than a
smooth conditional-PES contour. The actual seven-image variable-cell T→C
energy and strain stay in separate panels and are never read from this sheet.

## Pre-registered interpolation gate: failed, with a bounded model response

All three holdout calculations have now finished on `hfacnormal01`; the first
needed a cache-preserving stress-only continuation rather than a repeated SCF.
The independent raw-output auditor covered **6, 13, and 18** DFT evaluations,
respectively. The frozen-plan comparison is reproducible with
`scripts/audit_bto_qy0_preregistered_holdouts.py` and recorded in
`benchmarks/numerical_integrity/bto_qy0_three_holdout_gate_20260928.json`.
The first two energy errors exceed the predeclared `2 meV/BTO` figure gate;
all three open-subspace gradients, raw stresses, and full atom-plus-strain
coordinate errors pass their respective gates:

| `(Q_z,Q_x)` (`sqrt(amu) Å`) | Measured `E−E_C` (meV/BTO) | Linear predicted (meV/BTO) | Measured−predicted (meV/BTO) | Structure error (`sqrt(amu) Å`) |
| --- | ---: | ---: | ---: | ---: |
| `(0.15,0.15)` | −7.199 | −13.458 | +6.259 | 0.0733 |
| `(0.45,0.15)` | −31.295 | −34.956 | +3.661 | 0.0650 |
| `(1.05,0.20)` | −86.142 | −84.372 | −1.770 | 0.0607 |

This is an interpolation-model failure, not a failed electronic calculation
or a reason to alter `ecutwfc`, orbital radius, k sampling, or the physical
stress gate. In particular, a bilinear function of *amplitudes* misses the
quadratic energy behavior close to the cubic soft-mode origin.

We therefore froze an exploratory **even, axis-exchange-symmetric sixth-order
energy polynomial**, with its intercept fixed to `E_C`, fitted only to the
**nine previously measured grid nodes**. The script and coefficient/prediction
artifact are `scripts/plan_bto_qy0_even_mode_validation.py` and
`benchmarks/numerical_integrity/bto_qy0_even_mode_two_holdout_plan_20260928.json`.
Its training RMS is `0.213 meV/BTO`; errors at the four already-known interior
points are `0.455`, `−0.225`, `−0.133`, and `−1.040 meV/BTO`. **Those four are
retrospective checks, not independent validation of a model chosen after
seeing three of them.** Two genuinely new off-grid predictions are frozen at
`(0.45,0.25): −36.297 meV/BTO` and `(1.05,0.25): −86.706 meV/BTO`, with the
same `2 meV/BTO` energy and `0.10 sqrt(amu) Å` full-structure gates. Even if
they pass, representative restricted-subspace curvature and branch continuity
must still be checked before calling the chart a conditional-minimum sheet.
