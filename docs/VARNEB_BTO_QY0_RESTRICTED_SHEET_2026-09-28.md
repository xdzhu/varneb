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

The earlier independently audited open-\(Q_y\) pilot already contains the
`Q_y≈0` stationary branch at four off-center coordinates. Its energies are
recorded separately in
`paper/VARNEB_CPC/figures/bto_conditional_four_point_stage_2026-09-27_source_data.csv`.
They are not rerun merely to populate this sheet; their full raw cache,
\(Q_y\) amplitude, mode basis and 15-dimensional restricted gradients must be
cross-audited before promotion to restricted-sheet data. The previously
selected open-\(Q_y\) lower branches lie `28.308–57.467 meV/BTO` below those
stationary branches and must never be silently substituted into this sheet.

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
respectively. Their cache-preserving, `0.001`-gradient continuations are Slurm
array `27794708`, submitted only after independent audit. The `(1.2,0.3)`
initial element remained active at this writing; its outcome is not inferred
from the other three.

Before plotting an interpolated restricted surface, assemble a measured
path-covering grid from these points plus **audited, nonduplicated** existing
stationary branches; fill any missing T-side boundary point, measure interior
holdouts chosen before fitting, and compare energy **and full atom-plus-strain
coordinates**. Show every actual DFT node. If branch continuity or holdout
tests fail, keep a measured-point/triangulated stage figure rather than a
smooth conditional-PES contour. The actual seven-image variable-cell T→C
energy and strain stay in separate panels and are never read from this sheet.
