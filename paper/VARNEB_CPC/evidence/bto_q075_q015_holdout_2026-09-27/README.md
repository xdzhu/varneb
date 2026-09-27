# BTO conditional soft-mode plane: one internal holdout

This folder freezes the result and independent raw-output audit for the
pre-declared `(Q_z,Q_x)=(0.75,0.15) sqrt(amu) Å` holdout. It is evidence for
**one local interpolation test**, not a continuous 2D conditional PES or a
T→C activation barrier. The prediction and acceptance rules were written
before the result in
`docs/VARNEB_BTO_CONDITIONAL_HOLDOUT_PROTOCOL_2026-09-27.md`
(SHA-256 `5e85bf7c93812f32a93942a1f54ed569d7f1de0215e7f341e1d9d6355f1d5b7a`).

Slurm `27787487` ended `COMPLETED`, exit `0:0`, after `01:18:04` on hf.
The frozen ABACUS/PBE contract is `ecutwfc=100 Ry`, Ba/Ti/O full 10 au DZP,
electronic `4×4×4`; the Γ force-constant source uses the five-atom `1×1×1`
real-space cell. The independent audit re-read all 77 original DFT outputs:
three previously audited starting canaries (`27787471`) and 74 new
evaluations (`27787487`). It checked each SCF, 32-MPI record, input and log
hashes, energy, force, stress, fixed-Q projection, and geometry.

| Basin from unique audited stationary cache point | `E−E_C` (eV/BTO) | `Q_y` (sqrt(amu) Å) | Open-gradient norm (eV/(sqrt(amu) Å)) | Max raw stress (kbar) | Max atomic force (eV/Å) |
| --- | ---: | ---: | ---: | ---: | ---: |
| `Q_y=0` | −0.063047418038 | ≈0 | 0.002191801 | 0.274343 | 0.465964 |
| `+Q_y` | −0.106887865039 | +1.032331062 | 0.001507641 | 0.771759 | 0.103456 |
| `−Q_y` | −0.106887865032 | −1.032331062 | 0.001507641 | 0.771759 | 0.103456 |

All three satisfy the *conditional* gradient limit `0.003` and stress limit
`1 kbar`. The constrained `Q_y=0` point can have a nonzero raw atomic force
along the fixed mode coordinates; its open-coordinate gradient is the relevant
stationarity criterion here, not the ordinary NEB `0.10 eV/Å` threshold.

The locked four-corner prediction was `−0.106411735308 eV/BTO` for the low
branch and `−0.063858619978 eV/BTO` for `Q_y=0`. Observed minus predicted
is `−0.476130` and `+0.811202 meV/BTO`, respectively. The low-branch
error passes the pre-declared `2 meV/BTO` **single-holdout** line; the
observed low-versus-zero gap is `43.840447 meV/BTO`, and the signed pair
differs by only `7.73×10⁻⁹ meV/BTO`.

Important limit: the frozen runner's top-level JSON saved the selected branch
but **not explicit final records for all three starts**. A separate
deterministic replay of the frozen numerical optimizer, changing only its
branch-outcome recording, made 78 cache lookups; every one hit an already
raw-audited point and **zero new DFT** calls occurred. It reconstructed the
exact start-to-terminal mapping `frozen→eval-18`, `+Q_y→eval-47`, and
`−Q_y→eval-76`, including the original selected start. These terminals are
also the only three cache points meeting both stationarity limits. The
mapping is now reproducibly reconstructed, although it was not serialized by
the original job. A separate geometry/provenance preflight for the selected
`+Q_y` branch checked all 32 signed fixed-Q probes at a difference step of
`0.05 sqrt(amu) Å`: the smallest atomic separation is `1.809681819 Å`,
the smallest volume is `67.685878084 Å³`, and every probe passed the
predeclared minimum-distance guard. Its immutable output is
`curvature_0p05_preflight.json` (SHA-256
`f9c5f423696fdbc25cf2002a6c0605387f013e6dd0d0a7c078e85867397ce4d7`).
The corresponding 32-static-evaluation Slurm job `27787714` was submitted
on hf `hfacnormal01` and finished `COMPLETED/0:0` after `00:20:29`.
The independent auditor re-read **all 109** cache points, including the 32
new signed probes, from original ABACUS input/SCF/force/stress files and
reconstructed the same 16D Hessian. Its lowest eigenvalue is
`+0.00543824 eV/(amu Å²)` at this one step, and its antisymmetric relative
defect is `0.000507`. However, the maximum energy-versus-force diagonal
curvature mismatch is `0.01623171 eV/(amu Å²)`, **larger than the lowest
eigenvalue**; the maximum energy-versus-force first-derivative mismatch is
`0.00081562 eV/(sqrt(amu) Å)`. Thus the 0.05 result is a one-step screen,
not a numerical certification of positive curvature. The output, 109-point
raw audit and Hessian reconstruction are the three `curvature_0p05_*.json`
files in this folder. The full-domain branch continuity and T-endpoint
off-plane residual are also unverified.

After the 0.05 raw audit passed, all 32 signed `±0.10 sqrt(amu) Å` probes
passed an independent geometry preflight (minimum separation
`1.796364193 Å`, volume `67.630117027 Å³`; `curvature_0p10_preflight.json`,
SHA-256 `1589a26a00c5b5d775be8af267ae00775e26192966acff7ddbe84b6c635fc6f0`).
The second-step Slurm job `27787746` completed on hf `hfacnormal01` at
2026-09-27 15:27 CST (`COMPLETED/0:0`; elapsed `00:23:47`). Its independent
audit re-read all **141** cached original DFT points, including 32 new signed
probes. The two force-Hessian lowest eigenvalues are `+0.00543824` and
`+0.00546533 eV/(amu Å²)`; their common-basis lowest eigenvectors have
absolute overlap `0.999999657`. The 0.10-step antisymmetric relative defect
is `0.001571`; maximum energy-versus-force diagonal-curvature mismatch falls
to `0.00402230 eV/(amu Å²)`, but the 0.05-step mismatch remains larger than
the weakest eigenvalue. First-derivative mismatches are `0.000816` and
`0.003056 eV/(sqrt(amu) Å)` at the two steps. The hash-linked comparison is
`two_step_curvature_screen.json`, reproducible via
`scripts/audit_bto_selected_two_step_curvature.py`. These data show a
**reproducible positive force-Hessian sign**, not a certified conditional
minimum: direct energy/gradient probes along the mixed lowest eigenvector,
branch continuity and the off-plane T endpoint still need checking.
Therefore the manuscript may report this local predictive success only with
these limits; it may not promote the staged four-point figure into a smooth
conditional-PES contour yet.

The follow-up four signed probes along the common soft direction at
`±0.05/±0.10 sqrt(amu) Å` passed a separate no-DFT geometry preflight
(minimum atomic separation `1.82211245 Å`, minimum cell volume
`67.69342590 Å³`). Its frozen `soft_direction_preflight.json` has SHA-256
`da7f47f2cfc47f23cf497257999939d0069dd7a2a30d55ac35bc467ae2b40096`.
A complete read-only runner input check then passed. The latter
verified the same calculator/cache contract before any new static evaluation.
The guarded hf `hfacnormal01` Slurm job is `27791225`; its launch script is
`cluster/hf_bto_transverse_soft_q075_q015_direct_probe_20260927.slurm`
(SHA-256 `10c9fe21d27b3f098e0a22c5d55b39a3845c2702ed1f32b22e398d01a960480d`).
It completed `COMPLETED/0:0` at 2026-09-27 22:13:57 CST after `00:04:34`.
The independent auditor re-read all **145** original cached ABACUS points,
including the four new static evaluations, and checked each probe against its
preflight geometry and original INPUT/KPT/STRU, SCF, energy, force and stress.
The frozen records are `soft_direction_result.json`,
`all_points_145_raw_audit.json`, `soft_direction_raw_audit.json`, and
`soft_direction_center_result.json`. Their hash-linked line-integral check is
`soft_direction_work_integral.json`.

| Step (sqrt(amu) Å) | Energy curvature | Directional-gradient curvature | Absolute mismatch |
| ---: | ---: | ---: | ---: |
| 0.05 | +0.00512544 | +0.00543225 | 0.00030682 |
| 0.10 | +0.00504960 | +0.00543009 | 0.00038050 |

Curvature entries are in `eV/(amu Å²)`. The five-point energy-versus-
directional-gradient Simpson work residual is at most `0.006404 meV` over
the two half-lines. The mixed direction has `97.67%` of its metric weight in
strain, so this is particularly relevant to the variable-cell coupling. Both
energy and gradient independently give a positive sign **in this selected
direction**. This is not a spectral error bound for the full 16D Hessian:
the larger diagonal energy/force discrepancy in other directions, competing
branches, and the T-endpoint off-plane residual remain open. Do not promote
this one-direction result to a certified conditional minimum or a 2D PES.

An older offline projection of the archived seven-image T→C chain already
measures the T endpoint's atomic residual outside these two cubic Γ axes as
`0.28996193 sqrt(amu) Å` (cubic endpoint ≈0), with T-end symmetric strain
`ηxx=ηyy=−0.00940965`, `ηzz=+0.05252416`. Its report and reference/force-
constant/eigenpair SHA-256 values match the sources frozen in this holdout's
soft-direction preflight. That projection currently resides at
`outputs/batio3_t_to_c_pbe100_dzp10au/bto_transverse_soft_plane_offline_audit_2026-09-26.json`
in the local ignored result tree; it is **not** yet a portable raw-data audit
in this folder. A nonzero frozen-plane residual does not rule out a
*conditionally relaxed* surface containing T, because the orthogonal
coordinates may relax. The T-end branch correspondence and portable source
freeze remain open before any path-overlaid conditional contour claim.

`conditional_q075_q015_result.json` is the frozen job's output. The
`audit-*.json` files are the canary, 77-point raw-output, and cache-only
branch-replay audits; the
original ABACUS logs remain under
`hf:/public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030/run-q075_q015`.
`holdout_analysis.json` is reproducible from these files and the locked
protocol with `scripts/analyze_bto_conditional_holdout.py`. The replay
script is `scripts/replay_bto_conditional_branches.py`; the corresponding
regression is `tests/test_bto_conditional_holdout.py`.
