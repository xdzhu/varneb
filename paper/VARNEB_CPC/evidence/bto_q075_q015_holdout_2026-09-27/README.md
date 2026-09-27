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
the original job. The holdout Hessian, full-domain
branch continuity, and T-endpoint off-plane residual are also unverified.
Therefore the manuscript may report this local predictive success only with
these limits; it may not promote the staged four-point figure into a smooth
conditional-PES contour yet.

`conditional_q075_q015_result.json` is the frozen job's output. The
`audit-*.json` files are the canary, 77-point raw-output, and cache-only
branch-replay audits; the
original ABACUS logs remain under
`hf:/public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030/run-q075_q015`.
`holdout_analysis.json` is reproducible from these files and the locked
protocol with `scripts/analyze_bto_conditional_holdout.py`. The replay
script is `scripts/replay_bto_conditional_branches.py`; the corresponding
regression is `tests/test_bto_conditional_holdout.py`.
