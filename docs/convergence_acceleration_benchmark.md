# Convergence acceleration: design and benchmark protocol

## Why the current long-chain runs can stall

The default ASE FIRE optimizer represents the complete interior band as one
vector.  Its `maxstep` therefore limits the Euclidean norm of all image and
cell displacements together.  More images reduce the displacement available
to an individual image.  For `M` similarly moving interior images, a global
cap `Delta` implies an approximate per-image cap `Delta/sqrt(M)`.  An
image-count-neutral global control therefore uses `Delta = delta*sqrt(M)` for
a desired per-image scale `delta`.  In addition, the original candidate backtracking
multiplied the one global time step by every accepted shortening factor.  A
few difficult cells could reduce the effective time step by many orders of
magnitude while the manager continued launching new electronic calculations.

## Implemented experimental optimizers

`BlockFIRE` assigns independent FIRE state and a displacement cap to each
interior image.  `SplitFIRE` additionally separates the atomic and cell
generalized-coordinate blocks.  When the pre-calculation candidate validator
rejects a cell, only blocks belonging to the reported image are shortened and
reset.  Other images retain their momentum and time step.  A finite time-step
floor prevents irreversible numerical collapse; it does not alter any DFT
parameter or accept an invalid geometry.

`ImageScaledFIRE` keeps a global FIRE state but interprets `maxstep` as a
per-image scale and applies the `sqrt(M)` correction automatically.
`StagedFIRE` starts with that normalized coarse cap, then switches at a
configurable force (`0.12 eV/Angstrom` by default) to a `0.02` global cap and
resets momentum.  The stage check reuses the current force and costs no extra
calculator call.  The equivalent CLI controls are:

```bash
--optimizer ImageScaledFIRE --maxstep 0.02
--optimizer StagedFIRE --maxstep 0.02 --switch-fmax 0.12 --refine-maxstep 0.02
--optimizer BlockFIRE --maxstep 0.01
--optimizer SplitFIRE --maxstep 0.02 --cell-maxstep 0.01
```

Plain `FIRE` remains unchanged for backward compatibility and controlled
comparisons.

These names are intentionally experimental.  Image-local FIRE itself is not a
novel concept.  A publishable VARNEB contribution would require the full,
validated combination of variable-cell block metrics, feasibility-aware local
trust regions, distributed active-image evaluation, and auditable fallback
semantics, supported by first-principles benchmarks.

## First-principles BTO ablation

All rows start from the same archived T-to-C n=9 step-9 chain and use the same
ABACUS PBE/100 Ry/DZP settings.  Time starts at the common initial force
evaluation and ends at the first `fmax <= 0.10 eV/Angstrom` chain.

| Optimizer | Step cap | Steps | Time (s) | Estimated image evaluations | Final fmax |
|---|---:|---:|---:|---:|---:|
| global FIRE baseline | global 0.02 | 28 | 2669 | 205 | 0.096704 |
| BlockFIRE | per-image 0.02 | 10 | 963 | 79 | 0.086473 |
| SplitFIRE | atom 0.02, cell 0.01 | 9 | 863 | 72 | 0.099548 |
| global FIRE, image-count scaled | global 0.052915 | 8 | 772 | 65 | 0.078068 |
| global FIRE, matched 0.01 scale | global 0.026458 | 9 | 869 | 72 | 0.090979 |
| BlockFIRE, conservative | per-image 0.01 | 5 | 478 | 44 | 0.089923 |

Every completed accelerated chain passed the geometry audit and retained the
same monotonic endpoint-defined barrier, 0.0871289017 eV.  The conservative
BlockFIRE row reduced threshold time by 82.1% and estimated image evaluations
by 78.5%.  The matched global-cap control reached 0.108829 at step 4, rebounded
to 0.202851 at step 5, and converged only at step 9.  In contrast, conservative
BlockFIRE decreased monotonically and converged at step 5.  The comparison
therefore isolates a benefit from block-local FIRE state beyond the global
step-cap rescaling.  HfO2 remains the transfer test.

## HfO2 transfer result

The archived 7-image T-to-PO chain was used without retuning electronic
settings.  The historical global-FIRE baseline first crossed 0.10 at step 42.

| Strategy | Steps to threshold | Time basis (s) | Estimated image evaluations | Final fmax | Barrier (eV) |
|---|---:|---:|---:|---:|---:|
| global FIRE baseline | 42 | 5242 | 217 | 0.094173 at first crossing | 0.159677 final |
| one-stage global `sqrt(5)` cap | 30 | 3818 | 157 | 0.069988 | 0.162189 |
| coarse `sqrt(5)` through step 23, then 0.02 | 23 + 1 | 3513 | 134 | 0.092367 | 0.165150 |

The one-stage correction reduced outer iterations by 28.6%, threshold time by
27.2%, and image evaluations by 27.6%.  The staged trust-radius schedule
reduced estimated image evaluations by 38.2% and the reported time basis by
33.0%.  Its barrier differs from the tighter historical final chain by 5.5
meV, while the reaction energy is identical and the geometry audit is valid.
The staged chain also has a higher minimum adjacent-segment cosine (0.672
versus 0.581), so the saving was not obtained by accepting a folded path.

SplitFIRE was stopped after a fixed 12-step transfer window because its force
was 0.322028 versus 0.279179 for the scaled global control.  A later resume
reached only 0.241013 after two further steps and was stopped after the staged
and one-stage controls had completed.  Conservative BlockFIRE at 0.01 closely
matched the historical baseline for its first three steps and was likewise
stopped as dominated.  These negative transfer results remain in the raw log.

## Benchmark rules

- Identical initial chain, endpoints, calculator, pressure, spring, image
  count, climbing policy, and force threshold for every optimizer.
- Report total images and interior images explicitly.
- Count new SCF evaluations/cache misses, not only outer optimizer steps.
- Compare force histories, final image enthalpies, barrier, path geometry and
  any candidate rejections.
- Tune on BTO; transfer unchanged settings to HfO2 and one VASP case.
- Preserve regressions and failed hypotheses.  Do not select only the fastest
  successful seed.
- Do not call a method accelerated until at least two DFT cases improve in
  both electronic work and wall time without a material path change.

## Follow-up strategy matrix

| Strategy | Intended bottleneck | Required evidence |
|---|---|---|
| Per-image/block FIRE | global norm cap and global time-step collapse | BTO, HfO2, VASP safety case |
| ODE12r residual control | non-conservative NEB force and oscillation | same starts versus FIRE |
| Atomic/cell preconditioner | stiffness mismatch between modes | reduced iterations with stable barrier |
| Dynamic image activation | repeated SCF on locally converged images | fewer actual image evaluations |
| Coarse-to-fine insertion | high image count early in relaxation | equal resolved saddle at lower node-hours |
| Variable-cell IDPP/surrogate start | poor interpolated path | lower initial force and total DFT work |
