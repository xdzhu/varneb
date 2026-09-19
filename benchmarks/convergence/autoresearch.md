# VARNEB convergence-acceleration experiment log

## Objective and acceptance criteria

Reduce the number of *new electronic-structure image evaluations* and wall time
needed to reach the project default `fmax = 0.10 eV/Angstrom`, while preserving
the converged path, barrier and calculator settings.  An optimizer is promoted
only if it improves at least two real first-principles paths and does not
introduce a material barrier/path discrepancy.

Primary metrics are outer VCNEB steps, cache misses (new SCF calculations),
node-hours, and wall time to first reach the threshold.  Force decrease per
step and final barrier differences are secondary diagnostics.  A single lucky
force dip does not count as convergence; the final chain is audited normally.

## Controlled cases

1. BaTiO3 T-to-C, 9 total images including endpoints, ABACUS PBE/100 Ry/DZP
   10 au.  Every candidate starts from the archived serial step-9 chain.  The
   historical global-FIRE continuation is the baseline (about 28 further
   steps to cross 0.10 eV/Angstrom).
2. HfO2 T-to-PO, 7 total images, ABACUS.  This is the planned transfer case;
   it is not used for tuning.
3. A VASP path with native candidate validation is the safety case.  It tests
   whether a rejected cell shrinks only the responsible image/block instead of
   collapsing the global FIRE time step.

## Strategies

- `FIRE`: unchanged ASE global-band baseline.
- `BlockFIRE`: independent FIRE velocity, time step and displacement cap for
  every interior image; candidate rejection is localized by image.
- `SplitFIRE`: as above, with separate atomic and cell blocks and a smaller
  optional cell trust radius.
- Next experiments after the first DFT gate: ODE12r-like residual integration,
  variable-cell preconditioning in a six-component symmetric strain gauge,
  dynamic freezing/evaluation of locally converged images, and coarse-to-fine
  image insertion.  These must be benchmarked separately so savings are not
  conflated.

## Current evidence

The analytic scaling screen uses the same perturbed path at 5, 9 and 17 total
images.  At 9 images, global FIRE took 33 steps and both block variants took 16
steps.  At 17 images, however, global FIRE took 47 steps, BlockFIRE 72 and
SplitFIRE 123.  Therefore per-image FIRE is a hypothesis with a demonstrated
mid-size benefit, not yet a generally valid acceleration claim.  The negative
17-image result is retained and will guide trust-radius coupling work.

A local residual-rebound guard was also prototyped.  Although it matched the
16-step n=9 toy result, it stalled the n=17 case at 0.2400 eV/Angstrom after
300 steps.  It was therefore rejected at the analytic gate and removed from
the runtime API before any DFT allocation.

A second prototype made each image trust radius proportional to its current
force residual, `min(0.02, 0.05*f_i)`.  It converged, but was slower than global
FIRE at every analytic resolution (18/59/85 versus 11/33/47 steps for
5/9/17 images).  It was likewise rejected before DFT submission and removed
from the runtime API.

Real ABACUS benchmark jobs submitted on 2026-09-20:

- Slurm 27728996: `BlockFIRE`, BTO n=9, per-image maxstep 0.02.
- Slurm 27728997: `SplitFIRE`, BTO n=9, atomic maxstep 0.02 and cell maxstep
  0.01.
- Slurm 27729012: global `FIRE` control with its band-wide maxstep scaled by
  `sqrt(7)` to 0.052915.  This separates block-state gains from a simple
  image-count correction to the global displacement cap.

Both use an isolated source snapshot and do not modify running production
jobs.  Results are recorded in `autoresearch.jsonl` rather than overwritten.

Transfer jobs submitted without retuning the strategy logic:

- Slurm 27729066: HfO2 T-to-PO n=7, global FIRE maxstep scaled by `sqrt(5)`
  from 0.02 to 0.044721.
- Slurm 27729067: the same HfO2 initial chain with SplitFIRE, atomic cap 0.02
  and cell cap 0.01.
