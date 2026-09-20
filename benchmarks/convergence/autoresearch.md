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

## Production-chain continuations

On 2026-09-20 the successful controls were applied to the preserved lowest-force
complete snapshots of the unconverged VASP chains.  The calculator, endpoint,
pressure, spring, image-count and `fmax=0.10 eV/Angstrom` contracts remain
unchanged; only the optimizer state and strategy were restarted.

- CdSe cell mapping n=17: BlockFIRE 0.01 from baseline step 49
  (`fmax=0.111587`), Slurm 27732728.
- GaN hexagonal n=29: BlockFIRE 0.01 from step 35
  (`fmax=0.104281`), Slurm 27732740.  An earlier launch 27732729 failed before
  any DFT call because its requested cache namespace did not match the preserved
  cache metadata; the replacement deliberately uses a fresh isolated cache.
- HfO2 T-to-PO n=20: StagedFIRE first reduced `0.179442` to `0.151528`, then
  the native VASP lattice checker rejected the next candidate before DFT.
  Continuation 27732743 retains that accepted chain and uses candidate-aware
  BlockFIRE 0.01.
- HfO2 PO-to-M n=20: StagedFIRE from baseline step 103, Slurm 27732731; its
  first completed step reduced `0.278429` to `0.271560`.

## 2026-09-20 VASP transfer and feasibility results

The production-chain tests now separate optimizer acceleration from input
feasibility and transient calculator failures.

- GaN hexagonal n=29 converged in one conservative BlockFIRE update,
  `0.104281 -> 0.098719 eV/Angstrom` (Slurm 27732740).
- CdSe cell mapping n=17 converged in four BlockFIRE updates from the common
  step-49 chain, `0.111587 -> 0.109890 -> 0.106575 -> 0.101830 ->
  0.095892 eV/Angstrom` (Slurm 27732728, 2:06:28 wall time).  The historical
  CheckedFIRE branch from the same chain instead gave `0.124067`, `0.117012`
  and `0.161471` over its next three updates.  This is a successful rescue,
  but the baseline was censored before a later crossing, so no unsupported
  wall-time speedup ratio is assigned to CdSe.
- HfO2 T-to-PO n=20 exposed a distinct failure of conservative BlockFIRE: it
  plateaued near `0.188 eV/Angstrom` through step 28.  Candidate-aware
  StagedFIRE, resumed from the best preserved complete chain, converged after
  one accepted half-step, `0.101143 -> 0.099347 eV/Angstrom` (Slurm 27735411).
  The full proposal failed the native Bravais-consistency gate at image 10;
  no electronic job was launched for it, the half proposal passed, and the
  final chain passed the normal audit.  This demonstrates why feasibility
  backtracking and an optimizer are complementary rather than interchangeable.
  An exact no-backtracking control from the same chain, calculator cache and
  proposal failed at that full step in six seconds (Slurm 27735454), while the
  backtracking run accepted the half-step and converged.  This is a controlled
  algorithm ablation, not an inference from two unrelated trajectories.
- HfO2 PO-to-M n=20 subsequently completed with StagedFIRE and one exact-input
  image retry enabled (Slurm 27735441).  Its force history was
  `0.214936, 0.211147, 0.204578, 0.194083, 0.180491, 0.163464, 0.152925,
  0.142694, 0.148453, 0.152878, 0.153952, 0.154380, 0.149539, 0.138817,
  0.123297, 0.104358, 0.103834, 0.102725, 0.101297, 0.099462`
  eV/Angstrom.  The 20-image chain passed audit with barrier `0.33449254 eV`,
  reaction enthalpy `-0.33276529 eV`, and maximum deformation `0.1741`.

The HfO2 VASP lattice tests also show that `SYMPREC` is not a monotonic
"precision" knob.  Full-chain native probes gave different failure sets as
the threshold changed.  At `1e-5`, preserved rejected T-to-PO candidates also
failed real VASP before the first electronic step.  At `1e-4`, two such
candidates completed, but three images from an accepted T-to-PO chain failed
real VASP before electronic iteration.  A case-specific fixed `3e-3` contract
passed both preserved 20-image-chain scans and four representative full-SCF
canaries; for all
directly comparable structures, energies and maximum forces were identical
to the tighter-threshold calculations.  Runtime threshold switching remains
forbidden, and this evidence does not change the general `1e-4` default.
Endpoint statics were regenerated from the exact resumed
trajectories so the coordinate hashes, not only periodic structures, match.

One PO-to-M image then hit a ScaLAPACK `pdstebz/pzheevx` segmentation fault
after 12 electronic iterations.  Repeating the exact same structure and
input completed in 55 seconds.  This is classified as a transient calculator
failure, not a geometry rejection.  The Hefei launcher now permits an
explicit finite `IMAGE_RETRIES`; retries repeat the identical calculation and
never change the physical contract.  Successful images remain reusable from
the content-addressed cache.
