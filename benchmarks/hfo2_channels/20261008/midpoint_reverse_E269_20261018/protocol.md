# E269: upper-to-midpoint initial-path sensitivity

This is a finite mechanistic follow-up after seeing the E268 outcomes.
It is not an unseen-strain prediction, a new method win or a replacement
for the original E268 protocol, labels, frozen forecasts and scoring.

## Fixed scientific comparison

Compare two +0.5% continuations, preserving reversal and PO-to-T, to their
completed E268 zero-to-midpoint counterparts. Initialize only from the
corresponding audited ordinary-converged +1% polygons. Preserve ordered
atoms, channel lifts and chart; use the same E268 target endpoint geometries
and exact static caches. New interiors receive fresh SCFs, never old results.

First impose the target substrate plane. Distribute fractional endpoint
corrections and cell corrections linearly over the original nine indices,
using the already tested E268 geometric recipe in the reverse strain direction.
The target features are pinned by SHA256
0441a5e7e63ec86a68807068ea21db6901f4a762d740592473d91bf15bbd8f23.

## Unchanged calculation and limits

ABACUS/PBE, original six input/pseudo/full-10au-DZP bytes, ecutwfc 100 Ry,
Gamma 2x2x2, Hf4O8, P=0 and zero field. Original partially clamped plane,
third vector plus tilts free. Nine total images, seven interior; FIRE,
fmax 0.10 eV/A, k 0.2, maxstep .02 A, cell step .02, four geometry-step
retries, minimum pair distance 1.6 A and deformation .25. No climbing image.
At most 80 updates and 24 h per path. Do not stop for one force rebound.

Only two new Slurm allocations on hf/hfacnormal01. Each requests one node,
one controller task with 32 CPUs; each image launches native mpirun -np32,
one thread per rank. Images run serially, paths concurrently. No endpoint
recalculation. Other projects' queue entries remain untouched. No auto-retry
or source overwrite. Existing queue/intent/handle prevents duplicate submission.

Reuse the original immutable E268 source, verified against all 375 recorded
file hashes. New controls and seeds are separately pinned and deployed in
a fresh E269 namespace. Before release, check all eighteen native decks,
the four endpoint cache uses and exact held resource metadata, with fresh
DFT disabled. No CI, installation, publication or workflow release is involved.

## Report and decide

Retain both initializations, not only the lower peak. Compare full energy
profiles, ordinary residuals, updates, peak locations, signed oxygen patterns,
open-cell response and actual cost. Report T-minus-preserving for each
initialization and separate per-channel spreads. The original midpoint
pair improvement is 4.298861 meV/fu and preserving interpolation deviation
12.595430 meV/fu: a comparable observed start sensitivity prevents assigning
them uniquely to intrinsic continuous-barrier curvature. No artificial
statistical errors or barrier bounds follow from the ordinary residual.

Agreement only closes this start-direction check; it does not certify the
continuous MEP, a stationary bottleneck, endpoint stability, an exhaustive
network or a stable reduced surface. If paths disagree, resolve the local
structural difference before choosing any additional calculation. If a job
hits its bound or fails, retain its actual result and diagnose; no automatic
new submission or parameter change. The six-neighborhood follow-up remains
a scientific decision, not an automatic sampling queue.
