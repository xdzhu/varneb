# E273: does the released preserving branch form a connected lower path?

This finite follow-up is motivated by completed E272 zero-strain minus
native data, not by a new material, tighter DFT settings or another local
surface search. Original labels and frozen forecast scores remain intact.

## Physical question and prior evidence

At the original zero-strain preserving peak, holding Qx/Qz/Qy and the polar
structural proxy, the minus branch lowers A by13.599258 and J by13.617082
meV/HfO2. Inserting J as **free image4** in the original nine-image polygon
leaves local peaks at3 and5. Its sampled maximum is85.499889meV/HfO2,
not96.828183-13.617082. All-image same-boundary NEB replay gives0.125611
eV/A, above0.10, with the largest forces on the neighbors. The center has
a local turn (adjacent-segment cosine-0.4560). Therefore this is neither a
smooth new MEP nor a corrected barrier. A connection test is needed.

## One predeclared continuation, unchanged physics

Only the zero-strain minus J polygon is admitted. Minus is the first
predeclared signed branch, not the lower of a post-hoc energy competition.
The original two endpoints, ordered atoms and continuous periodic lift
are unchanged. All seven interiors, including image4, evolve freely under
the original partially clamped boundary; the four local intervention
constraints are **not** imposed during NEB. No TS endpoint, climbing image,
midpoint branch, extra Hessian, fitted model or automatic sign restart.

ABACUS/PBE100Ry, all six original pseudo/orbital/input bytes, full10auDZP,
Gamma2x2x2, ordered Hf4O8, P=E=0. Cell scale5.12968067458423A, k0.2,
FIRE maxstep0.02A and cell step0.02, geometry distance1.6A/deformation0.25,
four geometry-step retries, ordinary NEB0.10eV/A. At most24updates,
168new interior SCFs and6h. The exact nine initial E/F/stress caches are
reused only at their exact ordered periodic geometries; normal movement
invalidates interior caches. Endpoints are never recalculated.

## Admission and stop rules

`prepare.py` is inert: it produces a geometry-and-cache proposal, not a
submission. E27228948425 must first finish with an unambiguous native
terminal result; a failed local optimization within its finite budget is
an outcome, not a reason to relaunch. Before execution, freeze a separate
runtime bundle, verify the original E268 source and all six assets, check
all nine exact caches and input roundtrips with fresh DFT prohibited,
then use a fresh HF namespace and a held, verified hfacnormal01 allocation:
one node,32CPUs, one manager, each image mpirun-np32, one thread/rank.
At most two study allocations including E269; no speculative third
dependent job. No running source may be edited or overwritten.

Stop on ordinary convergence, the24-update/6h bound, or a genuine native
failure. Do not stop for a single force rebound, silently alter physical
parameters, repeat the budget, or manufacture missing energies. Retain
the original and new full profiles, forces, geometry and actual cost.

## Conclusions allowed

If converged, compare the connected sampled maximum, peak locations,
oxygen-pattern history and open-cell response with the original preserving
path and the original T route at zero strain. A lower converged polygon
would demonstrate that the original sampled bottleneck was not robust
to the omitted accommodation direction. It would not establish a unique
MEP, a first-order stationary TS, continuous sampling convergence or a
predictive improvement. If the branch returns or remains unconverged,
report that result. No additional strain job is automatically admitted.
