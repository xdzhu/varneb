# E274: finite atomic/cell crossover at the completed midpoint peaks

The completed E269 control changed preserving's sampled midpoint peak by
10.820844meV/HfO2 but T's by only-0.131935. Both internal patterns and open
cell rows differ. This batch measures the previously planned four-state
contrast, not a new prediction test, pathway search or relaxation.

A is the original E268 terminal peak; B is the E269 terminal peak at the
same +0.5% boundary. Reuse native AA=(s_A,H_A) and BB=(s_B,H_B) results.
Only AB=(s_A,H_B) and BA=(s_B,H_A) are new, for preserving and T: exactly
four frozen SCFs. Holding ordered unwrapped fractional positions while
changing H is affine strain, not fixed Cartesian positions. No atom
permutation, phase remapping, endpoint recalculation or relaxation occurs.
The native STRU writer canonically wraps periodic coordinates. The continuous
analysis lift is not mistaken for that serialized representation. Preflight
checks same-order periodic equivalence and registers the exact STRU bytes
produced by the immutable calculator's writer on HF; runtime must match them.

Use the original ABACUS/PBE100Ry, complete10auDZP, six physical asset bytes,
Gamma2x2x2, Hf4O8, P=E=0. Ordinary NEB0.10 is unchanged; this static batch
does not impose a force stopping threshold or enable CI. The immutable E268
runtime and binary are reused read-only. Original E269 handles must be
COMPLETED0:0. Native geometry/input and AA/BB cache checks precede submission.

One Slurm hfacnormal01 node, one manager task with32CPUs, four serial
mpirun-np32 calls, one thread/rank, one hour. At most two committed study
allocations; unrelated HZO jobs are untouched. Submit held, check resources,
record the exact accepted handle and then release. A saved intent or handle
prevents duplicate submission. No automatic retry or larger matrix.

The measured interaction is E_BB-E_BA-E_AB+E_AA. Report the internal
contrast at each cell and the affine-cell contrast at each internal pattern,
using250meV/HfO2 per eV/cell. Both intervention orders must sum to E_BB-E_AA.
Nonadditivity prevents a unique elastic/phonon energy partition. A missing
corner remains missing; no fit, imputation or replacement of frozen forecasts.
These are local frozen contrasts, not a new relaxed path or saddle certificate.

Retain every new native E/F/stress and cost. After the original C0 terminal,
archive portable inputs/logs and independently replay locally. Pseudo/orbital
bytes remain on HF. Stop this intervention after these four measurements,
including an uninformative result; do not expand merely to obtain a mechanism.
