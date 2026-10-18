# E271: one bounded numerical repair of incomplete E270 treatments

The E270 identity-initialized safeguarded-BFGS outcomes and costs remain
unchanged. Its zero-strain A exhausted12 attempts; its midpoint J did too.
Converged native SCFs and the measured first-step curvature diagnose a
conditioning problem, not a reason to change electronic inputs.

Only after the original E270 job completes normally may this batch run.
Reuse every successful A/J candidate, after exact-coordinate native replay;
repair only unavailable A or J. If A was unavailable, J may run once after
a repaired A succeeds. No A candidate means no J attempt. One new namespace,
one start per arm,12 new attempts per arm including failures, at most72
new attempts overall. No retry after this batch without a new scientific
decision; a numerical failure remains a result.

The sole optimizer change is its initial Hessian. Freeze the existing
zero-strain T raw-curvature probe at0.01A, using its atomic and open-cell
diagonal blocks as a block-Jacobi numerical seed. Transport this matrix to
each target chart by the already checked affine Jacobian. Do not fit,
clip eigenvalues, add damping, measure a new Hessian, or interpret the seed
as the strained physical Hessian. The0.02A source remains a sensitivity
diagnostic, not an alternative selected using the new outcomes. The actual
DFT mixed coupling is unchanged. Count earlier curvature preparation and
failed E270 attempts separately; this is not a demonstrated speedup.

Keep all four physical pattern values, three translation gauges, exact
fixed-cell A embedding and A-to-J sequence from E270. The definitions,
not their native values, match across0,+0.5,+1%. Retain its released-gradient
target0.03eV/A, initial trust0.05A, total norm0.5A, distance1.6A,
deformation0.25 and iteration cap200. Do not alter E270 source or outputs.

ABACUS/PBE100Ry/full10auDZP, the same six asset bytes, Gamma2x2x2,
P=E=0 and native mpirun32 are fixed. Ordinary NEB remains0.10eV/A, no CI.
Use only hf/hfacnormal01,1node32CPUs/6h and at most two committed study
allocations including E269. Preflight without DFT, commit/push the protocol,
submit held only after E270 is COMPLETED0:0, inspect exact resources, then
release. Preserve intent on ambiguous submission; never blindly resubmit.

Report original and repair outcomes, local energies, residuals, retained
variable defects and cumulative native cost. These are conditional gradient
candidates, not new path labels, certified minima/TSs, a unique energy
partition or an independent prediction. The remaining preserving-channel
question is not answered by this T-only batch.
