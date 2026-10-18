# E272 preparation: preserving-neighborhood signed branch controls

The measured zero-strain preserving curvature still has a negative atomic
direction after retaining the same four physical patterns used for T.
Its gradient component is about1e-6eV/A despite a released norm0.288eV/A.
Thus relaxation from the nearly symmetric origin alone can miss an omitted
branch. This is a nonstationary raw-matrix diagnostic, not a phonon frequency
or a certified saddle index. No new full Hessian is requested.

Prepare only the original0 and+1% preserving sampled peaks now. Keep their
ordinary path labels, raw input/logs and exact F energy/forces/stress. The
midpoint source choice remains pending the live E269 initialization control;
do not generate or admit a new midpoint result through this preparation.

Hold exactly the E270 physical definitions Qx,Qz,Qy and structural polar z
proxy, plus three translations. The seed amplitudes remain native to each
peak; this is not a fixed-Q strain derivative. Use the most negative A-space
eigenvector of the existing0.01A P raw matrix. Freeze its sign by making the
largest absolute oxygen component positive. The0.02A matrix checks the
diagnostic, but never replaces the preselected0.01A execution direction.
Transport the direction through the affine chart map and normalize in
target physical coordinates. Prepare just two starts,+/-0.10A, at each
peak's exact cell. No atom rematching, symmetry enforcement or random starts.
The measured curvature predicts only~0.022eV/A along this direction at0.05A,
already below the0.03 released-gradient target. Thus0.05A is not used as a
branch-discriminating start;0.10A is fixed before new DFT, not tuned to its
outcomes. The quadratic estimate is a rationale, not a force prediction.

The raw atomic block is indefinite, so T's whole-atomic-block numerical
initializer cannot be used unchanged. A prospective standard on-site Jacobi
initializer retains each atomic3x3 diagonal block and the open-cell3x3
block from the same0.01A source, removing only numerical off-site/mixed
blocks. Check positive definiteness without clipping/damping, then transport
it covariantly. All actual DFT couplings remain unchanged. This is numerical
initialization, not a fitted physical model or a new acceleration claim.

If admitted later, each signed A arm would have at most12 fresh attempts;
only its successful A may seed one J arm, also capped at12. Two source
neighborhoods x two signs x A/J gives at most96 optimization attempts. One start per arm,
gradient0.03eV/A, trust0.05A, full norm0.5A, distance1.6A, deformation0.25,
no retries. Report both signs and explicit failures, not just the lower one.
Candidate energies are not revised barriers or certified minima. A passed
gradient alone can still hide a negative direction: reserve exactly two
additional static probes,+/-0.01A along the transported original negative
direction, at each successful A. At most8 such attempts bring the total
cap to104. Geometrically forbidden probes remain unavailable, not repaired.
Their directional gradient secant checks this known direction only; even a
positive result does not certify the other28 released directions or a full
conditional minimum. No extra displacement amplitudes or retries follow
a negative result, and no whole Hessian is added.

This package is **inert preparation only**: no worker, deployment, submission
or new SCF is admitted by creating it. Both current study slots remain with
E269/E271. A later execution gate must preserve ABACUS/PBE100Ry/full10auDZP,
the original six asset bytes, Gamma2x2x2,P=E=0,32 native MPI,hf/hfacnormal01,
ordinary NEB0.10/noCI and all existing sources/costs. It must also verify
one free study allocation, a fresh namespace and exact held resources.
