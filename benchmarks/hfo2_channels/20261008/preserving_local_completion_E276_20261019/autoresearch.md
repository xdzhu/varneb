# Bounded experiment record

Metric: native released-gradient norm≤0.03eV/A for A and J, with the
original four physical observables retained and all native EFS replayed.
This is not a predictor-score or speedup contest. The slash-workflow tool
is unavailable in this session; existing shell/native readers implement
the bounded hypothesis→measurement→decision workflow.

H1: SCF failure caused the missing control. Rejected by24native successful
SCFs. H2: the old neighborhood/backtracking and finite budget obstructed
the local search. Supported by boundary approach, tiny accepted-trial
increments and large tangential residual; not proof of an exterior minimum.
H3: an explicit larger numerical window from the exact failed starts can
supply the original matched A/J. Supported by both original completed
handles: A 8/9 fresh calls, two probes each, J 3 each, 27 total. All four A/J
meet 0.03 eV/A after independent native replay. A lowers 17.172/17.175;
J adds 0.485/0.484 meV/fu. This is local finite convergence in the revised
window, not a global-minimum or full-stability claim. No further DFT
is required by this finite target.

One revision only in the initial protocol:1.0A guard,32A/16J+2checks per
original sign. Original electronic inputs, BFGS trust and numerical seed
are unchanged. Keep both signs, failed attempts, actual cost and original
labels. A successful local candidate is not a full-stability certificate.
