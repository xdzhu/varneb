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
supply the original matched A/J. Prepared, not yet tested by new DFT.

One revision only in the initial protocol:1.0A guard,32A/16J+2checks per
original sign. Original electronic inputs, BFGS trust and numerical seed
are unchanged. Keep both signs, failed attempts, actual cost and original
labels. A successful local candidate is not a full-stability certificate.
