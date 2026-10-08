# HfO₂ JCTC research design — not an evidence-complete manuscript

The execution plan and finite calculation matrix are in
[VARNEB_JCTC_HFO2_RESEARCH_PLAN.md](../../docs/VARNEB_JCTC_HFO2_RESEARCH_PLAN.md).
The existing software manuscript remains in `../VARNEB_CPC`; neither its title
nor its numerical claims are silently promoted to a new journal.

## Scientific thesis to test

Can mechanical boundary conditions decouple homogeneous polarization switching
from escape to a nonpolar phase in hafnia, and can a branch-aware mode–strain
reduction predict the channel-selective response rather than only fit it?

This is a hypothesis, not a result. Compare barriers from the *same* PO+
initial state. Treat T→PO formation separately (reverse of the measured
PO→T edge). A uniform-cell barrier is not a device retention time, coercive
field or a domain-wall nucleation barrier.

## Claim–evidence ladder

| Candidate claim | Required measurement | Null/negative outcome | Status |
|---|---|---|---|
| Boundary conditions change switching/escape selectivity | Same PO+, matched Hamiltonian, several channels and error-resolved barrier differences | All barriers move together, or differences below uncertainty | Unmeasured |
| The change is not explained by well shifts alone | Endpoint and transition-region responses, common reference, force–energy consistency | Endpoint-only model predicts equally well | Unmeasured |
| Stable orthogonal relaxation enables useful reduced predictions | Curvature resolution and stable eliminated block; independent DFT holdouts | Extra unstable directions or branch crossings invalidate one chart | Analytic foundation implemented; material validation pending |
| Adaptive mode/branch information improves efficiency | Equal seeds/contracts, DFT-call and core-hour accounting, mechanism agreement | Faster iteration changes pathway or fails held-out prediction | Unmeasured |

The Schur complement and envelope theorem are established mathematics. Neither
is claimed as a newly invented theory. Strain effects, trilinear coupling,
nonpolar switching variants and VCNEB also have substantial prior literature.
Novelty must be established against the specific 2019–2026 works named in the
execution plan, including Zhou–Zhang–Rappe 2022 and Qi–Singh–Rabe 2025.

## Proposed main-text sequence

1. A physical dilemma: switchability versus phase escape, not a backend survey.
2. NEB and joint atomic–cell paths; mechanical ensembles; modes and branches.
3. Local conditional reduction and its explicit instability failure boundary.
4. Audited hafnia variants and common-initial-state channel comparison.
5. Boundary response, coupling mechanisms and independently verified predictions.
6. Ablations, computational cost, limits, and reproducible VARNEB implementation.

Figures follow this causal order. Two-dimensional contours are measured and
validated local conditional branches, not decorative smoothing or assumed
global surfaces. Original backend-validation figures may support an appendix;
they cannot substitute for the missing hafnia predictions.

No abstract containing unmeasured scientific conclusions will be drafted yet.
Acceptance gates are evidence-based, not a fixed journal-acceptance promise.
