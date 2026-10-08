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
| Stable orthogonal relaxation enables useful reduced predictions | Curvature resolution and stable eliminated block; independent DFT holdouts | Extra unstable directions or branch crossings invalidate one chart | Analytic foundation and8-point fixed-cell T atomic pilot pass; joint-cell/channel prediction pending |
| Adaptive mode/branch information improves efficiency | Equal seeds/contracts, DFT-call and core-hour accounting, mechanism agreement | Faster iteration changes pathway or fails held-out prediction | Unmeasured |

The Schur complement and envelope theorem are established mathematics. Neither
is claimed as a newly invented theory. Strain effects, trilinear coupling,
nonpolar switching variants and VCNEB also have substantial prior literature.
Novelty must be established against the specific 2019–2026 works named in the
execution plan, including Zhou–Zhang–Rappe 2022 and Qi–Singh–Rabe 2025.

Qi–Rabe's distinct PRL135,046101 (2025) introduces a Cmma reference to organize
competing polymorphs through unstable phonon branches. It motivates a stronger
fixed-reference control: poor coverage in our selected T representation alone
does not establish that path adaptation is necessary. Reference choice and
orthogonal/branch response must be tested separately. The
[v2 prospective addendum](../../docs/HFO2_PREDICTION_PROTOCOL_V2_2026-10-09.md)
records this change after exploratory G1 analysis but before G2/holdout labels,
without replacing the original protocol or adding a production matrix.

Behara–Van der Ven (2022) additionally maps polymorph/variant paths in strain
coordinates, overlays a path on a two-shuffle energy landscape, and shows a
Pbcn-to-T intermediate change when the entire O cell is fixed. Neither path
projection nor a cell-constraint-induced mechanism change is our novelty claim.

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

The [connected main-text draft](MANUSCRIPT_DRAFT.md) now joins the physical
question, cited prior art, operative method definitions and dated pilot
results. Its missing boundary-response/independent-prediction evidence is
explicit; it contains no abstract or conclusion asserting H1/H2 success.
The [reading/build guide](README.md) provides the single-source LuaLaTeX
recipe and dated compile/visual receipt, without maintaining a duplicate
scientific body or claiming an evidence-complete PDF.
The newest [targeted neighbour check](../../outputs/HFO2_JCTC_RECENT_MODE_BOUNDARY_AUDIT_2026-10-09.md)
adds the 2026 phonon-pair domain-wall study and the 2025 Pbcn functional/boundary
benchmark. Neither multi-mode compensation nor boundary-dependent switching
alone is the proposed novelty. Their interface/full-cell quantities are not
equated with our homogeneous/partially clamped measurements.

## Prospective controls and prior-art boundary

The [prospective prediction protocol](../../docs/HFO2_PREDICTION_PROTOCOL.md)
is recorded before any G2 material labels. The +0.5% condition is an unseen
in-range test between 0/+1%, not extrapolation. A direct barrier-interpolation
baseline and both endpoint-response null models accompany frozen, atomic-only,
joint-cell and branch-aware reductions. Prediction files are frozen before
reading complete holdout path labels. Model abstention is reported as coverage,
not a successful numerical prediction. Two tested edges cannot certify the
whole candidate network minimum at the holdout condition.

Zhou (2022) already tests random starts in its constrained-mode subspace and
verifies mode-landscape switching predictions with NEB. Qi (2025)'s accepted
manuscript explicitly discusses an additional oxygen-crossing switching
category outside our current registered variants. These observations rule out
claims of first constrained-mode prediction/DFT validation or a globally
exhaustive four-edge network. The novelty test is quantitative advantage over
declared baselines for the matched-boundary switching/escape question, with
failure and coverage limits retained. See the primary-source
[audit and version limits](../../outputs/HFO2_JCTC_CLOSEST_PRIOR_ART_2026-10-08.md).

Definitions and the physical rationale for same-initial-state selectivity are
developed in [METHODS_DRAFT.md](METHODS_DRAFT.md). It is a theory/methods working
draft, not a completed Results section or an assertion of new theorems.

The first ordinary-residual-passed T--PO band now has an audited six-panel
energy/pattern/strain/reconstruction/Gamma-subspace/residual figure in
[PILOT_RESULTS_2026-10-08.md](PILOT_RESULTS_2026-10-08.md). It shows why good
local pattern coverage need not imply a complete global few-mode plane;
it is not a boundary-response prediction or a full-variable TS result, nor an
exclusion of more suitable fixed parent-mode representations.
