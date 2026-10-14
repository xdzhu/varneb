# E164 — atomic instability versus open-cell coupling at the two sampled peaks

Question: are the two negative directions of each existing local internal
matrix induced by releasing the open cell, or already present in the atomic
sector? This uses the complete E159 responses, not new DFT or post-relaxation
Hessians. Previous goal turn E163 was PROGRESS (native prefix audit, 218c913b).
This is a descriptive analysis chosen after seeing training spectra, not a
preregistered independent prediction or a new H1/H2 acceptance test.

## Observation in the declared raw quadratic

At both centers and both finite-difference steps, the 33-dimensional atomic
restriction already has two negative eigenvalues. The three-dimensional
open-cell restriction is raw positive. Coupling changes their magnitudes,
but does not add a negative direction in these measured matrices. The
counterexample tests explicitly allow coupling-induced extra instability;
the material observation is not hardcoded into the analysis.

| Center | Step, Angstrom | Lowest atomic-only | Lowest joint | Minimum open-cell | Negative-subspace cell trace |
|---|---:|---:|---:|---:|---:|
| T | 0.01 | -0.240169 | -0.276193 | 2.956214 | 0.0205773 |
| T | 0.02 | -0.237806 | -0.273596 | 2.956557 | 0.0206105 |
| Preserving | 0.01 | -1.271818 | -1.271819 | 3.090798 | 1.04022e-06 |
| Preserving | 0.02 | -1.281815 | -1.281821 | 3.092172 | 3.09248e-06 |

Curvatures are eV/Angstrom squared in the fixed joint length metric, not
phonon frequencies. The last column is Tr(P_cell P_negative), dimensionless
for the two-dimensional negative subspace, not an energy percentage. It is
invariant under orthogonal re-expression of the internal basis and under
rotation of vectors within a degenerate negative eigenspace. It depends on
the declared physical atomic/cell split and cell scale L=5.12968067458423
Angstrom, not an invariant under arbitrary rescaling of strain coordinates.

The T peak has an observable open-cell admixture and a more negative joint
lowest curvature than its atomic restriction. The preserving peak's raw
negative subspace is almost entirely atomic in this metric. Nonetheless,
the full atomic-cell coupling operator norms are approximately 7.72 and
8.10 eV/Angstrom squared: a large total coupling norm alone does not identify
which unstable direction is coupled. Nor does a small cell weight imply
insensitivity of an entire barrier to the prescribed substrate strain.

This narrows the mechanism analysis: the extra raw negative direction is
not merely an elastic instability introduced by the variable-cell degrees
of freedom. Multiple atomic directions must be considered. It does not
establish a physical stationary index, a lower complete-path barrier or a
cause for the observed strain response. The centers are nonstationary
ordinary sampled peaks, and the finite-step agreement does not bound DFT error.

## Algebra, nonzero gradients and limitations

The exact registered internal embedding excludes prescribed substrate control
and three translations. Its physical coordinate rows define the atomic
rank-33 and open-cell rank-3 subspaces; no fixed list of internal eigenvector
indices is assumed. Raw responses are retained at their original source;
the analysis explicitly uses the symmetrized internal block.

For each negative-space projector, the atomic, cell and cross Rayleigh
traces sum to the negative eigenvalue sum. These are directional curvature
identities, not contributions to the NEB barrier energy. Individual eigenvectors
are not exported or called phonon modes.

Only for raw-positive Hcc, the formal quadratic cell-stationary reduction
reports Haa-Hac Hcc^-1 Hca. Its nonzero linear term is retained:
c*(a)=-Hcc^-1(gc+Hca a), ga_eff=ga-Hac Hcc^-1 gc and
delta_E0=-gc^T Hcc^-1 gc/2. These are arithmetic properties of this local
quadratic, not additional SCF results, physical error bounds or certified
nonlinear conditional response. A zero/negative cell block is never clipped
or pseudoinverted. Full released-block stability at the newly relaxed
candidate structures remains unmeasured.

## Reproduction and scope

From the repository root, run:

```text
python -m scripts.analyze_hfo2_G3_atomic_cell --output a-new-output.json
python -m pytest -q tests/test_hfo2_G3_atomic_cell_signature.py
```

Existing outputs are refused. `analysis.json` pins the original manifest,
native archive and both losslessly archived native report hashes. The
source performs zero new DFT, scheduler action, raw SCF reparse, reference
selection, ASR or input change; it reads no private author arrays or physical
+0.005 holdout labels. All four center/step combinations remain in the result.

The original two 25-point conditional batches continue separately. This
analysis does not alter their axes, domain, trust radius, force target,
physical parameters, runtime or resource limits. Stable branches, independent
errors and the registered model comparison are still required for the full goal.
See [complete native curvature](../G3_complete_curvature_E159_20261014/README.md),
[actual local release](../G3_conditional_canary_E160_20261014/README.md),
[candidate-sheet audit](../G3_face_prefix_audit_E163_20261015/README.md) and
[analysis source](../../../../scripts/analyze_hfo2_G3_atomic_cell.py).

The related fifteen-module suite passes 278 checks locally and 278 in the
independent archive, with zero failures/errors/skips; 23 checks are new.
They include an explicit coupling-induced-extra-negative-direction
counterexample, nonzero cell-gradient bias, nonpositive-block refusal,
orthogonal-basis/degenerate-subspace invariance and all four actual rows.
These checks validate arithmetic and provenance, not a physical sign floor.
The independent command reproduces the complete saved analysis byte-for-byte
(SHA ef7798df3e71fbbc3e86bd468b608e4bf732cb13cd048e1518b8608d4753db3a).
Tested tree af3de69dab2f46a90024c893bdadaac5d34567cb and exact source/archive
hashes are retained in validation.json. This is not a new whole-repository
regression or raw SCF re-audit.
