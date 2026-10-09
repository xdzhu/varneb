# Fixed-control continuation of a restricted stationary branch

`vcneb.stationary_branch.stationary_quadratic_branch` supplies one missing
numerical step between the measured same-centre quadratic response and the
registered B2--B5 comparison. It is standard harmonic algebra, **not** a new
theorem, a branch-discovery algorithm or an independent material prediction.
It calls no calculator and does not authorize G2/G3 or held-out calculations.

## Internal freedom and an external control are different

Use one declared reference, metric, Hamiltonian and mechanical boundary for
E0, the physical gradient g and Hessian H. Q contains retained internal
directions, R contains only selected stable released directions, and c is one
prescribed control. Q/R/c are orthonormal. Q and R must lie in the explicit
`admissible_internal_basis`, and c must be orthogonal to that actual mechanical
space. For the biaxial chart, use its translation-free internal basis and
tagged control direction. Do not silently substitute free-cell strains.

The existing conditional quadratic first minimizes in stable R, retaining
its nonstationary offset and reference relaxation energy. Write the resulting
model as a function of retained internal q and prescribed t:

```text
E(q,t)-E0 = dE0 + gq.q + gt*t + (q.Kqq.q + 2*q.Kqt*t + Ktt*t*t)/2
q*(t) = -solve(Kqq, gq + Kqt*t)
dE_branch/dt = gt + Ktq.q*(t) + Ktt*t
d2E_branch/dt2 = Ktt - Ktq.solve(Kqq,Kqt)
```

All eliminated offsets follow the same model, not a stationary T Hessian
borrowed for a different centre. `expected_index=0` requires a resolved
minimum in the represented internal space; `1` requires exactly one resolved
negative retained curvature. Stable R adds no negative directions. The
controlled coordinate never enters that index, even if its curvature is
negative. Unresolved, wrong-index or unstable eliminated blocks are rejected;
no pseudoinverse or automatic branch repair is used.

Unlike elimination of stable variables, stationary continuation along a
negative retained direction can **increase** the control curvature. Therefore
stable-mode softening alone does not determine a bottleneck response. A
barrier response still requires a compatible initial-well model and the same
mechanical control; neither a reference-well shift nor one local branch is
by itself a material barrier prediction.

## Outputs and restrictions

The point keeps t prescribed and reports internal coordinates, selected
released coordinates, full displacement, energy change relative to E0 and
control work. It distinguishes the represented internal residual, the
unrepresented but physically admissible internal gradient, and mechanically
clamped reactions. A zero represented residual does not certify all other
degrees of freedom. Condensed Q eigenvalues are not phonon frequencies.

For the [biaxial coordinate contract](BIAXIAL_CONTROL_CURVATURE.md),
t=L*(epsilon-epsilon0), so dE/d epsilon=L*dE/dt and
d2E/d epsilon2=L**2*d2E/dt2. L is a registered geometric metric, not a tunable
cutoff or Hamiltonian. No mass/strain/unit conversion is implicit.

The algebra cannot establish probe-convex-hull coverage, anharmonic validity,
full transverse stability, continuity to a different branch, stationary DFT
geometry or barrier uncertainty. Check those separately using actual training
data. Ordinary NEB residuals and spring forces are not Hessian inputs. An
unconverged G1 peak must not be turned into a certified TS by this solver.

## Reproduction and status

From the repository root, choose a new output path:

```sh
python -m scripts.check_stationary_branch --output /new/analytic-check.json
python -m pytest -q tests/test_stationary_branch.py tests/test_quadratic_reduction.py tests/test_relaxed_curvature.py tests/test_biaxial_curvature.py
```

The bounded checker uses eight analytic combinations of index0/1,
frozen/stable release and rotated charts,24 prescribed points. It compares
with a separately assembled fixed-control solve and finite energy derivatives.
The tests additionally cover nonzero offsets, actual biaxial geometry,
external-scale chain rules, forbidden mechanical directions, unresolved
curvatures and input ownership. These are implementation checks, not HfO2
curvatures, B2--B5 prediction/freezing completion, model advantage or JCTC
readiness. The complete prospective controls and independent material labels
remain required. See the [E031 delivery](../benchmarks/hfo2_channels/20261008/stationary_branch/README.md).

The [paired-response interface](STATIONARY_GAP_RESPONSE.md) now combines an
initial minimum and bottleneck model under matching declarations, retaining
both anchor corrections and distinct external-coordinate scales. Its signed
local gap is not a certified activation barrier or completed material forecast.
