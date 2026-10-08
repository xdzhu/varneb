# HfO₂ mode–strain channel study: bounded experiment ledger

Execution plan: `docs/VARNEB_JCTC_HFO2_RESEARCH_PLAN.md`.
Date: 2026-10-08. No physical result is inferred from a unit-test fixture.

## E001 — Stable eliminated-subspace curvature

- Hypothesis: a stable orthogonal harmonic block gives a well-defined local
  conditional curvature; unstable/unresolved blocks invalidate this reduction.
- Change: new calculator-independent `vcneb/relaxed_curvature.py`.
- Measurement: analytic quadratic minimization, retained saddle, rotation
  covariance, clamped omitted directions, invalid/marginal blocks.
- Decision: retain analytic foundation: 23 focused tests and 684 full-worktree
  tests passed (one skipped); physical application remains pending. Added
  historical-summary normalization tests are reported separately below.
- Novelty boundary: the Schur complement is established, not a new theorem.

## E002 — Same-input replica of two historical HfO₂ peaks

- Hypothesis: historical ordinary and guided peak results are reproducible under
  genuine `mpirun -np 32`, without changing INPUT/KPT/STRU/pseudopotentials/orbitals.
- Change: new isolated Slurm array, two statics only; source images read-only.
- Measurement: source/replica hashes, one DSIZE=32, SCF evidence, all 12 forces and
  six ASE stress components, energy/force/stress differences, timing, affinity.
- Expected precision: 1e-5 eV/cell, 1e-4 eV/Å force, 0.02 kbar stress. These are
  reproducibility tolerances, not new endpoint/NEB convergence requirements.
- Decision: pending real SCF. A pass does not certify either peak as a TS or
  identify which nonpolar mode changed sign.

## Existing evidence and non-results

- ABACUS ordinary-7/9: 39.187725/39.0523 meV/f.u.; historical CI-7: 32.29305.
- Guided-7: 0.963925 meV/f.u., fmax=0.0866792767, passes 0.10 but not original 0.05;
  sub-meV barrier resolution and exact variant identity are not yet established.
- VASP 27727755/27727756: Slurm FAILED/exit1 after about 9.5 hours. Actual wrapper
  output reports `step_limit_reached`, fmax 0.2030435724/0.2784290788. This is a
  nonconvergence audit exit, **not a newly observed VASP crash/Bravais error**.
  Provisional barriers 0.01123034/0.44605894 eV/cell are not final production results.
- Source code of submitted jobs will be a clean, archived commit. Dirty user
  BTO docs/figures/scripts are excluded. Branch push does not run PyPI workflow.

## Reproduce locally

`python -m pytest -q tests/test_relaxed_curvature.py tests/test_joint_curvature.py tests/test_hfo2_static_replica_audit.py`

Submission metadata, real job handles and audited results are appended to this
directory at each evidence-bearing milestone; no repeated/duplicate submissions.
