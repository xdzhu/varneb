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
- Decision: ordinary peak passed real SCF in job `28237630_0` (node11), exact
  energy/force/stress agreement and SCF wall time 128.43 s. Job `28237630_1`
  failed at Hydra bootstrap on node26 before ABACUS began; no electronic result
  exists for that element. Preserve that failure, retry only the missing case.
  A pass does not certify either peak as a TS or identify which nonpolar mode
  changed sign.
- Launcher-only repair: within a one-node 32-CPU Slurm allocation use Intel MPI
  `I_MPI_HYDRA_BOOTSTRAP=fork`, then verify 32 ranks/affinity and the raw SCF
  again. The physical files stay byte-identical. This is documented by Intel:
  https://www.intel.com/content/www/us/en/docs/mpi-library/developer-reference-linux/2021-15/hydra-environment-variables.html
- Full local worktree: 690 passed, one skipped after normalization tests.
- Missing guided point completed as `28238926_1` (node11, single-node fork):
  genuine DSIZE=32, SCF 123.37 s, exact source energy/force/stress agreement.
  Both historical peak results now pass static reproducibility. This excludes
  replica launcher/input drift as the cause of their energy difference, but
  does not identify the channel modes or establish sub-meV barrier accuracy.

## E004 — Atomic/strain local work probes

- Eight new SCFs only: atomic component of the local chain secant and scaled
  symmetric xx+yy strain, both at ±0.01 and ±0.02 Å chart steps. P_ext remains 0.
- Geometry changes are the intended variable. INPUT/KPT/pseudopotential/orbital
  hashes remain identical to the historical 100-Ry/10-au baseline.
- Preparation requires actual HF ASE-ABACUS STRU roundtrip, ordered Hf4O8,
  positive volume and minimum distance >1.6 Å before any submission.
- Measurements: finite-energy derivative versus force/stress derivative, paired
  energy curvature versus gradient curvature, and two-step dependence. These
  are local nonstationary diagnostics, not Γ phonons or full-space TS certification.
- Budget: eight points, one node ×32 MPI each, max three simultaneously, 30-min
  wall-time cap per element. No endpoint optimization or full landscape rerun.
- Decision: pending true DFT. All sources must be immutable commit archives.
- Submitted array `28240257` from source commit `8923165`: first three complete
  with SCF/input audits at 12:25:48 Beijing; three running and two waiting for
  the declared array concurrency cap. First-wave element wall time 154 s.
  Estimated remaining compute completes around 12:30–12:34 if timing holds;
  this is an estimate, not a claim of TS/path convergence.
- Local full-worktree regression after the pilot code: 696 passed, one skipped.
- Final 12:31:38 check: all eight COMPLETED/exit0; scorer freshly rechecks every
  input hash and SCF/energy/force/stress, `work_probe_summary.json` is complete.
- Atomic energy derivative versus center force derivative differs by
  0.000573/0.001418 eV/Å (0.01/0.02 Å steps); atomic energy curvature
  4.77348/4.73788 eV/Å², force-gradient curvature 4.67123/4.67278.
- Biaxial scaled-strain derivative differs by -0.00000434/-0.000484 eV/Å;
  energy curvature 19.13423/19.05207 eV/Å², stress-gradient curvature
  19.07036/18.95555. Maximum curvature method spread is 0.10225 eV/Å².
- Decision: retain a numerically usable *directional* pilot; carry finite-step
  and energy/gradient spreads forward as resolution evidence, not zero error.
  Both probed curvatures are positive; unprobed mixed/orthogonal directions
  remain unknown and center gradient is nonzero. Therefore no TS index,
  Schur-complement material prediction or sub-meV barrier certification is claimed.

## E003 — Ordered geometry comparison

- Historical final snapshots: ordinary step30 and guided step300, seven total
  images each, downloaded read-only from the source runs.
- Both ordered initial/final structural hashes match exactly; nearest-link
  final periodic unwrap integers also match. Therefore this audit does not
  show an endpoint permutation or different final lattice winding.
- Maximum translation-free atomic RMS separation at equal normalized geometric
  arc is 0.1695786428 Å. The interpolation here aligns geometries only; it does
  not create DFT points, phonon labels, a topology proof or a new TS.
- Measurement: `geometry_comparison.json` contains complete per-image cells,
  displacement/translation records and source file hashes; together with the
  endpoint fractional coordinates the ordered structures can be reconstructed.
- Decision: retain distinct channel candidates; mode identity still pending.

## E005 — Parent operation registration and T Gamma/seed batch

- Geometric question: do two explicitly mapped inversion operations produce
  opposite-polar-distortion candidates without erasing nonpolar site identity?
- T/PO/M seed symmetry is stable at symprec 1e-4/1e-3/1e-2 Å, angle tolerance
  1 degree: SG137/29/14; ordered Hf4O8, positive cells and safe distances.
- Two candidates retain SG29 and reverse the geometric Gamma sector; parent
  site permutations are explicit, never inferred from nearest distorted sites.
- The PO projection on the historical T-pattern vector is only 5.7e-7 Å.
  Decision: reject using this scalar as an X2- sign or channel label. Register
  all Gamma/Xx/Xy/Xz translation sectors instead; these are not individual irreps.
- New declared budget (separate from completed ten-point G0): 16 central T-cell
  force calculations, two displacement amplitudes, 1x1x1; plus two inversion
  candidate statics and one literature M-seed static = 19 new SCFs.
- Parameters: exact historical INPUT/KPT/pseudo/orbital files. Every new STRU
  must roundtrip through the HF ASE-ABACUS reader before any launch.
- Measurements: complete point hashes/SCF/DSIZE32/energy/forces/stress; raw
  versus ASR/permutation-symmetrized force constants, two-step frequency/FC
  differences, translation-identified acoustics, numerical endpoint equivalence.
- Negative Gamma modes are fixed-cell local evidence, not a variable-cell TS
  certificate. Polarization branch and full irrep identification remain pending.
- Analytic scorer fixture and parent-operation tests are software checks only;
  no material force constants have yet been generated by this stage.
- Submitted `28246022` from immutable `eb9caa1` archive; all19 COMPLETED/exit0.
  Complete SCF/matrix audit is separate from this scheduler status.
- Rotated T-pattern triplet resolves a missing coordinate axis: T=(0,0,0.831575)
  Angstrom, PO=(-0.972368,0,0). The two inversion candidates retain/reverse the
  major x component. No phonon or X2- irrep label is inferred from this geometry.
- Fixed z reflection and one species assignment register the published PO->M
  seed to actual PO+ (RMS0.007456/max0.012612 Angstrom). Assignment is applied
  once to all20 images; corrected endpoint optimization remains pending.
- HF ASE3.23 optimizer API was inspected: it calls converged(forces), not the
  new gradient_converged hook. EndpointBFGS now enforces the same force/stress
  gates via both entry points; no physical parameter changed.
- Final fresh audit: all19 input/DSIZE32/SCF/full force/stress contracts pass;
  T fixed-cell Gamma has no negative optical mode at either step. Lowest
  frequencies1.020298/0.953709THz; max two-step spread0.066590THz.
- Both inverted candidates have energy-9783.24967581196eV/cell, maximum
  force0.000925eV/Angstrom and stress0.130872kbar. They need no repeat relaxation.
- M seed force0.136625/stress2.520136 requires BFGS. Its old published geometry
  energy is not promoted as a relaxed endpoint. Newly prepared optimizer uses
  physical force0.03/stress2kbar/maxstep0.02/max80steps, with bytewise original
  electronic INPUT/KPT and complete pseudo/orbital hashes on every SCF call.

## E006 — Sparse guided segment static audit

- Hypothesis: the historical sub-meV guided-chain barrier could be an
  under-resolved discrete maximum. Reject treating it as a precise barrier
  before filling the large pattern-space gap between images01/02.
- New budget: three points only, at fractions .25/.5/.75 of the declared
  linear fractional/cell interpolation, no atom remapping or rotation alignment.
- Inputs: exact historical INPUT/KPT/pseudo/orbitals, raw endpoint/snapshot
  geometry agreement required; genuine32MPI single-node fork, max3 simultaneous.
- Measurement: true sampled energies versus the same T baseline. A high
  interpolation energy would flag this straight reconstruction as unresolved;
  it cannot establish the globally optimal MEP or prove the two paths differ
  by a specific irreducible mode. Further path optimization is a separate step.
- All3 SCFs COMPLETED/exit0 in28251303 with fresh complete audits.
- Result:35.795717/59.386224/34.045090meV/f.u. above the common T baseline.
- Decision: reject the historical sub-meV discrete maximum as a resolved
  physical barrier or accelerator benchmark. Preserve raw history; insert the
  three audited geometries and optimize a10-total/8-interior ordinary chain.
  First allocation is a10-step health pilot, no CI or physical-input change.
  Sampled straight-interpolation energies do not determine the optimized MEP.

## E007 — Independent frozen/atomic-response curvature test

- Apply the two real T Gamma matrices to the3 rotated geometric patterns,
  release30 orthogonal nontranslational atomic directions, keep the cell fixed.
- Observed two-step operator spread0.0683446eV/Angstrom² is an operational
  stability gate, not a rigorous numerical-error bound. Eliminated minima
  0.100839/0.088140 pass with narrow margins; condition numbers490/561.
- Before any independent DFT: predict Qx frozen curvature4.813147 and
  linearly responded curvature1.979301eV/Angstrom² from the0.01-Angstrom matrix.
  The independent0.02 matrix predicts4.811954/1.976804 and is a step check.
- New budget8 SCFs: frozen/responded Qx at +/-0.05 and +/-0.10 Angstrom.
  Response must not be renormalized: the retained projection is the Q axis.
- Prespecified assessment: paired energy and force-gradient curvatures versus
  the fixed prediction; <=10% relative curvature error and energy/force
  disagreement at both amplitudes. Orthogonal force residual is reported too.
  Failure is recorded, not cured by changing cutoff or refitting these holdouts.
- This local harmonic response line is not fully minimized conditional DFT,
  joint-cell reduction, a barrier model or a new Schur-complement theorem.
- All8 new SCFs COMPLETED/exit0 in28257780; fresh complete audits pass.
  Frozen energy curvatures4.854014/4.836618 and force curvatures4.820535/4.846827;
  responded energy2.022935/2.000460 and force1.983617/2.004147eV/Angstrom².
  All4 prespecified pairs pass; max prediction error2.2045 percent, no refit.
- Decision: retain independently verified local atomic harmonic prediction,
  with orthogonal residuals and finite-step spread. Do not extrapolate it into
  full conditional relaxation, a joint-cell channel or an activation barrier.

## E008: endpoint electronic polarization (preregistered)

- Verify PO+ and two ordered inversion-mapped PO- with native ABACUS f7cb1d3 LCAO Berry.
- Budget3 output-only SCFs +9 fixed-charge NSCFs. Start only PO+ before validating other endpoints.
- Keep100Ry/full10auDZP/SCF2x2x2 and strict original physical contract. Output addition must reproduce E/F/stress.
- On the same charge use R3 NSCF222/224/228. These energies never enter barrier comparisons.
- Before DFT: sampled gap>.1eV;224-to228 modulo difference<=.01C/m²;
  POminus/Pplus inversion residual<=.01C/m². Reject self-inverse ambiguity.
- Record native2eR/V modulus and physical eR/V separately. No spontaneous-P or path-branch claim.
- Reuse hf installed binary and ASE. Three SCFs and all nine retained NSCFs
  now pass fresh raw-output audits (28267393_0,28267783_1/2).
  At228 the modern-SI branch-valued R3 components are+/-0.715814325C/m²,
  inversion residuals0 at print precision;224-to228 differences0.000526896.
  No absolute spontaneous polarization or continuous path branch is claimed.
- Two wrapper audits initially failed after successful DFT: missing NSCF band
  output, then a unit-period check that ignored source-specific SI constants.
  Corrected against the actual f7cb1d3 source, reusing completed outputs.
  Actual study cost3 SCFs+10 NSCFs, including one repeated222.
- Switching seed inspection caught a restored wrapped final endpoint causing
  a false long last segment in unwrapped VCNEB coordinates. Preserve the
  continuous ordered MIC lift and regress every initial segment length.
  No switching DFT was submitted with the defective seed. This repair changes
  only integer endpoint lattice translations, not physical geometry or inputs.

## Existing evidence and non-results

- ABACUS ordinary-7/9: 39.187725/39.0523 meV/f.u.; historical CI-7: 32.29305.
- Guided-7: historical discrete0.963925meV/f.u., fmax=0.0866792767 passes0.10;
  E006 now shows the straight reconstruction is under-resolved (59.386meV/f.u.
  internal sample). It is not a reliable low-barrier/acceleration benchmark.
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

## E009: exact clamped-plane boundary before the strain matrix

- Hypothesis: a fixed substrate requires an affine deformation subspace;
  reusing the slab mask would release the wrong variables. Holding two ASE
  cell rows fixes the plane exactly through deltaF=v outer n.
- Action: one explicit tilt option, pre-projection chain validation, a matching
  ASE BFGS endpoint filter, and stress-work derivatives with stated units.
  All atomic freedoms retained. No DFT parameter change or DFT submission.
- Evidence:22 new tests;29 focused including prior slab/joint tests.
  Clean staged-tree archive b8666a3:748passed/2skipped,80.94s.
  Module6034e42d byte-identical on HF; ASE3.23.1b1 analytic EMT checks give
  maximum derivative error9.42282336e-10eV/A, zero substrate drift and9 BFGS
  steps for both normal-only and tilt-released boundaries.
- Decision: retain as the G2 mechanical interface, select tilt-released
  constraints for the finite study. HfO2 clamped endpoint/chain DFT and
  channel predictions are pending; this is not an acceleration or novelty
  claim. Do not require zero substrate reaction stress.
