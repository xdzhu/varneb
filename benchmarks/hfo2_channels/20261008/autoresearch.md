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

## E010: common-substrate seeds and physical endpoint convergence

- Hypothesis: a single substrate and a rotation-invariant open-traction gate
  prevent independently strained phases or reaction stress from invalidating
  the finite G2 comparison.
- Action: proper cyclic orientation (new x/y/z=old z/x/y, unchanged atom
  order), five registered structures at epsilon=0/+0.01, and ten hashed
  geometry starters. Reserve+0.005 holdout. Keep original 100Ry/full10auDZP
  INPUT/KPT/PP/orbital checksums, drop all free-cell calculator results.
- The T-zero substrate imposes+3.12094%/-3.60936% PO in-plane length changes;
  explicitly distinguish this from free-cell P=0. No material conclusion.
- Matching BFGS tests atomic force<0.03eV/A and open-traction norm<2kbar.
  Traction formula cross-checks stress work in oblique cells and rotates
  covariantly. Dangerous geometry is rejected before atom mutation/DFT.
- Evidence so far:44 focused tests pass. Full clean-tree and installed-HF-ASE
  verification recorded below after execution. This preparation costs0 DFT;
  G2 SCFs/matrix remain gated on G1 and a reviewed bounded endpoint canary.
- Decision: retain the explicit geometry/mechanics contract, not a scientific
  claim of favourable strain selectivity or a completed HfO2 endpoint.

- Clean staged tree b0f7432, archive4e242396a8ea393fd18c4558b7135923c82892c21b34a88ee6817887ee7eed4d:
  766passed/2skipped/77.18s. First export51d0963 correctly failed one source
  checksum gate because Windows core.autocrlf converted hash-bound bytes;
  targeted .gitattributes preserves original bytes, without loosening guards.
- Actual HF ASE3.23.1b1 checks all10 archived hashes and independently regenerates
  all10 seeds. Both EMT endpoint variants converge in9 BFGS steps/10 evaluations,
  with0 plane drift; atomic force0.0154315eV/A and open traction0.678580kbar.
  In-plane/full stress14.2747kbar demonstrates why the reaction is excluded.
  Record:clamped_endpoint_preflight_hf.json. Total new DFT calls remains0.

## G1 switching dispatch R13 (not a strain experiment)

Two registered/electronically verified candidates are now queued as28288045
and28288063. Each is9total/7active,32MPI,ordinary0.10,10steps/4h. They reuse
immutable sourcecc536a2/R12 and validatedR11 seeds, with original physical
inputs and cached endpoints. `afterany:28274895:28275259` waits for both
current allocations to finish, preserving at most2 active study chains.
This is scheduling, not a claim that the predecessor paths have converged.
See submission_handles_r13.json; no switching results or G2 DFT claimed yet.

## E011: continuous-reference projections and actual force attribution

- Question: can the same rotated-T triplet represent both pilot bottlenecks,
  and is the gap-chain rebound primarily caused by springs or cell forces?
- Method: freeze seven complete exact-SCF observations (66 existing image
  records), replay E/F/stress, retain the initial integer gauge along each
  path, and separate atomic/cell/perpendicular/spring force diagnostics.
  Preserve original Gamma masses, record complete F/Green strain, no DFT.
- Measurement: gap step6 peak triplet displacement fraction0.939053 versus
  PO--M step10 fraction0.528366. At the rebounding gap image1, atom/cell
  residual0.344031/0.030582 and spring0.023657eV/A; true perpendicular
  force0.346935. PO--M step10 remains0.293325 despite Slurm COMPLETED.
- Decision: retain a descriptive three-pattern tracker, reject a complete
  common triplet landscape claim; continue ordinary relaxation, not parameter
  retuning or a rebound-triggered stop. Neither provisional barrier nor full
  T-Gamma reconstruction certifies a TS or predictive energy model.
- Evidence: chain_observations/analysis.json and the dedicated report;
  29 focused tests pass. Clean-tree/HF results appended after execution.

- Final clean tree001a7229/archivecc915dc6:791passed/2skipped/69.10s.
  The first archive lacked ignored.traj evidence (3 failures); explicit
  versioning fixes delivery without changing data. First HF comparison
  caught platform-dependent pathlib ordering; explicit casefold sorting
  followed by a second full clean regression and actual HF replay passes.
- Actual HF ASE3.23.1b1:66 records, maximum fmax difference5.55e-16eV/A,
  Gamma Q difference1.78e-15sqrt(amu)A, Green-strain difference0; DFT0.
  Record:chain_observations/hf_replay_check.json. Numerical reconstruction
  accuracy is not a material-convergence, harmonic-validity or novelty claim.

## E012: terminal-run exact-cache ordinary continuation

- Question: can the latest complete health-segment chain be continued without
  changing physical inputs, re-evaluating endpoints, or forking a live run?
- Guard: authoritative terminal sacct allocation; latest complete snapshot;
  freshly matched raw SCF inputs/geometry/logs/E/F/stress; ordinary replay;
  refuse existing output and already-converged residual. Fresh FIRE state is
  explicitly distinct from restored optimizer state or acceleration.
- Measured: both R12 allocations COMPLETED/0:0 at step10, but gap residual
  0.284995243 and PO--M 0.293325305eV/A remain above0.10. Actual immutable
  production-factory preflight on HF reproduces both using19 exact caches;
  all DFT calls forbidden in memory, no new SCF. Separate CLI checks pass.
- Decision: retain unchanged physics, enlarge execution segment to80steps/24h;
  submit once as28300425(gap) and28298794(PO--M), verifiedPending/Dependency.
  Wait for both R13 health allocations viaafterany, preserving2activechains.
  Do not cancel on a rebound or interpret allocation completion as convergence.
- Evidence: ordinary_continuation andsubmission_handles_r14.json. Source/template
  preparation archive100c0de9 passes804tests/2skipped/71.02s; with the two new
  delivered material-cache replay cases,15focused tests pass. Final clean
  deliverable regression follows. Two actual R13 first-imageSCFs take98.73/98.43s
  with32MPI, distinct from the zero-DFT preparation measurement.

- Final clean delivery treecdffc642, archivededec562:
  806passed/2skipped/77.27s, including both published material-cache replays.
  Helper and batch SHA agree bytewise between the clean delivery and actual
  HF preparation. Only documentation/receipts follow; no code/seed changes.
  Both switching candidates complete step0 at19:39 and continue normally;
  their initial fmax2.214746/1.300483 does not establish the optimized barriers.

## E013: explicit candidate coverage and common-initial-state comparison

- Question: can four observed channels be compared without mixing directional
  baselines, omitting a decay candidate or confusing relative selectivity with
  the stronger absolute noninferiority claim?
- Method: declare both flip variants plusPO--T/PO--M; verify ordered periodic
  commonPO+ geometry/raw energy/physical-input hashes/ensemble. Retain original
  force-replay metric when reversing the T--PO thermodynamic view. Numeric
  barrier/sampling-error evidence is separate from ordinary residual0.10.
- Measured:37existingSCF images, DFT0. PO-based provisional maxima116.277,
  95.199,145.779,435.169meV/fu; all unconverged and error/sampling bounds
  missing, soH1 conclusions are blocked. These are not final MEP rankings.
- Representation: peak T-triplet fractions94.92%,52.84%,54.29%,approximately0;
  both flip step1 residuals are cell-dominated and their central spring forces
  negligible. Retain joint cell/atomic freedoms and test later bottleneck
  bases, not a forced complete triplet plane or inferred modal energy partition.
- Software: required-role minima are withheld on missing coverage; deterministic
  minimum/difference intervals allow uncertainty in which channel is lowest.
  Same-family conditions only; explicit predeclared noninferiority margin,
  no Gaussian/independence assumptions or lifetime prediction.
- Verification:54focusedtests pass. Initialclean643a2af2/archive05c75b63
  gives844passed/2skipped/76.61s. An additionalNumPy-count JSON test then
  exposesnumpy.bool_ serialization; normalize the API input count, not the
  scientific data. Replay changes only the module SHA, no numeric result.
  Final clean regression and actualHF replay are recorded next.
- Unit audit: standard-library static checker finds0issues in2modules; an
  independentSI pressure-volume regression uses library electron charge.
  NoPint/environment/constant change; measured numerical bounds are not
  Gaussian standard errors or95%confidence intervals.

- Finalcleanb777ad49/archive70acd20c:845passed/2skipped/78.69s.
  ActualHF ASE3.23.1b1 replays4channels/37images with force/barrier differences0,
  patternQ difference3.68e-16A, GammaQ difference1.78e-15sqrt(amu)A,
  Green-strain difference0; allsource/module hashes agree. DFT0.
- At20:32, flips28288045/63 areconfirmedRunning atstep3/4, residuals1.802978/
  1.026364; bothR14 continuations remainPending/Dependency. Recent11.6--16.6min
  perstep projects the10step health segments to21:40--22:30, not convergence.
  No newallocation orphysical change in thisanalysis milestone.

## E014: close-prior-art audit and prospective prediction controls

- Question: does an apparent mode/NEB prediction increment actually go beyond
  the closest accepted work, and can the future holdout test distinguish it
  from interpolation or well shifts?
- Evidence: Qi2025 accepted manuscript scientific body/AppendicesA–E,
  visualFig4 per4fu-cell units and the additional oxygen-crossing category;
  Zhou2022 mainJATS/Methods and allSI, publisher-declared MD5-matched PDF,
  visualS7 and the500random-start restricted-subspace search.
- Decision: reject first-mode-landscape/NEB-validation or exhaustive-network
  novelty claims. Preserve current production matrix, disclose unexamined
  crossing mechanism; no speculative new variant added.
- Protocol: freeze+0.005 in-range holdout before materiallabels; retain linear
  barrier interpolation and both endpoint-response nulls plus frozen/atomic/
  joint/branch controls. H1strictprimary decay noninferioritymargin0; numerical
  ambiguity/abstention cannot manufacture a favorable result. Twoheldout edges
  do not certify allcandidate minima. Amendments retain first predictions.
- Limits: otherSI and Lee doubleDOI relation remain unresolved; no literature
  absence proof, newDFT, scientific H1/H2 verdict or finalJCTCreadiness.
- Live20:51: flips28288045/63Running,step5/6,fmax1.552917/0.890479;
  R14continuations28298794/28300425PendingDependency. Existing sourceunchanged.
- Follow-up21:05: Behara2022 accepted manuscript methods and polymorph/variant
  path sections, Fig9 visual audit. Already contains two-shuffle contours,
  strain-path mapping and a Pbcn-to-T intermediate change on fixing the Ocell.
  Retain restricted candidate coverage; reject first-mechanism-change or
  first-mode/path-overlay claims. SI remains unread, no matrix expansion.
- Validation:16localMarkdownlinks,39JSONLrecords parse before the validation
  receipt;43focusedchannel/material/atomic-holdout tests pass in8.87s. No new
  full-suite claim; source and scientific data unchanged. At21:15:55, flips
  remainRunning atstep7/8,residual1.314195/0.759024, bothR14PendingDependency.
  The10-step health segment ETA updates to21:40--22:00, not convergence.

## E015: joint probes use the endpoint/path mechanical subspace

- Question: is the free-cell symmetric-strain chart consistent with theG2
  tilt-released fixed substrate? It is not: the exact three rank-one allowed
  directions are generally nonsymmetric. Sampling six strains changes the
  experiment instead of measuring its open atom--cell coupling.
- Action: separate ActiveJointCurvatureCoordinates, shared exact boundary
  basis and raw atom/cell work; fixed-cell0, normal1, tilt3cell freedoms,
  translation removal and pre-evaluation geometry rejection. No DFT call.
- First prototype changed the legacy module and correctly failed one historical
  GaNsource checksum (66pass/1fail). Do not rewrite the old manifest. Preserve
  byte-identical legacy6c83d2ec and use an independent subclass; default
  geometry/gradient equality is additionally regressed.
- Evidence:70focusedtests/3.67s;10unrelaxedregisteredHfO2starters,50probe
  geometries,8CuEMTworkcases/232evaluations. Max gradient error1.0997911e-9eV/A,
  below prespecified3e-5. Clean delivery/HF verification follow, not yet claimed.
- Decision: retain consistentG3preparation; materialjointHessian/TS/barrier
  prediction remain pending. Current electronic inputs and active jobs unchanged.
- At21:47:26 both switch allocationsCOMPLETED0:0 atstep10, residual0.980657/
  0.700279; bothR14continuationsRunning since21:47:06. A health segment ending
  above0.10 is not scientific convergence. Exact-cache continuation is next.
- Delivery evidence: clean tree fed6a95f / archive dba8a5fe, 864 passed,
  2 skipped in 142.06s. Local/HF identical archive modules match all five
  hashes; ASE3.28.0/3.23.1b1, max work errors1.100e-9/1.840e-9eV/A,
  HF substrate drift8.0321e-20A. Windows worktree EOL-normalization differences
  are retained and disclosed, not confused with a changed historical module.
  Source code unchanged after testing; all checks made zero DFT calls.

## E016: terminal switching evidence and exact-cache continuation

- Both first10-step flip segments ended normally above0.10. Audit18 complete
  rawSCFs/geometry/input/fullE/F/stress; ordinary replay0.980656799/0.700278856.
  No newSCF, source modification, relifting, remapping or physics changes.
- Execute oldR12production factory with its DFT entry disabled: all18 cached
  evaluations pass and reproduce the seed residuals. Actual old production
  CLI geometry/periodic-lift preflight passes, not a replacement factory mock.
- Submit once each28319570/28319571, afterany28298794/28300425, 32MPI/OMP1,
  80steps/24h,9total/7active. Actualscontrol verifies PendingDependency;
  the maximum2active-chain cap and registeredcandidate matrix are unchanged.
- New frozen four-channel comparison retains the earlier step1 data. The
  preserving flip/PO--M discrete maxima differ by2.002meV/fu; all unconverged,
  sampling/error bounds missing, selectivityunset. No finalranking/H1 claim.
- Both flip peaks remain cell-dominated, but the reversing chain's largest
  residual moved to image2atomic motion. Peak-triplet capture57.39%/near0
  is descriptive, not a phonon/Hessian/energetic contribution orTScertificate.
- Regression caught the growing root (9observations) entering a historical
  7-observation baseline. Pin its original7named inputs, preserve old data,
  and add3independent newevidence tests.72focusedtests pass in33.58s; final
  clean delivered-tree regression is next. Verifier/raw bytes pinned by Git.
- Decision: retain same-physics ordinary continuations and use actual moving
  bottlenecks for subsequent joint-basis selection. No acceleration claimed.
- Final delivery:867passed/2skipped/84.72s from clean f575cd8f /2c00c537.
  All four ignoredtraj files are explicitly present; executed verifier byte
  hash agrees onHF/local/clean archive. Historical data retained. Only receipts
  and documentation follow the test; no extraDFT or sourcecode changes.

## E017: selected-direction joint probes retain transverse response

- Question: can the existing square-gradient assembler safely consume k<d
  selected directions? No. Discarding the full-chart gradient loses H*B
  coupling and cannot measure the missing complement's self-curvature.
- Action: calculator-free paired probe generator and a separate directional
  assembler with full H*B, raw/symmetric projected matrices and transverse
  action. No implicit normalization, atom remapping, boundary modification,
  calculator attachment or unsampled-block elimination. Physical enthalpy
  gradients are required; NEB/spring forces are not curvature measurements.
- Evidence:72focusedtests/4.32s;720inertHfgeometryprobes and180CuEMT
  evaluations with two amplitudes, normal/tilt release and global rotations.
  Local energy/gradient curvature difference9.17149e-5eV/A2, reciprocity
  defect1.04768e-5, two-step operator spread1.59885e-4eV/A2. These are
  implementation bounds, not material DFT uncertainty or TS certification.
- A misplaced assertion in an expanded test caused a NameError and was fixed
  in the test; the physical implementation cases did not fail. HF has no
  pytest; do not install it, use the standalone same-archive verifier.
- At22:47:58 twoR14jobsRunning,PO--Mstep5fmax0.230286 andT/POstep4
  0.136163; twoR17flipsPendingDependency. Single rebounds are not stop rules.
- Decision: retain backend-independent G3 preparation. No new DFT, production
  parameter change, extra active chain, material Hessian or prediction claim.
  Clean full-suite and same-archive HF delivery checks follow before publishing.
- Delivery: clean treea6dad26a/archive45dd4d67,895passed/2skipped/91.30s.
  Same archive verifies onHF ASE3.23.1b1 with six matching module byte hashes,
  720geometryprobes/180EMT; maxenergy-gradient difference9.17198e-5eV/A2,
  reciprocity1.04768e-5, two-step spread1.59885e-4eV/A2, substrate
  drift4.8736e-20A. Legacy module unchanged. No code changed after tests,
  only receipts/documentation. Actual22:58logs stillstep5/4above0.10.

## E018: extract real mixed atom--cell response from completed G0 probes

- Question: do the existing8 full-gradient G0 probes contain measurable mixed
  response, and is their two-direction slice sufficient for reduction?
- Action: reuse the original free-cell chart/atomic secant/biaxial strain,
  not G2's substrate chart. Assemble both42-coordinate Hessian actions with
  the new API; retain original numeric reports, pair tags and raw matrix.
- Preliminary archive evidence: mixedentries-3.31840/-3.31898eV/A2 at0.01A,
  -3.31790/-3.31942at0.02A. Positive projected eigenvalues3.94316/19.79843
  coexist with transverse norms5.86712/7.09425. Centergradient0.118695eV/A
  is nonzero; no full stability index, Schur elimination or TS is claimed.
- Decision: selected directions expose coupling but are not a closed conditional
  plane. This is exploratory reuse of measured DFT, not a new prediction test.
  Fresh HF raw-source/geometry/input/SCF/gradient audit and clean delivery follow.
  Initial46focusedtests pass/2.60s; metadata naming then clarified and rechecked.
- No extra DFT, production source overwrite, parameter retuning or G2 submission.
  At23:02R14Runningstep6/5,residual0.217313/0.101538; flipsPendingDependency.
- Actual delivery:c06169bc/archive250857fc,903passed/2skipped/88.61s.
  HF rechecks center+8rawSCFs, allgeometry/contract/DSIZE32/convergence/gradient
  and original log hashes. Full actions agree exactly with local archive;
  two script/API hashes match. Raw/normalized manifest hashes differ only
  through Git text normalization, content equals; DFT input hashes remain strict.

## E019: first ordinary G1 edge passes the residual gate

- Actual28300425COMPLETED0:0 at23:13:32,force_thresholdstep6, notstep_limit.
  Ten finalimageSCFs/geometry/input/E/F/stress audited; numeric-only replay
  0.05988161569025288HF/0.059881615690252875local. Original endpoints retained.
- DiscreteT→PO33.888971/PO→T115.210164meV/fu,deltaE-81.321193;peakimage3,
  largestresidualimage2atomic. No fullTS/sampling/error-bound orG1completeclaim.
- Actualcontinuation48completeinteriorSCFs,0endpointcalls,5186s/32ranks=
  46.09778allocationcorehours; earlierhealth/seedcostexcluded,notacceleration.
- Existingflip28319570capacitydependencycleared afterterminalpredecessor,
  noresubmissionorinputchange. At23:20Runningnode11;PO--MRunningnode26,
  flip28319571stillwaitsPO--M. Initialjournalretained,overrideinnewreceipt.
  This preserves atmost2activechains. ProductionarchiveSHA df4c12ee unchanged.
- 29focusedmaterial/source/replaytests pass/14.87s. Finalclean regression of
  the combinednewresultdataset andraw-hash packaging follows beforepush.
- Finalcombinedclean tree42994505/archive1ec261f4:905passed/2skipped/89.01s.
  Trajectory included andrawdatahashes preserved; analysis/API match actualHF
  validation bytes,legacyjointcurvature unchanged. Post-test changes only
  receipts/docs. At23:32:46PO--Mstep8fmax0.201984, preservingflipstep1
  0.900383, bothRunning; reversingflipPendingDependency. NoG1completionclaim.

## E020: descriptive coverage of the first converged ordinary band

- Question: does a compact reference-pattern description at the maximum also
  span the full path, and are the lowest T optical modes the dominant weights?
- Re-audit ten frozen SCFs and replay all archived descriptors before plotting.
  Triplet coverage95.54% at image3 versus31.43% at PO; T-reference5.22THz
  doublet/9.69THz singlet carry62.16%/28.53% of peak optical squared norm.
  The two lowest optical doublets sum3.54%. Bases/metrics differ explicitly.
- Keep complete degenerate groups, undefined zero-reference fractions and
  absent fixed-endpoint NEB residuals. No energypartition, localTSphonon,
  independentprediction, frozenconditionalplane or accelerationclaim.
- Six aligned panels, normal(a), four spines, framedtranslucentlegends.
  PNG reviewed;61SVGtext elements and446PDFcharacters/3embeddedTrueTypefonts.
  Script/sourceCSV/vector/raster/QA/caption provided, with no interpolation.
- Initial16focusedtests pass10.10s. Bundle-hash check and finalclean suite
  follow beforepush. NoDFT, calculatorchange, jobmutation or userfileoverwrite.
- Actualfinalclean treeb2663974/archiveb3c411ac:915passed/2skipped/100.09s.
  Samearchive regeneratedCSV andPNG are byteidentical. Vector timestamps/IDs
  are not claimed byte-deterministic. Post-test changes only receipts/docs and
  generated-SVG-specific whitespace metadata; no code/data/figure rewrite.
  Actual23:56:49PO--Mstep10fmax0.189975, preservingflipstep3fmax0.747910,
  bothRunning; reversingflip stillDependency. Twoactivechains unchanged.

## E021: stronger fixed-parent control before G2 labels

- Read the distinct Qi--Rabe PRL135,046101/arXiv2412.16792v2 scientific body
  and original20-page ancillary SI. Download signatures/hashes retained;
  Poppler pages3/5 check the literal coordinate issue and phonon folding.
- Existing G1 T-pattern coverage is exploratory evidence, not proof that all
  fixed references fail. Add a mapped Cmma control in prospective protocolv2,
  separating reference choice from frozen/atomic/joint/branch response.
  Keep v1 unchanged except its explicit dated pointer; no holdout labels exist.
- No copied PDFs/figures, new DFT, electronic change, job mutation, expanded
  channel matrix or claim that the stronger control is implemented/validated.
  The next preparation is inert ordered-cell/basis mapping, not a new phonon
  production matrix. Existing SI/DOI access gaps remain disclosed.
- Actual00:14:40CST: PO--Mstep11fmax0.186321, preservingflipstep4fmax0.677541,
  bothRunning; reversingflipPendingDependency. At most2activechains retained.
- Decision: retain the stronger test and reject novelty from a weak-reference
  comparison alone. Documentation/link/source-hash validation and existing
  material-figure regression precede ordinary milestone push.
- Actual delivery checks: both original PDFs match bytes/pages/SHA256;
  32local document links pass, all53ledger JSONL records parse, existing
  material-figure regression10passed/3.19s. This is a documentation-only
  milestone, not a fresh full-suite/material-prediction pass. User5tracked
  files and unrelated untracked artifacts remain excluded.

## E022: inert published-Cmma mapping and strict QE Gamma import

- Question: can the stronger fixed parent be defined from the author's data
  without silently correcting the SI, adopting LDA curvatures, or changing a
  production cell? Pin author commit a438e4ec and five original Git blobs.
- Implement a bounded first-Gamma flvec reader with explicit ordered masses,
  complex phases and measured normalization/Gram defects. Historical official
  QE6.3 source and current docs agree that flvec is displacement, not an
  orthonormal dynamical-matrix basis. No QR/ASR or imaginary-part truncation.
- Real-data checks initially rejected the repeated QE block separator and an
  incorrect proper-frame candidate. Support the actual separator after a
  complete3N block; explicitly compose the author's y/z exchange with Cmma
  inversion/translation. No tolerances or electronic parameters were retuned.
- Author QE->POSCAR site error8.29164e-7A; independent literal TableS2
  primitive->conventional reconstruction1.45467e-5A. The TableS1 x=.05000
  discrepancy is retained. Six->twelve only reconstructs an existing reference,
  not a production-cell expansion or new Gamma calculation.
- FirstGamma36modes: restoredGram4.00853e-6; three rigid translations overlap
  .999992296. Four resolved negative reference modes separate into two odd and
  two even centering characters; do not infer irreps or reclassify other q.
- Production-path gauge and a real/exact comparison subspace remain unprepared;
  the strong material baseline, PBEcurvatures, TS or prediction has not passed.
  Raw author files remain outsideGit, with no explicit redistribution license
  seen at repository root. Underlying flfrc/DFPT generation is not audited.
-34focusedtests pass. Clean stagedtree5d25f8d6/archivea4cab9a5:949passed,
  2skipped/108.99s; samearchive reference report is byteidentical. A separate
  working-tree diagnostic969passed/1skipped is not the authoritative delivery
  test because it includes unrelated untracked files. No code changed after
  clean tests; subsequent additions are receipts/docs only.
- Actual00:54:11CST:28298794Runningstep14fmax.176483,node26;
  28319570Runningstep7fmax.497080,node11;28319571PendingDependency.
  Original productionarchiveSHA df4c12ee unchanged; no job mutation/newDFT.
- Decision: retain audited parent-source preparation for a fair strong control,
  not a claim that the adaptive model beats it. Samearchive HF inert
  compatibility verification follows; G1/G2/G3 material gates remain open.
- ActualHFsamearchive inert audit succeeds with existing icu/ASE3.23.1b1;
  archive/twoAPIbytehashes exact, five raw source hashes exact. All133float
  fields agree within1e-12(max1.13687e-13); only library-version metadata
  and three capturedspglibdeprecation lists differ, retained in both reports.
  Nopytestinstalled/runonHF,DFTcallorproductionoverwrite. Successful source
  compatibility is not material prediction or TS validation.
- Packaging check rejected anLF/CRLF-only sourcehash drift; restored exactly
  the tested sourcebytes and added explicit raw-byte attributes. No semantic
  code/data edit after clean testing. Final delivery rechecks all six code/
  reporthashes, JSONL parsing and local document links beforeordinarypush.

## E023: register the strong reference in actual frozen production paths

- Question: can the published four-direction reference be compared without
  arbitrary axes, hidden atom relabeling, origin/metric changes or silent ASR?
  The naiveCmma->T registration has four equal-distance oxygen ambiguities;
  its arbitrary-axis low-coverage trial is not accepted as strong-control evidence.
- Verify an additional pinned authorT file via publicGitHubmetadata/blob/SHA.
  Enumerate192 proper signed-frame/quarter-origin choices, retain16 T-compatible
  frames, then all4 commonPO+ nearest-registration ties. Coverage/energy/forces
  do not enter selection. Global second-best assignment gaps1.47265A.
- Register reference rows once; leave every production image/lift untouched.
  Transport directions into the common originalT fractional/mass chart.
  Explicit real/imagSVD requiresrank4 andrecords complexprojection error2.08355e-6.
  Source translationadmixture norm.0036854 is retained, notASR-repaired.
- Ten convergedTPOrecords replayfmax.059881615690252875. CommonT-origin Cmma4
  peakcoverage.271810-.271901/PO.522401-.522402;Ttriplet3 peak.896905/PO.304151.
  Differentorigin/metric fractions from E020 are not mixed. Cmmaaffine-plane
  endpointresiduals6.43182/5.60197sqrt(amu)A block a global four-mode plane claim.
- Same unchanged registration analyzes9 earlier snapshots/84frozenimage records,
  including repeated endpoints, not84 new/uniqueSCFs. Reversingflip step10
  peakCmma~.4478 versusTtriplet~.2921, unlike formation. These are unconverged
  observations, not final competing barriers or prospective predictions.
- Synthetic phase/unitary/projector, explicitrank/error, full frequencydoublet,
  boundedframeenumeration/ambiguity and input-preservation regressions pass.
  Initialnoise test incorrectly inspected the last rather than the fourth
  singular value; correct the assertion, not the numerical policy.
-68focusedtests pass/1.20s. Clean treedc803124/archivec6658bc6:973passed,
  2skipped/94.46s. Samearchive convergedreport byteidentical. Subsequent
  changes are reports/docs only, no implementation/test edit after clean tests.
- Actual01:33:22CST:PO--M28298794Runningstep17fmax.166350,node26;
  preservingflip28319570Runningstep11fmax.319587,node11;28319571PendingDependency.
  No newDFT, electronicretuning, productionoverwrite, jobmutation or earlyG2.
- Decision: retain all registered strong-reference choices and channel-specific
  descriptive limits. G1complete/G2response/holdout prediction advantage remain
  unproven. Samearchive HF converged/historical replays and delivery checks follow.
- ActualHFsamearchive replays both reports successfully. Module/archive and
  referencesourcehashes exact;1451+3537float fields within1e-12,maximum1.13687e-13.
  Only originalsourceaudit library versions and capturedwarning lists differ,
  retained in both reports. NoDFT, installation, nodeallocation, production
  sourceoverwrite or source/test edit after clean regression.

## E024: nonstationary same-center response, not a merged-reference Hessian

- Hypothesis: the existing mixed G0 block can quantify a restricted local
  atomic-release correction without pretending its center is stationaryT.
  Reject merging the historical image02 block with the differentT Gamma matrix.
- Implement backend-independent nonstationary quadratic condensation with
  affine response, gradient and energy shift; omitted variables stay clamped.
  Stable-block gate, covariance, scalar/vector shapes and immutable arrays
  are regression-tested. Schur condensation is not claimed as new mathematics.
- Actual8existing G0SCFs: one-direction release softens strain curvature
  12.36%/12.43%; offset.021648A, energy lowering.00109457eV/12atomcell.
 35otheratomic+5cell variables fixed, no unmeasured stability assertion.
- Short-step gradients+center predict four previously seen long-axis energy
  changes with maximum.0000416931eV/cell; projected/full-gradient-action
  errors.0044273/.0050403eV/A. Retrospective, not blind/independentbarrier evidence.
- The entire requested response segment is outside the axis-probe convex
  hull. This is a real coverage failure: reject a releasedDFTbranch/TS claim,
  not the measured mixed response. No extraG0DFT is added to beautify a contour.
- Explicit qcell=sqrt(2)*L*epsilon units/Jacobian and offline static audit
  pass; no installedPint, fake statistical sigma or environment retuning.
-45focused pass/.80s; cleantree9038c05b/archive0c2c662b:997pass2skip/100.64s.
  Initialworktree vs archive differs only in three legacyLF/CRLF metadata
  hashes. Both reports retained; all other fields exactly identical, noDFT
  inputhash relaxation or source edit after clean regression.
- Actual01:59CST twoexistingchainsRunning:PO--Mstep19.159754 and
  preservingflipstep13.264869; reversingflipPendingDependency. Nojobmutation.
- Samearchive originalHFcenter+8SCF rawreplay succeeds; final numeric/source
  comparison and delivery receipt are archived before ordinary push.
- Decision: retain restricted response/API and explicit failure coverage;
  G1/G2/G3 materialprediction and JCTCcompletion remain unproven.

## E025: phase identity is not the zero of a truncated pattern coordinate

- Export complete PO--Mstep20 and preservingflipstep14 directly from original
  live HF SCFs:18 newly frozen image records, including cached endpoints, no
  new SCF. Combine19 reused records to make the dated37-image commonPO+ set.
  Original E/F/stress/INPUT/KPT/orbital/log hashes and32MPI provenance pass.
- Same PO+ energy reference: sampled maxima115.210/82.136/65.279/410.744meV/f.u.
  OnlyTPOordinary force criterion passed. No final ranking, sampling error,
  fullTS or switchability-window claim; discrete low profiles are not a
  continuous energy bound. G1 remains incomplete, no newG2/G0DFT.
- Independent immutable pymatgen structure check reports three symprec values
  with warnings, no standardisation/wrapping/reordering/oxidation guesses.
  PreservingcentrePbcn; reversingcentrePbca despite nearzeroTtriplet amplitude.
  OriginalT isP4_2/nmc. MpeakP1/P2_1 tolerance dependence remains explicit.
  These are structure labels, not stable-phase or TS-index certificates.
- Check Behara--Van der Ven's primary accepted Sec.III.D: stablePbcn switching
  intermediate and constraint-inducedT are prior art. Our unconvergedPBEpeak
  is not that SCAN intermediate, andPbcn occurrence is not our novelty claim.
- Four-panel energy/residual/shuffle/strain plot follows the user style,
  straight discrete connections only. Reversing410meV profile is fully in
  CSV/table, outside the disclosed low-energy(a)panel. Correct first-pass
  annotation overlap before final tests/actual PNG view.54visible texts stay
  in canvas;48editableSVG texts. Superseded own outputs moved recoverably
  to E:/TEMP, no user files removed.
-24focusedpass/29.71s. Clean tree8344357d/archive53aa5b90:1012pass,2skip,
  313warnings/113.95s. All clean report fields match; regeneratedCSV/PNG/TIFF
  bytes identical. Source/tests unchanged after clean regression.
- SamearchiveHF0DFTreplay passes all source/phase/force/energy fields:
  697material+2635network floats, maxdifferences2.84217e-14/1.42109e-14.
  Actual package versions and symmetry-warning lists differ and are retained
  on both sides; no rawDFT-input policy relaxed, installation or job mutation.
- Actual02:45:20CST:PO--M28298794Runningstep22.150184,node26;
  preserving28319570Runningstep17.199587,node11;reversing28319571PendingDependency.
- Decision: retain material-specific mechanism/evidence figure, not more
  generic API or an inflated claim. Finish G1, then the registeredG2/G3
  response and strong-control/independent-prediction tests. JCTCstillopen.
  Full receipt:network_update_20261009/validation_delivery.json.

## E026: registered null controls can execute, but G1 is not G2 training

- Preparation question: can the v1 B0/B1 definitions and training-only edge
  selection be executed without silently using target path labels, changing
  mechanical families, clipping impossible barriers or fabricating error bounds?
  The formulas and selection rule were preregistered before this implementation;
  they are conventional controls, not new theoretical results.
- Implement direct interpolation and both endpoint-response null models in a
  calculator-free module. Retain possible lowest channels under training bounds,
  lexical representative and upper-condition changes. Propagated training
  intervals explicitly do not bound unseen model error. Missing endpoint inputs
  stay unavailable; all supplied endpoint feature calls count in visible cost.
- Fresh HfO2 entry point records code/data/protocol hashes and caller unread-label
  attestation, not independent blinding proof. Original100Ry/full10au contract,
  tilt-open common substrate,0/+1% and unseen+0.5%, ordinary0.10 are fixed.
  It freezes B0/B1 only; B2--B5 forecasts and full holdout gate remain incomplete.
- Actual37-imageG1report is rejected by numeric and registered gates before
  output creation. Source SHA7bd2dbf7 retained. No material forecast is invented.
  Correcting G1 into fake strain/convergence labels is not used as a fixture.
- First clean run1044pass2skip/117.13s, then prepublication review identifies
  that rawE-only B1 features must not silently omitP*deltaV at finite pressure.
  Restrict that schema toP=0, retainP in output, and add a rejection regression;
  no HfO2 physical parameter changes. B0 still uses audited trainingH barriers.
-64focusedpass/1.68s. Finalcleantreee02e1417/archive7cd0277e:1045pass2skip,
  313warnings/114.03s; all three source/test bytes match the tested archive.
  Post-final-test changes are docs/receipts only.33new cases are mathematical workflow
  regressions, not independent HfO2 prediction successes.
- Actual samearchiveHFCLI rejects realG1 with the expected condition error and
  creates no output. Exact archived synthetic fixture functions independently
  reproduce B0(.275,.475),B1(.18,.19)eV/f.u. within1e-12 with existing libraries;
  source/archive hashes exact; finitepressureE-only misuse is rejected as well.
  These numbers are synthetic, not DFT labels.
- Actual03:09:48CST:PO--M28298794Runningstep23.146978,node26;
  preserving28319570Runningstep19.180742,node11;reversing28319571PendingDependency.
  No job mutation/newDFT/install/production overwrite or earlyG2/holdout.
- Decision: keep necessary executable controls and explicit limitations. Finish
  existing G1, then material boundary-response and independent prediction work;
  neither library tests nor hashes establish JCTC readiness. Receipt:
  prediction_controls/validation_delivery.json.
- Before delivery03:28:55CST recheck:PO--MRunningstep25.140413,
  preservingRunningstep21.162964,reversingPendingDependency. Same two active
  jobs, no resubmission or parameter/source change; both continue improving.

## Capacity-only scheduling update, 2026-10-09 03:40 CST

- Reversing continuation28319571 has an independent audited nine-record seed;
  its existing dependency waits for capacity, not PO--M physical outputs.
  Change only that queued job from `afterany:28298794` to
  `afterany:28298794?afterany:28319570`, allowing either active allocation to
  free the slot. This is not a new research experiment or faster-NEB claim.
- The first controller read was stale. Briefly hold only the pending job to
  verify the OR dependency, then release it at03:40:17CST. At03:40:53 it is
  PENDING Dependency with both parents unfulfilled; both parents are RUNNING,
  and the reverse continuation workdir has not been created. The hold is not
  left in place. At most two study chains remain active.
- Production archive/batch hashes are unchanged.32MPI,100Ry/full10auDZP,
  ordinary0.10, noCI, inputs and running sources remain unchanged; no new
  submission, DFT call caused by the adjustment, or earlyG2/holdout work.
  A failed parent may free capacity but still needs separate scientific audit.
  No actual wait-time saving has yet been measured. Full journal:
  `scheduling_update_0337CST_20261009.json`.

## E027 - Controlled biaxial input versus internal release

- Gap: a fixed-substrate chart samples internal freedom, but cannot measure
  mixed response to the imposed substrate strain. Add a tagged external input,
  not a new freely relaxed degree of freedom or a modified Hamiltonian.
- The new chart shares reference atomic coordinates and exact open-cell basis.
  Physical current-cell work supplies the external gradient, retaining substrate
  reactions. Its release basis excludes the control and translations. No
  augmented-control eigenvalue is counted as a fixed-epsilon TS instability.
-130focusedpass/6.42s;28new cases. Clean baseline407076a plus four unchanged
  code/test files:1073pass2skip313warnings/107.29s. User changes are excluded.
- Same tested source bundle c7d7d9a3 and checker a381619a verified on HF.
  There is no remote pytest; keep that failed invocation and use existing
  NumPy/ASE directly, without installing. Eight gradient configurations and
  two mixed-step checks pass,256EMT evaluations and10uncomputed training-seed
  geometries per host. Full-gradient error <=1.85e-9eV/A; it is an implementation
  residual, not a measured HfO2 error. The full suite ran locally only.
- Original production archive df4c12ee is unchanged. No ABACUS calls, new
  Slurm jobs, job mutations, calculator retuning, source overwrite or G2/holdout.
  The preserving chain's latest rebound is observed, not a stop trigger.
- At04:43:04CST squeue/sacct confirm M28298794Runningstep30/.123663,
  preserving28319570Runningstep27/.140772, reverse28319571PendingDependency.
  These log values are not promoted into a new frozen material result bundle.
- Decision: keep this necessary response-coordinate preparation. It does not
  complete B2--B5, material conditional branches or independent forecasts.
  Receipt: `biaxial_control/validation_delivery.json`.

## E028 - Residual-block crossover in complete production frames

- Freeze preserving steps25/27 and PO--Mstep30 in a new diagnostic namespace;
 27 numeric image records include6 reused endpoints. No new SCF, Slurm mutation,
 restart, parameter change or production-source overwrite.
- Exact raw E/F/stress, contract/geometry/log hashes and optimizer replay pass.
 Same unchanged analysis on HF/local matches1923float fields within2.85e-14
 and365nonfloat fields exactly.30existing focused tests pass/29.39s.
- Preserving sampled peak54.238773→51.683480meV/f.u.; atomic NEB residual
  .134176879→.140771896, cell block.124173502→.103344822eV/A. The dominant
  component is now atomic at image4, not the earlierstep14cell. Global spring
  contribution<=.000609eV/A cannot alone explain the residual magnitude.
- PO--Mstep30is atomic dominated atimage3/.123663133eV/A; its sampled peak's
 physical tangent-.365998838eV/A remains nonstationary. Ordinary NEB convergence
 must not be substituted for a full-variable stationary saddle certificate.
- Decision: keep measured diagnostics, continue both live G1 chains, do not
 stop from a rebound or open G2 from an unconverged provisional energy ranking.
 No acceleration performance, new phase identity or independent prediction
 is inferred. Receipt: `residual_update_20261009_0450/validation_delivery.json`.
- At04:58:23CST both parents are still Running:PO--Mstep31/.120496,
 preservingstep28/.143216; reversing continuation remains PendingDependency.
 These later log values are not substituted for the frozen raw-evidence frames.
- Staged treecea6266a/archive77001836 preserves all37 evidence files exactly;
 clean numeric replay is byte-identical and30focused tests pass/28.26s. User
 changes are excluded. Subsequent updates are README/receipt/journal metadata
 only, not source, tests or numeric inputs. Raw-byte Git attributes prevent
 line-ending conversion from invalidating the archived evidence checksums.

## E029 - Second ordinary pass and a changing preserving-flip profile

- PO--M28298794 ends06:29:53/exit0:0,step39; terminal summary and complete
 raw-SCF replay agree on ordinary residual.097904490eV/A. Sampled forward/
 reverse maxima71.582207/142.966807meV/f.u.; M endpoint P2_1/c at all three
 registered tolerances. Its sampled peak's physical tangent-.383495004eV/A
 prevents promoting an ordinary pass to a stationary-TS certificate.
- Preservingstep48has peaks3/5at39.069477/39.069483meV/f.u.; central4/Pbcn
 is lower17.440944but retains the dominantcell residual.204908776eV/A.
 The current bottleneck is not the old central Pbcn maximum; no stable
 intermediate, optimizer cause or acceleration advantage is certified.
- Reversing28319571starts automatically when M ends, not by resubmission.
 At08:47 both flips are Running, preserved48/.204909 and reverse12/.455054.
 Keep at most two active chains and the existing bounded continuations.
- New dated metadata adapter reuses the original strict physical analysis;
 legacy source and historical figure remain byte-unchanged/replayable. Two
 ordinary passes do not certify G1/TS/sampling/predictions.35focused tests
 pass/58.42s, five new cases. Full independent/clean checks recorded next.
-27newly frozen plus10reused image records are not new or independent DFT.
 Three-frame HF/local1923floats agree<=2.85e-14,365nonfloats exact. No new
 SCF, input change, job mutation, earlyG2 or holdout is used for the update.
 Evidence: `morning_update_20261009_0850/validation_delivery.json`.
- Clean staged treed12dd301 passes1078tests/2skip/313warnings in126.51s;
 code/test hashes remain unchanged after this check. All40 initial evidence
 files retain their bytes. Clean network replay is byte-identical to local.
 The exact tested archive b59a331c also runs the same CLI on HF in an isolated
 directory:697floating fields agree<=2.85e-14,882physical/source nonfloat
 fields exact. Versions and38warning fields differ and remain in both reports;
 no remote pytest or installation is claimed. Subsequent changes are delivery
 metadata and the retained HF report only, not implementation/numeric inputs.
-09:08:42CST Slurm confirms preserving50/.199713 and reversing13/.435211
 still Running. Later rounded logs do not replace the frozenstep48/12SCF data.
 No further DFT call, scheduler mutation, parameter change or production source
 overwrite is needed for this evidence delivery.

## E030 - Audited morning evidence integrated into the main manuscript

- The main four-channel table now uses E029 complete steps6/39/48/12 and
 their actual jobs, rather than leaving the second ordinary pass only in
 the pilot appendix. All four rows, rounding, pass labels and local links
 are checked against the frozen report. M's nonzero physical peak tangent,
 preserving double peaks and residual-block changes retain their limits.
- Both figure files and the earlier figure generator remain unchanged.
 Figure2 explicitly states its earlierstep6/20/14/10 source stage and is
 not represented as the latest table or a converged mechanism ranking.
- Existing LuaHBTeX/TeXLive compiles the single-source manuscript to8pages.
 All8Poppler renders actually reviewed; PDF retains10 checked numerical
 tokens and25URI links/12local targets, none missing. No problematic TeX
 log messages, install or new test-performance claim. Current and historical
 receipts are distinct; the old7page checksum is not the current PDF's hash.
- This is measured evidence integration, not G2 response or held-out prediction
 completion. No new DFT/job/parameter/source mutation, CI or broader matrix.
 Receipt: `paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.morning-update.validation.json`.
-09:22:11CST authoritative Slurm check confirms both existing switching jobs
 Running onhfacnormal01:preserving51/.197762,node11;reversing14/.415792,node26.
 These later rounded logs are not replaced into the frozenpaper table. M's
 terminal0:0is unchanged. Continue the same two handles, not a duplicate
 submission or a new matrix; remaining material gates are still incomplete.

## E031 - Fixed-control restricted stationary continuation

- Necessary gap: stable orthogonal condensation alone does not follow a
 minimum or first-order saddle as external strain changes. Retain its affine
 offset, solve the resolved internal stationary block and leave the control
 fixed. The supplied actual mechanical space contains Q/R but excludes c;
 unrepresented movable gradients and clamped reactions remain distinct.
- Standard harmonic continuation can harden a saddle's controlled curvature
 even when stable-mode release softens a block. This sign check is not new
 theory, a hafnia result or a barrier forecast. No pseudoinverse, automatic
 branch discovery or physically forbidden release is introduced.
-94focused pass/1.93s,38new cases. An initial test used a nonexistent boundary
 factory; corrected only that test to the canonical existing API. Clean tree
 63d310c4 excludes user changes:1116pass2skip313warnings/127.54s.
-8analytic groups/24points verify separate fixed-control linear solves and
 finite energy derivatives;0DFT/0calculator. Exact archived source also runs
 on HF:124floats agree<=1.51e-13,83nonfloats/source hashes exact; NumPy
 versions remain explicit. No remote pytest or installation is claimed.
- Original working report's raw hash check failed for one legacy module:
 workingLF versus archivedCRLF. Both raw hashes/reports are retained; exact
 newline-normalized byte equality is proven without editing that legacy file.
 Canonical clean-vs-HF source checks remain strict. New code/test bytes stay
 unchanged after the full suite; later updates are delivery documentation.
-09:55:19CST both G1 handles still Running:preserving54/.186039,node11;
 reversing17/.359954,node26. Production archive df4c12ee remains unchanged.
 No new SCF/job/mutation, calculator retuning, earlyG2 or holdout. Complete
 B2--B5 freezing and independent material advantage still require real data.
 Receipt: `stationary_branch/validation_delivery.json`.

## E032 - Pair the initial-well and bottleneck stationary response

- Added common physical/boundary/control/energy declarations and independent
  local parameter scales. Pair index0 and index1 at the same external parameter;
  retain raw-centre gap, stationary-anchor corrections and each energy response.
  Signed local gaps are not clipped or called certified activation barriers.
- Explicit caller-declared intervals/radii reject extrapolation and oversized
  offsets. Omitted admissible gradients and clamped reactions remain distinct;
  declarations/bounds alone do not prove full stability or probe-hull coverage.
- Initial test run: 109 pass / 1 fail, because a guessed omitted-gradient >0.1
  contradicted its direct value 0.071. Replaced that test assertion by the
  independent gradient calculation, not a revised physical threshold.
- Focused112pass/1.97s,46new cases. Clean staged tree ccfdb0b5 excludes user
  edits:1162pass2skip313warnings/133.93s,0errors/failures. New source/test
  bytes did not change after testing; later edits are delivery documentation.
- Eight analytic groups/24points and independent solves/finite derivatives,
  0DFT/0calculator. Same-source HF replay:196float fields agree<=2.09e-12,
  85nonfloats/source hashes exact, with NumPy versions retained. Original
  working and canonical reports preserve the unchanged legacy LF/CRLF hashes.
-10:11:33CST existing G1 handles still RUNNING:preserving55/.182393,node11;
  reversing18/.343628,node26. M remains COMPLETED0:0. Production archive
  df4c12ee unchanged; no job mutation, retuning, early G2 or holdout.
- Actual matched-boundary training, full B2--B5 selection/freezing, independent
  advantage and JCTC readiness remain unproved. Retain this necessary pairing
  step, not a substitute end state. Receipt:`stationary_gap/validation_delivery.json`.

## E033 - Put the controlled stationary response into the main manuscript

- Main Section2.3 now derives fixed-control internal stationarity, its indefinite
  saddle curvature correction and paired initial/bottleneck gap derivatives.
  Raw gap, anchor corrections, separate coordinate scales, omitted internal
  gradients, signed gaps and actual-TS limitations remain explicit. Standard
  quadratic mathematics is not presented as our new theorem or material result.
- Existing protocol permits visible audited target endpoint energies shared by
  every control. The text states this endpoint-assisted route, its cost and
  forecast lock; it does not silently add an initial-well Hessian centre to G3.
  Protocol hashes, two-bottleneck budget and original G1->G2->G3 order unchanged.
- One closest publisher extraction also compares strain-mode prediction with
  DFT(Fig.7), so this alone is not novelty. Added its known025DOI citation; new
  direct page opens fail, SI/dual-DOI relation still unresolved, no full-read
  claim or subagents. Primary-source limits are in the positioning sidecar.
- Existing LuaLaTeX/skill compiler succeeds with no new install or shell escape.
  Nine pages actually rendered/reviewed: equations, unchanged figures and table
  legible, no clipping/black glyphs/overlap. All10extracted tokens,28PDFURIs,
  14localPDFtargets and48Markdowntargets checked. Working layout retains large
  figure-preceding space and short appendix; not journal typesetting certified.
- Source/test bytes remain those of the1162pass E032 suite; no new pytest is
  claimed. Old material table/figures and old compilation receipts are retained.
-10:45:33CST same two handlesRUNNING:preserving58/.162157,node11;
  reversing21/.319623,node26. M remainsCOMPLETED0:0. No new SCF, retuning,
  production/source/job mutation, earlyG2 or holdout. Actual material and
  independent-prediction gates still incomplete. Keep the goal active.
- Receipt:`../../../paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.stationary-response.validation.json`.

## E034 - Lower central snapshot and tolerance-sensitive reversing maximum

- Freeze two complete existing switching observations, preserving59 and
  reversing22. Exact32MPI/SCF/input/geometry/raw-hash/log replay passes with
  zero new DFT. Reuse19 terminal T--PO/M records,37records including duplicates.
- Preserving peaks3/5 remain~35.045730meV/f.u.; central Pbcn is now
  -1.605722 below shared PO+, yet atomic/cell residuals .111271/.156442.
  It is not a stable basin or TS. Its original-T-x expansion may change the
  branch under the already registered plane: test this in existing G2, no
  extra Pbcn endpoint/path/Hessian or premature production.
- Reversing peak4 is Pbca/Pbca/Pa-3 at the original three symmetry tolerances.
  Whole chain .310501 is atomic-dominated at2 despite peak residual .050829;
  never select a tolerance or a small peak force as joint-TS certification.
- Unchanged clean E032 source passes35focused tests/58.74s and reproduces
  both new local reports byte-identically. HF/local compares1298floats/249
  nonfloats(two observations),697floats/882nonfloats(network), maximum
  difference2.842171e-14. Three package and38warning-field differences are
  explicitly retained, not normalized into agreement. No new full/remote
  pytest or install claim. Git attributes preserve all evidence bytes.
-11:17:53same handlesRUNNING:preserving60/.151932,node11;
  reversing23/.300858,node26. Recent11--15min steps and declining force
  support only conditional13--15h/15--18h estimates, not promised finish.
  Audit a terminal step cap before any further continuation. No job mutation.
- Bounded SI follow-up confirms Behara subscription requirement; public SI
  and Lee direct/PDF access still fail. No missing text read or novelty proof,
  no bypass or agents. G1/G2/G3 gates and protocols unchanged; goal active.
- Evidence:`switching_update_20261009_1105/validation_delivery.json`.

## E035 - Resolve the current Lee registry identity, not scientific novelty

- Official Crossref exact-DOI calls give one current025 record. A no-follow
  GET at11:37:08CST confirms026 HTTP301/Location ->025; two downloaded JSON
  payloads are23901bytes with the sameSHA. Four bounded native GETs total,
  no filters/pagination/authentication or new full-text/SI read.
- First026 exact-identity assertion rejected the returned025 despite a
  successful download. This was not an HTTP error; the explicit redirect
  resolves that guard failure. Do not infer why the identifier changed or
  equate indexed reference PDFs from equal metadata alone.
- Correct publication date3March2026;8April is metadata deposit/index, not
  a verified version-of-record date. Main manuscript already cites025;
  its code/tests/body/PDF, protocols, physical inputs and finite budget stay
  unchanged. No DFT/job mutation or earlyG2/holdout, no new pytest/compile.
- Bibliography/provenance correction only: no new material prediction or
  JCTC-ready claim. Raw responses outsideGit, reproducible URL/fields/hash
  in `../../../outputs/HFO2_LEE_DOI_METADATA_AUDIT_2026-10-09.json`.

## E036 - Preflight the real G2 delivery without evaluating a calculator

- Archive only committed 1e9259c source to a fresh local/HF namespace;
  22,725,666 bytes, SHA14420324. Runtime/test directories are unchanged from
  E032's recorded 1162-pass/2-skip suite; that suite is not rerun. New clean
  focused endpoint/epitaxial/input regression:48pass,0fail/error/skip,1.64s.
- Initial launcher was rejected before process creation because its workdir
  did not yet exist. Split archive/extract and clean pytest commands; no
  numerical or input-threshold repair. Preserve the failure in the receipt.
- Separately hashed case helper actually executes twice with existing HF
  ICU/ASE3.23.1b1/NumPy1.26.4. All10 seeds and6 physical-input hashes pass;
  constructor results stay empty, reserved directory stays absent, bash-n
  succeeds. No calculator evaluation, DFT, installed software or submission.
- This is delivery readiness only. Starters are still unrelaxed, phase/variant
  and G1/G2 material gates remain unpassed; holdout absent. The registered
  four-step PO+ BFGS canary still waits for G1, not auto-submitted by the check.
-12:04:16CST same two handles RUNNING:preserving64/.130327 on node11,
  reversing27/.257145 on node26. Production archive df4c12ee unchanged;
  ordinary0.10/100Ry/full10auDZP and at most two live study chains preserved.
- Keep finite G1->G2->G3 and actual independent evidence requirements. No
  material advantage/JCTC completion verdict. Evidence:
  `g2_delivery_preflight/validation_delivery.json`.

## E037 - Fix the actual clamped production entry before G2

- Found a delivery gap: API/endpoint clamping existed, but generic production
  full/fixed choices could not preserve the registered substrate. Do not
  reuse the G1 free-cell invocation as G2 or call the old E036 archive ready
  for this new path interface. Add explicit config/CLI selection, shared cell
  scale, raw-chain rejection and candidate-plane/geometry checks; no implicit
  mode-artifact composition, parameter retuning or boundary guessing.
- Clean final staged-source suite1197pass/2skip/0fail/error/313warnings,
  146.56s. Oblique normal/tilt tests include actual few EMT steps, both public
  prepare/run and direct runner, serial/threaded endpoint-once and unchanged
  factory options. Initial failures were new test roundoff/schema mistakes,
  fixed without relaxing mechanical or DFT tolerances. Earlier staged suite
  is recorded separately, not substituted for the final delivery.
- Final archive675bcce0 is transferred and actually checked once on HF.
  Ten valid unrelaxed seed paths and ten corrupted internal planes give the
  expected accept/reject results with zero calculator symbols loaded, zero
  DFT or submissions. Rejected internal planes are never projected into a
  claimed path. This proves entry integration only, not G2 material results.
-12:40:41same two32CPU handlesRUNNING:preserving68/.104515 on node11,
  reversing30/.220242 on node26; productiondf4c12ee unchanged. Preserve
  ordinary0.10/100Ry/full10auDZP, finiteG1->G2->G3 and max two live chains.
  No CI/cancellation/restart or earlyG2/holdout. Need terminal complete G1
  audit before the registered PO+ canary; physical forecasts remain missing.
- Evidence:`clamped_path_entry/validation_delivery.json`; no full-goal or
  JCTC-ready claim. Milestone receives one normal owned-files commit/push.

## E038 - Protect the reverse existing-directory boundary change

- Follow-up inspection: E037 refused a new clamped request against a
  different boundary, but omitting clamped flags could reuse that directory
  as free-cell. Guard both directions; fresh legacy runs are unchanged.
  Two additional regression cases preserve the old preflight bytes.
- Focused101pass/4.38s; final clean staged archiveac0532f9 gives
  1199pass/2skip/0fail/error/313warnings,137.82s. No older pass count is
  substituted for this changed source. Scientific/runtime solver and
  physical inputs remain unchanged; only a pre-calculator directory guard.
- The same archive actually tested clamped-to-free reuse on HF's prior
  geometry-only directory. Expected FileExistsError/entry exit1 is the
  safety success, shell check exits0, preflight SHA3a606216 identical before
  and after. Zero calculator loading/DFT/submissions; live chains untouched.
- Keep E037 data/source history, G1->G2->G3 gates and independent material
  evidence requirements. No new live poll, full-goal/JCTC verdict or new
  matrix. Evidence:`clamped_path_entry/boundary_guard_delivery.json`.

## E039 - Third ordinary pass, split peaks and physical-tangent limitation

- Actual preserving job28319570 ends COMPLETED0:0 at12:50:13, terminal reason
  force_threshold. Complete step69 replays0.099372458 at unchanged ordinary0.10;
  two Pca2_1 sampled peaks32.806meV/HfO2 flank a Pbcn centre-14.838meV/HfO2
  below common PO+. True side-peak tangents approximately+0.242486/-0.242483
  eV/A prohibit a stationary-TS claim. Lower centre has unmeasured stability;
  do not add a Pbcn endpoint/Hessian or call two certified saddle segments.
- Freeze reverting step32/.200001508 and retain its Pbca/Pa-3/Pa-3 tolerance
  sweep.13:23:51 same handle28319571 RUNNING on node26/hfacnormal01/32CPU,
  later rounded log step34/.186340. Conditional recent-trend estimate roughly
 15:30--18h today; not convergence assurance or reason to retune/restart.
- Both unchanged-code HF/local physical replays pass:1298 floats/249 exact
  two-frame fields(max1.42e-14),697 floats/882 exact network fields(max2.84e-14).
  Preserve all package-version/warning differences.18new frozen records+19
  reused, not37newSCFs. No DFT/job changes; productiondf4c12ee unchanged.
- New dated four-panel figure explicitly separates low/high-energy profiles,
  ordinary residuals and preserving physical tangents. Historical figure and
  data unchanged. Actual 9-page reading PDF, all pages reviewed, table/image
  and15local targets checked; no invented abstract/prediction result.
- Initial staged archive omitted two ignored raw trajectories; source-bound
  tests correctly report7errors. Add exactly those owned files; no weakening
  of assertions or change to source/DFT. Final extracted archivee9ef0982 gives
 1207pass/2skip/0fail/error/313warnings,139.37s;36numeric/figure/code files
  match the Git archive byte-for-byte. Keep initial failed archive/report.
- This third ordinary pass is not fullG1, TS/sampling, matchedG2 or independent
  G3 evidence. Keep finite gate sequence,100Ry/full10auDZP,noCI/max2chains;
  do not submit the canary/holdout before G1. Goal active; one normal scoped
  material-result commit/push. Evidence:
  `switching_converged_update_20261009_1255/validation_delivery.json`.

## E040 - Six actual sampling checks after an ordinary pass

- Cached step69 physical derivatives bracket maxima inside2->3 and5->6.
  Register exactly6fractional/cell reconstruction statics,3per segment.
  Hermite38.218meV/f.u. is selection evidence, not a DFT result or holdout.
  Preparation sourceff28bff/archivee0679852 passes1218clean tests/2skip;
  actual HFpreflight constructs the same calculator with0SCF and verifies
  all6physical hashes. Generated geometry tails differ at2.22e-16 inq,
  cells identical; both preparations retained, no electronic normalization.
- Actual job28374431 COMPLETED0:0,13:52:25--14:02:26,32CPUs/node1. Six
  real32-rank SCFs converge with full E/F/stress and identical physical
  INPUT/KPT/UPF/orbital bytes. New highest sample38.508448meV/HfO2 exceeds
  the old sampled32.806023 by5.702425. It belongs to the specified linear
  reconstruction, not an independently relaxed continuousMEP or TS.
- The15total/13moving inserted cached band still replays0.099372458
  ordinary residual with0optimizer steps/extraSCF. No new whole-band job
  is submitted. This is evidence that ordinary convergence and sampled
  peak resolution are distinct, not a claim of a new force threshold.
- Actual SCFsum542.979671s/4.826486coreh, allocation601s/5.342222coreh.
  Preserve all6raw logs/structures/EFS/physicalhashes and32unique rank
  probes. Large repeated orbital/pseudo payloads stay onHF; their bytes
  were verified there. Source/productionarchives remain unchanged.
- Initial offline figure tests3fail/4pass because stock localASE has no
  ABACUS I/O. Reuse the actual fixed production writer for offline token
  geometry checks; actualHF ASE round-trips unchanged.7focused tests pass.
  No install or production/DFT change. Final clean deliverytree6afd6a0f /
  archiveb7b790bf passes1225tests/2skip/0fail/error,314warnings,374.41s.
  93raw case/figure/runtime-code files match archive bytes; no source-bound
  assertions weakened. Receipts/ledgers added later, runtime source unchanged.
- New two-panel figure and10-page reading draft use actual samples only;
  all10pages/figurepixels,17local targets and numeric captions verified.
  Historical nine-image figures/data keep their scope. Paper/figure skills
  enforce no smooth-MEP/TS or blind-prediction label on these adaptive data.
-14:15:16 same reversing28319571 RUNNING, complete step37/.163922.
  Continue it, not another chain. FullG1/G2/holdout/TS claims remain false;
  keep100Ry/full10auDZP/noCI/max2, finite13new-chain matrix and goalactive.
  Keep the audit and make one scoped normal research commit/push. Evidence:
  `preserving_sampling_bridge_20261009/validation_delivery.json`.

## E041 - Five actual peak checks and required M refinement

- Five true32MPI statics28380672 COMPLETED0:0,544s/4.835556allocation coreh;
  SCF481.761171s/4.282322coreh. Same frozen lifts, P0/commonPO and all6
  physical bytes. Added peaks115.221524(T/PO),82.747071(M), increases
  .011360/11.164864meV/HfO2. Hermite selects points, not a DFT/holdout result.
- Twelve-image residuals.060037/.123408: T passes; M does not. Zero optimizer
  or extraSCF for replay. Three commonPO wells exactly match; endpoint force/
  stress screens pass, not Hessian stability. Static samples are not TS/MEP.
- Initial archive995239cc omitted15ignored logs/err:4failed/1240passed/2skip,
  204.55s; actualHF refuses raw-evidence gap, no SCF/job. Scoped ignore repair,
  no weaker assertions; retain failed source. Fixedtree512f21b/archive63dd3824
  passes1244/2skip/0fail/error in218.96s,112raw files byte-match. Later seed/
  receipts and Table2 editorial layout are not attributed to compute archive.
- ActualHF12cache/CLI/lift preflight pass with DFT forbidden. NecessaryM
  refinement28392675 submitted16:19:58, node2/32MPI,12total/10moving/2fixed,
  20steps/8h/at most200movingSCF, freshFIRE/.02/k.2. Actual step0.123408 with
  zero exportSCF; first newlogDSIZE32. Old9-image result preserved, no extra
  scientific edge, state-restoration or acceleration claim. Maximum2live.
- Ten-page/Table2 draft compiled and all final pages/numbers reviewed,18links
  valid. No paper success placeholders. Other reversing28319571 RUNNING at
  16:19,step44/.103449. No retuning, CI, G2/holdout or extra matrix. Cap14
  peakSCFs(11done+3max),13newchannels unchanged. Goalactive/incomplete.
  Evidence:`G1_peak_sampling_20261009/validation_delivery.json`.

## E042 - Two real terminal passes, still not full G1

R28319571 and refinedM28392675 COMPLETED0:0 at16:32:13/16:38:30.
ActualHF raw SCF/input/lift/EFS/log audit and local frozen replay pass:
R45/9images/.094022846/392.822905meV-HfO2; M1/12images/.096964013/
82.679957.21records reuse caches, export0SCF. Four candidates share identical
PO initial well and pass ordinary.10. No acceleration claim from one M step.
Archivedtreeef52921d/SHA0758a85c passes9focused tests20.87s; runtime unchanged
from1244full-pass E041. R/M allocations321.244444/9.884444coreh.
Complete remainingRsampling(≤3) and mechanism/variant/polarization gates;
Pbcn midpoint stability remains unresolved, M/T are not global escape coverage.
No extra path/Hessian, G2/holdout/CI/retuning; goalactive/incomplete.
Evidence:`terminal_G1_update_20261009_17/validation.json`.

## E043 - Two reversing peak statics and actual terminal mechanism/gap audit

28410062 COMPLETED0:0,204s/node61/32CPUs. Two actual energies392.822010685/
392.822010698meV-HfO2 below the old392.822905197 maximum; Hermite's tiny
increase is not observed. Cached11-image residual.094022846 passes; no
additional optimizer/SCF. Peak statics13/14, spare not automatically spent.
SCF147.275712s/1.309117coreh vs allocation1.813333coreh. Same six inputs.
Four actual terminal bands/40images replayed; Pbcn minimum−14.838201 and
Pa-3 reversing peak at three stated tolerances, not stability/TS evidence.
Forty original band files copied with raw hashes; sampled min gaps4.252037/
4.131676/4.589282/4.017321eV, no full-BZ/Berry-branch claim, analysis0DFT.
Early incomplete/failed transfers and scriptCRLF/source-hash mismatch were
preflight-only; noDFT then, retained failure sources, LF immutable runtime
archive2106f514/treefe559455 passes1260full/2skip261.06s. Actual delivery
data/tests subsequently checked separately. Goalactive/incomplete; no
G2/CI/retuning. Next concrete14-interior property audit reuses3PO endpoints,
bounded14output-onlySCF+42NSCF; no NSCF barrier energies or automatic branch
choice. Pbcn/global-escape coverage caveat preserved, no extra matrix.
Evidence:`reversing_peak_sampling_20261009/validation_delivery.json`.

Final actual-data archive4ff1f62a/SHA b532c60d passes26focused25.11s;
92casefiles exact-byte matched. Three initial local preparation files needed
targeted Git renormalization after attributes; snapshot-hash test added,
no physical geometry/DFT change. Receipt/ledger updates are subsequent metadata.

## E044 - Actual switching-path Berry data and manuscript figures

Array28418996/0–13 allCOMPLETED0:0;14 output-onlySCFs+42 fixed-chargeNSCFs,
3endpoint results reused, no peak/channel/endpoint repeats. Original source
treee1fff0ab/archive30acfa0e passed1285full/2skip283.65s before HF execution;
source remains immutable. Actual6physical bytes/geometry/EFS/charge audited.
EFS reproduction differences0 at raw-log precision, not a total-error bound.
Max224→228 sensitivity.000637680/.000838431C/m²; nearest sampled lifts both
unique, margins.729131/.611964, equal noninteger reduced increment1.623691.
Do not claim distinct winding sectors, quantized pump, absolute spontaneousP,
full-BZ insulation, transverseP, continuous-path/TS/fullG1 certification.
Integer nativeDFT3155s/28.044444coreh vs Slurm4029s/35.813333coreh; genuine
32MPI and observed disjoint concurrentCPU sets. Max2 allocations retained.

395 observable files hash-bound and local/HF raw replay identical, no
licensed pseudo/orbital/charge redistribution. Post-DFT NumPybool JSON failure
fixed in a separate analysis snapshotb08b9a09/archiveb5c8d036 (1290full/2skip);
no DFT repeated and original runtime untouched. Exact known CRLF/LF native
adapter pair recognized; arbitrary code hashes still rejected. Later generic
finite NumPy scalar JSON serialization and real-data figure regressions pass
27focused/27.86s. Failed figure export retained recoverably outside Git;
wrong-cwd partial test stopped only own exactPID, excluded from evidence.

New40record terminal network and18row Berry figures retain CSV/vector/raw
sources and actual pixel QA, regular(a) labels/aligned axes/top-right spines/
off-curve framed legends. Main manuscript updated without fabrication of
G2/G3/holdout results. Full final clean delivery regression and compile receipt
are recorded in `switching_path_polarization_20261009/validation_delivery.json`.
Keep measured data; winding interpretation inconclusive. Next formal G1review
then registered clampedPO+canary/G2, no unbounded extra property matrix.
Physical inputs/ordinary.10/noCI/13chains/14peakstatics remain fixed; goalactive.

Final clean treeb2b5905b/archive43d2993a passes1297full/2skip/0fail248.14s,
all395 observable bytes inventory-exact in the archive. Subsequent changes
are manuscript grammar/spacing and this delivery metadata only, not code,
Slurm, physical inputs or raw results. LuaLaTeX exit0,11pages,21local links,
no overfull/fatal log; actual page/figure inspection recorded in receipt.

## E045 - Whole-G1 candidate review and one clamped PO+ canary

Actual local replay and fresh HF40raw-record review both pass the registered
candidate-connectivity gate: common PO+, four terminal ordinary.10 chains,
free-endpoint force/stress screens, three symmetry tolerances, thirteen bounded
peak statics and E044 conditional sampled polarization lifts. Zero DFT for
this review. This is not continuous-MEP/TS/global-stability/escape-completeness
or distinct-winding certification; existing Pbcn coverage caveat remains.

Clean treefd738b1734/archive38c8c1da passes1316full/2skip/0fail269.19s and
63focused23.62s before execution. Actual HF read-only validator checks ten
seed geometries, shared substrates and the original six physical bytes without
evaluating a calculator. One registered zero-strain PO+ BFGS canary submitted:
28446324 at21:54:01CST,1node32CPU/hfacnormal01/directmpirun32,30min cap,
max4steps/5SCFs/.03atomic/2kbaropentraction/maxstep.02. No whole-G2/G3/holdout
submission, phase-restoring constraint, CI or physical parameter change.

At21:56:48 it is RUNNING166s/node6; native first-SCF DSIZE32 and ordinary
electronic iterations observed. Running sacct's0:0 is not a completion code.
Four-step cap with complete audited SCFs is canary health, not converged endpoint.
Actual final receipt is added only after observed completion. Preparation
command failures (Windows final-argumentCR and nonexistent ad hoc precheck
symbol) launched no DFT; existing read-only validator and final shell comment
resolved them without changing production code. Source remains immutable.

Evidence:`G1_review_clamped_canary_20261009/validation_delivery.json` and
`submission.json`; goalactive. The next decision is actual boundary/phase/SCF
health review, then only an explicitly bounded endpoint continuation if justified.

## E045 completion / E046 exact-cache continuation / E047 T endpoint / E048 pair

E04528446324 COMPLETED0:0,5SCFs/4steps,32trueMPI,original six physical bytes.
Terminal force.037333755/open6.948494kbar: healthy capped canary, not convergence.
All5 raw records exported; original HF audit retained. Local stock ASE lacks
general ABACUS STRU I/O; a strict case-writer Direct/Hf4O8 reader enables the
portable replay. No DFT repeated or original source overwritten.

Clean runtime tree2a5ca85c/archivea9e23aff passed1330full/2skip/0fail275.30s,
71focused1.81s before E046. Full HF raw/byte/geometry/log checks and prohibited
external-call preflight confirm zero-DFT cached initial E/F/stress. E04628453655
COMPLETED0:0; fresh BFGS Hessian,2steps/2newSCFs,force.022284634/open1.928434990,
Pca2_1at3tols.236.308852SCFs/2.100523coreh vs244s/2.168889allocationcoreh.
Cached initial SCF cost belongs to E045; not an acceleration benchmark.

E04728454923 uses the same immutable source and registered epsilon0T seed,
COMPLETED0:0.1SCF/0steps,force.000499915/open.009594677,P4_2/nmcat3tols;
125.001048SCFs/1.111120coreh vs130s/1.155556allocationcoreh. Same substrate
and six physical files. E(T)-E(PO+)=12.359975meVfu vs81.321193freecell.
This is a boundary-dependent well separation, not a barrier or H1/H2 test.

E04828456312M/28456313POminus-preserving submitted23:00:09CST,each1node32MPI,
20BFGSsteps/max21freshSCFs/2h,unchanged.03/2/.02/P0/E0/electroniccontract.
23:14bothRUNNING,5/2steps; no rebound stop/live repeat or unbounded restart.
Max2study allocations, unrelated user jobs untouched. No reversing-minus/+1%
matrix, G2bands/G3/holdout/CI/release/parameter tuning automatically submitted.

49+26+17observablefiles exported without pseudo/orb/charge/binaries. Actual-data
replay and common-substrate numerical tests added; first focused20pass5.00s.
ManuscriptSection3.6 added with endpoint screens and limitations, not fabricated
abstract/conclusion predictive gains. Final clean delivery and reading-build
checks are recorded separately; original runtime remains immutable. Goalactive.

Actual delivery tree7f5520e5/archive9b6bd7f9 passes1336full/2skip/0fail272.29s.
Only four metadata/portable-analysis reports added afterward: final tested
tree5cb82bfd/archive93c2ee35 passes33focused1.40s and all92exported observable
files are byte-exact. This final receipt/ledger update changes no production
code, Slurm script, geometry, physical input, raw result or manuscript.
Reading PDF12pages/23Markdown local targets/20PDF local-file actions all valid;
changed pages11/12 actually rendered and visually checked, no overfull/fatal.
Missing G2 barriers/G3/holdout results remain explicitly missing, not certified.

## E048 terminal actual-data update and E049 registered endpoint

Both E048 endpoints actually COMPLETED0:0. Preserving-minus28456313:7steps/
8SCFs,force.014888417/open1.692852091,Pca2_1at3tols,17.427534SCFcoreh vs
18.008889allocationcoreh. M28456312:17steps/18SCFs,force.018206569/open
.843275452,P2_1/cat3tols,22.707467SCFcoreh vs22.773333allocationcoreh.
Original six physical bytes/native32MPI/fullEFS reviewed on HF,224new raw
observable files exported without licensed basis/charge. Source/runtime never
overwritten. Zero-DFT invocation scripts retained; initial preserving audit
shell-final-blank-CR error occurred after successful audit, rerun corrected
read-only invocation exit0, no DFT repeated. Physical gate/space group is not
Hessian, ordered-pattern, electronic-polarity or G2barrier certification.

Nominalepsilon0four wells now screened. Actual M-minus-commonPO+=-92.540494
meVfu; preserving-minus=-.035259, not forced equal or called intrinsicbias.
Table3/Section3.6 updated, original two-endpoint build receipt historical.
Material replay tests35pass1.42s; final clean archive/build checks recorded
in `clamped_endpoint_matrix_20261009/validation_delivery.json`.

E04928464144 submitted23:50:34CST after actual E048terminals. Existing
registered reversing-minus seed only, same2a5ca85c/1330preDFTfull/source6bytes,
32MPI/20steps/max21SCFs/2h/.03/2/.02. One study jobRUNNING23:51:45, no
unknown/live repeat, CI, parameter change, new phase/strain or holdout matrix.
Earlier 23:29both-running receipt remains historical. 0133a96normalpush
verified localHEAD=GitHubmain, before these newly completed result updates.
Goalactive; next candidate endpoint/pattern audit then bounded registered
+1%training matrix, G3 and genuinely frozen independent predictions.

Final E048result tree848bd52b/archive84ef6c60 passes35focused1.75s; all316
owned observable files exact-byte matched. Production/helper source unchanged;
prior1336full regression does not falsely count the two subsequently added
actual-data tests. Four-well PDF12pages,23Markdown/20PDF local targets valid,
changed pages11/12 rendered/reviewed. Three space-group math labels corrected
after the test snapshot, no numeric/DFT/code change; updated compile receipt
preserves earlier two-endpoint history. Goal still lacks G2/G3/independent gain.

## E049 terminal and E050 registered training endpoints, 2026-10-10

E04928464144 completed00:11:23CST,7BFGS steps/8fresh SCFs, original six
physical bytes/native32MPI/full EFS/Pca2_1 three tolerances checked on HF.
Force.014888416801/open1.692852094521 pass. Cost11.061641SCF coreh and
11.102222allocation coreh;72observable files exported, no licensed basis,
charge or binary. First read-only audit passed then trailing PowerShell CR
caused shell exit1; trimmed retry exit0, no DFT repeated. First export stopped
before writing an archive on a mismatched filename; exact six endpoint filenames
were inspected and72explicit observable files exported, no broad scan/source edits.

New read-only ordered-gauge audit confirms parent-declared translation(.5,.5,0)
with permutation[2,3,0,1,9,8,11,10,5,4,7,6] relates the two minus seeds and
endpoints. Actual relaxed geometric residual1.85e-11A, E difference7.28e-12eV,
force covariance3e-10eV/A, stress covariance2.38e-12eV/A3. Same ordered
indices differ1.3615A: strict ordered cache remains ineligible; no production
remap/whole-path equivalence/topological distinction claimed. Dominant T-pattern
signs-.924790/+.924790 are fixed-reference geometric observations, not
origin-independent phase labels or electronic polarization. Five endpoint
representations are not five independent physical phases. Kept the mapping
review/tests because it prevents an incorrect scientific or cache identity claim.

After all earlier study allocations were verified COMPLETED in sacct, two
existing +1%training seeds T/PO+ passed original seedSHA/boundary/P0/E0/six
physical-byte preflight. E05028469832/28469833 submitted once and started
00:25:50CST, node6/node28,20steps/≤21SCFs/2h/1node32MPI/thread1/.03/2/.02.
Only two study calculations, unrelated pc jobs untouched. Immutable source
2a5ca85c/preDFT1330pass reused, not overwritten by new offline analysis.
Terminal squeue IDs were absent; sacct reconciled COMPLETED before submission,
not a new failure/restart. +1%M/two minus still unsubmitted;+.005holdout
uncreated/unread, no new material, G2band, G3, CI or physical parameter changes.

47focused actual-data/transport/variant checks pass1.51s in the working tree;
clean full delivery/build recorded separately, not inferred from this result.
Manuscript adds endpoint correspondence/cost37SCFs54.408286coreh excluding
canary, keeps G2/G3/independent predictions missing and goalactive.

## E050/E051 terminal review, 2026-10-10 10:16 CST

Both E050 jobs completed0:0: +1% PO+/T-derived each6BFGSsteps/7SCFs,
native32MPI and all original six physical bytes/raw EFS reviewed. Their
SCF costs6.358583/7.935382coreh are training evidence, not optimizer or
prediction benchmarks. T-derived is already Ccce before DFT because the
declared partial clamp breaks the two original short-axis equivalence;
all seven evaluated geometries retain the original single geometric T pattern.
An incorrect new replay assertion expected free-T P4_2/nmc: corrected the
assertion and added a before-DFT metric/pattern test, without modifying any
input, structure or production source. Actual well gap14.145883meVfu is
not a switching or escape maximum. Each64observable files exported.

E05128472493 +1%M completed0:0 at01:18:38,16steps/17SCFs, force.007755026/
open1.452557036 pass; all17native32MPI/raw full EFS/six physical bytes
and boundary replayed without new DFT. P2_1/c at all three tolerances,
M-minus-PO+ well gap-114.746973meVfu, not escape barrier. SCF19.369124
coreh vs19.466667allocation coreh,144observables exported without licensed
pseudo/basis/charge or binary. Eight of ten registered endpoint representations
screened, not eight distinct physical phases. Working-tree focused replay
51passed1.96s; this does not yet name a new clean full archive.

10:16:42squeue empty. Two+1%minus seeds unsubmitted; no surviving job,
live repeat, new material, CI or parameter change. G2barriers/G3/strong-control
independent predictions still missing;+.005holdout remains uncreated/unread.
Goal tool currently reports blocked, not an SCF failure. This status request
does not change the goal state or append new submissions; finish the tested
delivery of completed work, then obtain user direction to resume computation.

Clean staged tree4bd9e808/archivef8389c76 completed the full regression:
1354passed/2skipped/0failures/0errors,314compatibility warnings,298.55s.
All660owned observable files exact-byte matched with the archive. Current
13page PDF compiled exit0,25Markdown/23PDF local targets present,37link
actions, no overfull/fatal log entries; actual changed pages11-13 rendered
and visually checked. The delivery/build receipts preserve earlier history
and distinguish the verified executable/data tree from later editorial records.
Five unrelated user-owned tracked changes remain excluded from staging.

## E052/E053 — resumed bounded G2 research, 2026-10-10

The user requested continued research; the existing goal is active. E052
submitted only the two existing +1% minus endpoints. Both completed0:0
after6BFGSsteps/7freshSCFs; full raw original-six-byte/DSIZE32/EFS reviews
passed, Pca2_1 at all three tolerances, .013274962eV/A/.99580664kbar.
Total14freshSCFs14.250463coreh;140observable exports without licensed assets.
All ten training endpoint representations now pass, not ten distinct phases.
No forced symmetry, CI, parameter changes or unknown/live repeats.

E053 prepares zero-strain preserving flip and PO→M as the first matched
G2 comparison. Nine total images = seven fresh internal images plus two
fixed, hash-pinned *clamped* wells. Frozen G1 geometry contributes the
starting mechanism only, never free-cell E/F/stress. The ABACUS STRU writer
wraps Direct positions: record an explicit short-step lift of the complete
raw endpoint history (including the PO canary before its restart), anchor
the original G1 endpoint lift to the registered unrelaxed seed, and smoothly
apply the BFGS endpoint corrections. No atom reordering/nearest endpoint
fitting or claim of distinct winding is made. The continuous lift and full
geometry/boundary preflight must pass before any DFT.

Each first G2 segment is limited to10FIREsteps/4h/oneHFnode/32MPI/thread1;
no more than two study allocations concurrently. A healthy step-limit is
not an SCF failure and remains the same chain for any later reviewed
continuation. Four channels × two training conditions remains the <=8chain
cap. +0.5%holdout/G3 surfaces are ungenerated. G2/G3/strong controls and
independent predictions remain missing, so no JCTC-ready claim is made.

## E054 — same-boundary G2 review and exact-cache continuation, 2026-10-10

Both E053first segments completed0:0 at14:12:02/13:45:24CST, healthy10step
caps, not convergence. Preserving flip fmax3.541469→1.466699; PO→M
.994612→.615276. Exact same-clamped-subspace replay confirms both logs;
the free-cell exporter cannot stand in for this mechanical ensemble.
All154fresh internal SCFs retain six original physical bytes, nativeDSIZE32
and complete E/F/stress; no new endpoint or analysis SCF. SCF transport
cost154.451016coreh, scheduler allocation155.822222coreh. These unfinished
discrete peaks287.182317/114.427812meV/f.u. cannot rank channels or label H1.
Worst residuals are real transverse forces: flip image4 has1.466700 versus
spring.000085eV/A; M image7 has.619056 versus spring.020765. No evidence
for changing cutoff, spring or symmetry, and no stop for an isolated rebound.

New narrow exporter requires the pinned production script, runtime metric,
substrate reference and exact ordered raw SCF cache. Fresh geometry
continuation prepares nine caches, reuses the complete current frame and
keeps endpoints fixed. Interior caches invalidate after motion. New FIRE
state is explicit, not falsely described as velocity/time-step recovery.
Continue at most20steps/12h per segment, same two independent G2chains.

Two zero-DFT delivery failures are retained, not hidden: final PowerShell
CRLF makes the here-doc delimiter a trailing Python`PY`token *after* full
receipts were written; separate18-cache verification passed without
re-export/overwrite. First clean regression1392passed/7failed/2skipped
correctly caught Git archive's LF→CRLF drift of unprotected E053substrate
(984→1004bytes). Protect both G2artifact namespaces bytewise; fixed clean
archive's26new audit tests pass, no relaxed gate or physical change.
Portable4snapshot/36image numeric replay needs no licensed assets and
matches HFfmax to≤1.12e-16. Full regression and actual continuation handles
are separate receipts; submit only after the clean full suite passes.

Fixed clean archive2b03239full regression passed1401/2skipped/0failed/errors
(666.41s), JUnit9bae00ea…. Independent HF verification exited0; identical
executable files compared between old analysis and fresh fixed archive.
Actual continuations28661019/28661020started20:07:55/56CST, bothRUNNING
at20:09:08, node6/node25, hfacnormal01/32CPU. New step0exactly reproduces
old step10cached force/energy; first new interior logs showDSIZE32.
Expected20update segment duration~4.79/3.99h from actual preceding SCFs,
roughOct11midnight–02:00, not a promise to reach.10. No other jobs changed.

## E055 — Same-data nested response controls and full visible-call cost

Hypothesis: measured physical Hessian actions can support explicit frozen,
atomic, joint and training-instability-promotion controls without inventing
missing stable blocks, silently discarding affine gradients or charging
shared/failed/unused preparation as free. No material labels are added.

New calculator-independent interfaces preserve measured action coverage,
internal index, prescribed-control scale/domain, well and bottleneck response,
signed gap and omission diagnostics. B5 is only local training promotion,
not an anharmonic branch search. Stable B4=B5 must agree; equivalence is not
an accuracy improvement or a novel Schur theorem. Strong T/Cmma/gauge matched
references and complete future forecast freeze remain required.

175focused tests pass; immutable86d2637clean full suite1459pass/2skip/0fail
in318.96s.48independent paired analytic points replay on HF1.26.4NumPy;
local/HF response residuals~1e-17eV/fu are implementation checks, not DFT
error estimates. HF has no pytest: retain that zero-DFT harness failure,
separate no-install numeric/cost check, do not claim a remote pytest pass.
The actual completed E053154interior calls independently re-audited retain
all original six input bytes, full native EFS and DSIZE32. Cost recomputes
154.451016coreh versus155.822222scheduler allocation; endpoint/preparation/
ongoing continuations are explicitly outside this partial cost scope.

Actual unfinished clamped observations are refused as complete training;
no predictions emitted, no+.005structure/label read, no new DFT/independent
chain/allocation, CI, parameter change, release or automation. Retain the
delivery's import/cwd/incorrect failure-filename checks and LF-versus-executed
CRLF source-hash failure. Exact execution-source bytes are archived instead
of pretending a later checkout is the measured runtime. Six delivery replay
tests pass with unchanged physical gates. See nested_response_E055_20261010.

20:59:23CST actual two E054 jobs remainRUNNING, no actualvcneb_failure.json:
flip step2/logfmax1.148299, Mstep3/.538259. These live log forces are not new
raw-snapshot/TS certificates. Mean fresh SCFs151.29/109.40s indicate another
~5.30/~3.62h to this20step cap (~Oct11 02:20/00:40); no convergence promise.
Goal active: PROGRESS (new tested code, independently checked cost and gates),
not a blocked wait. No further unchanged-state polling in this delivery.

## E056 - Remaining finite training starts and manuscript integration

The six unstarted registered G2 channels were actually prepared and
preflighted on HF at21:18:05CST: the other two zero-strain starts and all
four +1% starts. All54geometries and12fixed raw-EFS caches pass.338runtime
Python files match the immutable86d2637tested archive bytewise. G1 supplies
only ordered geometry; no free-cell labels contaminate clamped endpoints.
No DFT/endpoint rerun, job submission, source replacement or environment
installation occurred. The existing two continuations are not duplicated.

Eight registered training starts are ready, **not eight optimized barriers**.
Eight new material-delivery test cases plus21existing chain tests give
29passed/0fail/0errors/0skip in104.664s; HFpytest is not claimed.
The manuscript now records measured-column coverage/strong-control scope
and partial actual SCF cost. The detailed methods explain symmetric action
extension and local instability promotion, including required abstentions.
Existing LuaLaTeX builds the13page preview; changed pages4/12/13were
rendered and visually checked. No new material accuracy claim is made.

PROGRESS: reproducible remaining starts and evidence-bound paper methods,
not a blocked wait. The next gate stays finite G2 matched paths and errors,
then original-budget G3 and a complete prospective prediction freeze.
No+.005holdout structure/label, extra material, CI, release or new monitor.
Evidence: clamped_G2_remaining_E056_20261010 and the manuscript's separate
training-matrix validation receipt; historical receipts remain unchanged.

21:38:01CST actual E054 jobs remain RUNNING/32CPU, no true failure report:
flip step4/log force0.863395 and M step6/0.455959. Mean fresh SCFs163.825/
114.497seconds suggest segment caps Oct11~02:45/~00:45, not guaranteed
ordinary convergence. The first additional94-case run used working-tree
cwd despite preparing an archive; preserve its report and do not label it
clean. Delivery testing must use the explicit independent archive cwd.

Independent archive-cwd replay then passed94/0fail/0errors/0skip in107.75s,
tested tree2487220b4a696bbdcdf0022588229f9e793a5141, tarSHAfa98b50d....
All76then-present case files and both paper source bytes match;28local
paper targets exist. Both worktree and actual clean JUnit reports remain.
Post-test additions are documentation/receipts, not executable/data changes.

## E057 - Full-joint harmonic partition null and real-chart design counts

Hypothesis: an exact fully measured harmonic stationary response is invariant
to the Q/R retained subspace and dimension, not merely its coordinate gauge.
Five partitions, two gauges and index0/index1 centres give114stationary and
54paired checks against independent original-full-B solves. Two invalid
Schur releases correctly abstain; promotion restores the same full solution.
Frozen and joint energies differ. Keep the direct full-B solve as an
equal-information harmonic reference even when a chosen R block fails.
Pure partition/promotion equivalence cannot supply material accuracy novelty.

Ten actual registered starter charts count33atomic/3open-cell/1controlDOFs.
Full central actions require37columns,75uncached evaluations at one step or
149at two, centre included. This explicit-stencil count is neither a universal
algorithm minimum nor a new allowance. Reuse same-centre data across gauges;
do not silently charge Hessian preparation to81surface nodes or omit it.

Working/independent clean80tests pass(2.29/3.01s); HF original environment
replays the168analytic points and10chart counts with exact eight source bytes,
no pytest install/claim and zeroDFT/external executable calls. Preserve the
first checksum rejection of a still-live90MB upload (namespace absent), wait
for that actualSCP handle, then execute only after matching archiveSHA.
No production/source/input or held-out condition is changed. Manuscript
records the null and keeps material claims pending;13page PDF passes visual
checks on pages4/5 and all29local links resolve. New fixture checks expand
the delivery suite to82, separate from original80source tests.

22:10actual continuations remain RUNNING/32CPU, no true failure report:
flip step5/.725003, Mstep8/.403503. New means176.643/116.739s imply remaining
5.15/2.72h to segment cap, Oct11~03:20/~00:55, not convergence or TS proof.
PROGRESS: a verified falsification boundary changes the method comparison
and prevents duplicate curvature jobs/unsupported claims. Goal remainsactive;
next finite G2 labels, then explicitly budgeted G3 and independent prediction.

Final delivery archive21283742d8b6fd9f82fc4bdf944ac5bfb8ee8e18passes82related
cases in3.26s, zero failure/error/skip.20case files and script/test/paper source
bytes match. This is not an additional full-project or HFpytest claim.
Only result receipts/documentation are added after testing; user dirty files
stay un-staged and physical inputs/running source are unchanged.

## E058 - Second zero-strain G2 wave with execution-success dependencies

Advance two remaining registered starts without new concurrency: PO_to_T
28692775 and PO_flip_T_pattern_reversing28692776 accepted22:38:42/43;
both actuallyPENDING/Dependency at22:41. Both wait for successful completion
of28661019and28661020. Each unchanged ten-step/four-hour pilot is one32CPU
allocation, seven fresh serial32MPI images and two exact cached endpoints.
Only two study allocations run. This is a bounded independent channel wave,
not a monitor, retry, duplicated parent or new trial outside the eight G2 starts.

Actual HF preflight rechecks338immutable source files, original archive/pilot/
binary/preparation pins, all18geometries and four native32MPI/rawEFS endpoint
caches with no newDFT launch. The+1% four starts remain unsubmitted; the+.005
holdout is neither generated nor read. No CI, cutoff/input change, extra
Hessian/mode/material, publication or running-source overwrite.

Initial offline Windows checks exposed path-rendering/mock issues before any
remote execution; retained11pass/3fail, fixed14pass local and independent
two-file extracted source. The first checksum command's PowerShell quoting
failed without starting preflight; root absence and actual checksum reconciled
before corrected invocation. Existing chain/material tests plus14guards
pass43cases locally. Actual accepted handles and queued scheduler proof are
added as two separate offline fixture checks. No HFpytest installation/claim.

Parent observation22:37flip step7/.467214, Mstep10/.355683, bothRUNNING32CPU,
no actual failure report. Fresh SCF means177.17/117.95s give approximately
4.48/2.29h remaining to segment caps (Oct11~03:06/00:55), not .10 promises.
afterok is success, not convergence; timeout/real failure requires dependency
reconciliation, not automatic replacement. Journal exclusive creation and
per-attempt/accepted-handle flush protect partial/ambiguous submissions.

PROGRESS: real finite material jobs queued, goalactive/incomplete. Next matched
G2 training paths, then explicitly budgeted G3 and frozen independent forecast.
No material barrier/response advantage or TS/JCTC-ready claim from acceptance.

Final actual archive-cwd delivery passes45selected cases/0failure/error/skip
in45.82s, treece16b41102e9644ef075ef18287fc4564910c112. All19then-present
tracked case files and test bytes match. Initial comparison included ignored
bytecode and stopped before tests; correctedgit-ls-files comparison used the
same fresh archive. Later additions onlyreports/docs, not source/data. Existing
five user dirty files remain un-staged and retain their original numstat.
