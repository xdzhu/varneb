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
