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

## E059 - Actual well response and explicit endpoint-only B1 nulls

Hypothesis/decision: a measured strain effect on shared well energies must
be separated from a material bottleneck response before claiming a new mode
mechanism. Use ten existing0/+1% screened endpoint representations, no extra
DFT, strain, mode matrix, path selection or+.005holdout access. Add the
calculator-free physical row-cell directional work and separate imposed-plane
partials from residual atomic/open-cell work on the ordered endpoint chord.

PO+ shifts+6.806502meV/fu, M-15.399977, hence relative well change-22.206479.
Registered B1_fixed training increments are all-6.806502; B1_follow gives
T+1.785908/M-22.206479/both flips+.035259. These are seen training null
implications without absolute anchor barriers, not H1/H2 or forecasts.
Signed trapezoid defects0.019--0.347meV/fu and chord defects-.001--.204
combine unmeasured curvature/numerical/residual effects, never error bounds.

Local17new tests pass; related75pass10.82s and actual archive-cwd75pass9.13s,
tree76dd34f6812bb6cb4dda8cef3f323837b4fade4b. Eleven scientific/source/case
files byte-match. HF original Python/ASE/NumPy replays all points at23:15,
max numericdifference2.84e-14, exact eight runtime/archive bytes and all
ten original terminal PP/orbital/INPUT/KPT/STRU/log/native32MPI/fullEFS checks.
No environment install or HFpytest/DFT/source/job mutation is claimed.

First actual-HF fixture wrongly compares checkout raw bytes to Git archive;
75pass/1fail retained. Second LF-only fix still wrongly assumes archiveLF;
75pass/1fail retained. All three legacy files have identical line contents,
and all eight executed HF bytes already match the tested archive. Export
raw and LF-canonical checks separately for both representations; require
new science/paper bytes exact, allow only legacy newline normalization.
No calculator/input normalization or physical change follows these harness issues.

PaperTable4 avoids collision with earlier Tables1--3. Existing LuaLaTeX
compiles14pages, updated12--14visually checked; text/API/evidence links
resolve. The work function itself is not claimed as new theory or a predictor.
23:16parentsRUNNING/noactualfailure, flip8/.363645 andM13/.289108; current
segment cap estimates Oct11~03:28/00:53, not convergence deadlines. Existing
two dependency starts remain unchanged. PROGRESS: material null controls,
useful backend-independent gradient interface, raw HF replay and paper
integration. Goalactive/incomplete; matchedG2/G3/independent forecasts stillneeded.

Final corrected fixture passes76/0fail/error/skip locally10.510s and actual
clean archive cwd7.305s. Full clean-project regression1515pass/2skip/0fail/
0error in383.231s. Optional ase.io.abacus and licensed VASP POTCAR are the
two explicit skips. Tested tree88b9214d8185dba6ae92d761132459ebf2be5e73,
archive1b1b7c35b4a9dee19eb309a64f5985b2ed97a03cf79f653eca13f21d5625065a.
Science and manuscript exact bytes unchanged after tests; later reports/docs
only. Parent12hlimits end08:07Oct11, not an unrecorded4hdeadline. 23:30
four-handle snapshot remains twoRUNNING/twoDependency, no new job mutation.
README/manual expose the new physical work API and distinguish partials,
periodic lifts, endpoint wells and actual bottlenecks. Goalnotcomplete.

## E060 - Queue the remaining registered +1% matrix, no scope expansion

Previous goal turnE059 isPROGRESS: tested strain work/material well controls/
14-page paper sources pushedb8fd08b. Current evidence confirms healthy
four existing handles, preserving five user dirty files. The hypothesis gate
cannot advance from endpoint wells to barrier mechanisms without the full
matched eight-path training set; queue only the four already prepared starts.

HF00:04:13preflight verifies338immutable runtime code files, original pilot/
binary,36geometries and8raw/native32MPI/full-EFS/six-physical-file endpoint
caches;28interiors remain fresh. No new endpoint SCF, electronic setting,
condition, Hessian or+.005holdout is created/read. Ongoing parent calls
are not zero-cost. The exact helper and actual scheduler/source graph are
checked before the single-use submission, not inferred from manifests alone.

Actual00:05:58accepted+1%T/M28709788/89dependsafterok28692775:28692776;
flips28709790/91dependsafterok28709788:28709789. Independent00:08:31Slurm
proof shows2RUNNING+6Dependency and1node/32CPU/1controller/4h/10steps each.
All8registered uniqueG2chains now have handles, only2currentlyrun; submission
is not convergence, H1/H2, TS, material advantage or independent prediction.
Old jobs/sources untouched, original100Ry/full10auDZP/.10/noCI retained.

Initial70selected cases pass; after causal dependency/resource guards77pass
locally42.31s and in actual independent archive cwd47.40s (tree1409207c7bdcac5a044acd66828b13a18645da03,
archive7e4e65ba7832ae2b6708c3245206e5d95722134d772be0153af20e1879a5ee67).
Actual-delivery fixtures then79pass locally43.25s. No fresh full-project or
HFpytest claim; runtime86d2637previousfull1459pass/2skip remains historical.
Exclusive flushed journal preserves partial acceptance; failed/ambiguous
reply, unexpected handle/dependency or existing namespace blocks blind retry.

Restricted13-file export includes only receipt/source/cache JSON, no licensed
PP/orbital/charge/wavefunction. Evidence:clamped_G2_training_queue_E060_20261010.
PROGRESS: full authorized training matrix queued; goalactive/incomplete.
Next actions remain complete-chain convergence/input/mechanism audit, finite
two-bottleneck G3, registered strong controls, freeze and unopened holdout.

Final actual-delivery clean archive cwd79pass/0failure/error/skip49.464s;
local79pass43.247s. All21then-present tracked case/test bytes match,
tree5410b33ead6a3be1722e107f2fd5cb4abb94fa26 and archive
4adca7418798ebe27960aeb681a6ea3abb015b2e316846c9f8b16a0e77fb0cee.
Later records/docs only, executed helper unchanged. 00:20actual eight-handle
accounting remains2RUNNING+6PENDING, not failed or converged. Five user dirty
files retain original numstat and are not staged. Goalactive/incomplete.

## E061 - Exact historical sources without locking future paper prose

PreviousE060isPROGRESS, pushed84b88f8; actual00:27eight handles2RUNNING/
6Dependency, not terminal or restarted. Review the finiteG3/gates instead
of generating more material/condition/Hessian tasks. Identify a real E059
fixture defect: mutable current paper/library byte checks falsely reject
prose-only changes despite unchanged science. Preserve the exact historical
eight sources separately and strengthen current numerical replay coverage.

Actual00:34HFsource-only export rechecks original E059 receipt and eight
raw executed source hashes, including legacyCRLF;147582bytes become55195
byte zip5178718515dfdfae3a03b98b67d82b7c30270f587f1dea1c8de284619b8f24c5.
NoDFT/input/holdout read, job/source/parameter mutation, installation or new
HF numerical replay/pytest. Existing receipts/raw outputs stay immutable.

All3869current material numeric entries/key sets now match historical HF,
maxdelta2.84e-14;1e-9is computational tolerance, not aDFTerror bound. Six
work/well mutations and re-signed corrupted source fail as intended; simulated
prose changes pass without editing actual paper. Local84related tests pass
6.82s. Scientific runtime/input unchanged, source/archive check stays exact.
Evidence:replay_delivery_E061_20261011. Not a new barrier/G3/TS/forecast;
goalactive/incomplete and original finite matrix/+.005unopened retained.

Actual clean archive cwd84pass/0fail/error/skip8.146s; six tracked case/test
files byte-match. Tested treef5209fde9ae8626ecc3aceb8f6b0dc509c86e5cc,
archive7c9a427401be4fb599dcf3425332e4760e039f15b61d31cfa32c358023bcbfe8.
Later record/docs additions only. 00:43actual eight-handle accounting remains
2RUNNING/6PENDING, not terminal; no new job or payload revision is performed.
PROGRESS: remove a false future drafting failure and broaden actual-material
numeric oracle, not a new physics result. Five user dirty files preserved.

## E062 — completed PO→M segment audited; same-chain continuation actually running

E061 milestone pushed as3b209a5. Parent28661020 COMPLETED0:0 at00:49:09,
max-step exit not0.10 convergence. Actual01:01HF audit replays both complete
9-image frames and all140new interior nativeSCFs with original6physicalbytes,
DSIZE32/fullEFS. fmax.615276→.161869; discrete forward104.109701meV/FU is
an unfinished observation, not H1 label/channel rank. Exact9caches seed the
same chain, new FIRE state explicitly disclosed, no extra independent chain.

Only28722320 submitted, initiallyheld. Scheduler read-after-write visibility
errors failclosed; actualacceptedhandle/partialwrites corroborated separately,
both existing waiters updated once, originaljob released01:17:31 and RUNNING
32CPU start01:17:32. No resubmit or duplicatedependencyupdate, no oldsource or
DFTparameter mutation.01:20step0matchesaudited.161869; firstfreshnativeDSIZE32;
otherflipstep14/.220938.2active6pending; fixedfinitebudget/holdoutunchanged.

Registered-resume audit now separates exactinitialcaches from freshmoving
interiorSCFs. Actualguard/comment and Gitignoredframe/log packaging defects
recorded and repaired, not gates weakened. Final82local/82independentarchive
passes;149source/material/receipt files byte-match treebdf47db73ac2d3e094747403efc32d1319eae83b,
archive8c2d64dc4ec2d0271d0581b0f745ab8f11150fabc3087d6476ffb3814b040624.
EvolvingJUnitoutputsnotbyteproof; laterdocs/recordsonly. No freshfull/HFpytest,
G3/TS/forecast/newmechanism claim. PROGRESS: realterminalrawaudit, verified
samechainexecution and completeportableevidence. Overallgoalactive/incomplete.

## E063 — zero-DFT diagnosis excludes spring dominance of a short plateau

Actual 01:48:31 HF export of complete flip steps13/14/15, original six inputs,
native E/F/stress and DSIZE32, unchanged clamped metric. All334 archived
runtime Python files byte-checked. The limiting image5/O index8 remains
.221736→.220818; spring-only fmax .005341→.005117, whereas fixed-geometry
perpendicular-only fmax .220034→.219289. Open-cell-block max grows
.101105→.204300: atomic plateau plus significant evolving cell residual,
not a proven Hessian stiffness mechanism or a validated acceleration gain.

Decision: do not tune k/electronic inputs, stop on these three frames, or
add a duplicate trial. Keep the registered live segment, inspect its
complete terminal chain. No DFT/submission/source/job/holdout/G3 mutation.
Local and independent Git-archive related suites each59pass/0fail/0skip;
181 new source/material files byte-equal. Initial byte check included
normalized root Git metadata, then corrected its scope; no material loss.
01:53 two32CPU RUNNING/six PENDING, M step2/.150404. PROGRESS: measured
residual diagnosis and portable replay tool/evidence, not new barrier/TS
labels. Overall goal remains active/incomplete.

## E064 — complete forecast panel and correlated-reference scoring implemented

Bounded zero-DFT arithmetic experiment, not actual HfO2 forecast evidence.
Require every registered model/case, explicit abstention, one frozen model
artifact and shared comparison contract. Report coverage/maximum errors,
pairwise gain and strongest-available baseline envelope under the SAME
reference interval. Missing strong references are unavailable, not bad
accuracy. Identical B4/B5 and already-accurate B0 reject manufactured gain.
This does not fit/select/freeze B2--B5 or prove chronological label blinding.

48random intervals versus10001-point independent grids plus an adversarial
off-grid cusp verified; cusp distance7.4e-6 is grid discretization, not DFT
uncertainty. Actual02:10HF isolated NumPy kernel verified only; historical
source bytes and current full numeric replay are independent contracts.
Local122/independent-Git-archive122pass,0fail/0skip;14new source/synthetic
evidence files byte-matched. All actual materials/holdout/DFT untouched.
02:13twoRUNNING32CPU/sixPENDING; Mstep4/.134093, flipstep16/.238148.
PROGRESS: missing prospective comparison interface implemented and verified,
not a physics/new-theory result. Full goal active; real G2/G3/models and
independent material forecasts remain required.

## E065 — full clean-checkout software regression, declared skips retained

Unfiltered complete tests on independent archive of pushed753c122, not user
dirty worktree. OMP/MKL/OpenBLAS1 and unrelated pytest autoload disabled;
actual library/scoring imports verified inside archive.1627collected,
1625pass/2skip/0fail/0error in326.04s;314warnings retained, not repaired.
Skip optional ASEabacusreader absent and licensed VASPPOTCAR not distributed:
these exact two tests remain unverified, not silently counted as passes.
Complete individual JUnit outcomes, environment/source/archive proof saved;
no claim of full raw console, new HF full pytest or every dependency version.

No code/driver/test/physical input/job changes, no DFT/submission/holdout.
Five pre-existing user tracked edits and unrelated untracked files preserved.
02:32twoRUNNING32CPU/sixPENDING; Mstep5/.129359, flipstep17/.260776.
PROGRESS: whole-repository software scope actually validated after E063/E064,
not another narrow green check or material/JCTC completeness. Full goal
active; G2/G3/strong models/frozen forecasts/independent labels still missing.

## E066 — protect pending wave before terminal material audit

afterok also releases on a successful step-cap exit; it is not physical
convergence.22mock/local and22independent-archive guard tests passed before
one actual HF hold of each registered pending28692775/76 at02:54.02:55
separate read-only snapshot confirms both JobHeldUser/Priority0, exact original
ANDdependencies unchanged, two parents RUNNING32CPU/hfacnormal01. Four+1%
successors/other projects untouched.0new jobs/DFT/chains/input/source changes.
Release only after actual terminal audit and verified same-chain handoff, or
proof no continuation needed; ambiguous writes leave journal and no repeats.

Related82local/82independent-archive pass,0fail/skip;11new source/receipt files
byte-match testedtree05f35938. Evidence/clamped_terminal_guard_E066_20261011.
02:57latest complete Mstep7/.120235 and flipstep18/.272906. PROGRESS is a
real safe queue-state change/verified guard, not new material barriers,
strong-model forecasts or G3/holdout success. Full goal remains active.

## E067 — terminal flip material audit and actual same-chain continuation

28661019completed20updates at cap, fmax.270056 notordinary-converged.
Actual03:21HF audits140fresh SCFs original6physical bytes/DSIZE32/nativeEFS,
334immutable runtimePython files, complete0/15/20nine-frame observations,
and9exact latest caches with0newDFT. Discrete peaks287.182→224.982→216.508
meV/fu remain unfinished, not H1labels. Residual bottleneck shifts atomimage5
to open-cellimage4; spring-only.004511 is notdominant but no Hessian cure or
acceleration efficacy is inferred. Continue latest20, not minimum-force15.

One heldjob28722810; two dependency writes onceeach. Second readackstale>15s
failsclosed, separate reads provebothwrites, existingjob releasedwithout
resubmission/rewrite. Actual03:28:19start32CPU, seed.270056 and nativeDSIZE32
confirmed03:30; firstnewSCF completionnotclaimed. Successors remainheldfor
materialaudit; at most2allocations. No source/input/holdout/budget changes.

Fixture failures/partialjournals/execution sources retained. Corrected105
preparation tests pass locally/cleanarchive; final actual material/receipt
scope110local/110cleanarchive pass,0fail/skip,219new source/material/receipt
files byte-match treea320001b. Evidence/clamped_flip_continue_E067_20261011.
03:36two32CPU RUNNING; Mstep10/.107920,newflip0/.270056. PROGRESS: actual
terminal native evidence and safe same-chain resumed production, not full
G2/G3/strong-model/independent prediction or JCTC completeness. Goal active.

## E068 — actual ordinary M convergence and one rolling-slot admission

Native accounting proves28722320 COMPLETED0:0/step13/.099214088, cumulative43.
91fresh SCFs/334immutable runtime files/nine clamped frames audited on HF;
zero newDFT in audit. Forward98.458122/reverse190.998616/FS-IS-92.540494
meV/HfO2 are discrete-chain results, not TS/error/H1 labels. Final-segment
90.730841coreh excludes earlier preparation/segments. Manuscript Table5 and
discussion updated;14-pageLuaLaTeX compile passes, newpages13/14 visually
checked. Complete scientific goal remains active/incomplete.

Flip28722810 capped20/.197683, not converged. Both old successors' native
dependencies expired, holds retained. Persistent sacct source tuple handles
expired scontrol records without weakening scientific identity. Revised
source-r2 tested separately; never-executed source-r1 preserved. Actual
10:17:57 one release of existing28692775,0sbatch/0dependency writes; actual
10:17:58start32CPU and nativeDSIZE32.28692776held/four+1%pending unchanged;
next spare slot reserved for latest complete flip20 same-chain continuation.

117local/117independent-archive actual material and refusal tests pass;
74source/material/receipt files byte-match testedtree44752edd. Earlier106/115
preparation scopes are not actual convergence evidence. No new independent
chain/physical input/running source/holdout/G3/phonon budget changes. PROGRESS
is one actual material ordinary pass and existing-channel startup, not a
replacement software-only goal or method-advantage assertion.

## E069 — real latest-complete flip continuation with finite longer segment

Native10:40:53audit proves28722810successfulcap20/.197683006, not ordinary
convergence.140fresh SCFs/334immutable runtime files/nine exact caches audited;
0DFT audit/preflight, original six bytes/native32/full EFS. Unfinished
180.008887/180.044146meV/HfO2peaks are not converged switching/competition
labels. Latest20/cumulative50 geometry retained, not historical force minimum.

One heldsubmission+verifiedrelease28723655,0dependency writes; other28692776
staysheld. Actual10:41:39startup32CPU/node25; step0replays.197683;10:43:25
fresh nativeDSIZE32/stderr0. Two allocations with28692775, no unrelated-job
change. Extend finite new segment20→40updates/unchanged12h, all physical,
boundary, metric, spring, FIREstep, guards and.10/noCI unchanged. Reset FIRE
history explicitly recorded; no acceleration-efficiency claim. Previous
3:30:45/20updates suggests roughly7h/40if speed persists, not convergenceETA.

139local/archive preparation then141local/archive actual native replay and
single-use mutation tests pass,0fail/skip;91source/material/receipt files
byte-match testedtree4bed5496. Portable108entries/source capsules retained;
no new independent chain/G3/Hessian/Γ/holdout budget. PROGRESS actual material
audit and resumed production, not H1/H2/method gain/full G2/G3/forecast or
goal completion. Continue original finite scientific goal, not a software-only
replacement; five pre-existing user modifications remain preserved.

## E070 — primary phase-proxy applicability and manuscript validation

Targeted primary-source reading separates a phase-resolved proxy from the
initial/final-well response null. Manuscript introduction and stale G2 heading
corrected; no new baseline, predictor, material label or superiority claim.
Missing applicable evidence is not a poor model result. Access/version limits,
failed PDF screenshots and incomplete SI review are explicit in the report's
sidecar, not silently treated as evidence of novelty.

Existing LuaLaTeX compilation succeeds:14pages, final pages1/2/14 visually
checked, ordinary M Table5 values preserved. Editorial validation/source
digests are in phase_proxy_audit_E070_20261011/validation.json; no new software
test or material convergence is claimed.11:00:40HF read-only snapshot shows
two32CPU live bands, latest complete PO-to-T1/.323053 and flip1/.197669.
No DFT/job mutation, input/source/budget change, model freeze or holdout read.
PROGRESS is a primary-evidence-backed manuscript correction; full finite
G2/G3/strong-control/prospective-material goal remains active and incomplete.

## E071 — real clamped-chain reference transport, not extra DFT

Reparse9portable native terminal frames and replay the original clamped
residual. An initial trial omitted the prescribed G2 cyclic Cartesian/cell
frame; its misleading coverage is rejected and retained with executed bytes.
The corrected geometry-only inverse uses the original seed provenance,
identity atom order and all4E023 reference registrations. No production
geometry/EFS/input/optimizer or reference selection changes.

Actual Mchain image3 versus FS fractions: Ttriplet57.9237→24.4400%, Cmma
rank4span36.4746–36.4865→56.3537–56.3667%. These are atomic squared-displacement
fractions, not energy partitions or predictive gains; cell strains separate.
No G3 bottleneck chosen, new curvature/Γ/chain budget or holdout access.

Preparation97local/97clean-archive pass and actual independent material replay
match all parsed values except the explicit CRLF/LF source-SHA metadata.
Final99local/99clean complete JUnit pass,0fail/error/skip. Finalclean exec
session ID was not retained; no direct exit-code claim/restart. Executed source
capsules preserve raw bytes. Manuscript Table6,15-page compile and newpages14/15
visual check complete.11:15two32CPU bands alive, Tstep2/.211448,flip3/.197376.
Evidence/clamped_reference_E071_20261011. PROGRESS is an actual material
analysis and tested source-level frame contract, not G2/G3/H1/H2/forecast
completion. Full scientific goal active; five user edits preserved.

## E072 — joint clamped escape-path figure, not new material labels

Build the manuscript's Figure5 from frozen E068/E071 scalar reports. Nine
discrete energies, atomic squared-norm coverage and original-T finite Green
strain share the unmodified generalized arc. All four Cmma frames retained;
no smoothing, energy partition, local phonon, TS certificate, G3 selection
or holdout access. Open-vector Exy reaches -7.5209% at M, despite zero imposed
substrate strain; this is geometry, not a measured shear energy benefit.

New figure source and same-report/frame checks pass28 locally and28 from the
actual independent archive root, direct exit0. A preliminary mislocated
"clean" command is excluded. Source CSV byte-reproduces all9rows; same-host
libraries, not a cross-environment rendering guarantee. Standalone PNG and
final manuscript15/16 visually checked; embedded PDF fonts, editable SVG,
16-page LuaLaTeX success without layout warnings. Evidence is the dated
paper figure bundle; no new DFT/scheduler writes/physical or budget change.
11:32:54two32CPU running: Tstep4/.127622,flip4/.197105; held/dependencies
unchanged. PROGRESS is visualization, tested delivery and manuscript
integration, not full G2/G3/H1/H2/forecast or JCTC readiness. Goal active.

## E073 — permitted physical gradient versus prescribed-coordinate reaction

Add explicit allowed/reaction force diagnostics while keeping the legacy
saddle raw-active fields and all force/optimizer equations. Reparse9native
terminal frames: Mpeak rawmax2.0651,excluded2.0631,allowedmax.0900,physical
tangent-.1243; ordinary residual.099214 unchanged byte-for-byte by diagnostics.
No Hessian/index/energy-error claim or stricter NEB target. Support reactions
are excluded from internal stationarity, not from external strain work.

39local/39actualarchive tests pass, native analysis byte-identical. Full local
suite1842pass/1optional ase.io.abacus unavailable skip/312deprecation warnings,
0failure,363.45s. Standalone force/conjugacy checks exit0 on synthetic models;
no material CI/DFT/HFpytest.16-page paper compile and2/14/15/16numeric/layout
inspection. Evidence/clamped_stationarity_E073_20261011. Runtime edits are
local; immutable live HFsource/physical contract/G3 budget/holdout unchanged.
11:47two32CPU alive:T5/.103118,flip6/.196321. Read-only native priority/dependency
audit identifies downstream capacity risk after replacement jobs: four+1%
starters need a separate held-slot guard before releasing held reversing job.
No scheduler write belongs to this diagnostic milestone. PROGRESS is tested
core interpretation and actual cached gradient evidence, not completion.

## E074 — replacement-segment capacity guard before next admission

The old E060 afterok DAG does not count preserving continuation28723655.
Actual exclusive11:54–11:55transaction holds exactly four pending+1%starters,
one native write/readback each, preserving all original dependencies/32CPU/
source/inputs. Two running study allocations and held zero-strain reversing
job are untouched;0new submissions/releases/DFT/dependency changes/holdout.
Next admission must verify terminal material and an actual spare slot before
releasing one registered waiter, never all four. No cap breach occurred.
Offline safety and separate real native-receipt tests distinguish preparation
from execution. Evidence/clamped_capacity_E074_20261011; goal remains full
HfO2G2/G3/strong-model/independent-forecast active, not scheduler-only.

## E075 — correct new-starter step-zero cost before a future admission

Audit now requires the registered recipe/runtime/cache policy to distinguish
7fresh interior SCFs at a starter's step0 from0for an exact nine-cache resume.
Old M91calls and receipts unchanged; extra/missing calls still reject. Real
native Tpreflight validates this policy.144local/144actual full-archive related
tests directexit0;4fresh deployed analysis/admission source hashes plus native
preflight match. No production/input modification, DFT, terminal audit,
submission/release/dependency write/G3selection/holdout. Prepared one-slot
existing reversing admission remains explicitly unexecuted.
12:06:46native bothT28692775/flip28723655RUNNING;Tstep6rebound.143524after
.103118does not stop the chain. All5waiters held. PROGRESS is corrected actual
workflow cost/safety, not another completed material label or scope reduction.
Full finite HfO2G2/G3/strong-controls/independent-prediction goal remains active.

## E076 — resolve exact public Ma-v3 mapping evidence, not new DFT

Anonymous public v3 PDF actually downloaded/hash-identified. Selected main/SI
sections and renderedS2/S12 clarify competing pure M-PO mappings, a prose/figure
ordering conflict and the vacancy-only SCAN comparison. Not a full48page read
or publisher-SI identity claim. Introduction now separates prior mapping work
from our candidate identity/coverage, without adding a second M chain, changing
fixed inputs or G3 selection. Strong T/Cmma controls and prospective prediction
remain required.17-page existing LuaLaTeX compile directexit0; actual pages
1/2/3/16/17 inspected and retained M numbers checked in Markdown/PDF. Initial
preview clipping is not a source-PDF layout defect; wrapper unchanged.

12:23:32two native32CPU jobs live:Tstep6/.143524,flipstep9/.194546; five waiters
remainheld.0scheduler writes/newDFT/holdout/model freeze. Evidence is the
public-source report/provenance and manuscript Ma-mapping validation receipt.
PROGRESS is a resolved primary-source gap and manuscript boundary, not a
completed material response or reduced scientific objective. Goal active.

## E077 — native substrate work and the held-internal-coordinate distinction

Reparse9portable native M frames, reuse existing biaxial chart/raw stress and
actual cell Jacobian. Peak fixed-open-row partial-599.3443meV/fu/unit strain,
affine-carried-row-594.0228; difference is measured open-row work, not an
uncertainty or new mechanical ensemble. Relative-PO-370.9297 is not a relaxed
barrier slope or blind forecast. Stationary released-branch agreement and
full finite-strain G2/strong-controls remain unproved.0DFT/scheduler writes,
physical/core/optimizer or G3/holdout changes. Current M receipt admitted
before label access; unknown receipts reject first.

40final local/40actual archive tests directexit0, nine-frame material fields
exactly match. Sole source-hash metadata difference is existing biaxial
CRLF/LF; normalized text identical. Static units skill audit0findings is
heuristic, not native DFT energy/stress consistency or an uncertainty bound.
17-page manuscript directcompile0, final16/17 body/layout and PDF header
coordinates checked. No all-page final submission audit. Four own preliminary
outputs moved recoverably to E:/TEMP, final bundle kept compact.12:37:01both
32CPU live:T7/.174570,flip10/.193797; five waiters held. PROGRESS is actual
raw-stress evidence, tested reusable diagnostic and main-text integration,
not complete G2/G3/H1/H2/prediction/JCTC readiness. Full goal remains active.

## E078 — actual dual-gradient representation keeps Cmma a serious control

Reparse the same9native clamped M frames and reuse all E023 reference choices.
Correct virtual-work pullback, original ordered G2 frame and explicit
translation-free coordinate restriction preserve ranks3/4, all four Cmma
frames and raw translation diagnostics without source ASR/force repair.
At sampled peak3, retained atomic-gradient squared fractions are9.0282%/
32.7096%/93.5085--93.5087% for T-triplet/T-optical/Cmma. Low Cmma displacement
coverage36.48% is not local-gradient failure or a predictive-accuracy result.
Omitted movable modes are not support reactions; open-cell work stays separate.
No energy partition, Hessian, conditional branch, G3 selection or blind gain
is claimed. No new DFT/physical parameter/production/optimizer/holdout change.

68local/68actual full-Git-archive related tests directexit0; both actual native
analysis reports byte-identical. Static units audit0findings is heuristic;
finite differences use an analytic pair energy, not extra native DFT.
17-page manuscript compile direct0, rendered15/16/17 inspected; no final
all-page submission audit. Own preliminary outputs retained in E:/TEMP.
13:01:37native two32CPU jobs live:T9/.144511, preservingflip13/.191214,
five waiters held and untouched. PROGRESS is tested reusable dual projection,
actual strong-reference evidence and main-text discussion, not completion
of G2/G3/H1/H2/independent forecasting. Full scientific goal remains active.

## E079 — actual PO-to-T convergence and one existing reversing admission

Native28692775 COMPLETED0:0at13:11:49, step10/.0935897657 ordinary residual.
The complete77-call audit includes seven initial fresh interiors, original
six physical bytes/native32rank/fullEFS and338immutable runtime files.
Nine portable terminal native frames replay the same clamped ensemble.
Forward/reverse peaks84.416624/72.056649, well difference12.359975meV/fu;
the14.041498difference below the M discrete peak is not certified continuous
escape ordering. Ordinary convergence is not a TS/error certificate. G2
now has2/8 ordinary-converged training bands, not a completed network.

One accepted release of existing28692776 at13:21:40, no sbatch/dependency
write/new independent chain. Separate native proof:RUNNING32CPU/node6from
13:21:41 and first fresh SCF DSIZE32; two allocations, four+1% user holds.
0physical/source/optimizer/G3/holdout changes. A mistaken archive path in
the first read-only audit invocation exited1 before output creation; the
corrected native audit exited0. Admission exits0 and must never be replayed.
Evidence, tests and manuscript checks are in clamped_T_terminal_E079_20261011.
This is material PROGRESS, not a wait-only commit or full goal completion.

## E080 — paired escape channels separate geometry and energy-gradient coverage

Replay9actual T native frames with the unchanged E023 references, original
G2 frame/dual pullback and all4Cmma registrations; reuse M scalar reports.
At T sampledpeak4 geometrictriplet displacement99.0957% versus gradient14.3013%;
Cmma gradient54.0718--54.0765%, versus M93.5085--93.5087%. This is descriptive
channel dependence, not energy partition/model accuracy, TS certification,
G3selection or independent predictive advantage. Same cached PO+ native
inputs/log/EFS and ordered initial geometry; last-digit POSCAR serialization
and boundary-source provenance differences recorded without repair.

63local/63actual clean-archive related tests direct0/noerrors/failures/skips;
actual clean gradient report byte-identical. Displacement differs only in
existing reference-script CRLF/LF byte hash; all other fields/normalized code
identical. Tested sources/data byte-identical. Static units0findings heuristic.
Six-panel PDF/SVG/PNG/TIFF+CSV, undefined T-displacement endpoint omitted,
normal labels/aligned axes/framed external legends tested and visuallychecked.
Main-text paired table/discussion/Figure5 integrated,18-page compile0 and
actual15--18pages inspected. Preliminary own overlapping figure retained
recoverably outsideGit. No final all-page/submission-ready claim.

13:48:16native both32CPU live:preserving16/.188172,reversing0/.221799;
four+1%startersheld.0schedulerwrites/newDFT/physical/runtime/optimizer/G3/
holdoutchange. PROGRESS is actual paired evidence, tested reusable plot,
case and main-text integration. Full finite HfO2scientific goal staysactive.
Evidence:clamped_T_reference_E080_20261011/validation.json.

## E081 — use the exact clamped chart in conditional calculator evaluations

The historical six-symmetric-strain evaluator cannot impose the actual G2
two-vector substrate/open-row boundary. New public ClampedModeCoordinates
reuses the unchanged active geometry/physical gradient and existing ASE
factory/raw validator, with explicit Angstrom-scaled cell coordinates and
unwrapped common-reference registration. It does not release substrate strain,
select modes or replace the old chart/resume contract. Raw stress is checked
before ASE can average an asymmetric tensor into Voigt form.

232 independent EMT work evaluations and 18 actual portable T/M frames verify
geometry, units, cache recovery and native-gradient replay. No new native DFT
finite differences, material surface/Hessian/TS or prediction is claimed.
70 related local/70 actual clean-archive tests pass; full local regression is
1951 pass/1 optional ase.io.abacus skip, direct exit0. Five changed Python
sources are byte-identical in the tested archive. Clean verifier values match;
only existing active_curvature CRLF/LF source-hash metadata differs, with
normalized source text identical. Static units audit0findings is heuristic.
The existing whole-grid interpolation gate is not loosened for the eight-point
G3 validation budget. User changes and untracked work remain untouched.

14:13:17 native32CPU jobs still RUNNING: preserving18/.185892 and reversing
2/.171046 latest complete rounded-log residuals; exact stderr files0bytes.
Four+1%startersheld, 0scheduler writes/newDFT/production/G3/holdout changes.
G2 remains2/8 ordinary-converged candidates. PROGRESS is a necessary tested
calculation interface, not material G3 or complete H1/H2/JCTC readiness.
Full scientific goal staysactive; evidence:clamped_mode_flow_E081_20261011.

## E082 — partial Lee SI text clarifies the novelty and response units

Publisher extraction supplies four-page Parts1--3 and caption text, but no
successful source PDF binary or graphics inspection. The short attributed
reading and explicit access failures are retained; coefficients, vectors and
plots are not imported. Restricted release lowering a coupled interpolation
is prior work, not our method innovation. Internal length-valued mixed response
and elastic stiffness have different dimensions; the registered formula and
same-boundary material tests, not the acoustic argument, must supply evidence.

No code or numerical tables/figures changed. Related127tests pass without
failure/error/skip; static two-file unit scan0findings is only heuristic.
The enclosing initial check failed on GBK decoding after successful pytest;
corrected explicitUTF8/final PDF checks exit0. Actual final compile exits0,
18-page contact layout plus changed pages2/3 reviewed. Existing34table rows
and5figure references equal HEAD. This is not an external-figure/full-source
version audit, fresh full regression, material surface or submission-ready claim.

14:49:31native two32CPU jobsRUNNING; subsequently read actual complete logs:
preserving21/.182030,reversing4/.113133. Four+1%starters JobHeldUser;
0schedulerwrites/newsubmissions/physical or production-source changes.
Full finite G2/G3/strong-controls/frozen-independent-test goal remainsactive.
Evidence:lee_SI_text_E082_20261011/validation.json.

## E083 — third matched-substrate ordinary band and actual +1% startup

Reversing28692776 COMPLETED0:0 at14:51:28CST, step5/.09218242631.
Native audit checks42fresh SCFs including seven step0 calls, original six
physical hashes, DSIZE32/fullEFS and338 immutable runtime files. Nine portable
terminal frames separately match native input/log hashes, ordered geometry,
EFS and ordinary replay. Forward/reverse/well discrete energies338.968965/
339.004224/-.035259meV/HfO2; image4 highest. SCF transport47.667216corehours
excludes endpoint preparation. This is3/8 G2 ordinary candidates, not a
continuousTS, sampling-error certificate, complete network or forecast gain.

117related local/actual independent full-archive tests direct0, nofailures/
errors/skips. Two new source files differ only CRLF/LF in the initial archive,
normalized code equal; native helper byte-matches local tested source.
Initial source import diagnostic selected an old exporter through cwd;
corrected explicit import proves three registeredE075 modules before the
successful native audit. Initial archive byte-equality diagnostic fails on
line endings, not behavior; later raw-material/replay checks exit0.

Actual one-release journal admits existing+1%PO-to-T28709788, no sbatch or
dependency rewrite, no new independent chain. Native RUNNING32CPU/node6 from
15:06:04;15:08:11 first SCF header DSIZE32/stderr0. Three others remainheld,
28723655 occupies the other slot. Never replay accepted/ambiguous intent.
Main text/Table5 updated, final18-page compile0, pages1/15/16 individually
and late14--18 contact reviewed. Shortened prose keepsTable6 caption/data
together, original34rows/5figures unchanged except one added material row.
0physical/production-source/optimizer/G3/holdout changes; full scientific
goal remainsactive. Evidence:clamped_reversing_terminal_E083_20261011.

Archive hygiene correction: force-adding the otherwise ignored native logs/
trajectory also included one generated local pyc. Its tracked entry is removed
in a normal follow-up commit, with local cache retained. All71 final portable
case files byte-match the actual evidence; tested code, physical data, native
journal/handles and scientific conclusions are unchanged. No history rewrite.

## E084 — finite independent-probe coverage without loosening whole-grid gates

The E081 finite-reporting gap is closed: registration binds training-only
geometry/metric/basis/grid/evidence, probe positions, declared starts/rule,
evaluation IDs and thresholds. It reports signed energy/metric coordinate
errors, failed calculations and pending points, with unchanged per-patch
and whole-grid screens.8 passing points in8 of64cells never become a whole
surface or continuous-region bound. Refits, declared reuse and unregistered
observations are rejected. Declarations/hashes do not prove source independence
or prospective blinding; publish the exact manifest before real labels.

Synthetic verifier0calculator calls. Related local/actual Git-archive90tests
each, adjacent88, full local2021passes/one optionalase.io.abacus skip, all
direct0 and final full process exit observed. Five reviewed Python sources
have equal normalized text; fourCRLF/LF raw hashes differ, not a code mismatch.
Registration/both numerical reports byte-identical; verifier receipt differs
only source hashes. Two-file static unit audit0findings is heuristic.

15:44:11native28723655 step25/.175897 and28709788 step1/.224071, bothRUNNING
32CPU/stderr0; three+1%startersheld.0scheduler writes/newDFT/production changes.
Own+1%firstupdate16:54 suggests10update segment18:10--18:40,wall19:06;
preserving40update segmentaround19:00,neither a convergence guarantee.
G2 remains3/8; noG3 selection/DFT/Hessian/+.005 labels or paper/PDF changes.
PROGRESS is a required finite-budget reporting implementation, not H1/H2
material validation or an innovation claim. Full scientific goal staysactive.
Evidence:conditional_probe_coverage_E084_20261011/validation.json.

## E085 — completed reversing band's unchanged reference and frozen-work traces

Existing native terminal data supply the third channel's common-reference
description, without selecting or rotating a reference from its result.
At peak4, Ttriplet squared displacement/gradient fractions26.2915%/.0320%
contrast with Cmma46.5637--46.5639%/65.1046%. Absolute/omitted dual-gradient
norms and all four registrations remain. Fixed-open-row external work is
-27.18358meV/HfO2/unitstrain, not a relaxed-barrier derivative or TS test.
Three completed bands' cached initial native inputs/logs/EFS are identical;
4.44e-16A ordered geometry serialization difference is retained, not repaired.

61 related tests local/actual Git archive each pass/direct0. Initial exact
floating geometry assertion1fail/60pass is retained with corrected explicit
serialization bound/exact SCF checks. Eight reviewed Python sources have
equal normalized text; two original sources differ CRLF/LF, so two report
source-hash fields differ, while all numerical fields and gradient report agree.
Actual final compile0,19pages; final contact sheet and16/17/18 individually
reviewed. Original35table lines/fivefigures retained, three actual Table6 rows
added. Table caption/no-break hook and shortened Figure5 caption preserve
scientific interpretation without orphaned captions. Initial table guard had
an overbroad regex, corrected and verified; no physical calculation failure.

16:10:01both native32CPU allocations remainRUNNING: preserving27/.172483,
+1%PO-to-T2/.156615; three other training startersheld. Wrong guessed stderr
paths were rejected and resolved against actual read-only Slurm job fields.
No scheduler writes/newDFT/production physical or optimizer changes/G3 choice/
holdout access. G2 remains3/8, full scientific goalactive. PROGRESS is new
material description and reproducibility, not an independent-prediction gain.
Evidence:clamped_reversing_reference_E085_20261011/validation.json.

## E086 — tested next +1% rolling handoff, no premature admission

The E083 one-use release assumes four original held starters; after its
accepted T release the real queue has three. A separate private M handoff
now reuses the pinned read-only seed/runtime checks but requires an actual
ordinary-converged +1% T terminal receipt and native spare slot before one
M release. Both flipping starters stay held with their original M dependency.
Exclusive durable intent blocks accepted/ambiguous replay; no sbatch or
dependency updates. Pending-input preflight is separate from admission and
does not infer convergence while T is live.

172related local/actual Git-archive tests each pass/direct0, including46new
mock transaction tests. Two new sources differ CRLF/LF only, normalized
text matches; deployed helper byte-matches locally tested source. Actual HF
read-only preflight0 checks9images/two native cached endpoint SCFs, original
six bytes/DSIZE32/fullEFS,338 runtime files/binary/pilot and three held jobs.
Local receipt byte-matches native. No native failure injection or whole-project
regression claimed. No admission journal, DFT, release or dependency write.

16:31:39native T28709788 remainsRUNNING step4/.111567, preserving28723655
step28/.170654; both32CPU/stderr0. G2 remains3/8, with three+1%startersheld.
Finite segment ETAs18:10--18:40 and around19:10, not convergence promises.
PROGRESS is required operational preparation, not new material/innovation/
prediction evidence. No source/physical changes, G3 selection or+.005 access;
paper/PDF unchanged and full scientific goalactive. Actual terminal audit
must precede any reviewed once-only M admission.
Evidence:clamped_plus1_T_handoff_E086_20261011/validation.json.

## E087 — fourth actual G2 terminal, first matched-strain T peak response

28709788ends16:43:06COMPLETED0:0; lateststep5replays.0949479822ordinary residual.
Actual HF audit42freshSCFs/338runtime files/originalsix bytes/nativeDSIZE32/EFS
passes/direct0. Nine portable terminal frames reparse independently; not every
earlier raw SCF or licensed orbital is redistributed/re-audited locally.
Forward/reverse/well96.763650/82.617767/14.145883meV/HfO2, samplepeak4.
Relative to zero-strain T, discrete peak grows12.347027 while well separation
grows1.785908. This seen training contrast is not an independent prediction,
continuous TS/error certificate or full competing-network H1 conclusion.

Existing M28709789released once after actual audit/capacity/originalseed checks,
exclusive intent and second native recheck.0sbatch/0dependency writes/0new chains.
Do not replay accepted journal7cca3198. Separate startup16:52:43RUNNING32CPU/node3
from16:50:44, rawDSIZE32/stderr0. Registered production DFT now runs; audit/
analysis/admission itself made0extraSCFs, not a claim of zero DFT for the turn.

185related local/actual Git archive tests each pass/direct0;13actual-case tests.
First184pass/1fail local replay used Windows backslashes for native Linux path;
retained failed JUnit, corrected test uses original PurePosixPath semantics,
native receipt/executed source unchanged. Initial inline SSH command was
rejected by PowerShell parser before remote execution; file invocation passed.
71initial case files byte-match archive; force reports byte-identical. New test
normalized text matches, rawCRLF/LF differs. No full-project/HF pytest claim.

19-page manuscript compiled/direct0 and contact+1/15/16/17/18 individually
inspected. All38oldtable lines/fivefigures retained; one actualTable5row added.
At17:05:23preserving31/.164759 andM0/.983175 stillRUNNING32CPU/stderr0;
two+1%flipsheld. G2now4/8. SegmentETAs18:35--19:10 and19:00--19:40 are not
convergence guarantees. No G3 selection/extraGamma/Hessian/+.005access or
physical/production changes. Full scientific goalactive; PROGRESS is real
same-boundary training response and authorized operational handoff, not H2gain.
Evidence:clamped_plus1_T_terminal_E087_20261011/validation.json.

## E088 — shared visible initial response, no extra initial Hessian

Implements the existing optional-feature route: one immutable, matched
P=0 initial energy increment and recorded preparation cost with any controlled
bottleneck arm. Keep anchor correction, separate scales/domains, signed gap
and omission/reaction diagnostics. No repeated correction of measured initial
energy or derivatives/curvature inferred from two endpoint energies. No
material feature-choice/forecast freeze or physical initial Hessian request.

39new tests; local198 then203related pass. First actual clean full run2102pass/
1fail/2skip; previous HEAD4d6e05d archive reproduces the same old T/M registry
assertion. E085 admitted the third pinned reversing receipt but the old test
still required only two. Exact three-hash equality plus actual receipt byte
check fixes the test, not production policy. Final independent Git archive
2103pass/2skip/314existingwarnings/direct0. Retain both failures, all source
raw/normalized comparisons and final JUnit. The two skips are optional ABACUS
reader and intentionally absent licensed POTCAR, not backend passes.

NumPy-only12point direct replay completes locally/HF/direct0, max identity and
cross-environment energy difference2.78e-17; seven source/checker byte hashes
match. HF has no pytest; no installation/native-pytest success claim. Known
synthetic initial Hessian supplies independent truth only. Numeric output
protection actually rejects existing file/exit1 with original bytes preserved.
These are implementation residuals, not material accuracy or independent H2gain.

17:41:17two study allocations stillRUNNING32CPU/stderr0: preserving34/.158213,
plus1M3/.849185. Two flipping starters held; G2remains4/8. SegmentETAs18:50--19:05
and18:55--19:15 from actual recent timings, not convergence guarantees.
0newphysicalDFT in validation/analysis, while registered production continues
normalDFT;0schedulerwrites/admissionretry/extraHessian/G3choice/+.005access.
Material MD/PDF and five user-dirty tracked changes preserved. Full goalactive.
Evidence:shared_initial_response_E088_20261011/validation.json.

## E089 — unchanged reference spaces under seen +1% T strain

Nine actual +1% T terminalSCFs reparse/direct0 in original displacement,
dual-gradient and raw-stress analyzers. Original T/Cmma masses/chart/ranks,
all four registrations, integer gauge and electronic input hashes remain.
Peak triplet displacement99.0957->99.0045%, gradient14.3013->69.7245%; Cmma
gradient54.0718--54.0765->10.5750--10.5766%. Absolute internal norms.010196606
and.038221369 remain explicit; percentages are not energetic accuracy or TS
certificates. No choice based on partialG2; different relaxed PO caches under
the two strains are not falsely made identical. Each profile has its own
PO zero and individual normalized arc. Seen peak increment12.347027meV/HfO2
is not an independent forecast, relaxed-barrier derivative or H1/H2 gain.

Figure6 contract predates code; all frames, undefined zero-T fraction and
defined strained-T endpoint retained. SourceCSV18records, no replicates/
confidence band/fitting. Actual21-page PDF/direct0/contact+19/20/21reviewed,
39oldtable lines/fivefigures and38local links retained, new main-text section.
Compiled build stays ignored by existing policy; editable figure/source stored.

Local60 and actual final archive60related tests eachpass/direct0. Actual first
clean full2114pass/2existing skips/314warnings/direct0,342.778s. Six changed
code/test sources exactly equal between archives. Final16case/figure files
byte-match archive; every export andCSVhash matches after scoped raw-byte
attributes. Three independent clean analyzers preserve all numerical fields;
gradient/CSVbyte-identical, two existing CRLF/LF source hashes explicitly differ
with normalized text equal. Initial strict comparison, missing filename,
unit-normal exact comparison, metadata and PDF/runtime/encoding corrections
are retained in receipt; none is a physical DFT failure or input adjustment.

18:11:43native preserving36/.153532, +1%M6/.725166 remainRUNNING32CPU/stderr0;
two +1%flipsheld. Finite segmentETA18:50--19:05 not convergence promise. G2still
4/8. No scheduler writes/admission replay/new chains/extraDFT/G3choice/+.005
access. Registered production continues normalDFT. Original five user-dirty
tracked changes untouched; full scientific goalactive. PROGRESS is descriptive
material strain evidence and its audited figure, not prospective superiority.
Evidence:clamped_T_strain_reference_E089_20261011/validation.json.

## E090 — bounded latest-chain continuation prepared, not admitted

One inert 40-update/12-hour/32-MPI wrapper accepts only the two registered
parents after an actual successful unconverged cap audit. Original E054/E055
runtimes, physical inputs, clamped boundary, spring/geometry guards and .10
remain unchanged. Full 280/77 parent SCFs, nine exact current caches, native
terminal identity and every original runtime file must pass; failed/live/
converged or different-parent evidence rejects. FIRE history resets explicitly,
not a full checkpoint or acceleration result. Future channel-named logs and
zero-initial-SCF/280-update-call policy are tested against the existing auditor.

23 new synthetic tests; local related 134 pass. Initial 18 fixture failures
come from monkeypatching generic Path's platform dispatch, not production.
Concrete platform Path fixes the fixture while leaving every production gate;
the original failed JUnit and corrected intermediate 132 pass are retained.
Actual clean Git full suite 2137 pass/2 existing skips/314 warnings/direct0,
352.014 seconds. Three changed source bytes exactly equal the full archive.
Actual HF bash -n/direct0 is syntax only; no new native pytest or DFT/admission
success is invented. A first SSH quoting probe fails before diagnostic execution;
the stdin repeat records native scheduler/log state and syntax correctly.

18:35:03 two allocations still RUNNING32CPU/stderr0: preserving38/.148593,
+1% M8/.645983. Median last-five updates793/641s imply finite-cap ETA18:54:50/
18:52:59, not convergence estimates. Two original +1% flips remain held; G2still
4/8. Preparation contributes no material result, new independent chain, deploy,
scheduler write, accepted-journal replay, extra DFT/Gamma/Hessian/G3 choice or
+.005 access. Registered production continues normal DFT. MD/PDF and five
user-dirty tracked files remain untouched. Full scientific goal remains active.
Evidence:clamped_rolling_continue_E090_20261011/validation.json.

## E091 — actual terminal audits and two admitted same-chain continuations

Both registered parents finish successfully at caps, not ordinary convergence:
preserving28723655 step40/cumulative90/.1434371761; +1%M28709789 step10/.5701098698.
All280/77freshSCFs, unchanged six physical bytes, actualDSIZE32, completeE/F/stress,
SCF convergence, original334/338runtimefiles and9latestexactcaches pass onHF.
Portable18rawframes and both same-boundary residuals replay locally. Latest
sampledpeaks remain provisional; G2ordinary4/8, no complete channelranking,
TS certification, model advantage or frozen independent forecast is claimed.

Actual first preserving audit fails because frozenE054 seed preparation admits
only the pilot script, not registeredE069. Original code/partial namespace and
actual rejection preserved. Existing tested continuation-aware preparation is
pinned in audit-only resume_support, not the worker source. Full native checks
rerun in fresh terminal-r2 before continuation. M audit used retained v1 driver.
Two exclusive once-only journals admit then release28727368/28727382: 2held
submissions+2releases,0dependency/old-waiterwrites,0newindependentchains. Same
latest geometry/caches, original electronic and clamped boundaries,40updates/
12hours/32MPI each. FIRE history resets explicitly, not a checkpoint/gain claim.

15newtests; preparation clean v1/v2 each95pass; local material suite75pass.
Final actual independent archive98pass/0fail/0skip/direct0,11.009seconds.
161case artifacts and2newtests+5existing support/wrapper files byte-identical.
Two real delivery attempts97pass/1fail each reject line-ending-converted old
source bytes. Exact snapshot -text/-eol plus scoped index renormalization fix
delivery; original digest and test unchanged, both failedJUnit retained.
No new full-suite/nativepytest result is invented; E0902137pass is historical.

19:40:32CST bothchildrenRUNNING32CPU/stderr0/no failure. Eachnewstep1:
M.552680 and preserving.142522. Actual first updates861/820s plusparenttimings
give bounded40-updatecapETAOctober12about03--06CST, not convergenceETA; wallends
07:15:33/07:21:30. Original +1%flips28709790/91remainJobHeldUser/nullDependency.
No replay/extraHessian/G3choice/+.005holdoutaccess,CI/release/installation or new
independentchain. Audits/tests add0DFT, admitted production runs normalDFT.
MD/PDF and original5userdirtytrackedfiles unchanged. Scientificgoalactive.
PROGRESS is actual validated cap material and same-chain production admission,
not final convergence or JCTC completion.
Evidence:clamped_rolling_native_E091_20261011/validation.json.

## E092 — actual signed path-work check, not mode energy partition

Reusable calculator-independent physical-covector trapezoids and a strict
four-completed-G2-band adapter preserve the declared atomic mass chart,
dual pullback, explicit translation gauge and original permitted open-cell
space. Keep rank3T triplet, rank4completeT doublets and all4Cmma registrations.
Signed retained/omitted/translation/cell work adds to the full contraction;
raw native force/stress contraction independently agrees within1.12e-16eV/cell.
No source-mode/forceASR, fitted frame, geometry/input change, reference selection
or spring-force integration. Author vectors remain external and byte-pinned.

Actual36portableframe records from0M/T/reversing and+1%T contain34distinct
rawlog/STRU identities; common initial caches are not new independent calls.
Energy-minus-work at highest samples is20.897743/8.460687/35.275323/9.418257
meV/HfO2, while endpoint closure is3.006016/.728018/.115373/1.997185. Coarse
endpoint closure can hide interior discrepancy. All strong-reference signed
work and omitted contributions remain; no model winner, separate mode energy,
physical/barrier error bound, conditional branch or certifiedTS is inferred.
Coarse nonlinearity and native force/stress error are not separated. This is
real training-data consistency evidence for the required conditional/strong
control validation, not a new harmonic theorem or prospective advantage.

35newtests; initiallocal63pass, finallocal67pass. Actual clean full2187pass/
2existing skips/314warnings/direct0,361.051JUnitseconds. Four material analyses
execute independently from that archive and each scalar report is byte-identical.
Eight initialcasefiles+fournewcode/testsources equal the full archive. A later
combined metadata patch fails context before writes; actual absence checked,
then corrected. No DFT or scientific-data failure, no weakened tests. Source/
rawmaterial/test logic remains unchanged after full regression.

20:02:30CST jobs28727368/28727382 remainRUNNING32CPU/stderr0/no failure; newstep2
forces.525706/.140697. Originaltwo+1%flipstartersheld/nullDependencies;G2still4/8.
0schedulerwrites,admissionreplays,newchains or analysisDFT. Liveworkers continue
normalDFT. NoG3choice/extraHessian/+.005physicalholdoutaccess,CI/release/install/
automation/submission. MD/PDF andoriginal5userdirtytrackedfiles untouched.
Scientificgoalactive. PROGRESS is tested execution plus actual36frame work
diagnostic; incomplete G2, conditional branches and frozen independent tests
still prevent a final JCTC claim.
Evidence:clamped_path_work_E092_20261011/validation.json.

## E093 — pinned compatible child audit support, not a new physical result

The E075 exporter predates the accepted E090 continuation wrapper. Bind the
current auditor/exporter/preparer/residual together in a new audit-only HF
namespace; preserve both original E055/E054 production runtimes. Only actual
children28727368/28727382 are admitted, with their exact native submitted
parent receipt, seed, source, work/log namespace and wrapper. No native RUNNING
row is relabeled COMPLETED. Preflight and terminal claims remain separate.

Actual native compatibility checks both pass, importing all4 pinned support
modules and checking338/334 original runtime files against their archives.
Actual terminal attempts on the two live jobs intentionally exit1 before
material reads/imports or terminal output creation. These are negative guard
tests, not production failures. Terminal mode still delegates all SCF/EFS/six
physical bytes, cache/geometry and ordinary-force checks to the existing auditor.
Only a successful finite cap may prepare a seed; this never submits a new job.

39new tests; first local135pass, final137pass. Actual independent Git archive
137pass/0fail/0skip/direct0,13.998JUnitseconds. Driver/test/README/two priorJUnit
files equal archived bytes. Initial root.gitattributes byte comparison rejects
LF/CRLF normalization; actual decoded content is identical and all five code/
case/test files subsequently compare byte-equal. No scientific digest or test
is weakened. No new full-suite/HFpytest result is inferred; E0922187pass remains
historical. Only reports/validation/ledger are added after archive tests.

20:29:50CST bothjobsRUNNING32CPU/stderr0/no failure, eachnewstep4:
+1%M.470113 and zero-strain preserving.135313. Measured recent update times
M861--1127s, preserving803--820s imply nominal40-update caps around04:18--04:28
for preserving and04:58--07:38 forM onOctober12. M's wall limit07:15:33 can
arrive first; this is not a convergence ETA. Held+1%flipstarters remain
JobHeldUser/nullDependency. G2still4/8; noG3selection/extraHessian/+.005holdout.

0schedulerwrites/admissionreplays/newchains/auditDFT; registered workers
continue normalDFT. NoCI/release/install/automation/submission. Original5user
dirty tracked files and manuscript remain untouched. Goalactive/incomplete.
PROGRESS here is bounded, actually verified audit engineering needed for
forthcoming material results, not a new scientific prediction or JCTC completion.
Evidence:clamped_child_audit_E093_20261011/validation.json.

## E094 — existing four-band work evidence written into the main manuscript

Sections2.5/3.8 andTable7 connect physical-covector trapezoids to the actual
four completed training bands. Endpoint closure .115373 at the reversing
flip can conceal peak discrepancy35.275323meV/HfO2. Retain all signed
retained/omitted/translation/cell work and every originalT/Cmma representation.
No state-function mode energies, physical error bound or forecast winner is
inferred. No extra sampling/Hessian budget is opened by this diagnostic.

The initial combined patch context failure is verified before writes using
the unchanged original manuscript hash, then reapplied separately. A manual
hyphen/linebreak artifact is corrected before the final build/test archive.
Duplicate3.6 headings become3.6/3.7. Actual final compileexit0,22pages/
1049375bytes. All22pages rendered/reviewed at contact-sheet scale; changed
pages5/21/22 inspected separately130dpi. All16Table7 values occur in actual
PDF text, with no final overfull/underfull/missing-glyph/LaTeX/package warnings.

11newtests; actuallocal46pass/0fail/0skip,2.174JUnitseconds. IndependentGit
archive46pass/0fail/0skip/direct0,2.225seconds. Main source/newtest/firstJUnit
byte-equal; README decoded-content equal before validation-link-only update.
No material report, figure or runtime source change, no new full-suite/native
pytest result invented. Existing strong-reference ranks and frames stayfixed.

20:47:49CST jobs28727368/28727382 RUNNING32CPU/stderr0/no failure, newsteps5/6
and forces.441986/.129644. Both original+1%flipstarters remainJobHeldUser.
G2still4/8; conditional/strong-control and frozen independent prediction
evidence remain missing. 0schedulerwrites/replays/newDFT/chains/Hessians,
G3selection or+.005holdoutaccess; productionworkers continue normalDFT.
NoCI/release/install/automation/submission. Original5userdirtytrackedfiles
preserved. Goalactive;PROGRESS is verified substantive main-text delivery,
not final JCTC scientific completion.
Evidence:paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.path-work-E094.validation.json.

## E095 — nested controls reuse one visible initial-well feature

The existing nested gap route required both well and bottleneck Hessians,
whereas the existing visible-initial-energy route accepted only one control
at a time. Connect these actual kernels: all B2--B5 arms receive the same
immutable measured well increment, bottleneck columns and declared domain.
No initial Hessian, minimum index, derivative or additional anchor correction
is fabricated. Unavailable arms retain their original assumption/measurement
reason, feature and cost reference. Full source-ID cost deduplication remains
external; one shared endpoint record is not incurred anew for four models.
This is standard-kernel integration, not a material accuracy/novelty claim.

27newtests use separately solved original physical blocks, not the returned
partition. First local178pass/1fail is a wording assertion about the inherited
missing-curvature reason, corrected to exact reason and existing status.
Expandedlocal215pass; independentarchive215pass/0fail/0skip. Actual full
archive2264pass/0fail/2existing skips/314warnings/direct0,355.117JUnitseconds.
Core sources, nested docs and firstJUnit are byte-equal; test, shared docs,
initialREADME and snapshot undergo Git LF/CRLF normalization with equal
decoded content and both hashes recorded. No scientific hash is weakened.
One PowerShell comparison parser rejection occurs before execution; corrected
comparison succeeds. An initial append guard incorrectly assumes Tail3
returns one JSON record; actual E094 last record is then verified and kept.
After archive tests only delivery/verification metadata and generated JUnit
records change, not model/test logic or physical data.

21:24:54CST both children remainRUNNING32CPU/stderr0/no failure/summary:
+1%Mnewstep7/.386236 and zero-strain preservingnewstep9/.120476eV/Angstrom.
Two original+1%flipstarters remainJobHeldUser;G2ordinary still4/8. Runtime
parameters stay fixed. 0schedulerwrites/replays/implementationDFT/newchains/
Hessians; registeredworkers continue normalDFT. NoG3choice/physicalholdout,
CI/install/release/automation/submission. Maintext/materialreports and the
original5userdirtytrackedfiles remain unchanged. Goalactive/incomplete.

PROGRESS is the actually tested integration needed to avoid unnecessary
initial curvature preparation if the optional feature is preregistered.
Material response columns, all strong references and complete prospective
forecast freeze are still missing; no independent path is released.
Evidence:nested_shared_initial_E095_20261011/validation.json.

## E099 — three existing matched-substrate candidates in the main figure

E096/E097/E098 are verified waits on the same live native children, with
no empty code/ledger updates. Current Figure5 still lacks the already audited
clamped reversing profile. Extend it using the original scalar records and
unchanged pinned T/M plotting adapters; preserve the older dated figure.
One PO+ native zero and one energy scale compare T/M/reversing peaks
84.416624/98.458122/338.968965meV/HfO2. All original T ranks3/4 and four Cmma
rank4 registrations remain; displacement and gradient coverage are distinct.
Individual normalized arcs are not common absolute coordinates. Nine discrete
points/chain with straight connectors are not a validated smooth MEP orTS.

27records include25distinct manifestlog/STRU identities and a reused PO+
cache; none is claimed new or independently reaudited. No source eigenvector
array is redistributed. New183x207mm PDF/SVG/PNG/600dpi TIFF and scalar CSV
have source hashes and layout tests. Nonbold labels(a)--(i), all four spines,
aligned grids, outside semi-transparent framed legends and missing undefined
fractions honor user styling. Updated Figure5/caption stay descriptive.

11newtests; firstlocal27pass, expandedlocal50pass. Actual independent Git
archive50pass/0fail/0skip/direct0,4.475JUnitseconds; standalone plot CLIexit0
and complete scalar CSVbyte-identical. All14initialcase/figure/sourcefiles
compare byte-equal. Root attributes/paperREADME decode identically despite
Git newline normalization. After archive testing only verification metadata,
ledger and copied JUnit change, not tested source, maintext or plotteddata.
No new full-suite or HFpytest claim; E0952264pass remains historical.

Actual manuscriptcompileexit0,22pages/1054618bytes. All22pages rendered and
contact-reviewed; new figure/captionpage18 separately inspected130dpi.
No final overfull/underfull/missing-glyph/LaTeX/package warnings. Optional
fitz probe fails before graphics; existing Python PDFium performs QA, with
no dependency install or backend switch. Figure PNG also actually inspected.
Compile/layout checks do not prove the remaining material/prediction claims.

22:01:39CST children28727368/28727382 remainRUNNING32CPU/stderr0/no failure:
newsteps10/11 and forces.306576/.113740eV/Angstrom. Originaltwo+1%flipstarters
remainJobHeldUser. G2ordinary still4/8;0schedulerwrites/replays/plotDFT/
newchains/Hessians. Existingworkers continue normalDFT unchanged. NoG3choice,
physicalholdout,CI/install/release/automation/submission. Five original dirty
user files and all old material/figure sources preserved. Goalactive and
incomplete. PROGRESS is an actually tested figure/manuscript delivery using
already admitted physical evidence, not a novel-method or prediction win.
Evidence:clamped_three_band_figure_E099_20261011/validation.json.

## E100 — ordinary preserving result and one existing +1% rolling handoff

Zero-strain preserving child28727382 finishes22:39:34CST/COMPLETED0:0/32CPU.
Actual pinned E093 native audit checks15newupdates/105fresh SCFs, exact seed,
334originalruntime files and sixphysicalbytes/fullEFS/nativeDSIZE32. Ordinary
same-boundary residual0.0988168085874911eV/Angstrom passes0.10. Nine portable
terminal raw frames and complete receipt are exported without licensed PP/
orbital or author vector arrays. Forward/reverse sampled peaks atimage4 are
96.828183003/96.863441806meV/HfO2, not a continuous TS/error certificate.
Segment transportcost11800.945525seconds/104.89729432rankcorehourproxy does
not include earlier continuations/endpoints or prove an optimizer speedup.

Ordinary G2 becomes5/8. The1.63meV/HfO2 sampled difference from M is not a
robust continuous-barrier ranking; all prior strong references/gates stay.
The original held+1%flip inputs are rechecked using pinned E068 seed logic,
E060preflight/338-file runtime and original pilot/binary. Each has9checked
geometries and2nativeendpoint caches; zero preflightDFT, not newreplicates.

New helper admits only28709790/28709791, never sbatch/dependency changes.
Exclusive per-target intent journals, one shared lock, committed-slot counting
including runnablepending, and native post-ack corroboration prevent duplicate
or racing releases. Both0strainflips must actually auditordinaryconverged;
live/capped/failed results and changed source/seed stop. Original read-onlyr1
helper is preserved; finalr2 is uploaded separately. No active source rewrite.
One existing28709790 release succeeds and separately starts22:44:35/node1/
32CPU. At22:45:55 its first SCF hasoneDSIZE32 header but is not complete.
Other28709791 staysheld; +1%M retains its firstslot. Noextra independentchain.

Firstlocal63pass/4fail finds Windows fixture quoting, corrected without a
native guard weakening. Expanded75pass, post-ack76pass, then actualmaterial/
admission80pass locally and in independentarchive/direct0/JUnit3.112seconds.
All74initial casefiles and testsource byte-equal; attributes LF-normalized
content equal. All9portableSCFs are locally reparsed. No freshfullsuite or
HFpytest claim; E0952264pass/2skip remains historical. Testedhelper/test/raw
material do not change afterarchive; only reports/status/README/ledger follow.

22:51:52CST two32CPU live handles28727368/28709790, stderr0/no reportedfailure;
Mstep13/.234106 and starter3completedcallaudits/notcomplete0step. One observer
stderrpath typo fails before snapshot and is corrected, without calculation
or scheduler mutation. Ordinary5/8, G3/bottleneckchoice/extraHessians/physical
holdout remain closed. Manuscript/figures and five original dirtyuserfiles
stayunchanged. Goalactive/incomplete. PROGRESS is a new actually audited
training candidate plus a tested existing-handle handoff, not a forecastwin.
Evidence:clamped_held_flips_E100_20261011/validation.json.

## E101 — admitted preserving reference and ordered flip comparison (PROGRESS)

Continue E100's actual terminal result, without a new DFT or scheduler write.
One fifth completed receipt is admitted to the strict training whitelist;
unknown/reserved labels still reject before access. Existing four analyzers
produce the zero-strain preserving displacement/dualgradient/pathwork/
externalwork reports. T ranks3/4 and all four rank4Cmma frames are unchanged.
At peak4, T3displacement82.6301% contrasts with gradient12.9583%; Cmma
gradient52.7203%. Internal/omitted T3norms.070196452/.065490597 have units
eV/(sqrt(amu)Angstrom), not mode energy. No ASR repair or reference winner.

Figure7 compares actual preserving96.828183/reversing338.968965meV/HfO2
profiles on common PO+ zero/energy scale.18records/17nativeidentities, not
replicates. All45workrecords have42identities. Preserving peak/endpoint
energy-minus-work14.523948/.032398 illustrates local discrepancy despite
endpoint closure, not a physical error bound. Tables5--7 and new3.9 add this
evidence;1.629939meV sampled separation from M is not a robust ranking.

Actual six-panel PNG and final24-page PDF contacts are inspected; changed
pages15/16/21--24 reviewed at130dpi. Number-range paragraph overflow is fixed
without changing values; finalTeX diagnosticcount0. Reading page22 has large
space before unbreakablefigure23, not final journal typesetting. PDF/source
pins and visual/figure CSV receipts are retained; no new curvature/G3/holdout.

Initial plot7pass/4fail is a canonical-vs-raw observation pin mistake,
corrected to committed raw bytes. Related93pass/1fail is the stale4receipt
test, updated to actual5receipt membership. First independent108pass/1fail
is only oldcurvatureLF/CRLF source metadata; exact numerical fields agree.
An intermediate reversed LFpin fails locally108/1 and is corrected with
both normalized-byte comparisons. Original reports/data/archives are kept.
Final109pass locally and independentr2/direct0/JUnit8.170s,1existingfontTools
warning each.27initialscience/figure/test/evidence files are byte-equal;
later changes are receipts/status/README/ledger only. No new fullsuite or
HFpytest claim; production/physical data byte hashes are not normalized.

At23:40:25CST M28727368step16/.168904 and preserving28709790step4/2.467425
are live on two32CPUhfacnormal01allocations; stderr0. Reversing28709791 stays
held/AllocCPUs0. G2ordinary5/8; no scheduler mutation, new DFT, CI, publication,
new phase, Hessian, conditional-face selection or physicalholdout. Five
original dirtyuserfiles are not staged or edited. Goalactive/incomplete.
PROGRESS is actual completed material analysis plus validated/main-text
evidence, not a status-only ledger commit or a prospective prediction claim.
Evidence:clamped_preserving_reference_E101_20261011/validation.json.

## E102 — same group, different sampled boundary role (PROGRESS)

Bounded zero-DFT comparison reuses actual ordinary free-cell preserving69
and clamped preserving15 observations, original exact readers, six physical
hashes, source/terminal receipts and fixed symmetry sweep. Both image4 centres
arePbcn at.001/.01/.05Angstrom/1degree. Own-well energies-14.838201/+96.828183
meV/HfO2 and centre-minus-adjacent-mean-47.644192/+11.480619 distinguish
sampled trough/high point. T-chartFxx1.122022/1.0 provides boundary geometry
evidence. No stable-intermediate/TS/all-clamped-branches-absent assertion.
Different mechanical boundaries, initial wells and original numeric cell
metrics are explicit; this is not the matched-boundary strain derivative.
Free38.508448sampling-union maximum remains distinct from nine-point32.806023.

New2panelFigure8/Section3.10 adds actual material mechanism evidence, not a
prediction-model win. Existing structure audit validates12orderedHf4O8,
finite/right-handed periodic cells/distances before symmetry, no standardize/
wrap/reorder/oxidation guess.18records, no independent replicate/newSCF.
First render wrongly assumes common cellscale and rejects; corrected by
retaining bothoriginalmetrics. One patch context rejects before writes.
Initial PDF caption split/plain underscore formatting corrected. Final25page
PDF compiles0, TeXdiagnosticcount0, allpagecontacts/final24/25detail inspected.
External Perl localewarnings remain. Figurecaption complete on24; old22blank
region remains readinglayoutlimitation, not journalready. Ddisk fills at r2
PNG creation; only own regenerable r1/r2 reviews move recoverably toE, finalr3
review runs onE. No unrelated source/science/user data deleted or moved.

Actual117focused tests locally/independentarchive pass,0fail/skip,1existing
fontTools warningeach.8science files byteequal to tested archive. No freshfull
suite/HFpytest claimed. G2still5/8; latest liveM28727368step18/.129055 and
preserving28709790step8/1.680009, reversing28709791held. No schedulerwrite,
DFT, G3/reference selection, curvature, physicalholdout, CI/publication or
physical/running-source change. Five original dirtyuserfiles remainunstaged.
Goalactive/incomplete. PROGRESS is measured material comparison plus tested
scientific code/maintext, not another status-only ledger commit.
Evidence:preserving_boundary_role_E102_20261012/validation.json.

## E103 — Audited same-chain +1% preserving continuation

Previous goal turn:E102 is PROGRESS (material mechanism evidence, tested
scientific code/maintext). Current finite question: can the original10-step
starter's exact latest chain continue without recomputation/input drift or
overcommitting the two study slots? This is production completion work, not
a claimed new acceleration method or another independent channel.

28709790 completes natively0:0 at00:34:32, step10/1.3257040697. HF audits all
77SCFs/338runtime files/six fixed physical bytes/DSIZE32/EFS; all9latest caches
pass. Latest geometry-only seed, not historical best; fresh FIRE history is
explicit. New40-update/12h wrapper preserves the E090 driver byte-for-byte.
Legacy flip_continue filename's historical log-prefix rule is caught during
review, avoided and covered by the future native-accounting regression.

147predeployment local+clean-archive checks pass; actual HF audit/cache pins,
native terminal status, fresh work and committed-slot policy admit precisely
one held submission/one release as28737465. Shared E100 lock remains until
native released-slot corroboration. No new independent chain/dependency
write or changed running source/physics. Runnable pending jobs reserve slots;
ambiguous replies retain durable intents and cannot be retried.

Child starts00:45:12/node25/32CPU, exact step0 parent force at00:46:02,
9caches/0new initial SCFs and first new runDSIZE32. Initial00:47snapshot has
0completed fresh calls; later00:54snapshot has4, first109.8438846seconds.
Both snapshots retain their temporal meanings. Authorized production runs
new DFT; audit/admission/capture themselves do not. M28727368step20/.106688
remains live;28709791reversing remains user-held. G2still5/8, no G3/mode/
bottleneck/curvature/physicalholdout promotion or unconverged paper labels.
OldE100helper correctly rejects new unknown child until a reviewed adapter
recognizes this accepted receipt/journal; do not bypass it or reset intent.

152final focused tests bothlocal/independent archive pass0fail/error/skip,
26.831/29.520seconds;83scientific code/native data/test files byteequal (five
delivery metadata excluded after correcting the comparison list). No freshfullsuite or
HFpytest claimed. Only1MB portable raw evidence added, no licensedPP/orbitals.
Large archives/temp onE; five original dirty userfiles untouched/un-staged.
One documentation-context patch rejects before writes, then corrected.
Paper's25-page scientific text/figures unchanged; only honest status pointers.
Goalactive/incomplete. PROGRESS is verified production audit/cache handoff,
guarded actual continuation and tested delivery, not status-only bookkeeping.
Evidence:clamped_flip_continue_E103_20261012/validation.json.

## E104 — Actual +1% M completion, last original training start released

Previous turn:E103 is PROGRESS. Finite question: can actual M ordinary
completion free one committed slot for the ORIGINAL reversing starter,
without weakening the old input/seed guards or resetting release intent?

M28727368 finishes01:14:54CST, step22/.09858846968401362. Original pinned
E093 audits154fresh SCFs/338runtime files/six physical inputs/nativeDSIZE32/
full EFS and all9terminal images. Portable62files/1,075,312bytes preserve
native raw evidence; nine terminal raw SCFs and the same boundary force are
reparsed locally. Native full154audit is not called a154local re-audit.
Sampled forward/reverse97.2277185374/211.9746919066meV/HfO2, four formula
units/highest image3, are ordinary training labels. Neither stationary TS,
barrier error nor significant~1.23meV strain trend is established. G2now6/8.

New E104 adapter imports exact unchanged E100/E103 helper bytes, recognizes
only the actual accepted preserving child via receipt/journal/native source
and explicitly substitutes only the in-process slot counter. It requires
actual M ordinary audit and keeps old seed/input/338source/binary guards,
original exclusive journal/shared lock and all zero-strain terminal gates.
Two live allocations/runnable pending still reserve the same ceiling;
ambiguous reply leaves intent and prohibits re-release. Initial158local/
cleanarchive tests and actual HF read-only preflight pass before release.

Exactly one release of original28709791 is accepted/corroborated01:27:25,
zero new sbatch/dependency/independent chain. It starts01:27:26/node6/32CPU,
original4h/10-update cap/2cached endpoints/7fresh interiors. At01:29:51,
one fresh SCF125.5713876s/nativeDSIZE32, stderr0/no failure/no terminal
summary; no full initial chain/fmax claim. Preserving28737465 stillRUNNING,
step3/.867654. Admission/audit/capture zeroDFT is not entire-turn zero
productionDFT. Both original source/seed/helper bytes and live code unchanged.

Final163focused tests local/independent Git archive pass,0fail/error/skip,
15.499/15.641s.70case code/native data/README/test files byte-equal to that
archive;3new post-archive delivery metadata files excluded. No freshfullsuite
or HFpytest claim. Large temps onE, five original user changes preserved.
No CI, G3/reference/bottleneck/curvature or physical+.005holdout, new monitor,
software install or publication. Scientific figures/tables/25-page PDF are
unchanged, with honest current readiness pointers; no unconverged flip peak.
Goalactive/incomplete. PROGRESS is actual completed material and guarded
production start, not status-only bookkeeping or a new acceleration claim.
Evidence:clamped_reversing_release_E104_20261012/validation.json.

## E105 — Actual paired M strain response and strong-reference controls

Previous turn:E104 is PROGRESS. Finite question: does the sampled M escape
maximum follow its final-well change, and do the unchanged T/Cmma spaces
represent the new strained band's atomic gradient? This uses only the actual
audited +1% M terminal28727368 and unchanged zero-M28722320, not a new DFT
trial, G3 bottleneck choice or physical+.005holdout.

Four original scalar analyzers independently reparse the same nine portable
terminal SCFs. The historical E104/E093 full154SCF/338runtime/six-input/
nativeDSIZE32 audit remains; it is not repeated locally. The strained sampled
forward/reverse peaks97.2277185374/211.9746919066meV/HfO2 enter Tables5--7.
The forward peak change-1.230403416 contrasts with final-well change
-22.206478717 and the B1 fixed-peak/final-following nulls. It is descriptive
training evidence, not independent prediction or significant small-barrier
response. Both escape pairs are retained; unfinished switching prevents
the complete selectivity conclusion.

At strained highest image3, T triplet displacement/gradient coverage is
54.9251%/3.9008%, versus91.5873--91.5896% gradient coverage in ALL four
original Cmma frames. The original T doublets/ranks/mass chart/integer lift
stay fixed. Strong controls cannot be replaced by a weak geometric control
to manufacture an innovation. Fractions exclude cell response and are not
mode-energy allocations, curvature or transverse stability certificates.
Native/dual signed work agrees to5.55e-17eV/cell; peak/endpoint coarse
energy-minus-work21.305139/2.658882meV/HfO2 are not physical error bounds.

The existing guarded pair builder and renderer are reused. A thin M entry
provides adequate negative-energy margin; legacy T defaults are unchanged.
Figure9/Section3.11 preserve both matched PO+ references and individual arcs,
all four strong registrations and absolute-gradient context. The figure
contract precedes code. Actual PNG and all27reading-PDF page contacts are
inspected; detail pages15/16/25/26/27 and same-page Figure9caption are clear.
Actual LuaLaTeX compile exits0, no TeX warnings/missing glyphs/overflow,
52local manuscript links resolve. Reading draft is not submission-ready.

Initial full dirty-tree run2436passed/19failed/1skip identifies exact old
receipt expectations, deliberately pinned renderer dependencies and an old
analyzer metadata hash, not bad SCFs. Strict known code upgrades are reviewed
and re-pinned; original data/figures/report bytes and material/holdout guards
stay unchanged. Failed JUnit retained. Initial README percentage transcription
is corrected from the actual CSV before delivery. Corrected focused local
93tests pass; final93pass/0fail/error/skip in10.291s after all guard changes.

Final clean Git tree53bf3c2/archive1d2c3209 passes2436fulltests/2skips/0fail/
error in470.924s. Skips are absent optional local ASE-ABACUS reader and
nonredistributed licensed VASP POTCAR, not HF validation passes.314existing
deprecation warnings remain. Independent four scalar replays preserve every
scientific field; only two checked LF/CRLF code metadata fields differ. New
figure CSV/PNG/TIFF are byte-identical.29owned scientific/README files match
the tested archive:25byteequal, four exact normalized text; six later delivery
metadata files excluded. Native DFT bytes are never normalized.

02:12:33native slots remain28737465/28709791 RUNNING32/hfacnormal01, latest
complete steps7/.337210 and2/.150233; both stderr0. No stop for a rebound,
source overwrite, submission/release/dependency change or physical parameter
drift. G2still6/8, no G3/reference/curvature/holdout access. Finish/audit both
training flips, then follow original registered gates. No new matrix,
install/monitor/CI/release/PyPI/submission. Original five user changes remain
unstaged/untouched. Goalactive/incomplete. PROGRESS is actual material
analysis, tested figure/reproducibility and scientific manuscript integration.
Evidence:clamped_M_strain_reference_E105_20261012/validation.json.

## E106 — Verified wait and existing terminal-audit readiness

Previous turn:E105 is PROGRESS. This turn is VERIFIED_WAIT, not a new
material result, method implementation or failed experiment. At actual
2026-10-12T02:33:10+08:00, native squeue and sacct corroborate two RUNNING
hfacnormal01 allocations, 32 CPUs each: preserving child28737465/node25
and original reversing28709791/node6. Latest complete updates remain
8/.278775 at02:22:56 and3/.130212 at02:25:47 eV/Angstrom. Both stderr
files have zero bytes. Recent cadence is approximately12/15 minutes per
update, not a guaranteed convergence time. G2 remains6/8; neither current
force nor a live accounting ExitCode0:0 is a successful terminal result.

Inspection confirms the E103 continuation helper already implements
audit-child --parent28709790 and audit-parent --parent28709791. No second
adapter is needed. Actual HF read-only readiness probe completes
02:28:07.497112CST: exact helper SHA6c64c0155800082e28cbbb280511ac761929bdb6ccc5757fa8d68d156b32eac7,
accepted child28737465 from the existing receipt/journal, exact E103 wrapper,
all338 original runtime files and four pinned E093 support hashes match.
Both live jobs are rejected at the native successful-terminal identity
gate with 'exact native successful terminal identity required'. Zero
terminal outputs are created. This is an audit dependency/live-rejection
check, not a full fresh-SCF terminal audit or new HF test-suite run.

No new DFT calls are initiated by this turn, scheduler writes, production
source changes, new Gamma/Hessian work, G3 selection or physical+.005
holdout access. Existing production DFT continues in the two live jobs.
After actual native terminal status, reuse the existing helper to audit
SCFs, raw EFS, source/input hashes and ordinary convergence versus a step
cap. Do not repeat a live audit-refusal probe near termination or confuse
walltime/cap with convergence. No duplicate job, stop for one rebound,
new monitor, installation, CI, release or publication.

Five original tracked user edits and old untracked files remain untouched
and unstaged. HEAD/origin main milestone is c6a7cff9a6ca195b5005c224ce6468db52fab641.
Only the research ledger is updated; no status-only commit or push. Goal
stays active/incomplete. The verified wait is not a blocker and does not
constitute completion of the scientific objective.

## E107 — Verified wait; registered next-gate scope rechecked

Previous turn:E106 is VERIFIED_WAIT. Native02:36:24CST corroborates the same
two live32CPU hfacnormal01 jobs. Preserving28737465 now has complete update
9/.250953 at02:35:06; reversing28709791 remains3/.130212 at02:25:47. Both
stderr files are empty. This is VERIFIED_WAIT, G2 remains6/8. Re-read the
original v1/v2 prediction protocols and main-plan scientific/gate sections:
complete both training flips before G3 selection; retain B0/B1 and all four
strong Cmma registrations, two-bottleneck/two-face limits, prospective
prediction freeze and only two physical+.005 holdout paths. No new model,
DFT/scheduler write, source change, extra probe or holdout access. Original
five user edits remain unstaged; only existing ledgers are dirty fromE106.
No status-only commit, push, pause or completion claim.

## E108 — New publisher SI visual evidence and manuscript qualification

Previous turn:E107 is VERIFIED_WAIT. E108 is PROGRESS: an actual normal
publisher download succeedsHTTP200/7471830bytes/four pages/rawSHA
c4f05b4b3da81f9cc6eb4d31f159cd6822750b73ab1fed1d154209daa1515f09.
All four pages are extracted/rendered/read, including S1--S3 andTableS1.
This closes E082's binary/visual access gap; its previous failures remain.
The new source receipt distinguishes interpolation/amplitude normalization,
per-cell/per-formula-unit labels and each initial energy reference. Restricted
relaxation is prior work, not our prediction/TS innovation. The current
publisher April8 Version of record label also qualifies E035's too-broad
'only deposit' statement; historical reference-PDF identity is not inferred.
No publisher raster, coefficient, curve or vector enters the repository.

Original six electronic assets/optimizer/live source/registered references
and B0--B5/G3/holdout policy remain unchanged. The manuscript gains only the
source-grounded comparison. Original55 numeric table lines and9 figure
references match c6a7cff exactly;105 local links resolve. Existing
TeXLive2023/latexmk/LuaLaTeX compiles27pages/1374783bytes/exit0, actual TeX
warnings/overflow/missing glyphs0. All27 page contacts and2/3/26 full-size
pages inspected; Figure9/fullcaption stillpage26. No code change, new full
suite/cleanarchive claim or new material result. Reading draft remains
unready for submission until complete scientific/independent evidence.

Actual02:45:07CST native jobs28737465/28709791 bothRUNNING32/hfacnormal01,
complete9/.250953 and4/.120019. G2 remains6/8. No new submission/scheduler
write, stop for a rebound, new DFT initiated by this turn, G3 selection or
physical+.005 label access; authorized original production DFT continues.
Five original user tracked edits and olduntracked files preserved/unstaged.
The meaningful source/manuscript milestone, not E106/E107 waiting alone,
is prepared for ordinary explicit-stage commit/push. Goalactive/incomplete.
Evidence:[manuscript receipt](../../../paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.Lee-SI-visual-E108.validation.json).

## E110 — Fully prescribed conditional points obey the geometry guard

Previous turn:E109 is VERIFIED_WAIT. E110 is PROGRESS: a reproduced local
software defect allowed `safeguarded_bfgs` to bypass the supplied geometry
guard when no orthogonal optimization direction remained. Four new
parametrized cases fail before the correction, then pass. The existing
validator now runs before the static call in that closed-subspace branch.
Default L-BFGS-B and nonempty orthogonal optimization remain unchanged;
this is a safety correction, not a material prediction or acceleration win.
An empty projected gradient is not an unmeasured-stability certificate.

Corrected focused tests pass69, expanded tests pass116. A prematurely
started full-suite attempt is stopped by verified own Python PID66200,
because archive extraction completion had not been observed before launch;
it contributes no validation claim. Once extraction has ended and all5766
tracked files are present, a distinct final retry imports from the exact
clean archive and passes2440/skips2/failures0/errors0 in466.907JUnit seconds
(466.92console seconds),312warnings. Skips are the optional ASE ABACUS
reader and nonredistributed licensed VASP POTCAR. Tested tree386b47f17abe9a45d55e65aed2d319b82c5cfa97;
archiveSHAe90a923ccc0824035a1f92e8a1c223ca0fdd18601d54a88221f14fe27bb39d30.
The three changed code/test/doc files match the archive after text-only LF
normalization; delivery metadata is outside that tested tree. Failed and
passing JUnit records and raw hashes are retained in the bounded evidence
bundle. No source archive, pseudopotential or orbital is redistributed.

Native03:24:24CST still shows28737465/28709791 RUNNING32/hfacnormal01,
complete13/.178737 and6/.112694, stderr0. Preserving's previous complete
step12/.164542 has one rebound; this is not a stop condition. G2 stays6/8.
The local correction is not deployed over the immutable live E055 runtime.
No new DFT initiation, scheduler mutation, changed physical parameters,
G3 selection or physical+.005 holdout access. Authorized existing DFT
continues. Original five tracked user edits and olduntracked files remain
untouched/unstaged; the actual tested code milestone, not waiting alone,
is prepared for ordinary explicit-stage commit/push. Goalactive/incomplete.
Evidence:[receipt](closed_conditional_guard_E110_20261012/validation.json).

## E112 — Byte-bound complete-panel forecast storage, no material freeze yet

Previous turn:E111 is VERIFIED_WAIT. E112 is PROGRESS: the calculator-free
generic freezer closes a missing artifact boundary between the implemented
local B2--B5 response kernels and full-panel scorer. It snapshots registered
models/cases, signed abstentions and the same visible input bytes for every
model; loading requires an externally pinned manifest and complete hash graph.
It neither computes/certifies material models nor authorizes holdouts. Smaller
generic panels do not satisfy the HfO2 v1/v2 protocol. All four strong Cmma
registrations and original scientific gates remain necessary for the real batch.

Focused checks first pass102, then231. The initial same-platform Git transport
exercise passes1but masks index normalization. The strengthened CRLF-staging/
LF-export experiment actually fails1atinputs/0003; both raw receipts remain.
The freezer now emits and hashes its own fixed namespace-local byte-preserving
Git rule, with no global configuration change. Final focused tests pass234,
including actual temporary-repository transport and an externally pinned
manifest check. The simulated symlink tests are not native Windows experiments.
Equal synthetic predictions have exactly zero shared-reference gain.

The earlier full regression is stopped by verified own Python PID46228 after
the real transport failure; no complete JUnit/full-suite claim is made for it.
Final tested tree71608b7f94608598e79892c7bcd066559c14a2b3, archiveSHA
769086cc9551dd4f4aa7117fdd7bb22041e74cbd3920c11974319742054e501f.
All5776tracked files are present before final launch; imports come from that
exact archive. Full regression passes2488/skips2/failures0/errors0 in460.100
JUnit seconds(460.12console),314warnings. Skips are the optional ASE ABACUS
reader and undistributed licensed VASP POTCAR. Four source/test/docs files
match after text-only LF normalization; delivery metadata is outside that
tested tree. All six raw JUnit artifacts and their hashes are retained.

Native04:07:34CST jobs28737465/28709791 remainRUNNING32/hfacnormal01,
complete16/.255895 and9/.124033, stderr0. G2 stays6/8. Rebounds do not stop
the registered finite windows or trigger parameter/optimizer changes.
No new DFT initiation, scheduler mutation, live source/input overwrite,
G3 bottleneck/reference choice, physical+.005 holdout access or paper-result
change. Existing authorized DFT continues. Original five user edits remain
untouched/unstaged. This actual tested software milestone, not waiting alone,
is prepared for ordinary explicit-stage commit/push. Goalactive/incomplete;
no material accuracy, novel theory, TS or JCTC readiness claim.
Evidence:[receipt](forecast_panel_freeze_E112_20261012/validation.json).

## E113 — Audit the original reversing cap and continue the same chain

Previous turn E112 is PROGRESS; E113 is PROGRESS. Original28709791 ends
COMPLETED/0:0 at04:19:38CST, but latest10/0.12692249290222424 is an
unconverged cap, not an ordinary label or a failed job. All77fresh SCFs,
338original runtime files, six physical bytes, genuine DSIZE32 and9latest
caches pass the pinned native audit. Local portable replay reparses nine
retained raw terminal SCFs only; the full77-call audit remains native HF
evidence. No historical-best frame substitutes for the latest chain.

After156pre-admission tests, the unchanged E103 helper journals one held
submission and one release of28764704; only28737465/28764704 commit slots.
The new segment starts04:27:20 onnode6,32MPI/12h/40updates, nine exact
initial caches/zero fresh initial SCFs. At04:29:28 its complete step0
reproduces the parent force and first native log declares DSIZE32, but no
fresh call audit has completed. This is startup evidence, not fresh SCF
completion or convergence. FIRE history resets explicitly; physical inputs,
optimizer arguments and frozen production source remain byte-unchanged.

Ten pre-capture material/synthetic checks pass. Final related regression
passes164 locally/164 independently with no failures/errors/skips; not a
fresh full-suite or HFpytest run. Before independent tests, all5874tracked
files are extracted and88case files have identical bytes in the archive.
Archive preparation catches thirteen logs/trajectories omitted by generic
ignore rules before any independent run; a finite case-local exception fixes
delivery. The superseded first snapshot is kept and not counted as tested.
Only README/validation/JUnit delivery records and research ledgers change
after the tested source/data snapshot. No production source is redeployed.

Native04:36:03CST both jobs remain RUNNING32/hfacnormal01/stderr0:
preserving18/.244708, reversing initial0/.126922. G2 remains6/8. No G3
selection/curvature or physical+.005holdout access; no manuscript-result
change, new independent chain, periodic monitor, CI/release/install or
submission. The original five user edits remain untouched/unstaged. This
is a real material computation handoff, not a status-only commit. Goal stays
active and the full scientific objective remains incomplete.
Evidence:[receipt](clamped_reversing_continue_E113_20261012/validation.json).

## E114 — Verified wait for the two accepted training continuations

Previous turn E113 is PROGRESS. E114 is VERIFIED_WAIT: native04:42:57
corroborates both exact accepted handles RUNNING32/hfacnormal01, followed
by a bounded50.0179second wait and native04:45:12 recheck. Preserving
28737465 completes19/.225112 at04:36:32; reversing28764704 completes its
first new update1/.125715 at04:43:44. Both stderr files remain zero bytes.
An observation expiry is not terminal; neither handle is restarted.

At the measured roughly12minute preserving cadence, the remaining21
updates imply an08:40--09:15CST cap window. Reversing's first new update
takes15m34s; allowing the original chain's recent cadence, remaining39
updates imply14:30--16:20CST, bounded by its16:27:20 native wall limit.
These estimate segment caps, not ordinary convergence or goal completion.
G2 stays6/8; v1/v2 selection, strong controls, error and holdout rules are
rechecked without changing them. No new DFT initiation, scheduler mutation,
source/input/model/paper change, G3 selection or physical holdout access.
Only this local wait ledger changes; no status-only commit/push or new
monitor. Goal remains active/incomplete and the five user edits are preserved.

## E115 — Reconcile the actual six-band signed-work manuscript cohort

Previous turn E114 is VERIFIED_WAIT. E115 is PROGRESS: a real manuscript/data
consistency error is reproduced by one failed pre-fix regression. Table7
already had six bands, but section3.8 still said five/45records/42identities
and listed only five signed-work groups. Replay of54portable native EFS
records from six pinned completed cases establishes50raw-log/STRU identities.
Four zero-strain paths share one initial cache; two +1% paths share another.
Counts and the sixth work contributions are corrected; all55numerical table
lines and9figure references stay unchanged. No new DFT or model labels.

Every retained reference's full/open-cell work matches direct native physical
contractions within1.1102230246251565e-16eV/cell. This is algebraic consistency,
not physical uncertainty or saddle certification. All four Cmma registrations
and signed omitted/translation terms remain. Full remote runtime/assets are
not reaudited locally and original author arrays are not regenerated.

Final scoped regressions pass56 locally/56 independently, no failures/errors/
skips; no full-suite/HFpytest claim. All5885tracked archive files are present
and five changed science/code/test sources compare normalized-equal before
the verified run. Authentic XML/JSON bytes use finite case-local attributes.
The initial invalid-directory process launch executes no command; an earlier
archive run before completion verification is not the final test receipt.
Only delivery metadata/ledgers change after the tested snapshot.

Existing LuaLaTeX compiles27pages, no TeX warning/box/glyph problem. All-page
contacts and changed page21 are actually inspected. Native05:12:44CST both
accepted jobs RUNNING32/hfacnormal01/stderr0: preserving22/.164352,
reversing2/.125607. G2 stays6/8; no G3/holdout access, parameter/source change,
new agent/monitor/install/CI/release or submission. Five user edits remain
untouched/unstaged. A tested manuscript correction, not status-only delivery.
Goal active and full scientific objective incomplete.
Evidence:[case](manuscript_work_cohort_E115_20261012/validation.json).

## E116 — Verified wait on the accepted live training handles

Previous turn E115 is PROGRESS. E116 is VERIFIED_WAIT: native05:16:24
confirms28737465/28764704 RUNNING32/hfacnormal01; an actual50.0116second
bounded wait and native05:18:00 recheck corroborate both remain live.
Complete latest steps are22/.164352 and2/.125607, stderr0. No next complete
iteration appears in that observation window; expiry is not job termination.
No restart, cancellation or repeated submission occurs. G2 remains6/8.

Recent preserving updates take about12minutes; reversing's first two take
15m34s/16m52s. Remaining18/38updates estimate segment40caps near08:50/15:20,
not ordinary convergence or completion of the scientific goal. Original
native wall limits remain12:45:12/16:27:20. Full v1/v2 gates are reread and
unchanged: complete G2 before training bottleneck selection, retain all strong
controls, freeze real predictions before physical holdout access, no false
error/TS certificate from an ordinary residual. No new DFT initiation,
scheduler/source/input/model/paper write, G3/holdout access, new monitor or
status-only commit/push. Only this local wait ledger changes. Original user
edits stay untouched/unstaged; goal active and scientifically incomplete.

## E117 — Verified wait and inspection of the existing training gate

Previous E116 and current E117 are VERIFIED_WAIT. Native05:19:42 confirms
both accepted handles RUNNING32/hfacnormal01; actual50.012second wait and
native05:21:43 recheck corroborate them live, stderr0. Latest preserving
step22/.164352 and reversing3/.126334; the slight reversing rebound does
not trigger stopping or restart. G2 stays6/8. Read the existing E055 actual
incomplete-pilot gate, offline terminal-guard tests and control API docs;
no redundant gate is created and no new test/performance result is claimed.
Complete training, measured error and frozen full comparisons remain required.
No new physical/scheduler/source/model/paper writes, G3/holdout access,
agent/monitor or status-only commit/push. Only local wait ledger changes;
original user edits preserved, full objective active and incomplete.

## E118 — Guarded grid integration and all-start callback costs

Previous E117 is VERIFIED_WAIT. E118 is PROGRESS: two actual wrapper defects
are reproduced by failing tests. Existing guarded single-point options cannot
reach the grid; a4point/3start prescribed grid makes12callbacks but reports4.
The grid now forwards optimizer/validator/trust/start policy unchanged and
counts every start, including unconverged ones. Selected single-point counts
remain compatible; missing legacy multi-start outcomes give unknown totals.
Exact-cache callbacks are not new SCFs/core-hours. This is integration/cost
accounting, not a new physical method or HfO2 prediction advantage.

62focused tests pass, including nonlinear signed starts, rejection before
static evaluation, invalid options, unsuccessful-start costs and the exact
clamped ASE/EMT evaluator/cache flow. The independent full suite passes2518,
2skips,0failures/errors,314warnings, JUnit486.829s. Skips are optional
ase.io.abacus unavailable and licensed VASP POTCAR not redistributed. No
install or false pass. All5892archive files are present and three changed
code/test/doc sources compare normalized-equal before full tests. Only delivery
evidence/ledgers change afterward; no fresh local full-suite or HFpytest claim.
Two rejected patch contexts make no file changes before successful editing.

Native05:40:19CST both accepted production handles RUNNING32/hfacnormal01,
stderr0: preserving24/.159496, reversing4/.127876. The reversing rebound
does not trigger cancellation. G2 stays6/8. No new physical DFT, scheduler
mutation, running-source deployment, parameter change, G3/reference selection
or physical holdout access. Manuscript results/figures and five user edits
remain unchanged; no new agent/monitor/CI/release/install/submission. This
tested software milestone is eligible for a normal Git push, not status-only
delivery. Goal remains active and full scientific objective incomplete.
Evidence:[case](conditional_grid_E118_20261012/validation.json).

## E119 — Verified wait for the unchanged two production continuations

Previous E118 is PROGRESS; E119 is VERIFIED_WAIT. Native05:43:46 confirms
both accepted jobs RUNNING32/hfacnormal01, actual50.0114second bounded wait
and native05:45:36 recheck corroborate them live, stderr0. Latest complete
steps24/.159496 and4/.127876 are unchanged within that window. Mild reversing
rebound does not stop/restart a live job. G2 remains6/8. Recent preserving
cadence about11.5minutes implies remaining16updates reach the40cap near08:40;
recent reversing cadence about17.5minutes implies remaining36near15:45--16:15,
bounded by original16:27:20wall limit. These are cap estimates, not ordinary
convergence or goal completion. No new DFT initiation, scheduler/source/input/
model/paper change, G3/holdout access, new agent/monitor or status-only push.
Only local wait ledger changes; original user files preserved unstaged and
full objective remains active/incomplete.

## E120 — Verified wait, preserving continuation reaches step25

Previous E119/current E120 are VERIFIED_WAIT. Native05:49:16 confirms both
accepted handles RUNNING32/hfacnormal01, then actual50.0054second bounded
wait and native05:50:53 recheck corroborate them live, stderr0. Preserving
completes25/.156340; reversing latest complete4/.127876 with its next round
still running. An expired observation is not terminal, and rebound is not a
stop rule. G2 stays6/8; no material bottleneck/reference or holdout selection.
No new physical/scheduler/source/input/model/paper writes, agent/monitor or
status-only commit/push. Only local wait ledger changes; original user edits
remain untouched/unstaged and full scientific objective active/incomplete.

## E121 — Verified wait; optional pause preference requested, not assumed

Previous E120/current E121 are VERIFIED_WAIT. Native05:53:30 confirms both
accepted jobs RUNNING32/hfacnormal01, actual50.0128second bounded wait and
native05:55:28 recheck corroborate them live, stderr0. Preserving latest25/
.156340; reversing completes5/.129821 and continues, not stopped for rebound.
G2 remains6/8. An optional async question asks whether the user wants to pause
automatic goal turns while HF jobs continue, to avoid repeated waiting cost.
No explicit answer is received and no pause/status change is executed; a
preselected UI option is not authorization. This is not a material blocker.
No new physical/scheduler/source/input/model/paper writes, G3/holdout access,
agent/monitor or status-only commit/push. Only local wait ledger changes;
original user edits remain preserved unstaged and full objective active.

## E122 — Verified wait, no terminal or source change

Previous E121/current E122 are VERIFIED_WAIT. Native05:58:25 confirms both
accepted handles RUNNING32/hfacnormal01; actual50.0078second bounded wait
and native05:59:44 recheck corroborate them live, stderr0, latest25/.156340
and5/.129821 unchanged. G2 remains6/8. No explicit pause confirmation is
received and no goal-state change occurs. No timeout-based restart, rebound
stop, new physical/scheduler/source/input/model/paper write, G3/holdout access,
new monitor or status-only push. Only local wait ledger changes; original
user edits remain untouched/unstaged and full scientific objective active.

## E123 — Verified wait, preserving step26

Previous E122/current E123 are VERIFIED_WAIT. Native06:03:19 confirms both
accepted handles RUNNING32/hfacnormal01, actual50.008second wait and native
06:05:46 recheck corroborate them live, stderr0. Latest26/.152847 and5/
.129821; neither is an ordinary pass or terminal. G2 remains6/8. No explicit
pause confirmation, status change, new physical/scheduler/source/input/model/
paper write, G3/holdout access, monitor or status-only push. Only local wait
ledger changes; original user edits preserved unstaged, full objective active.

## E124 — Verified wait, preserving27/reversing6

Previous E123/current E124 are VERIFIED_WAIT. Native06:10:00 confirms both
accepted handles RUNNING32/hfacnormal01; actual50.0123second bounded wait and
native06:12:24 recheck corroborate them live, stderr0. Latest27/.149167 and
6/.131919; neither meets ordinary.10 or is terminal. Reversing's recent rise
does not trigger a premature stop. G2 remains6/8. No explicit pause approval,
goal-state change, physical/scheduler/source/input/model/paper write, G3/
holdout access, monitor or status-only push. Local wait ledger only; original
user edits preserved unstaged and full scientific objective remains active.

## E125 — Verified wait, no terminal or mutation

Previous E124/current E125 are VERIFIED_WAIT. Native06:17:37 confirms both
accepted handles RUNNING32/hfacnormal01; actual50.0113second wait and native
06:19:09 recheck corroborate them live, stderr0, latest27/.149167 and6/
.131919 unchanged in that window. G2 remains6/8; observation expiry is not
terminal. No explicit pause confirmation, goal-state change, new physical/
scheduler/source/input/model/paper write, G3/holdout access, monitor or
status-only push. Local wait ledger only; user edits preserved unstaged,
full scientific objective active/incomplete.

## E126 — Verified wait, preserving28/reversing7

Previous E125/current E126 are VERIFIED_WAIT. Native06:26:34 confirms both
accepted handles RUNNING32/hfacnormal01; actual50.0167second wait and native
06:30:53 recheck corroborate them live, stderr0, latest28/.145328 and7/
.133973. Both remain above ordinary.10; G2 stays6/8. Neither a timeout nor
rebound causes a restart/stop. No explicit pause confirmation, goal-state
change, physical/scheduler/source/input/model/paper write, G3/holdout access,
monitor or status-only push. Local wait ledger only; original user edits
remain preserved unstaged and full scientific objective active/incomplete.

## E127 — Verified wait on both live handles

Previous E126/current E127 are VERIFIED_WAIT. Native06:33:32 confirms both
accepted handles RUNNING32/hfacnormal01; actual50.0158second wait and native
06:35:38 recheck corroborate them live, stderr0, latest28/.145328 and7/
.133973 unchanged in that observation window. G2 remains6/8. No explicit
pause confirmation, goal-state change, new physical/scheduler/source/input/
model/paper write, G3/holdout access, monitor or status-only push. Local wait
ledger only; original user edits preserved unstaged, full objective active.

## E128 — Verified wait, preserving29/reversing8

Previous E127/current E128 are VERIFIED_WAIT. Native06:40:32 confirms both
accepted handles RUNNING32/hfacnormal01; actual50.0114second wait and native
06:43:32 recheck corroborate them live, stderr0, latest29/.141436 and8/
.135846. G2 remains6/8. No explicit pause confirmation, goal-state change,
new physical/scheduler/source/input/model/paper write, G3/holdout access,
monitor or status-only push. Local wait ledger only; original user edits
preserved unstaged and full scientific objective remains active/incomplete.

## E129 — Verified wait, unchanged live handles

Previous E128/current E129 are VERIFIED_WAIT. Native06:45:50 confirms both
accepted handles RUNNING32/hfacnormal01; actual50.018second wait and native
06:47:13 recheck corroborate them live, stderr0, latest29/.141436 and8/
.135846 unchanged. G2 stays6/8; no terminal assumption or timeout restart.
No explicit pause confirmation, goal-state change, new physical/scheduler/
source/input/model/paper write, G3/holdout access, monitor or status-only push.
Local wait ledger only; original user edits preserved unstaged, goal active.

## E130 — Verified segmented wait; fewer native scheduler polls

Previous E129/current E130 are VERIFIED_WAIT. Native06:51:04 confirms both
accepted handles RUNNING32/hfacnormal01. Five interruptible50second sleeps
actually total250.0559s; no single blocked wait exceeds60s. Native07:00:42
recheck corroborates both live, stderr0, latest30/.137561 and9/.137528. The
observed native window is578s including intervening work, not just sleep.
Only two native scheduler observations occur in that window; no new monitor
or unattended process is created. G2 remains6/8. No explicit pause approval,
goal-state change, new physical/scheduler/source/input/model/paper write,
G3/holdout access or status-only push. Local wait ledger only; original user
edits preserved unstaged and full scientific objective active/incomplete.

## E131 — Verified wait; authoritative goal status remains active

Previous E130/current E131 are VERIFIED_WAIT. get_goal confirms active with
the full unchanged objective. Native07:10:38 confirms both accepted handles
RUNNING32/hfacnormal01; actual50.0055second wait and native07:14:52 recheck
corroborate them live, stderr0, latest32/.130113 and10/.138981. G2 remains
6/8. No explicit pause approval, goal-state change, new physical/scheduler/
source/input/model/paper write, G3/holdout access, monitor or status-only push.
Only local wait ledger changes; original user edits preserved unstaged.

## E132 — Verified wait, preserving33/reversing11

Previous E131/current E132 are VERIFIED_WAIT. Native07:25:00 confirms both
accepted handles RUNNING32/hfacnormal01; actual50.0205second wait and native
07:26:30 recheck corroborate them live, stderr0, latest33/.126614 and11/
.140244. G2 remains6/8. No pause confirmation, goal-state change, physical/
scheduler/source/input/model/paper write, G3/holdout access, monitor or
status-only push. Local wait ledger only; original user edits preserved.

## E133 — Verified wait, preserving34/reversing12

Previous E132/current E133 are VERIFIED_WAIT. Native07:37:28 confirms both
accepted handles RUNNING32/hfacnormal01; actual50.0159second wait and native
07:39:31 recheck corroborate them live, stderr0, latest34/.123646 and12/
.141321. G2 remains6/8. No pause confirmation, goal-state change, physical/
scheduler/source/input/model/paper write, G3/holdout access, monitor or
status-only push. Local wait ledger only; original user edits preserved.

## E134 — Verified wait, unchanged live handles

Previous E133/current E134 are VERIFIED_WAIT. Native07:48:52 confirms both
accepted handles RUNNING32/hfacnormal01; actual50.0123second wait and native
07:50:15 recheck corroborate them live, stderr0, latest34/.123646 and12/
.141321 unchanged. G2 remains6/8. No pause confirmation, goal-state change,
physical/scheduler/source/input/model/paper write, G3/holdout access, monitor
or status-only push. Local wait ledger only; original user edits preserved.

## E135 — Verified wait, preserving35/reversing13

Previous E134/current E135 are VERIFIED_WAIT. Native07:53:58 confirms both
accepted handles RUNNING32/hfacnormal01; actual50.0061second wait and native
07:55:34 recheck corroborate them live, stderr0, latest35/.121285 and13/
.142231. G2 remains6/8. No pause confirmation, goal-state change, physical/
scheduler/source/input/model/paper write, G3/holdout access, monitor or
status-only push. Local wait ledger only; original user edits preserved.

## E136 — Verified wait; correct ordering of local wait notes

Previous E135/current E136 are VERIFIED_WAIT. Native07:59:29 confirms both
accepted handles RUNNING32/hfacnormal01; actual50.0146second wait and native
08:01:03 recheck corroborate them live, stderr0, latest35/.121285 and13/
.142231. G2 stays6/8. Own E133--E135 MD notes are reordered chronologically;
their contents and prior JSONL records are unchanged. No scientific/source/
input/model/paper/scheduler change, G3/holdout access or status-only push.
No explicit pause confirmation or goal-state change; original user edits
preserved unstaged. Only local wait/evidence metadata changes.

## E137 — Verified wait, preserving36/reversing14

Previous E136/current E137 are VERIFIED_WAIT. Native08:05:02 confirms both
accepted handles RUNNING32/hfacnormal01; actual50.0093second wait and native
08:07:02 recheck corroborate them live, stderr0, latest36/.119267 and14/
.143033. G2 remains6/8. No pause confirmation, goal-state change, physical/
scheduler/source/input/model/paper write, G3/holdout access, monitor or
status-only push. Local wait ledger only; original user edits preserved.

## E138 — Verified wait; resume yielded read sessions, no restart

Previous E137/current E138 are VERIFIED_WAIT. Local read sessions34231/62049
yield and then complete exit0 when resumed, not restarted. Native08:13:18
confirms both accepted handles RUNNING32/hfacnormal01; actual50.0205second
wait and native08:16:56 recheck corroborate them live, stderr0, latest37/
.117605 and14/.143033. G2 remains6/8. No pause confirmation, goal-state
change, physical/scheduler/source/input/model/paper write, G3/holdout access,
monitor or status-only push. Local wait ledger only; user edits preserved.

## E139 — Native-cache diagnosis, not another status-only milestone

User explicitly resumes research. Both original accepted children remain
RUNNING32/hfacnormal01 at08:28:37: latest38/.116223 and15/.143674. G2 remains
6/8, full goal active. No pause or completion inferred from automatic turns.

Existing pinned E103/E093 analysis support verifies338original runtime files
and exports preserving37/reversing15 in a separate new analysis namespace.
All nine exact raw SCFs per snapshot pass six physical hashes/DSIZE32/EFS/
ordered geometry and same-boundary residual replay. No new DFT is initiated
by analysis; the two already-authorized production jobs continue normal DFT.
No scheduler/production/input/optimizer/geometry/endpoint changes.

Both are limited by central image4's generalized atomic block. fmax values
.117604576798/.143674141082; same-vector spring projections2.8811661e-7/
3.6067163e-12. No spring dominance at these geometries. This does not prove
the cause of the force rise, an optimizer defect, faster alternatives, TS/
barrier accuracy, final mechanism or G3/holdout admission. No paper labels
updated. Segments run to their existing finite limits without stopping on
rebound, relaxing.10 or changing100Ry/full10auDZP/Gamma2^3/P0/E0/noCI.

Portable114files/1,423,411bytes retain18raw SCF records, not the undistributed
physical assets or all historical calls. Six new case checks and relevant
replay/admission regression total123pass locally and123in actual Git archive,
0failure/error/skip. Before independent tests all6013tracked files exist,
115raw files byte-equal and4source/metadata texts newline-normalized equal.
Initial comparison rejects metadata CRLF/LF before tests; no physical bytes
changed. Not a new full-suite/HFpytest run or physical prediction evidence.

Evidence is in clamped_flip_residual_E139_20261012; protected five user files
and older untracked artifacts remain unmodified/unstaged. One inline tooling
script initially rejects a duplicate JS variable declaration before execution
or any writes; corrected command runs once. The scientific milestone is this
native diagnostic/replay package, not the preceding repeated wait records.

## E140 — Successful child40 cap audited; ordinary G2 remains6/8

Previous E139 is PROGRESS, not an unexecuted plan or status-only turn. Five
actual bounded50s waits total250.0408s between native polls; no new monitoring
automation.28737465 ends nativelyCOMPLETED0:0 at08:50:08, final40/.114146631939.
This is not failure or ordinary convergence. The pinned existing E103 child
audit runs once only after actual terminal state, auditing280fresh native
SCFs/338original runtime files/six physical bytes/DSIZE32/fullEFS;9latest
cache evaluations initiate0SCFs. Receiptdf34573384f76650e123691dbb508a52d7d8ee529e6e50dca0fce5e4e6835c02,
seed1a222cc63b89344816f0829bbf079793fa06306d3d4d13f5cc3c2b401f5fbf3f.

Local77files/1,250,797bytes retain9raw terminal SCFs and the exact seed, not
the undistributed physical assets/all280calls/HF runtime. Four new real-case
checks and related regression total127pass locally and127in actual clean
archive,0fail/error/skip. Before independent tests all6097tracked files are
present;78raw files byte-equal/4source-metadata texts newline-normalized equal.
No physical evidence bytes changed. Not a new full-suite/HFpytest run.

Healthy segment declines1.325704→.114147, including .137561→.114147 over last
10updates; no SCF/geometry fault, E139 excludes instantaneous spring dominance.
Review supports one more finite40-update12h same-chain segment under unchanged
physics/optimizer/threshold/noCI, at most2study slots; not an automatic rolling
dispatch or new independent chain. Original E103 admission only knows starter
parents: accepted journal must not be reset or re-used to masquerade this
child40 as a ten-step starter. Next action is tested explicit child40 ancestry
and real spare-slot admission, then one held submission/release if all pass.
No such new submission in E140.280calls/29,003.679s≈257.810477call-time proxy
core-hours versus actual29,096s/258.631111allocated CPU-hours; future40updates
roughly8.1h at measured rate,12hwall ceiling, not a convergence promise.

09:06:33 only reversing28764704 RUNNING32/hfacnormal01, complete17/.144601,
stderr0. OrdinaryG2 remains6/8. Current turn is PROGRESS_NATIVE_TERMINAL_AUDIT,
not a blocked or completed goal. No scheduler writes, new DFT initiated by
audit, physical/optimizer/source/paper changes, G3 selection or holdout access.
Protected five user files/old untracked work remain intact and unstaged. One
metadata patch rejects its MD context before any write; actual unchanged
files are checked before the corrected patch, not silently assumed written.

## E141 — Execute the reviewed child40 continuation, not an extra trial

Previous E140 is PROGRESS_NATIVE_TERMINAL_AUDIT. Actual parent28737465 remains
COMPLETED0 with unchanged df3457...receipt and1a222c...latest seed. Only
reversing28764704 occupies one study slot. New explicit child40 ancestry
wrapper/helper are developed in a fresh namespace; original code, audits,
journals and physical settings unchanged. Runtime command after module purge
is byte-identical to E103;100Ry/full10auDZP/Gamma2^3/P0/E0/.10/noCI/32MPI.
Two development string escapes are fixed before first tests/deployment;
only new code LF-formatted. Initial138tests precede8new parent-verifier
checks. Final pre-admission146local/146clean archive,0fail/error/skip,
all6106files present/4live-code-JUnit byte-equal/3source-metadata normalized.

HF exact-byte delivery and zero-DFT actual parent/startup preflight pass.
One held submission/one release admit28812952 with five fsynced journal
events and independently corroborated two slots28764704/28812952. No new
independent chain/dependency write; old intent is not reset or replayed.
09:49:29 starts node75/32CPU.09:52:20 startup: step0/.114147 reproduces parent,
9exact caches/9cache audits/0missing, first logDSIZE32, one fresh-call audit
then completed, stderr0/no final summary/failure. Fresh first raw SCF separately
passes original six bytes/DSIZE32/EFS replay at97.363037s. Readers initiate0DFT/
scheduler writes; newly authorized production does DFT, not an all-zeroDFT turn.
FIRE history resets explicitly; no optimizer settings or physical parameters
change.40updates/12h finite cap, no automatic rolling dispatch.

Three real admission/start/first-SCF checks raise final regression to149pass
locally and149in independently extracted material Git snapshot,0f/e/s.
All6112tracked files present before tests,10live-code/evidence byte-equal,
3source-metadata normalized. Final treebf5ff40ad266ad6dfa638d09c7185277cf4ad5cd;
case clamped_preserving_continue_E141_20261012 retains all actual receipts.
No full-suite/HFpytest/install claimed. Only delivery metadata/independent
JUnit/ledger added after final tested snapshot; admitted live code unchanged.

10:09:21 both RUNNING32/hfacnormal01/stderr0: preserving28812952 step1/.114339,
reversing28764704 step22/.144428. Small preserving rebound is not a stopping
rule. G2 remains6/8; no G3/curvature/reference selection, physical+.005holdout,
paper numerical label, CI/release/monitor or original user-file change.
Full scientific objective stays active/incomplete. Next actual terminal
28812952 uses E141 helper audit-child; the reversing handle28764704 still
uses unchanged E103 audit-child --parent 28709791. Later R admission
must recognize the accepted E141 child explicitly, never loosen unknown-job
guards. Protected five tracked changes/old untracked work remain unstaged.

## E142 — Distinguish ideal nested-path ordering from relaxed-well barriers

Previous E141 is PROGRESS_ACTUAL_BOUNDED_PRODUCTION_HANDOFF. Native10:52:35CST:
28812952 RUNNING32/node75 step5/.115421,28764704 RUNNING32/node6 step25/.142540;
both hfacnormal01/stderr0. No failure or convergence claim from live Exit0:0,
and no stopping for the preserving rebound. G2 remains6/8.

The known Behara subscription-SI gap is boundedly rechecked, not closed.
The existing public accepted PDF SHA5f12ceb8... is rendered with existing
Poppler and visually read at printed pp.2/7/8: cubic-axis strain coordinates,
fully fixed orthorhombic switching routes and a separate planar-boundary
thermodynamic construction cannot be silently substituted for our registered
homogeneous T-substrate ensemble. No author coordinates/barriers imported.

Manuscript2.2 now states ideal minimax ordering only for identical endpoints
and potential with nested continuous path sets. Local finite-image candidate
bands do not inherit a guaranteed ordering. Appendix B's dimensionless
elementary potential shows higher absolute bottleneck energy with a smaller
barrier when newly relaxed wells rise more. This is standard interpretation,
not hafnia evidence, an independent predictor or theoretical novelty. Existing
B0/B1/H1 and the finite scientific matrix remain unchanged; no G3 or holdout.

Ten analytic checks plus the actual six-band manuscript regressions give30
local and30real independent-archive passes,0fail/error/skip. Tested tree
8ea0f78d3744b81b56015f102d0762566172b5ad, all6119tracked files exist before
tests,5new source texts newline-normalized equal and1retained JUnit byte-equal.
All55numeric-table lines and9figure-reference lines retain their old hashes.
Existing LuaLaTeX/latexmk -norc builds the actual28-page reading PDF with
no warning/overfull/underfull/missing glyph. All pages are rendered/contact
reviewed and changed pages3/28 inspected at full reading size; ten prior
numeric sentinels and Figure9 caption remain present. PDF SHA19de0a18...
under E:/TEMP/varneb-E142-manuscript-build; not submission-ready.

No physical inputs, production code, optimizer, material records, scheduler
state or live allocation source changes. No new DFT initiated by this
clarification; already authorized running production still performs DFT.
Only metadata/copied independent JUnit/ledger follow the tested snapshot.
Protected five user hashes and older untracked files remain unchanged and
unstaged. This is PROGRESS_MANUSCRIPT_BOUNDARY_CLARIFICATION, not a G2/G3 gate
advance or complete goal. Keep the full goal active without extra allocations,
CI, install, monitor, release or unsupported scientific labels.

## E143 — Keep all starts in the returned conditional grid and its registration

Previous E142 is PROGRESS_MANUSCRIPT_BOUNDARY_CLARIFICATION. Code inspection
finds a real dataflow gap: the grid correctly totals all-start costs but drops
the point objects holding alternative/unconverged final outcomes. A before-
feature regression actually fails1case. The bounded remedy adds complete owned
read-only point_records to sampled grids, with selected-node/count/shape/type
consistency gates. Optimizer, energy/gradient calls and winner selection do not
change; an all-failed node still raises, and legacy missing records stay unknown.

Prospective snapshots recursively serialize the complete new records and bind
even unselected/unconverged outcomes. Removing/changing them fails registration.
Patch screens preserve the corresponding record slice. Old empty-record
snapshot/registration signatures remain exactly0c1c8d.../20fa5e..., measured
before the change and checked afterwards. Curvature/whole-grid/TS gates are
not relaxed. Final outcomes are not initial seeds/trajectories or raw audits.

28new tests;132local and132actual independent related passes,0f/e/s. Full
tracked-source archive regression finishes2600pass/2skip/0fail/0error,
314warnings,646.275s JUnit. Actual skips are optional ase.io.abacus and a
licensed nonredistributed VASP POTCAR; no package or environment installation.
All6126tracked files exist before tests,7sources newline-normalized equal and
2original JUnit bytes equal. Tested tree dbc9e92d59b3abb5dae164a238050411812655b4,
archive SHAaa42b04e...; only delivery metadata/README/copied JUnits/ledger follow.
Toy evidence retains4nodes/4points/12starts;24actual/all-start callbacks versus
8selected-only. Not an acceleration or material H2 benchmark. The initial
analytic seed fixture is refined to exact wells, not a timing comparison; one
doc patch rejects its context before any write and actual state is checked.

11:26:23 both original handles RUNNING32/hfacnormal01/stderr0:28812952 step8/
.113387,28764704 step28/.139222. G2 remains6/8. No scheduler writes/new material
DFT initiated by this milestone, source deployment, optimizer/physical changes,
G3 bottleneck/mode/Hessian/physical+.005holdout, CI/install/monitor/release or
paper/PDF change. Already authorized production continues normally. Protected
five user hashes and old untracked work are preserved unstaged. Outcome is
PROGRESS_IMPLEMENTED_BRANCH_PROVENANCE, not completion of G2/G3/full goal.

## E144 — Cumulative attempts are finite without changing physical inputs

Previous E143 is PROGRESS_IMPLEMENTED_BRANCH_PROVENANCE. The evaluator's
completed-call counter resets after restart and cannot bound failed/retried
fresh attempts. Its optional admission API is first reproduced missing1case.
A separate immutable operational policy now reserves exclusive indexed slots
before factory invocation, consumes failures/interruption conservatively and
allows exact cache hits without new reservations. Uncapping/enlarging ordinary
resume, missing/corrupt history and stale writers are refused, not repaired.
Legacy physical cache contract1308f346... stays identical; no electronic,
optimizer or geometry tolerance is retuned. These are attempts, not raw SCFs,
electronic steps or core-hours; native resource audits remain mandatory.

37new checks;172local/172actual independent related pass,0f/e/s. Full tested
archive2637pass/2skip/0fail/error,314warnings,413.513s JUnit. Actual skips are
optional ase.io.abacus and nonredistributed licensed VASP POTCAR; no installs.
All6140tracked files exist before independent tests,6source texts normalized
equal/6original receipt bytes equal. Tree b43c1ce0d76c6dbc1291fc8ae096e944a07f7f51,
archive97a7db1d...; only delivery metadata/README/copied reports/ledger follow.

Local, actual-archive and actualHF analytic replay agree on15fields except
host path:2factory invocations/1intentional failure/2reservations, exact cache
resume creates0new completions and two policy changes are refused. All three
are0DFT/0scheduler writes, not material/G3 or prediction/acceleration evidence.
HF isolated test namespace receives57Python source files only, tar f8beafda...,
core B1bf5c28... and script6f527d... byte-matched to the actual tested archive.
Working script's28e643... newline difference is explicitly retained; normalized
equality was proved before tests. No fullHFpytest or 338-file production update.

12:19:07 both study handles RUNNING32/hfacnormal01/stderr0:28812952 step13/
.105734,28764704 step31/.134434. G2 remains6/8; no terminal label promotion,
new allocation/material DFT initiated, selected production budget, G3 mode/
centre/Hessian or physical+.005holdout. No CI/install/monitor/release or paper/
PDF/optimizer/physical-input changes. Original five protected user hashes and
old untracked work stay unchanged unstaged. Already authorized production
continues. Outcome is PROGRESS_IMPLEMENTED_FINITE_ATTEMPT_ADMISSION, not full
goal completion or a narrower substitute for the independent material evidence.

## E145 — Seventh ordinary-converged G2 band, not a favorable-window claim

28812952 completes +1% preserving at step17/.099283376. Original native
119SCFs/338runtime/sixinputs/DSIZE32 audit passes; nine portable E/F/stress
frames independently reparse. The sampled peak125.424431meV/fu increases
28.596248 over zero strain, not the hypothesized favorable switching drop.
All four original-reference scalar analyses independently replay in a clean
archive; no new modes, G3 or physical+.005label is selected/read.

The manuscript adds actual Table5--7 rows/section3.12 and reconciles seven
bands/63records/58different native identities. All42old numeric lines and
nine figure references remain. Full28-page PDF compiles, all contacts and
changed details are visually inspected; no compile/layout warning/error.
55local/55independent related pass; first full run's2stale assertion failures
are retained, then final actual archive2642pass/2skip/0fail/error,314warnings,
600.660s JUnit. Skips: optional ASE ABACUS reader and licensed POTCAR.
Tree6dad7459db95e24dd13d0b554e69ca822c0ec271, archive3010861c5bdd55cf0927a3bbf8a4afa0b3ae62f2c604ae6bced0262c39528e9b;
all6220files present before tests,62native files byte-equal/9sources normalized
equal. Only delivery metadata, copied final reports, plan/ledger follow.

G2 is7/8; full network response, conditional stability and independent material
prediction remain incomplete. Original physical/optimizer contracts and five
protected user files stay unchanged. See clamped_plus1_preserving_terminal_E145_20261012.

## E146 — Audit the reversing cap, then continue the original chain once

28764704 ends successfully at step40/.109959903, not ordinary convergence.
The original E103 audit checks280fresh SCFs,338runtime files, six physical
hashes, DSIZE32 and nine exact latest caches. Portable77files/9raw SCFs are
retained unchanged. Last five residuals decrease; no physical fault or
single-rebound stop is inferred.

The reviewed continuation changes ancestry/namespace only. All commands
from module purge onward match E103 exactly.78local/78independent related
checks pass/0f/e/s, including exact startup, cap/seed replay, unsafe/ambiguous
admission refusal and actual one-handle receipt/journal. Treeeb710856a0b951afe16b375ad41d7f9750f1ca4b,
archive53405e98000fb940fdf5ee5f958ba27824b3713adbe3c083651157d78ebc6df7;
all6308tracked files present before tests,77native and4admission/source files
byte-equal. No whole-new-suite/HFpytest or new scientific method claim.

After actualHF preflight, one held28852667 is submitted and released once.
It starts16:16:01CST/node75/32MPI/hfacnormal01. At16:23:52 the initial nine
caches remain exact, first fresh native SCF reparses six hashes/EFS/DSIZE32,
108.154542seconds, stderr0. New segment has40updates/12h maximum; FIRE history
is reset, not restored. Only one committed study slot, no new independent
chain, source replacement, physics/optimizer-argument drift, CI or monitor.
G2 remains7/8; no G3 mode/centre/curvature or reserved physical+.005access.
See clamped_reversing_cap_continue_E146_20261014. Full scientific goal active.

## E147 — Complete production histories, not just the final segments

The seven completed G2 chains span16successful32CPU/hfacnormal01 allocations,
1918fresh interior SCFs/1966.956848call-time-proxy core-hours, versus
1976.364444allocation core-hours. All capped ancestors stay included.
+1% preserving is476calls/415.495917core-hours, not just119/99.555915.
The original response_cost API is reused unchanged; original receipt pins,
six-input declarations and exact7-image fresh-call grids are checked.
Different actual source calls are not merged by common geometry/log hashes.

24new material/fault plus18existing API checks pass42local/42independent,
0f/e/s. Tree506390df64fb82f51d6094c3acdd9062224833a0,
archive76404675f3317561b4b97ad658eb43c462ef52961da8adfd42dec76d6b92c27c;
all6323files exist before tests,16receipt/accounting bytes identical and2source
texts normalized-equal. Actual independent CLI reproduces every cost/other
field exactly; first byte comparison refuses only the actual script LF/CRLF
source metadata difference. Both original exports remain unmodified.

Not complete endpoint/preparation/G1/phonon/failed-exploration/G3 or whole-study
cost, nor acceleration/prediction advantage. No newDFT/scheduler writes,
fullsuite/HFpytest, paper/PDF or physical/optimizer/runtime changes. The advanced
manual distinguishes attempts, callbacks, nativeSCFs and allocation cost.
At16:48:52:28852667RUNNING32/node75/hfacnormal01, step2/.106872/stderr0.
G2 remains7/8, no G3/physical+.005access and full goal active.
See clamped_chain_costs_E147_20261014.

## E148 — Complete eight-band sampled response, retaining the negative result

28852667 ends17:03:53CST,4updates/.099568630/28freshSCFs. Native audit passes
338original runtime files/sixphysical inputs/DSIZE32;62portable files unchanged.
The same reversing chain totals385freshSCFs/432.318217SCF-proxy core-hours,
not just the last28/24.609091. G2 ordinary-audited is8/8. Reversing falls
44.024331meV/fu but lower preserving rises28.596248; selectivity proxy falls
16.249222. Favorable raw H1 signature is absent; unknown barrier errors prevent
error-certified rejection. The+1%T/M escape gap is only.464068. H2 untested.

Actual independent scalar/network/cohort replay,68local/68independent related,
full2700pass/2expected skip/0failure/error,314warnings/630.333JUnit seconds.
Tree56429c7d110f67eb13b002d8e988ef2d6556adff,
archive00601efe797931c1b5ad3da68792ad7848aea82b023cac3483663fe3e324d410;
all6403files present before tests,62native byte-equal/11sources normalized-equal.
Exactly three actual executed-source SHA fields differ only by verified source
LF/CRLF; no physical evidence normalized. Five subsequent prose-only cost-term
corrections change no numbers/code; final68related pass again and current cohort
export regenerated. No post-prose full-suite rerun claimed.

Tables5--8/section3.13 include full response and unknown-error boundary. Current
30-page PDF compiles without found errors/layout warnings; all pages rendered,
four contacts and16/17/22/29details actually viewed. Existing9figures unchanged;
page6 whitespace/chronological narrative still need consolidation. No G3
basis/centre/curvature, physical+.005access, newDFT/scheduler/CI/install/release/
monitor or physical/optimizer/source changes. Not submission-ready/TS-certified
or independently predictive success. Full goal active; next is the explicit
controlled-strain conditional bridge, not a substitute easier fixed-cell face.
See clamped_plus1_reversing_terminal_E148_20261014.

## E149 — Prescribed strain survives the conditional calculator interface

New explicit biaxial mode chart reuses original geometry/physical gradients,
with mechanical cache identity including anchor/scale/control map. A covector
guard on the actual released complement refuses substrate drift before any
calculator call; it covers both conditional optimization and orthogonal
curvature. Existing reference/fixed-substrate cache contracts remain unchanged.

200local/200actual independent related pass,0failure/error/skip. Includes
independent EMT finite energy gradients, both conditional optimizers, recovery/
attempt budgets/invalid input and18native HfO2zero/+1% common-chart replays.
No new HfO2DFT or numerical-error/TS/model-success claim. Tree6884cf4ca3152e19b3da5be32fc6b90ddf1284b8,
archive5526ab7b335544ccf6be16d50b5d47129356ca4685782911c6f7fb42b9e29b1f,
all6412files present and6source texts normalized-equal before tests. Physical
evidence not normalized. At18:09:36 the independent full suite remainsRUNNING
in exec17076; no final result is claimed or substituted from E148.

No G3 mode/centre/face/Hessian/holdout selection or calculation, production
source/electronic/optimizer change, CI/install/release/monitor. HF18:03:10
shows no remaining iai806hfacnormal01 jobs; final reversing isCOMPLETED0:0.
Implementation prerequisite delivered, not full material validation. Full goal
active; next remains the original finite G3 branch/response experiment.
See controlled_strain_bridge_E149_20261014.

## E150 — First finite native G3 local-curvature batch

Prior goal turn is PROGRESS: actual E148 eight-band matrix and tested E149
interface pushed. Its still-live independent full run is resumed, not restarted;
it actually finishes2742pass/2expected skips/0f/e,1171.095JUnit seconds. This
result is E149, not later E150. The original RUNNING observation is retained.

Current progress is two exploratory raw-minimum zero-strain candidate centres,
preserving/T atimage4. Full36internal plus1control, paired.01/.02Angstrom,
148new static attempts each/296total. No invented barrier-error interval or
certified-lowest/TS/final-predictive selection. Original uncertainty/tie policy
and two-centre cap remain; a changed representative does not buy a third.
Physical+.005stays sealed. No conditional face/mode choice/holdout forecast
is completed; the derivative probe basis itself is explicitly frozen.

99local/99independent related pass and actual independent preparation yields
identical318476-byte manifest5f047f54ca83cb2dceb8cb71687eb9206a85c04205cdf26431565bc8a9b98a5a.
Tree0222d4204439ebc450f6c14e97d724e3845ee6bf,
archive191fea45ead2448abc851655cf68315ec4b543287755b04373c62b81ca20f8b4;
all6424tracked files present before tests, three source texts/raw manifest
match. No new whole E150 suite claimed. Analytic mock callbacks are not DFT.

Native3706880-byte payload a4e405d0bcf3be1345ecc5f3a7e91f92acb2f39325c0de64601393835e0255f5
verifies343Python byte hashes/15critical source texts/original six assets and
binary/native centres; local/HF mechanical contracts exactly equal. Initial
two submissions28855614/28855615 are accepted but overly literal NumNodes=1
guard refuses Slurm's valid pending1-1. Original helper/journal preserved.
Separate reviewed recovery only releases these exact handles,0new submissions,
no production source/input/wrapper changes. Final13local/13independent helper
checks pass on isolated4file treebca7020158e298d92d85194fa33cbc9e65fe67c7,
archivead2a2355449ce35cea2dfdc5988f044c32aeb18691d0273601c1e6690cad6b2f.
Original fixture failure is retained, not relabeled as production DFT failure.

Actual starts18:39:04CST/node10,node27,32CPU each/hfacnormal01.18:41:58both
RUNNING/stderr0, first fresh native DSIZE32 each; no completedSCF claim at
that snapshot. At most2committed allocations, exactly2submissions/2releases,
no new NEB/material/orientation, physical/optimizer/old-runtime changes, CI/
install/release/PyPI/monitor. Derivative SCFs are separately charged from later
two81-node/eight-test faces. Not material gate/full-goal completion. Goal active.
See G3_curvature_E150_20261014 for original raw evidence and scope.

18:44:02native follow-up:2/148probe evaluations complete in each job. Both
first fresh SCFs pass electronic convergence, original six byte hashes,
uniqueDSIZE32 and full E/F/stress/40-coordinate gradient; call times
101.678342/104.632985s. Real material data, not whole Hessians or TS/face/
prediction certification. Same unchanged allocations continue; roughly23:00
to01:00CST only if later cost is comparable, not a guaranteed completion time.

## E151 — actual existing-band cross-strain work, no new DFT (2026-10-14)

Outcome: PROGRESS; previous E150 outcome remains partial native launch, not
material certification. Actual admission/pairing uses8bands,72records,
66distinct log/STRU identities,36pairs in a single ordered unwrapped chart.
Two72-record reader passes are not144SCFs. Endpoint trapezoidal work splits
atomic/explicit translation/open-row/prescribed-control contributions;
independent physical contractions agree to6.662e-16eV/cell. All actual peak
indices match at0/+1%. Keep connection dependence and nonzero residuals:
preserving control-12.859/atomic+36.895/observed+28.596/residual+4.639meV/fu;
reversing residual-11.495meV/fu. These are neither causal mode energies nor
an envelope derivative, frozen-slice measurement or physical error interval.

Six analytic complete-harmonic-release checks tie under retained rank1/3/4
and rotated references, nonzero gradients/internal index1. No reference
rename alone improves this exact limit; strong baseline and preregistered H2
are retained, not weakened. Real branch/finite-basis/nonlinear/cost evidence
is still required. Initial source caches are exact; binary64 serialization
lift4.441e-16Angstrom/work3.469e-18eV exposed an overly strict equality.
Explicit1e-12 algebraic guards are not physical accuracy bounds or a repair;
an integer initial relift still fails. No source E/F/stress/ASR change.

72local/72independent related checks pass with0failure/error/skip; actual
independent CLI reproduces all scientific fields. Tested tree
646b7cfdc2e80448ddaa98bd6855d429b6424079/archive
a0b88df333ea1377c475640978dbc8aea6ce22a0a094452ce2abdf48ac4c9017;
all6449tracked files present before tests. Verified LF/CRLF code SHA metadata
only differs in replay. No whole E151 suite or paper/figure change claimed.
At19:20:33CST both original curvature jobs remain live,23/148complete each,
stderr0,receipt means102.655/104.513s; raw outputs not reparsed by this
status check. No DFT/scheduler/live-source/reserved-geometry change. Current
G3/full goal active and incomplete. See cross_strain_work_E151_20261014.

## E152 — actual native prefix raw replay, not complete curvature (2026-10-14)

Previous E151 classified PROGRESS. New read-only native auditor reuses the
unchanged E150 runtime, native result validator and existing curvature/work
tools. It snapshots only atomic complete prefixes; no calculators, scheduler
writes or live-source changes. Actual34+34native probe results reparse with
original six byte hashes, oneDSIZE32/converged SCF, geometry, raw E/F/stress
and exact full-cache gradient. Two parser reads/log and a cached-center
reparse are not extra DFT. Each measured block is40x17action/17x17projection,
not a fabricated full Hessian; no full internal index, plane or TS selection.
Partial reciprocity1.45486e-4/1.19830e-4; max energy/gradient diagonal mismatch
.0584274/.1890415eV/Angstrom². Preserve these diagnostics; no invented error
interval or acceptance threshold. Two-amplitude comparison remains pending.

95local/95actual independent related checks pass,0failure/error/skip;
29.475/32.278JUnit seconds, tree09f4aa45e3c823f999ecfc31cea4446c402a0acc,
archive800aef31eed8801869efa1c440ad7f8e23ac28d0fdd8df2b12b34ed0618626e0,
all6455tracked files present. First53-test attempt1failure was a test-only
work-field name error, corrected without production/DFT change. Initial
snapshot shell mistakenly tested original cwd: retain its95pass output as
repeat local, not independent. Separate actual fresh-cwd imports/test run
prove the independent result. Raw native report copies SHA-match HF; all
source text normalization is restricted to code, never physical evidence.

19:42:32CST exact two original jobs28855614/28855615 remain RUNNING,32CPU
each/node10,node27. Call-time sums3509.200/3603.502s are existing E150 subset
cost, not new SCFs or allocation cost. No full E152 regression, complete
curvature/branch/prediction gate, paper/figure or reserved-geometry change.
Full goal active/incomplete. See G3_curvature_audit_E152_20261014.

## E153 — evidence-first manuscript and verified reading build (2026-10-14)

Previous E152:PROGRESS. Main Results13->5scientific jobs; E151 cross-strain
work and its chart/units enter the body. Old74table lines survive with
multiplicity; all9original figure bytes/13external citation targets remain.
Free-cell/sampling/details move to article appendicesC/D, no separate SI.
Public paper README is consolidated rather than retaining conflicting
historical incomplete counts. Missing physical error/conditional/independent
prediction evidence stays explicit; no novelty based on harmonic renaming.

96local/96actual fresh-checkout checks pass with0failure/error/skip;
24.702/26.674JUnit seconds, one existing fontTools warning each. Initial
95pass/1failure was an old editorial section anchor, corrected with the same
numeric/scope checks. Treeb17d8d566ea73aa263859a61bf37d112274fbf1d/archive
38fc8fc84b48fbcdeb7f6fc1fb846ff9517db090629542f9221263299462189d;
all6467tracked files present and actual import root verified. Nine tested
scientific/code/test/entry texts unchanged; actual8band/72record/66identity
CLI replay exports byte-identical in both roots, including manuscript SHA.

Final26page reading PDF5af5c9f64eb009d33704e859053edf05f9e2ebd0e3a544c221a85726bb53b73b
uses existing TeX Live/LuaLaTeX, no installs/rc/shell escape. All pages
rendered; four final contacts and details7/10/15/26actually inspected. All9
whole captions stay on one page; original image size and dependency bytes
unchanged. First26page split-caption and28page samepage/white-page trials
retained/rejected. Post-archive TeX-only layout change is separately compiled
and rendered, not folded into Python archive claims. Not submission-ready.

20:08:16CST originals28855614/28855615RUNNING,32CPU each,node10/node27,
50/148and49/148complete,stderr0; status receipts only, no fresh raw reparse.
No DFT/scheduler/live-source/physical-input/reference/holdout change, no H1/H2/
TS or whole E153 regression claim. Full goal active/incomplete. See
manuscript_narrative_E153_20261014 for original reports and evidence mapping.

## E154 — actual strong-reference/local-chart adapter (2026-10-14)

Previous E153:PROGRESS. Backend-independent tangent/unit/gauge transport now
connects registered mass-weighted original-T spans to local length charts,
with explicit integer basisP and immutable source directions/forces. Rank
loss is rejected; translations are normalized analysis directions, not an
ASR repair. Source/local amplitudes and the excluded-work term remain
explicit. No source frequency/old LDA Hessian becomes a PBE local curvature.

Actual2frozen E150centers x7spaces=14maps preserve ranks3/4/33 and all four
separate Cmma rank4frames. CompleteT optical33 comes from the existing matrix,
not new DFT or free preparation. Max reconstruction2.160e-16inverse-sqrt-amu,
independent atomic-work covector6.939e-18eV/(sqrt(amu)Angstrom), measured-span
defect2.341e-15. Frozen force records, not fresh raw-SCF audits. Only scalars
and source hashes exported; private source/transformed arrays not committed.

Initial helper overchecked the open third reference row as fixedT; refusal
produced no output/DFT. Corrected only to two prescribed rows; full original
manifest/source remains pinned. Final author-output-directory guard changes
only source SHA relative to valid v2; all scalar/input fields identical.
Preliminary exports retained outside Git. A zero excluded-translation norm
remains valid; no fake positive translation signal is required.

77local/77actual fresh checkout checks pass,0failure/error/skip,3.905/5.966s;
tree47a1923b9e6c918fe04c97b1d4f4dc6e16383d88/archive
094f9ab21d03a54157c6cc46a09c911575ad32e33a78de73c4ed6e28a5b0066e,
all6481tracked files present/import roots checked. Six source/test/manual
texts normalized-equal. Actual14map scalar CLI replay fields exact; only
two source SHAs differ by verified code-only LF/CRLF, never physical data.

20:43:41CST same28855614/28855615RUNNING,32CPU/node10,node27,70/148and69/148;
receipt read, not raw EFS reparse. E153 remote main now natively confirmed.
No newDFT/scheduler/live-source/input/reference/conditional/holdout change;
no fullE154suite, Hessian, H1/H2, TS or paper-readiness claim. Full goal active.
See reference_transport_E154_20261014 for actual source, scalar data and scope.

## E155 — complete first amplitude, negative transverse complements (2026-10-14)

Previous E154:PROGRESS. The unchanged E152 native auditor replayed76T/77
preserving completed probe SCFs from the same E150 allocations. All original
six-byte hashes, DSIZE32, convergence, geometry and full E/F/stress/gradient
checks passed; no failure artifacts. Reports are matching native byte copies,
not normalized data. These153prefix records overlap the earlier audits and
are not153new calculations. Both have37complete signed pairs at.01Angstrom;
each has only one complete pair at.02. No full two-step result yet.

New offline arithmetic replay separates36internal directions from prescribed
control. Raw first two internal eigenvalues: T-.276193429/-.164648089,
preserving-1.271818849/-.447449785eV/Angstrom². The normalized measured path
tangent's35-dimensional complement still has one negative direction in each
case, minima-.176131009/-.447450308. Internal gradients.108189564/.305934312
are nonzero. This is not a certified index/TS or a change to ordinary.10NEB.
Two raw negative directions forbid a positive codimension-one complement in
this matrix; no arbitrary positive Schur elimination. A34Dpositive spectral
complement is only an algebraic diagnostic, not final axes or a relaxed branch.

Initial analysis guard wrongly demanded zero translation in the native path
tangent; it refused before any output. Corrected to explicitly report the
excluded translation norm(.101688097/.007027557) and normalize the measured
internal tangent. Original forces/path/source/SCF unchanged. No physical
threshold or author-mode array added. Second-step, numerical, branch and
strong-reference independent-prediction gates remain. No automatic third
center/3Dface/cutoff change or physical+.005geometry/label. See
G3_first_amplitude_E155_20261014 for native evidence and actual test/replay scope.

## E156 — strongest references, not a weak-baseline victory (2026-10-14)

Previous E155:PROGRESS. On exactly the same E155 first-step measured internal
H/g, all seven registered reference spans at both centers now have explicit
atomic-only and joint-release restrictions. All four separate Cmma frames,
Ttriplet3, two complete optical doublets4 and fullToptical33 remain. No raw
or transformed private author arrays exported, no new DFT/SCF/holdout input.

All small-reference atomic/joint complements retain raw negative directions.
The strongest T33 reference instead leaves three raw-positive open-cell
directions(minimum2.956213824T/3.090797687preserving eV/Angstrom²), but its
retained Schur matrix still has two negative directions. Do not conflate
partition failure with failure of every fixed reference, or a positive R
with an index1TS/physical certificate. Complete harmonic stationarity remains
reference-partition invariant; no adaptive accuracy/cost advantage measured.

Behavior-preserving source-bank extraction separately committeda5bd18e:
22existing tests pass unchanged; actual old14map scalar/source fields exact,
only source helper SHA changes. New diagnostic:97local/97actual independent
checks pass; private-input14row CLI replay has exact non-code fields. First
56check attempt had55pass/1test-fixture alias failure, retained; fixture
returned a view, so its mutation also changed the measured B. Corrected only
to a copied test span. No production/input/physical-sign-floor change.

Second amplitudes still incomplete. Ordinary.10,100Ry/full10auDZP, two-center/
two-face finite scope and sealed physical+.005 unchanged. Final axes/branches,
domains and B0--B5 forecasts not selected. See
G3_reference_complements_E156_20261014 for actual evidence and scope.

## E157 — optional conditional-curvature seed and clean full regression (2026-10-14)

Previous E156:PROGRESS. Backend-independent safeguarded BFGS can initialize
from a full physical-coordinate approximate Hessian, restricted by U.T H0 U
in the actual metric-orthonormal open space. Reject nonpositive/unrepresentable
seeds before callback; no clipping, input/force changes or physical stability
certificate. None retains the identity default. Grid forwarding, geometry,
trust radius, final tolerance and existing conservative evaluator caps remain.

The exact quadratic implementation replay uses the pinned E155 H/g, not new
DFT. Callback counts129to5T/126to7preserving reach the same analytic solution.
The quadratic's two spectral axes are not final physical sampling axes.
This expected exact-H benefit is neither novel nor a measured material speedup.
Actual conditional branches and fair total-cost DFT tests remain incomplete.

131local/131independent related checks pass. Final clean tracked-code suite:
2904pass/2expected skip/0failure/error,504.106JUnit seconds,314existing warnings.
First clean archive:2902pass/2fail/2skip; two E155 assertions demanded bitwise
derived floats across BLAS threading. Separate19check replay reproduced only
those2failures; max difference3.907985046680551e-14. Derived arithmetic tests
now use absolute2e-12/relative0 while exact native bytes, metadata/types and
overwrite/drift guards remain. Original first report retained; full reports
losslessly zipped with decompressed-byte verification. Initial untracked-
workspace run was interrupted/reclassified, never counted as a passed suite.

Final tree8888fe35cf6367406e99043fde8268ea8002d9b7/archive
387fb2873074349d2cf48d775d95879ab8c285f09f1a0a1d018206acd2c72ebd,
6514tracked files present/import roots checked. Six code/test/manual texts
normalized-equal, physical evidence untouched. Independent surrogate integer
counts/statuses/input pins exact; residual floats differ at most9.2733e-15,
not byte-identical scalar reports or physical error estimates.

No newDFT/scheduler/live-source/physical-input/conditional/holdout change.
100Ry/full10auDZP/ordinary.10 and finite matrix preserved; goal active/incomplete.
See conditional_preconditioning_E157_20261014 for code, evidence and scope.

## E158 — native bounded conditional-canary worker, not yet submitted (2026-10-14)

Previous E157:PROGRESS. Prepare/preflight/run now links complete native
two-step audits, the exact original center/cache, the36Dinternal map and
the original six-input32MPI evaluator to two independently capped optimizer
arms. One origin modeling node per existing center,12reserved attempts/arm,
24total, no refund/retry/dispatch. Actual partial E155input is refused before
configuration output. Synthetic fixtures are explicitly not material data.

The case policy fixes Q0, strain/P0, released-norm target.03, step.05,
total internal norm bound.5Angstrom/minimum distance1.6Angstrom and2h/one
node32MPI. This is not a changed ordinary.10NEB criterion, a physical
stability/TS certificate, full face, forecast or measured DFT speed gain.
Parent296curvature/source preparation costs remain additional/shared.

First19checks:17pass/2failure from forwarding the identically zero control
pullback as an active36Dconstraint. Fix only the new adapter: retain40Dcontrol,
prove its exact exclusion, expose remaining nonzero pullbacks. Generic
control checks stay strict; missing/mixed/complex controls refuse. Original
report retained. Final112local/112actual independent checks pass,0f/e/s.

Treec70f94dc1a6ba537350d257e5079338c0e9fc920/archive
3143bc151627379fb5e3f6052bfe11154b72949262d9be02761e84668d4546a8,
all6527tracked files present/import roots checked; four code/test/wrapper/
manual texts normalized-equal. Native copied wrapper has identical SHA,
bash-nexit0 andRUN_DFT0exit2. First shell status capture failed due to
PowerShell expansion; read-only native Python rechecked actual exit codes,
no script change/recopy and no DFT. No new allocation or live E150 change,
no final axes/face/holdout. Goal remains active/incomplete; see
conditional_canary_worker_E158_20261014 for evidence and scope.

## E159 — completed two-step native data and real DFT canaries (2026-10-14)

Previous E158:PROGRESS. Both original jobs COMPLETED/0:0; the unchanged
native auditor reparsed296/296SCFs, original binary/six assets/ordered
geometry/EFS/gradient/SCF/genuine32MPI passed, no failures. Lossless archived
native reports pin original bytes; earlier prefixes overlap, not extra data.
Both measured steps have two raw negative internal directions per center;
the tangent complement still has one. This is not a stationary/TS/error
certificate. All seven references retained at both steps yield28scalar rows;
the rank33optical reference's three-cell release remains raw-positive, so
no claim that all fixed references fail or that unequal ranks are equivalent.

The tested015719b378Python-file runtime is deployed fresh, not over E150.
Both complete native configs pass zero-DFT preflight; the reviewed canary-only
two-negative plane leaves a raw-positive complement at both steps. Actual
same-input DFT jobs28859804preserving/28859983T use1node/32MPI/2h each,
at most two committed allocations, separate12attempt optimizer caps, Q0only.
Original intent/accepted/held validation/release receipts retained. Both first
SCFs independently reparsed:95.46s/115.13s, original EFS/DSIZE32/hashes and
gradient replay passed; all378runtime files still pinned. No DFT speed,
equal-branch, stability, TS, completed face or independent-prediction claim.

154local/154independent checks pass, including8new measured-data checks;
all6538snapshot files present, actual import roots verified. Complete-summary
CLI byte replay exact; the actual private registered reference-bank replay
has exact scientific fields, only actual source SHA/line-endings differ.
No author arrays exported. No new full-suite result, physics/cutoff/NEB.10
change, CI, release or physical+.005access. Goal active/incomplete; next
audit actual two-arm results, then proceed in the existing two-face budget.
See G3_complete_curvature_E159_20261014 for exact data/pins/receipts/tests.

## E160 — actual conditional DFT benchmark and native replay (2026-10-15)

Previous E159:PROGRESS. Both original jobs COMPLETED/0:0. Identity exhausted
12SCFs at each center; curvature initialization reached the common.03target
in2(T)/5(preserving)SCFs, residuals.02630/.01959.31new native SCFs reparsed:
original six inputs/SCF/32MPI/EFS/ordered geometry/gradient passed. No new
submission, source/cutoff change, failure refund or repeated center DFT.

Read-only strict cached replay in the pinned native environment reproduces
all four original optimizer arms' trial order/cache hits/guards/stop reasons/
final coordinates; missing/reordered points refuse without a calculator.
Last completed SCFs are not accepted iterates. Original identity convergence
costs are censored, so no12/2or12/5speed factor or full-NEB speedup.296curvature
preparation SCFs and prior G2/centers remain additional, not free.

Preserving Q0release lowers.04948493eV/cell. Its common translation-free
chart distance to the original T center decreases.36375to.15906Angstrom,
but remains.19178from the separately relaxed T candidate. Their small
.00081345eV/cell energy difference has no resolved physical error bound.
No channel-merging, phase, new barrier, stable release block or TS claim.
Proceed to the existing two local faces, not another material/Hessian center.

19new checks;173local/173independent pass,0f/e/s. Snapshot6555tracked files,
actual import roots verified; three edited source texts normalized-equal,
native files/archives byte-identical, actual structural CLI values exact.
First16synthetic replay checks had one fixture omission of final verification;
only expected callback/cache counts1→2changed, original report retained.
No new full suite, CI/release/PyPI or physical+.005access; goal active/incomplete.
See G3_conditional_canary_E160_20261014 for complete evidence and boundaries.

## E161 — bounded first5x5conditional-candidate sheets launched (2026-10-15)

Previous E160:PROGRESS. Two existing centers, no new Hessian center/material/
full band. Prospective axes[-.2,-.1,0,.1,.2]Angstrom,25nodes including reused
E160origin;24new nodes×12reserved attempts=288maximum/face.56remaining
training/8independent points untouched. Same physics and ordinaryNEB.10,
same conditional.03/radius.05/domain.5/minimum distance1.6. No+.005access.

Fixed Manhattan/lexicographic order and nearest-safe-successful-parent warm
continuation, one branch only. Budget/nonconvergence stays a hole, last trial
not promoted; actual backend failure stops remaining nodes. Recoverable
exclusive point JSONs and separate per-node conservative namespaces. A
temporarily incomplete metadata JSON is not evidence of a native failure.

16new checks,189local/189independent pass; actual import roots checked,
6565snapshot files present. Both independent prepare CLIs reproduce exact
frozen configurations. Runtime from testede3811a09bd1835f72fa00549a6fd06a200792e6e,
380Python files/10inputs exact; shell syntax/zero-run guard and native old
center/conditional-origin EFS preflight pass without new DFT. No old runtime
overwrite, pseudo/orbital/ecut/SCF drift or external-control relaxation.

Actual one-shot held admission/release yields28862319T and28862441preserving,
hfacnormal01/one node32MPI/12h each, at most two committed study allocations.
No repeated zero-point SCF or new monitor/CI/release. This is a training
conditional candidate sheet, not certified stability/global PES/TS/contour,
full-cost acceleration or frozen prediction. Goal active/incomplete; see
G3_conditional_faces_E161_20261015 for exact input/source/test/native evidence.

## E162 — completed local evidence enters the manuscript (2026-10-15)

PROGRESS after E161; commit caa159da. Complete E159 curvature, all registered
reference complements and the bounded E160 release comparison were integrated
into the connected manuscript, with fifteen source-derived rows and a terminology
ledger. Original eighty table lines, nine figures and thirteen external citations
remain. 37 local/37 independent related checks pass; the 29-page reading PDF was
fully rendered and inspected. No new DFT, source-runtime change or holdout access.
Stable branches, independent accuracy and the JCTC methodology claim remain open.

## E163 — completed-node native evidence and optimizer replay (2026-10-15)

Previous E162:PROGRESS. Keep the standalone read-only face auditor after
255 local/255 independent related checks and actual HF verification. Native
T prefix:9 points including cached origin/29 new SCFs reparsed; preserving:
6 points including cached origin/19 new SCFs. All thirteen new nonzero-Q
candidate optimizer traces replay exactly from native E/g caches. Source
identity, original six assets, DSIZE32, SCF, EFS, full40D gradients, fixed Q,
substrate rows, nearest-safe parent, reservations and final identities pass.
Original jobs are still RUNNING at those observations; current unfinished
nodes are excluded, not promoted or called failures. No new calculator call,
submission, source modification, held-out label or physical parameter change.

New report hashes and native source are in G3_face_prefix_audit_E163_20261015.
Runtime380 files are unchanged; deploy tree61c7c2e, final evidence tree46243b2a.
The two candidate faces are not yet complete, released-stable, unique branches,
valid smooth contours, TSs or independent predictions. Keep the full goal active.

## E164 — atomic versus open-cell negative curvature (2026-10-15)

Previous E163:PROGRESS (218c913b). Descriptive analysis of the original complete
training quadratics, not a new physical experiment or independent forecast.
In all four center/step blocks, the atomic restriction already has two raw
negative directions, the open-cell block is raw positive and coupling adds no
new negative direction. T's negative subspace has cell projector trace about
.0206, versus1.04e-6/3.09e-6 for preserving. T's lowest raw curvature changes
more between atomic restriction and joint matrix; small cell weight does not
imply an entire path is insensitive to prescribed substrate strain. Total
coupling norms alone do not identify the coupled unstable subspace.

Retain all four rows, the fixed length metric and measured nonzero cell-gradient
term in formal quadratic stationarity. No phonon identity, ASR/force edit,
barrier-energy partition, physical sign bound, nonlinear conditional minimum,
post-relaxation stability, TS or JCTC improvement is certified. No new DFT,
scheduler action, center, material, surface node, author array or+.005 label.

23 new checks/278 local/278 independent pass, including a counterexample where
coupling really does create an extra negative direction and invariant traces
for degenerate negative subspaces. Independent CLI reproduces the complete
analysis bytes. Tested treeaf3de69d; evidence G3_atomic_cell_signature_E164_20261015.
Updated the stale P4/P5 status board. Continue original faces/stability/independent
prediction work; do not replace that objective by this descriptive analysis.

## E165 — finite-domain agreement of original local quadratics (2026-10-15)

Previous E164:PROGRESS, e9c1de67. Compare original two-step formal conditioned
quadratics with all audited T9/preserving6 candidates, including cached Q0.
The common6 numeric-Q cohort is their complete intersection, not error-selected
or matched physical configurations. Retained eigenaxes make first-step Q-R
coupling/Schur curvature correction numerically zero by construction. The
nonzero gR energy offset is kept, not erased.

At Q=(-.2,0), even seen-origin-aligned shape disagreement is+8.5241meV/cell
for T and-20.5218for preserving (second step+8.5878/-20.7206). These are
retrospective candidate-level errors, not independent accuracy, barrier errors
or physical error bounds. Finite candidate-energy/partial-gradient secants
are not certified Hessians/envelope derivatives. Do not uniquely attribute
the disagreement to nonlinearity; branch/step/residual effects remain possible.

Keep the unmodified raw quadratic as a limited comparator, not a whole-face
predictor by assertion. Next discriminator, not yet generated/submitted:
four axial frozen-R controls per existing face within remaining56 training
slots (if admitted:52 left), original Q0 cached, independent8 unchanged.
Respect the two-allocation gate; no new center/material/strain/full path or
physics retuning. 307local/307independent related checks pass,14new case
checks; independent full analysis bytes exact. Tested treefa0fa344. Evidence:
G3_quadratic_agreement_E165_20261015. Full goal active; no new DFT this turn.

## E166 — eight matched frozen-R statics prepared and native-preflighted (2026-10-15)

Previous E165:PROGRESS,7d0c8895. Prepare original-axis four axial Q controls
per face, R=0 about each original center. Q can include cell motion; no
optimizer or fully fixed-cell relabeling. Reuse original Q0 and audited
released counterparts. Unique E166 reservation is4 slots/face, total29/81,
52uncommitted and8independent untouched. Do not rewrite old E161 policy56.

Frozen config5ffb1fe0 reproduces byte-exact in the independent CLI. Deploy
first tested tree56f11176 to a fresh E166 root,383Python/14inputs exact.
Eight real inert native decks pass unchanged six hashes and ASE roundtrip;
shell syntax/zero-run exit2 pass, no calculator directory or new SCF.
Exact-coordinate cache check rehashes28 registered matched trial JSONs,
no identical frozen coordinate found; not a global cache claim.

Native admission observes original two RUNNING jobs,ready=false. No third
allocation is submitted/held/queued. Future single compound job:1node32MPI/
1h, preserving then T, one attempt/point; existing parent normal completion
and at-most-two allocation guards precede DFT. Missing handles/non-normal
parents require review, not auto-resubmit. Negative signed differences kept.

331local/331independent pass before deployment; final333/333 include26new
checks and native/config evidence, final treefc2aa5b7, worker unchanged.
Toy call counts are not material runs. Full goal remains active; completed
controls, complete stable faces, independent errors and predictions still open.
Evidence:G3_frozen_R_controls_E166_20261015.

## E167 — completed local matrix/model tests enter the connected paper (2026-10-15)

Previous E166:PROGRESS,2936fcb0. Integrate E164/E165 into Methods2.7,
Results3.6/3.7, main Tables5/6, Discussion and the terminology ledger.
Eight rows retain both original centers and probe steps. Raw atomic
restrictions already have two negative directions; positive cell
restrictions and metric-dependent negative-space cell traces do not
certify TSs or substrate responses. Fixed-original-quadratic shape
disagreement is8.52–8.59/20.52–20.72meV/cell on the historical audited
six-Q cohort, not independent accuracy or whole-path barrier errors.
No affine gradient dropped, model fit, new reference or unique causal
attribution. Optimizer initialization and fixed-model accuracy are separate.

Preserve101 old table lines,9 image targets/bytes,13 external citations
and TeX wrapper. Final102local/102independent checks pass,13new; retain
initial single five-section expectation failure and its strict seven-section
update without deleting the original preservation checks. Tested tree5e3660bc.
Compile/render30-page internal PDF; all contacts and new-evidence pages
reviewed, nooverfull/underfull/missing characters. Not submission-ready.

RO03:03 confirmsT28862319/preserving28862441 live. T19/25candidates;
preserving14processed/13candidates, point13 caps12completeSCFs and
35/48geometry rejections, no observed SCF crash. Keep hole and continue
the original batch; no automatic retry or parameter change. Zero new
DFT/queue mutation/runtime overwrite/physicalholdout access in this turn.
E166 eight statics remain prepared only. Stable conditional branches,
native complete-face audits and frozen independent comparisons remain
required; full goal active. Evidence:manuscript_curvature_validation_E167_20261015.

## E168 — native replay identifies the guard, not an unreachable minimum (2026-10-15)

Previous E167:PROGRESS,b39db373. Reuse existing E163 auditor and immutable
E161 runtime. Audit14preserving nodes/60existing SCFs/13nonzero trajectories,
including actual point13nonconvergence; original380Python hashes and
six-input/DSIZE32/EFS/fixed-Q/substrate/budget identities pass.

Observer returns original guard decisions and verifies exact cache-only
trace/outcome. All35rejections are original total-radius short circuits;
13guard passes/12completeSCFs, no observed backend crash. Last radius
.4998853735A, g_Rnorm.0735963713eV/A, radial-.0136123401 and constant-radius
tangent.0723265516. Minimum observed norm.0337907153still exceeds.03;
continue through rebound until original12-attempt cap, no last-trial promotion.
Rejected radius steps were not distance-tested. No boundary-KKT, impossible
domain, stable-branch or unique-cause claim follows. Specific next hypothesis:
feasible-direction step handling at identical start/domain/target and measured
cost, not enlarging the live radius or changing physical settings.

41local/41independent pass before native execution; final50/50 pass with21new
checks and native/arithmetic/table evidence. Independent full analysis bytes
exact,a4adb825; native report7c9581eb. Final tree3bef044d. No new DFT, scheduler
write, source repair, physicalholdout access or extra training reservation;
paper bytes unchanged. Original faces continue; E166 controls retain resource
priority. Full stable-face/independent-prediction goal remains active.
Evidence:G3_bound_diagnostic_E168_20261015.

## E169 — opt-in ball-feasible conditional solver, analytic controls only (2026-10-15)

Previous E168:PROGRESS,0d65557b. Add a separate backend-independent single-start
entry using existing coordinates/metric/control checks/result types. Project
quasi-Newton proposals into the declared released ball; non-descent proposals
fall back to projected steepest descent, with actual-displacement Armijo and
the original physical guard. Success still uses full released gradient; no
boundary-KKT promotion, start repair, negative seed clipping or new default.
Released radius is explicitly distinct from total original displacement.

Four frozen analytic controls, both arms with common start/seed/.5radius/
.03target/.05trust/12unique-call cap:4/4 interior,7/6 near-wall,9/9 nonlinear;
exterior11/8both no candidate. Keep all results/no-gain controls, no speedup
ratio for failures and no DFT/HfO2/new-theorem claim. Analytic ASE E/F/stress
budget test verifies no extra factory call after the existing cap.

218local/218independent related checks pass,19new,10modules not full regression.
Independent full output bytes exact20c5d743; tested tree4b1a61ca. Existing
mode_surface/default and19E166scientific code pins remain unchanged. Add API
doc/manual link, not a production promotion. Zero new DFT/native reparse/
scheduler write/training reservation/physicalholdout access/runtime mutation.
Actual matched bounded material experiment and final stable-face/independent
strong-control evidence still required. E166 controls retain resource priority.
Full goal active. Evidence:feasible_steps_E169_20261015.

## E170 — real matched arm frozen; native start/cache/novel proposal checked without DFT (2026-10-15)

Previous E169:PROGRESS,f0c9cce1. Fix original preserving point13, exact requested
start/actual first SCF,parent5,first-step36Dseed,.05trust,.03full-gradient target,
.5total-radius/1.6distance guard and original six inputs. Released radius
.4472135954999579 derives from the same fixed Q/metric; no domain expansion.
Hash-pin13lossless original JSONs,cacheZIP02cf4f9b. Exact physical-cache hits
only; same-contract delegate keeps separate conservative budget/error handling.

Reserve one preserving training slot:30/81committed,51left;T29/81,52left;
all8independent slots/face unchanged. Future12new-attempt/1node32MPI/1h arm
must include old12SCF/cache/curvature/path preparation in net cost. Entry has
prepare/probe only, not a material runner or submission. Execution/allocation/
cost manifest still needs freezing; E166 statics keep priority.

Configf937a2e3 independently reproduces exact bytes. Native old380source hashes,
source/origin EFS and12original pointSCFs reparse; old baseline strict replay
passes. New solver reuses7old callbacks including exact original start, then
requests a new coordinate:totalnorm.49999999999999994A,mindistance2.0505811584A,
fixedQ/control/substrate and actual STRU/sixhashes pass. Inert deck only;
no invented E/g, no DFT or benefit claim. Native report2c22c66c.

233/233checks pass before native;final234/234pass,16new/11related modules,
not full regression. First transfer's POSIX/Windows path mismatch failed
before SCF/queue state; fresh-r2 validates all payload before archive creation.
Prototype coordinate-hash label corrected before final freeze; drafts retained,
physical inputs unchanged. Actual E166 admission observes twoRUNNING originals,
readyfalse/no third job. RO03:57T24/25,P17processed/15candidates.
Full stable-face/independent-prediction goal active. Evidence:G3_feasible_material_E170_20261015.

## E171 — original T batch audited; statics admitted and matched-arm execution frozen (2026-10-15)

Previous E170:PROGRESS,a12abd90. T28862319 completes normally:25/25gradient
candidates including cachedQ0,118newSCFs. Original native audit reparses all118
and exactly replays24new optimizer traces, preserving original380source bytes,
six inputs,DSIZE32,EFS,Q/control and per-node budget identities. This remains
candidate convergence, not released-stability, global surface or TS evidence.

After the real slot opens, the existing8static E166 controls are submitted
once as28865656, held/resource-validated and released. At04:27CST it is pending
underAssocGrpCpuLimit; preserving28862441 is live at19processed/17candidates,
103completeSCFs and no observed backend failure. Preserve the two original
holes and all original source. No third committed allocation is created.

Add the E171 execution supplement for the already reserved E170 point13arm:
same original start/seed/Q/target/.5totaldomain/.05trust/sixinputs, fresh12new-
attempt cap,1node32MPI/1h. E170's reservation is not charged again. Exact first
callback reuse is checked before any DFT. Separate old hits, new backend hits,
new SCFs, failed/reserved and denied requests; capped baseline is censored,
not a speedup denominator. All preparation costs remain additional.

Native external-bundle preflight passes:380oldsource files,12oldpointSCFs,
strict old replay and7cache hits before novel query; no new optimizer DFT.
E166 completion/all8native statics remain required before admission and run;
that true branch has not yet been exercised. Manual admission correctly
returns not-ready, no intent/job/output and zero scheduler writes. No timer.

250local/250independent checks pass,16new/12related modules, not full regression;
independent prepared execution bytes exactb80dc015. Retain initial portability/
synthetic-test failures, incomplete E-drive extraction and raw-JSON archive
conversion failure. Fail before a mismatched starting SCF, never approximate
physical cache coordinates. Preserve exact native auditor bytes; case-specific
Git attributes keep raw JSON/XML hashes. Final independent archive onC passes
without weakening raw-hash checks. No local file deletion, physical-input
change, live runtime overwrite, holdout access or manuscript promotion.

Final tested treea22be6ec. Execution preflight1ffb76f6; completeT auditb8e78367.
Full stable-face/independent-strong-control prediction goal remains active.
Evidence:G3_feasible_execution_E171_20261015. Next:process the original8statics
when genuinely complete, audit them, then admit the bounded matched arm.

## E172 — completed-sheet consistency analyzed and integrated, no new sampling (2026-10-15)

Previous E171:PROGRESS,6472844b. Reuse the complete T audit, all25candidate
points/118sourceSCFs and cached released origin. Analyze40nearest edges and
16elementary loops using the existing sampled_path_work. Direct retained/
released formulas and independent full40D atomic/open-cell/control work agree;
fixed-control work is zero. No new physical result, fit, interpolation or point.

Final T radii range.04579–.46378A, inside the old.5A domain; released gradients
range.006525–.029784eV/A. Energies span-38.530500to0meV/12atomcell relative
to cached origin. Largest retained-only endpoint energy/work difference1.209125
meV/cell contains signed released work1.093874, leaving.115251full discrepancy.
The maximum full discrepancy across40edges is1.039406. Largest absolute retained/
full loop circulations2.000384/.578259occur on different loops; raw loop energy
changes telescope to zero. No claim of nonconservative DFT forces, a branch jump,
or a measured physical/interpolation/barrier error follows.

The endpoint Cauchy inequality bounds only released endpoint work. An explicit
smooth cubic with zero endpoint derivatives and finite energy change prevents
it from becoming a false energy-error bound. Finite residuals remain visible;
the original .03conditional/.10ordinaryNEB targets and all electronic inputs
are unchanged. No released-stability or exact-envelope derivative certificate.

Add51lines to the main manuscript's existing Methods/Results/Discussion and
update detailed methods/terminology. Preserve all old nonblank lines in order,
113table lines,9figure bytes,13external citations and TeX wrapper. No new
abstract/conclusion, optimizer advantage or independent prediction claim.
Existing LuaLaTeX compiles31pages; all pages rendered/contact-reviewed, new
methods/results/discussion pages7/18/19 examined in detail. Nooverfull,
underfull or missing characters. This remains an internal incomplete draft.

83local/83independent related checks pass,23new/7modules, not full regression;
full independent analysis bytes exact7696aeade3. Tested tree055ec328. No new
DFT/nativeSCFreparse/scheduler write/training reservation/holdout access.
RO04:50confirms preserving28862441live,20processed/17candidates,115SCFs/no
observed backend failure. The admitted8static controls28865656now RUNNING,
3/8complete onnode6; no duplicate submission or third allocation. Full stable-
branch/independent-strong-control prediction goal remains active. Next native
audit and analyze these controls, then bounded matched-arm admission only.
Evidence:G3_sheet_consistency_E172_20261015.

## E173 — eight frozen statics audited; matched arm truly running (2026-10-15)

Previous E172:PROGRESS,30b7da58. Controls28865656 normally complete8/8onnode6,
1021allocation seconds/32CPUs. Existing E171 native admission reparses all8
E/F/stress/sixasset/DSIZE32/one-attempt records plus12old point13cache SCFs,
retains original380/383source snapshots and strictly replays the old capped
baseline. No new parameter or training reservation. Static call proxy8.136814
core-hours; allocation9.075556;296curvature and prior path/candidate costs extra.

After normal completion, the already frozen one-shot driver submits28866132
once, validates held1node32CPU/1h resource identity, then releases it. At most
two committed study allocations, original point13start/Q/seed/.5domain/.03full
released target/.05trust and12new-attempt cap. Old cached data not recomputed.
This is an actual running material test, not merely an inert input package.

Matched four-axial-point ablation keeps both centers, both raw probe steps and
all16model-point comparisons. Maximum frozen/released shape errors are
P:.075880/4.671404and.116535/4.720758meV/12atomcell;T:.416264/1.327177and
.403160/1.346071. Exact signed pointwise decomposition separates frozen error
from missed release-shape change. Pvariation4.606959,T1.168761; no coefficients
fit, sign clipping or independent-accuracy/microscopic-cause claim. Maxima
are not subtracted across unmatched points. This is not the full strong-
reference/control matrix or a whole-path barrier error measurement.

First new SCF independently audited. Energy/gradient/stress/Q/domain match,
but initial read-only audit falsely rejects validation metadata because
`bundle/../arm` and `arm` spellings differ. A diagnostic also encounters NumPy bool
JSON serialization. Resolve only the exact expected real native directory
inside permitted arm before reproducing validation through original spelling;
all fields still compare exactly. Regression rejects other directories and
root escapes. No raw result, physical tolerance, live worker or job is repaired,
stopped or restarted. First geometry matches the E170 novel proposal exactly;
radius.49999999999999994passes the original guard, no roundoff loosening.

RO05:11shows28866132RUNNING,4complete/5reservednew attempts, no backend failure,
no final arm summary. First released norm.069855is not convergence; continue
under original cap, do not stop on a rebound. Preserving28862441also live.
57local/57independent related checks pass,26new/3modules, not full regression.
Full independent analysis bytes exact7d886069; tested tree4cf5a3a7.
No ablation DFT or extra reservation; physical+.005and all8independent slots/
face untouched. Main manuscript f23742d8unchanged this turn. Full stable-branch/
independent-strong-control prediction goal remains active. Next inspect the
bounded arm's actual terminal outcome and finish/audit the original preserving
face; do not promote this launch/first-SCF check to acceleration evidence.
Evidence:G3_frozen_release_E173_20261015.

## E174 — manuscript ablation and terminal matched-arm negative result (PROGRESS)

E173completed frozen/released ablation now enters connected Methods2.7 and
Results3.7/Table7. Both centers, both probes, all16model-point records and
signed pointwise identity retained; no subtraction of unmatched maxima,
fit, clipping or independent-accuracy/unique-cause claim. All113old table
lines/9figures/13external citations/wrapper preserved;6table lines added.
Old preparation-only prose/test superseded by actual completed-native proof.
32page LuaLaTeX PDF contact-reviewed throughout, newMethods/Table7at higher
resolution. One long coordinate math-span overflow fixed without data changes;
final0overfull/underfull/missing-character warnings. Existing5user edits
unchanged. Nature-writing evidence logic separates frozen and release-response
discrepancies, rather than treating a favorable cancellation as model accuracy.

28866132COMPLETED0:0,1417s,node6/32CPUs/hfacnormal01. Exactly12new SCFs,
no backend failure, no gradient candidate before original cap;7old cached
hits precede new calls. Best new trial norm.04585524above.03, old12-trial
minimum.03379072. Trial minima are not accepted iterates or an equal-cold-
cost benchmark. New call-proxy/allocation12.049803/12.595556core-hours;
all old12/296curvature/other preparation remain additional. Censored old
baseline cannot supply a speedup denominator. No default promotion,
extended budget, restart, new allocation or new training reservation.

Reparse12new+12old native E/F/stress/assets/DSIZE32/contract/reservation;
both original-order cache-only optimizer replays pass. Preserve initial audit
failures: wrong external helper source path; NumPy bool diagnostic print;
one strict basis-roundtrip norm.5000000000000001at originalfullnorm.5;
unset-login-thread exact replay mismatch; absent optional threadpoolctl;
standalone solver relative-import name. No calculator/result/physics repair.
Auditor instead requires replayed original internal proposal's unchanged
strict guard and bitwise full embedding, with no epsilon/inward move. All
native fields stay exact; roundtrip failure retained, one-ULP/outside-proposal
regressions reject. Original OMP/MKL/OpenBLAS1 and real package import restore
both replay sequences. Original380runtime/9bundle files rechecked.

112local/112independent related checks,27new/8modules, no failures/errors/skips;
not full regression. Archive tree63365093, raw material proof and new-cache
bytes identical in independent extraction. Final PDF29f008cb/mainc6b33de6.
RO05:37:54preserving28862441RUNNING,22/25processed,17gradientcandidates,
holes13/14/19/20/21,139completeSCFs; no final-face/stability claim. Keep
original job alive. T29/81,P30/81 and8independent slots each unchanged,
physical+.005sealed. Next complete/audit preserving face, choose training-only
valid domain/representations and test frozen strong-control predictions;
no unregistered third center/new material/full path. Full goal active.
Evidence:G3_release_manuscript_E174_20261015.

## E175 — actual B3 atomic-only release treatment submitted (PROGRESS)

Prior turn PROGRESS_E174_639b1847. Current tree, finite plan, v1/v2prediction
protocol and exact28862441handle rechecked before action. No /jobs or
/autoresearch slash expansion is available; use real Slurm and durable logs.
ABACUS native audit workflow applied with HF-only restriction, not235defaults.

Within the original two centers/faces, reserve5atomic-release training-node
treatments per face: origin and the four existing axial .1A Q points. No
new Hessian center/material/strain/full path/independent label. T29→34/81,
P30→35/81, remaining47/46;8independent each and physical+.005sealed.
Same coupled Q/frozen seeds;31released atomic directions have exact zero
cell rows and are Q-orthogonal. Cell is fixed to each Qseed, not one common
cell across Q. Original full36D Hessian seed/.05trust/.5totalradius/1.6minsep,
.03Euclidean atomic-released target; not a fulljoint-gradient or stability
certificate. Original safeguarded BFGS, not failed E171feasible optimizer.
Eachnode12new attempts includingfailures,60perface/120total. Cachecenter/
8frozen seeds exactly; firstcallbackmiss refuses beforeDFT. Unknown/backend
failures propagate, no lasttrialpromotion, single-rebound stop or signclip.

R1native preflight rejects regenerated frozen embedding bitwise before
DFT/submission. Preserve it. In fresh r2bundle, native map explicitly uses
literal original full seed plus B*(u-u_seed). This removes cross-machine
matrix-product roundoff from cache identity without relaxing cache matching,
Q arithmetic threshold, geometry domain or changing any old seed/input.
1ULPbasis-variation regression keeps seed exact. Development configs moved
to Cscratch for provenance, no raw calculation or live source overwritten.

Corrected native preflight reparses all8original frozen EFS/assets/DSIZE32,
oldcenters/origins and380runtime files;10inert STRU/input decks pass exact
fixedcell/control,17digit geometry roundtrip and6physical hashes. One external
immutable5file bundle; old E161runtime unchanged. Local/independent73related
checks each across4modules,16new,0failures/errors/skips; not full regression.
Finaltestedtree79338ebb/archive9c0d6500, material/config/source bytes exact
in independent extraction. No private author arrays exported.

One-shot job28867177submittedheld/validated/released:1node32CPUs/1task,
hfacnormal01,mpirun32,OMP/MKL/OpenBLAS1,5hcap, two centers serial. Intent is
saved before sbatch; ambiguous intent never retried. At most2committed
allocations includingheld/pending. RO06:00:19confirms PENDING due
AssocGrpCpuLimit,0newSCFs/points, no backend failure; this is a verified wait,
not terminal or a reason to duplicate/stop other jobs. Original preserving
28862441RUNNING,23/25processed/17candidates; full face not yet present.
No new B3result/convergence/acceleration/independent advantage is claimed.

JCTC mainc6b33de6unchanged, finite plan updates P4/P5 only. This is the
registered release-treatment factor, not completed T/Cmma/B0--B5 matrix.
Next audit actual atomic candidates and compare signed EF→EA→EJ at all
five Qs/centers with costs; complete original preserving face without
interpolating holes. Stablebranches/finaltraining-onlyrepresentatives,
strongcontrols/independentpredictions and numericalresolution remain.
No install/otherhosts/CI/newmonitor/release/PyPI/journal submission.
Full goal active; evidence:G3_atomic_release_E175_20261015.

## E176 — pre-label B3 expectations and exact read-only audit (PROGRESS)

Previous goal turn PROGRESS_E175_785b93bd. HEAD/tree and exact28862441/
28867177handles rechecked. Jobs/autoresearch slash commands unavailable;
use real Slurm and durable artifacts, ABACUS native-audit workflow with
HF/hfacnormal01only override. No newDFT/submission/restart/reservation.

Freeze original-center B3quadratic expectations before any new atomic-only
result:2centers×2originalsteps×5existingQ=20model-point rows, no coefficients
fit, no favorable step, raw affine offsets retained. Formal atomic/joint
origin releases(meV/12atomcell):T.01−.389345/−.866958,.02−.397438/−.874429;
P.01−50.387439/−54.521222,.02−50.515678/−54.705828. These are model values,
not actual newDFT or the protocol's independent validation. Frozen/joint
labels were already seen; all atomic-treatment labels remain unread.
HF06:11:35freeze binds code/forecast hashes and confirms0newSCFs/points.

Signed EF→EA→EJ identity preserves positive release and separate origin
offsets; missing atomic origin abstains only for shape, not known absolute
differences. No clipping, unique mode-energy or global nested-minimum claim.
No source author arrays or physical+.005geometry/label accessed.

Read-only B3auditor binds original handle/config/5filebundle. Native branch
checks EFS/stress/6assets/DSIZE32, exact seed cell/Q/control, contract and
12attempt reservations; recorded path alias must resolve to the exact
expected calculation. Cache-only unchanged-worker/optimizer replay requires
exact request sequence, candidate/counters/source identity. Original affine
anchor/.5guard/.03atomic norm retained; joint norm separately reported.
Partial pointJSON is pending, no last-trial promotion or restart. Complete
mode additionally requires normal original completion/all10points/summary.

Actually exercised HFprefix has0completed nodes/0newSCFs: newjob PENDING
AssocGrpCpuLimit. Thus native completed-node/replay branches are prepared
but NOT certified on new material data. Synthetic quadratic replay is
explicitly synthetic. First18check run17pass/1fail because fixture backend
double-counted final repeated coordinate unlike real evaluator cache; add
exact duplicate caching in fixture, no production gate weakened.

75local/75independent related checks across4modules,18new,0fail/errors/skips;
not full regression. Full original20expectation bytes independently exact
321ba700. Testedtreeffdaf77d/archive32afa8c7; freeze941412d4, prefixb7027bd8.
5protected user edits and mainc6b33de6 unchanged; no running source overwrite.
RO06:16:19PjobRUNNING,23/25processed/17candidates, holes13/14/19/20/21/22;
newB3PENDINGCPUquota/0newSCFs/0points/no failure, bothrealhandleslive.
T34/81,P35/81 and8independent each unchanged. Next audit actual B3and
finish original preserving face, compare all signed treatment differences
and unchanged expectations with cost/coverage; no new material/Hessian/
fullpath/CI/install/release/PyPI/submission/monitor. Strong-reference/
independent advantage, releasedstability and final domain/representatives
remain missing; full goal active. Evidence:G3_atomic_audit_E176_20261015.

## E177 — close Fan publication visual-access gap (PROGRESS)

Previous goal turn PROGRESS_E176_d5ff34e5. Exact HF handles unchanged:
RO06:27:00CST original preserving28862441RUNNING/24processed/17candidates,
no face.json; B3atomic28867177PENDINGCPUquota/0newSCFs/points/failures.
No scheduler write, new computation, reservation or holdout access.

E070's primary-source visual gap is closed: normal publisher HTML/PDF
fetch and author-hosted published PDF, exact same PDFbytes/d28d955d,
11pages rendered; scientific main pp.1–9 read, details pp.4/6/7/8 inspected.
Explicit comparison matrix distinguishes force-training parity/energy MAE,
PO-based prescribed lattice/our partial Tclamp, intermediate-phase proxy/
endpoint null and interface-energy/activation units. Existing E108 Lee SI
work reused, not counted again. Publisher page lists no verified SI;
this is not proof of universal SI absence. No private author array export.

Main introduction already covers the pertinent boundaries, so unchanged;
no code edit, new pytest/compile claim or material-accuracy/TS claim.
Documentation/source hashes and local evidence links checked. Five user
edits preserved,100Ry/full10auDZP/.10NEB unchanged. T34/81,P35/81 and8
independent each unchanged. Next actual original-face/B3native audit and
signed treatment/cost comparison; strong controls, stable branches and
independent material predictions still needed. Full goal active/incomplete.
Evidence:outputs/HFO2_FAN_VISUAL_BOUNDARY_AUDIT_2026-10-15.md.

## E178 — executed signed frozen/atomic/joint comparison (PROGRESS)

Previous goal turn PROGRESS_E177_3a345809; original HFhandles and5user
tracked edits rechecked. Jobs/autoresearch slash commands unavailable;
use actual Slurm and bounded shell/test workflow. No DFT/scheduler write,
new training reservation, geometry generation or physical+.005access.

Implement and execute a hash-bound three-treatment postprocessor with
all2centers×2originalsteps×5registeredQ=20rows, no fitted coefficient or
winner selection. Reuse sealed E176expectations and native-audited E173
frozen/joint observations; atomic supplied prefix is explicitly empty.
Missing A retains known J-F; missing only Aorigin retains absolute A
changes and J-Fshape. Preserve signed positive differences and separate
origin bias, enforce pointwise treatment/error identities. Statistics list
actual axial IDs/denominator4; empty errors null, not perfect0; no partial
cohort ranking. Atomic norm is separate from joint norm/stability; processed
hole differs from pending. Source accounting/time and audited-prefix-only
cost explicit,296curvature and original joint/preparation additional.

Actual output reproduces old four axial F/J training comparisons; no new
atomic material evidence. For preserving.01, max frozen/joint shape
discrepancy .075880/4.671404meVcell, missed release term4.606959; not a unique
mechanism, physical uncertainty or independent accuracy. Bothsteps retained.

First related104checks103pass/1fail; old75alone74pass/1fail. Old E176formal
replay wrongly required bitwise float equality. Diagnostic124floatdiffs,
max2.842170943040401e-14; nonfloatmetadata/sourcepins identical; exact
low-level cause unproved. Sealed original forecastSHA321ba700unchanged.
Test now pins its bytes and compares only model floats at1e-12abs/0rel,
exact types/metadata; rejects1e-9andmetadatachanges. This is not physical
resolution; raw input/native optimizer/cache/gradient gates unchanged.
Preserve both failed XMLs and diagnostic; no silent golden-file update.

Final106local/106independent related checks over5modules,29newcomparison
plus2newnegative formal-replay checks,0fail/errors/skips, not full suite.
Testedtree5b710adf/archive d1d6ac49; all20report rows independently exact
e467df51. RO06:39:44CST originalP28862441RUNNING/24processed/17candidates,
no face.json; B3atomic28867177PENDINGCPUquota/0SCFs/points/failures. Main
manuscript and6physical assets/settings unchanged. T34/81,P35/81,8independent
each unchanged. Next actual complete original-face/B3native audits, signed
A-F/J-A comparison and cost/coverage; full strong controls, releasedstability,
final training-onlyrepresentatives/domain and independent predictions remain.
Full goal active/incomplete. Evidence:G3_treatment_comparison_E178_20261015.

## E179 — original preserving batch complete and all-hole diagnosis (PROGRESS)

Previous goal turn PROGRESS_E178_bc185c1d. Exact original handles/worktree
rechecked. Original28862441normalCOMPLETED0:0,21951s,node321/32CPUs/
hfacnormal01;25points processed,17gradientcandidates incl cachedorigin,
8holes13/14/19/20/21/22/23/24. No resubmission or source repair.

Use ABACUS native-audit skill with explicit HF-only override, not its235
examples. Original E163complete audit verifies380sourcefiles, original
config/job/cache, all172newSCF E/F/stress/sixassets/genuine32MPI/reservations/
fixedQ/substrate and24exact original optimizer replays. Cachedorigin not
recomputed; failure paths empty. Nativeproof982d2657, facecb97225d. Cost
193.859072782completed-call proxy/195.12allocation core-hours; original
curvature/canary/G2 additional. Not complete stablePES or TS certification.

Read-only newdriver rechecks exact nativeproof and every used point/result
hash, then reuses unchanged E168observer on all8holes. No duplicate native
EFSparse or DFT. Original callbacks/stopreasons/guard outcomes restored and
unchanged. Sevenholes12calls,point23ten,94trials total,nonebelow.03. All436
recorded rejects occur at old.5A radius short circuit; rejected proposals
notdistance-tested. Lasttrial radii.49971--.50000A establish domain contact,
not nonexistence,KKT or physicalbranch termination. Completeholeproof
d7c89c51; newdrivera3c277b5. No fillingholes, enlargedradius or capped-arm retry.

Writing skill used for bounded Results/Discussion evidence, existingterms
retained. Mainpreface/3.7 now report coverage/cost/limit; no successful
prediction or causalmechanism paragraph.119oldtablelines,9figures,13external
citations preserved. Existing LuaLaTeX32page build passes0overfull/underfull/
missingchars; all32contacts and pages1/19detail inspected. Main38769b36,
PDF12dbe776 in Cscratch; no figurecreation or new publication action.

Initial local63pass/1fail was incorrect newtest assumption allholes12calls;
correct expectedactualpoint23=10, do not mutate records. Firstindependent
65pass/1fail from Git CRLFconversion of unchangedTeXwrapper. Exact normalized
source identical; add one explicitLFattribute, no raw-hash test weakening or
wrappercontent edit. Both failureXMLs retained. Final66local/66fresharchive
related checks over3modules,16new,0fail/errors/skips; not fullregression.
Testedtree88510e04/archive4249b066; both nativeproofs/driver/main/wrapperbytes
survive exactly. Independentcodecheck is not independentmaterialvalidation.

Fiveusertracked edits remain exact; physical100Ry/full10auDZP/.10ordinary
NEB, .03conditional target, original two rawHessian centers unchanged.
T34/81,P35/81,8independent each and physical+.005sealed. No newDFT,
SCFschedulerwrites/reservations/otherhosts/CI/release/PyPI/monitor/submission.
RO06:59:50CST atomic28867177PENDINGCPUquota/0SCFs/points/failures. Next
actualatomic treatment evidence, stable training-onlydomain/strongcontrols,
prospectively frozen independentpredictions; full goal active/incomplete.
Evidence:G3_preserving_terminal_E179_20261015.

## E180 — complete-reference equality and shared-Q control limits (PROGRESS)

Previous goal turn PROGRESS_E179_4fb42d95. Actual28867177and protected
worktree rechecked; jobs/autoresearch unavailable as slash commands, use
real Slurm and bounded shell/benchmark evidence. No new physical task.

Add backend-independent analysis-only affine_quadratic KKT utility:
minimize rawg.x+x.H.x/2 within orthonormal allowedU, subjectto A.T x=q.
Nonzero affine gradient and observable scale retained; full constraint
rank and strictly positive allowed-null block required. No unstable
pseudoinverse/clipping/fallback, no caller-array alias or default change.
Generic analytic tests compare independent Schur/KKT solutions, basis
invariance, gradients, clamped omission, emptyfree and invalid blocks.

Actual posthoc experiment:2oldcenters×2steps×8formulations×25originalQ=
800scalar model-point records. Reuse full E159H/g, E171/E179complete25/17
candidate labels and all7registered sourcebanks including Cmma103/106/114/127.
CompleteT33+3cell spans same36space as adaptive2+34complement, identical
Qconstraints. FullT/KKT vs adaptive/Schur maxenergy2.6368e-16eVcell,
coordinates1.5588e-14A. Standard same-problem identity, not new theory or
prediction gain. Both retain T28.777--28.872/P20.522--20.721meVcell seen-
origin shape discrepancy on24/16nonorigin candidates; newTmax uses all25
seenpoints, not earlier6pointcohort. Original biases retained; no fit.

Small-reference control clamps all omitted directions and releases only
allowed-null combinations; not E159's completeoutside-reference release.
Allrank4T/Cmma and rank33strongcase kept, no dimension/frame/step truncation.
Allreference common nonorigin geometryeligible setempty. Raw extrapolated
values remain diagnostics, not poor material predictions; null common
statistics, no subset ranking. This is compatibility with earlier coupled
Qaxes, not universal failure of native Cmma/Tcoordinate models. Complete
native-chart/independent strong comparison remains, not declared done.

81local/81independent related checks over6modules,29new,0fail/errors/skips;
notfullsuite. Sameauthorized privateinput CLI reproduces complete800row
reportbytes39e31d79. Testedtree40bf6241/archiveae30bf74. Authorfiles checked
before/after, no raw/transformed basis/response arrays or geometries exported.
296curvature and118/172conditionalSCFs/Gamma/canary/G2/other preparation
additional, no freepreparation/equalcost/speedupclaim. Exploratory metadata
key corrected beforefinalreceipt; earlierdevelopmentreports stayC, original
prospectiveE176expectations unchanged. Main38769b36unchanged/no newPDFclaim.

Fiveuseredits preserved, originalproductionruntime/100Ry/full10auDZP/
.10ordinaryNEB/.03conditional/.5radius unchanged. T34/81,P35/81,8independent
each and physical+.005sealed. No newDFT/rawSCFreparse/schedulerwrites/
reservations/otherhost/Hessiancenter/material/wholepath/CI/install/release/
PyPI/monitor/submission. RO07:21:40CST atomic28867177PENDINGCPUquota,
0points/SCFs/failures. Next actualatomic treatment and nonlinearstablebranch/
strongnativechart controls followedby frozenindependent predictions, not
another linearbasis advantage claim. Full goal active/incomplete.
Evidence:G3_shared_Q_controls_E180_20261015.

### E181 2026-10-15 — PROGRESS: explicitly budgeted released-base stencils admitted

Previous E180 PROGRESS/3dd0f095. Reported the missing actual post-release34D
curvature gate, rather than treating raw-peak Hessian/optimizer secants as
its proof. Frozen explicit budget amendment: existing twoE160Q=0anchors,
136new derivative attempts each272total, both signs/.01/.02steps, original
296andallpreparation additional. Twoadditionalverificationbases disclosed;
no thirdchannel/material/referencephonon/wholepath/control condition.
Not12optimizercalls/point, no newQslots, T34/81 P35/81 and eightindependent
each unchanged. Literalnative anchors reused; full literal probes/frame
frozen, no cacheidentity weakening or .5radius/domain repair.

All272safegeometries/nativepreflight with original380runtime and exact
EFS/sixassets/binary passed,0newSCF. Largestprobe norm.3187221853A,
minseps2.041856089/2.065504072A.86local/86independent relatedchecks over5modules,
23new,0fail/errors/skips, notfullsuite. Actualfrozenconfig, twoamplitudes/
negative signs/guard refusal/resource policy and all272syntheticcalls plus
firstfailurecharged/no retry tested; mocks notmaterialDFT. Testedtree4db6c165,
archive62e780c5. Cscratch retains initial83check/devconfig, no failedtest hidden.

Nativepreflight9630d041, config71a2e1b8, immutablebundle2e61cfc1.
Once-only held/validated/released actualSlurm28868827 at07:43:56CST:
hf/hfacnormal01,1node32CPUs/genuine mpirun32/threads1/12h. Original28867177
counts pending, two committed total; no otherjobs cancelled or sourceoverwritten.
RO07:45:51CST bothPENDING/AssocGrpCpuLimit; no new
completedDFT.9.5h approximate after allocation only, queue delay unknown.
Native32rank behavior in newrun still awaits first SCF, not implied by tests.

Both raw34x34and40x34actions, reciprocity/symmetrization/twostepdifference
required; no clipping/pseudoinverse/ASR/signcertificate. Positivefinite-step
spectra and .03residual notphysicalerrorbound/globalbranch/TS. OtherQnodes
and eight preservingholes remain uncertified. Fixed100Ry/full10auDZP/
Gamma2^3/P0/E0/.10NEB/.03conditional/.5radius unchanged. Fiveuseredits,
E176forecasts and manuscript38769b36 preserved. No+.005/privateauthorarray
access, CI/install/release/PyPI/monitor/otherhost/submission. ActualB3,
strongnativechart/independentmaterialtests and fullgoal remain incomplete.
Evidence:G3_released_curvature_E181_20261015.

### E182 2026-10-15 — PROGRESS: manuscript fair-reference null and finite-domain evidence

PreviousE181PROGRESS/1f50ddc8. Native-complement checks already existedE159,
not repeated. Integrated executedE180 same-observable controls into2.7/new3.8/
discussion/D.7, distinguishing complete-reference KKT/Schur equality from
prediction gain. Standard constrained quadratic, finite affineg and raw
null-stability criterion explicit.2.64e-16eVcell/1.56e-14A numerical identity
notDFTerror. Full-space models share training discrepancies T28.776664/
28.872042,P20.521839/20.720640meV12atomcell over24/16nonorigin candidates.
Sevenrefs/fourCmma retained, no common nonorigin eligible cohort, no ranking
ownsubsets/extrapolation or universal nativeCmma failure. Strong native-chart
and independent material tests not declared complete.

12new quantitative rows pinned to39e31d79source, original119rows/9figbytes/
13externalcitations/TeXwrapper preserved.113local/113independent relatedchecks
over9modules,20new,0fail/errors/skips, notfullsuite. Threeeditorialtests updated
only for eighthResult/two legitimate tables; source checks not removed.
FinalsourcecheckCLIexactb3fb3b65. Finaltestedtreeb745699e/archive5dbe91fc.
Mainraw92f15239. ExistingTeXLive/LuaLaTeX34pagePDF9d52dba9, allpagesrendered/
inspected, method8/results20–21/appendix34detail. Table6/8caption splits identified/
repaired withoutnumber/wrapperchanges; firstPDF/metadata/113successfulchecks
and earlierarchives retainedC. Final22pagePNGidentical, all12changedpages
re-inspected;0overfull/underfull/missingchars. Not submission-ready.

RO07:59:54CST28867177/28868827bothPENDING/AssocGrpCpuLimit,0newSCFs/labels.
No duplicate/cancel/sourceoverwrite/DFT/reservation/schedulerwrite/private
authorarray/physical+.005access/CI/install/release/PyPI/monitor/otherhost/
submission. Fixedphysicalinputs/.10NEB/.03conditional/.5radius unchanged,
fiveuseredits preserved. ActualB3/postreleasecurvature/stablebranch/resolution/
strongcontrols/frozenindependentforecasts remain. Fullgoalactive/incomplete.
Evidence:manuscript_shared_controls_E182_20261015.

### E183 2026-10-15 — PROGRESS: released-stencil native reader and complete-batch refusal

PreviousE182PROGRESS/f9bd3497. Added read-only E181prefix/complete/smoke
auditor in fresh external HFnamespace, frozen worker untouched. Bind exact
E181bundle/config/accepted/held/preflight, old380runtime/sixassets/binary.
Exact ordered literal reservations,136cap, native EFS/gradient/stress/volume/
distance/force and source hashes. Complete mode requires normal Slurm32CPU/
one-node completion, both136recordchannels and native-rebuilt two-step
actions/blocks/spectra/reciprocity/cost summaries. No clipping/ASR/forced
positive/stability/TS/prediction certificate. Start inventory explicit,
pending/inflight JSON/failure/completed distinct; never retry/submit.

HFactualprefix08:14:33PENDING0records/reservations/blocks cb1aa8b8.
ActualoneoldnativeE160T smoke425539d9 re-read7cf9ab9b/log09cd83bd,
genuine32MPI and allnativechecks passed.0newDFT, not newstencil data.
Actualcomplete136/272new native branches still awaitdata; syntheticfull
inventory is not physical validation. AuditorSHA d477fc05 matchesHF/local.
73local/73independent related checks4modules19new,0finalfail/errors/skips,
notfullsuite. First68pass/1fail caught missing-directory complete mode
returningnotstarted. Audithelperfixed and outer136guardadded; no production/
physicschange. RawfailureXML retainedunchanged, later69/70successfulreports
remainC. First2tests inheritedEtemp; allsubsequentfixture/archiveworkC,
no unrelateddeletion. Testedtree8535c272/archive4c0a1ea7.

RO08:19:36both28867177/28868827PENDING/AssocGrpCpuLimit,0newSCFs/labels,
bundleunchanged. Fiveuseredits/main92f15239 preserved. T34/81,P35/81,
eightindependentfacepointseach/physical+.005sealed. No newDFT/reservation/
schedulerwrite/privateauthorarray/sourceoverwrite/parameter/radius/threshold/
material/wholepath/CI/install/release/PyPI/monitor/otherhost/submission.
ActualB3/releasedcurvature/stablebranch/resolution/strongcontrols/frozen
independentmaterialforecasts remain. Fullgoalactive/incomplete.
Evidence:G3_released_audit_E183_20261015.

### E184 2026-10-15 — PROGRESS: actual post-hoc nonlinear energy-control experiment

PreviousE183PROGRESS/28bc8b29. Used original25/17candidate labels and fixed
grid-parity9coarse/16fine split, no energy/errorselection. All labels already
seen: exploratory coarsening, not blindvalidation. Bilinear/tensorquadratic/
refittedtotalquadratic/twooriginalH models all kept. T common16: maxseen-
originshape errors6.900894/2.003468/7.667090/12.901228/12.975249meV12atomcell,
RMSE4.219064/1.043743/4.006016/5.177484/5.219402. Refittedquadratic originbias
+1.546418, rawmax6.120672 retained; no beststep/metricselection. Standard
tensorinterpolation includesmixed totaldegree4, not a new harmonic theory.

Preserving5coarse/12fineobserved,4coarse/4fineholes remain. Quadraticrank5
refusesglobalfit, allmodelcommonempty/nullnotzero. Otherarms support only
fouraxispoints1–4; tensor0.470404error is 1Daxisdiagnostic, not2Dsurface
validation. Negative/nonzero weights retained; missingrequiredlabels abstain.
No fictional hole values, endpointforce-as-envelope assumption or errorbound.

Crucial leakage/cost: all8nonorigincoarse starts in BOTHcenters inheritfine
parents. Fineenergywithholding is not wholeworkflowblindness. Historical
coarse/fineSCFs47/71(T),58/114(P),includingPcoarsehole46/finehole48,not9SCFs.
Original296/Gamma/canary/G2/seed preparation additional; differentinformation
for interpolants vs unfittedH, no equalcost/speedup/DFTsaving/B5noveltyclaim.
Decisionretaintheseelementary controls for future independent comparisons,
not replace B0/B1barrierstrain or strongT/Cmma. No forecast freeze fromthisfit.

37local/37independentrelatedchecks3modules12new0finalfail/errors/skips,
notfullsuite. Independentfullreportbytes377f2c3e; testedtree3bbd688b/
archiveca97e2dc. First36pass1fail testassumed48coarseholecalls; actualsource
arithmeticalready46becauseonecornerused10calls. Correctedexpectation from
nativeproof, notdata. Rawfailurepreserved, initialreportC, allscratchC.
RO08:36:18bothjobsPENDING/AssocGrpCpuLimit0newSCF/labels. Fiveuseredits/main
92f15239 preserved. No newDFT/reservation/schedulerwrite/rawSCFreparse/private
authorarray/+.005/independentpointaccess/parameter/radius/threshold/material/
wholepath/CI/install/release/PyPI/monitor/otherhost/submission. ActualB3/
releasedcurvature/stablebranch/resolution/independentforecasts remain,
fullgoalactive/incomplete. Evidence:G3_coarsening_E184_20261015.

### E185 2026-10-15 — PROGRESS: prospective elementary local-face predictions

PreviousE184PROGRESS/c250e79c. Fixed central9 training nodes at each original
center; eight novel half-grid targets each, no test labels. Five standard
arms/80 scalar forecasts, same origin calibration, both H steps retained.
Single training-only starts all selectQ0 with exact integer-distance ties;
no test-to-test continuation. .03target/.5radius/1.6distance and physics
unchanged; max192newattempts, no new allocation. Independent slots now
reserved8/face/computed0; modeling T34/81,P35/81 unchanged. Full B0--B5,
strong native-reference/domain and physical+.005 freeze NOT complete.

Engineering reportingtarget1meV/fu=4meV/cell chosen after seen training,
before labels, not DFTerrorbound. Tensor-H01 disagreements only T.502421,
P1.465334meV/cell. Could be unresolved, no winner guaranteed or target
movement to force gain.9trainingtargets historically28/36SCF; fullfaces
118/172, raw296/Gamma/canary/G2/branch preparation additional. No equal
budget/DFTsaving/acceleration/new harmonic-theory claim. Original B3 and
released-Q0 native audits and final candidate predictions precede test.

38local/38independentrelatedchecks3modules11new,0finalfail/errors/skips,
notfullsuite. CLIbytes aa77cdd7, testedtree04f584c7/archive9ab42037.
First independent35pass1metadatafailure: WindowsCRLF vs GitLF rawsource
hashes. Originalforecast283f2b9d and failureXML retained. V2normalizes only
text-code hashes, all numerical rows/starts/policies checked unchanged;
nativeDFT/JSON/input/binary pins still exact. Allscratch/snapshots C.

RO08:55:09 both28867177/28868827quotaPENDING,0newSCF/labels. No newDFT,
schedulerwrite, live-source overwrite, privateauthorarray/physical+.005
access, parameter/radius/threshold/material/wholepath/CI/install/release/
PyPI/monitor/otherhost/submission. Fiveuseredits/main92f15239 unchanged.
Fullgoalactive/incomplete. Evidence:G3_local_forecast_E185_20261015.

### E186 2026-10-15 — PROGRESS:16 actual inert native initial decks

PreviousE185PROGRESS/7efb2b56. Original16targets/80predictions unchanged.
Fresh HFnamespace, old380runtime/sixassets/binary checked, two old released
native origins re-read. All16literal40vector starts written as inert
INPUT/KPT/STRU and roundtripped; fixed substrate/control/order/Q/radius
guards pass. Minimum distance2.052545887A; volumes134.568448--134.666503A3.
No new SCF/allocation/schedulerwrite.192attempt cap/reserved8/computed0each
are inherited, not additional. TrainingT34/81,P35/81 unchanged.

Helper only prepare/preflight, not a completed production test executor;
actual B3/released-Q0 native evidence and final pre-label model forecasts
still precede test execution. Literal anchor policy retained for future
evaluation. No parameter/template substitution, CI or branch repair.
Native report5bc192b8 and64file ZIP6015de45 downloaded;48STRU/INPUT/KPT
hashes and16manifests verified. No pseudo/orbital/binary redistribution.

32local/32independentrelatedchecks3modules9new,zero fail/error/skip,
notfullsuite; independentCLIconfig27ff592f exact. Testedtreefa9af21d/
archive8b6930f4. Allscratch C. Read-only09:06:53both original handles
quotaPENDING, notfailure/restart. Fiveuseredits/main92f15239 andlive source
unchanged. No privateauthorarrays/+.005/independentlabels/newmaterial/
wholepath/CI/install/release/PyPI/monitor/otherhost/submission.
Fullgoalactive/incomplete. Evidence:G3_independent_inputs_E186_20261015.

### E188 2026-10-15 — PROGRESS: registered strong-reference forecasts on locked targets

Previousgoalturn E187verifiedwait, no artificial status commit. All seven
original spans/allfourCmmaframes/twoHsteps at E185eighttargets eachcenter:
224formalrows,140acceptedpredictions/84geometryabstentions,0newDFT.
Ttriplet/Tdoublets/fullT/eachCmma coverage T6/8/8/6, preserving0/8/8/2,
bothsteps. FullT andtwo-doublets cover all8includingoff-axis; cannot omit
them to claim weak-reference failure. Allreference commonT1--4axesonly,
Pempty. P Cmma1/4only; shared-Q restricted models, not universalnative
Cmma failure. Omitted modesclamped, not full outside-reference release.

FullTvsE185H max1.6653e-16eV/cell standardidentity, notgain/independent
algorithm. Two-doublet/fullT disagreementsT474.975/474.958,P542.960/
542.668meV/cell mathematical only, no new DFTbarrier/accuracy conclusion.
Private source arrays read-only/in-memory, no raw/mapped/response export;
LDAfrequencies not PBEcurvature. No framesdeleted/worstorientation picked.

49local/49independentrelatedchecks4modules9new,zero fail/error/skip,
notfullsuite. All224 formalenergies/norm2 independently checked from old
E180scalarpolynomials; actualprivate-bank CLIbytes97692cb0. Testedtree
fd0871fa/archive84b6704b. Allscratch C. Costsoriginal296/118/172/Gamma/
G2/canary/mapping/branchadditional; no equalprepcost/DFTsaving/newtheory.
E185targets/starts/scoring/192attemptcap unchanged; trainingT34/P35,
independentreserved8/computed0each, physical+.005sealed. FullB0--B5/native
chart/finaldomain/admission/stability/TS/independentgain NOTcomplete.
RO09:20:20both oldhandlesPENDINGquota0SCF, no newschedulerwrite/runtime/
physics/material/wholepath/CI/install/release/PyPI/monitor/otherhost/submission.
Fiveuseredits/main92f15239 unchanged. Goalactive/incomplete.
Evidence:G3_reference_forecast_E188_20261015.

### E189 2026-10-15 — PROGRESS: nonlinear control integrated into the real manuscript

PreviousE188PROGRESS/d3406d26. E184evidence now inMethods2.7/Results3.9/
Table9/Discussion. All5same16Tfine targets retained: tensor2.003468 vs
H12.901228/12.975249meVcell, notnewlabels/independentgain. Cohortdiffers
from Table8, so smallerHmax notmodelimprovement. Rawfitmax6.120672/
originbias1.546418/P rank5holes/axis-only diagnostic retained. Coarse
fineparent leakage explicit, historical47/71,58/114SCFs+prep, no saving.

Fivequantrows sourcechecked377f2c3e; all135oldrows/9figbytes/13externallinks/
Texwrapper retained. Twoledgerterms.126local/126independentrelatedchecks,
10modules13new,zero finalfail/error/skip, notfullsuite. First124pass2fail
stale8sections/135rows updatedto9/142, all oldsourcechecks retained;
rawXML06ea5db4. CLI4500752b byteexact. Testedtree144c799d/archivef70f911c.

ExistingTeXLive2023/LuaLaTeX actualMarkdown9figures compiled35pageworkPDF,
allworkC. Firstphrase-textlookup linewrap fixed in verifier only. Visual
reviewall35contacts/detail8/21/22 caught splitA9caption34/35; notesmoved
toprose, shortcaption+table together, noquantchange. OriginalPDF/reports
retained. Final1--33pixelidentical, changed34/35detailsrechecked;0overfull/
underfull/missingchars. MD25738d69/PDF72a105d6, notsubmissionready.

RO09:37:19botholdhandlesPENDINGquota0newSCF. No newDFT/allocation/scheduler
write/privatebank/+.005/physics/radius/threshold/material/wholepath/CI/
install/release/PyPI/monitor/otherhost/submission. Fiveuseredits unchanged.
TrueB3/releasedcurvature/stablebranch/independentforecasts/tests stillmissing;
fullgoalactive/incomplete. Evidence:manuscript_nonlinear_controls_E189_20261015.

### E190 2026-10-15 — PROGRESS: real startup and proven provenance replay correction

PreviousE189PROGRESS/04a1ab77. Both committedjobs started09:42node6 and
09:44node9, genuine32MPI/threads1 confirmed. Original E183 reader falsely
rejected first real derivative: every physical scalar/gradient/stress
matches exactly, validation pathspellings differ bundle/../ vs canonical.
Strict actual resolved directory identity same. Initial direct-replay
diagnostics used storedparent and thus passed; tracing originalfunction
localized differingvalidation. No numerical tolerance or SCF problem.

E183/production/nativefiles preserved. NewE190reader only executable-line
change validator(raw.parent), AFTER unchanged namespace/resolvedidentity
checks; all physical/input/byte/reservation tolerances stay. Oldfailure
reproduced; wrongrealpath/physicalfield/MPI mutations stillrefused.
ActualcorrectedHFprefix13probe records passed b80d7571, nextpending and
Pnotstarted explicit. RealfirstSCF sevenfiles+manifest zipd7bd98e9,
resultc35035e2, no pseudo/orbital/binary. NativeEFS/gradient recheckedlocally.
71local/71independentrelatedchecks4modules13new,zero fail/error/skip,
notfullsuite. Testedtree74357434/archive10808dee, allscratch C.

Original B3native reader firstTQ0 passed1newSCF/cache-onlyreplay: atomic
norm.02147688 vs joint.11042982, not jointstationarity/stability. Complete
batch/pairedaxis errors stillawait. RO10:16:44bothRUNNING: atomic6/10
candidates15SCF, curvature17/272,0SCFfailure. Remainingatomic47attempt
cap excludes finishedpointunused budgets, no rollover; mean131s implies
fullceilingcost~12:00. Hmean108s gives~17:54. Costscenarios, notguarantees;
walllimits14:42/21:44. No additionalDFT/scheduleroperation/sourceoverwrite/
parameter/threshold/budget/privatebank/+.005/independentpoint/newmaterial/
CI/install/release/PyPI/monitor/otherhost/submission. Main25738d69/fiveuser
edits preserved. Goalactive/incomplete. Evidence:G3_released_audit_E190_20261015.

### E192 2026-10-15 — PROGRESS: measured atomic versus joint treatment comparison

PreviousE191VERIFIED_WAIT on the same two live handles. E192native prefix
audited7finished nodes/20newSCFs/exact original optimizer replay: T0--4,
preserving0--1, preserving2pending. Immutable E176audit673714cb consumed
by unchangedE178, all20pre-label model-point rows/two steps retained.
Reproduction16b9dde6 byte-exact; no coefficients fitted or newDFTinitiated.

PreservingQ0 A-F=-48.348691, J-A=-1.136239meVcell: large release is already
atomic there, not primarily extra cell release. Candidate-treatment
difference, not unique mode-energy partition/global nesting/fullpath cause.
T same4axial shape discrepancies F.416264/.403160, A.501879/.515704,
J1.327177/1.346071meVcell. Pointwise error identities, not subtraction of
unmatched maxima. Preservingonly1/4axial A labels; no partial ranking.
Atomicorigins jointnorm.110430/.194104 remain above.03, notjointstationary.

75local/75independent relatedchecks4modules10new,0fail/error/skip,
notfullsuite. Native20callproxy22.331631coreh excludes in-flightwork;
8frozenstatics/296originalcurvature/joint/preparationadditional. Actual
releasedstability/fullB0--B5/independentpredictivegain stillmissing.
Main25738d69/fiveuseredits/production/forecast unchanged. No newDFT,
schedulerwrite/parameter/budget/domain/+005/independentlabel/material/
CI/install/release/PyPI/monitor/otherhost/submission. Fullgoalactive.
Evidence:G3_atomic_treatment_E192_20261015.

### E194 2026-10-15 — PROGRESS: completed native B3 and full matched treatment cohort

PreviousE193VERIFIED_WAIT, no pending-job restart.28867177normal COMPLETED
0:0/4577s/node6/32CPU.10/10atomic candidates/36newSCFs. Original E176
complete reader checks allnativeEFS/input/MPI/order/cachedoptimizer replay
andterminal summary, proof8a2fdf2a. Unchanged E178complete analysis3cc150ec,
all20prefrozen predictions/two steps/5Qs each retained, no coefficient fit.

Preservingall5 A-F=-44.06to-48.35,J-A=-.80to-1.14meVcell: large release
alreadyinatomicstep for this local cohort, not entirepath/uniquemodepartition.
Its4axialmodelshape errors A4.349824/4.399770,J4.671404/4.720758, so
omittedcellalonecannotexplain mismatch. Missedatomicrelease shape4.284867/
4.284892 versusextraJ.332117/.331406, signedpointwiseidentity notsubtract
differentmaxima. Atomicall10jointnorms>.03, notjointstationary/stable/TS.
Tsame4axial values retained; old7nodes/20native source identities unchanged.

83local/83independent relatedchecks5modules8new,0fail/error/skip, notfullsuite.
IndependentCLIbyteexact. Atomic9/27SCFs T/P, callproxy39.749558coreh,
allocation40.684444, oldfrozen8/curvature296/G2/joint/prepadditional,
notbudget120versus36speedup. FullB3 prerequisite nowdone, releasedstability/
finalcandidateforecast/B0--B5/independentprediction remains. RO11:08:05
28868827RUNNINGnode9/32CPU/41of272/noSCFfailure, mean119s scenario18:45,
notguarantee. No newscheduler/DFTinitiation/runtime/physics/domain/budget/
independentpoint/+.005/privatebank/newmaterial/CI/install/release/PyPI/
monitor/otherhost/submission. Main25738d69/fiveuseredits preserved.
Fullgoalactive/incomplete. Evidence:G3_atomic_complete_E194_20261015.

### E195 2026-10-15 — PROGRESS: author-approved four-part manuscript and real treatment integration

PreviousE194PROGRESS/e3c599d0. Author explicitly chose four-part Results.
Existing9sections→4: channelcompetition/mechanicalboundary/modedescriptors
andatomicrelease/predictioncontrols. Sixdetailed diagnostics movedAppendixE,
no data/qualification deleted. Tables1--3/A11--A16referencesupdated.
E194actual10point/36SCFtreatment nowMethods2.7/Results3.3/Discussion/D.8A10;
4newquantrows sourcepinned3cc150ec/8a2fdf2a, samecenters/bothsteps/allarms.
All142oldtablelines/countsmultisetretained+6=148, ninefigureorder/bytes,
13externalcitations/Texwrapperunchanged. Twoledgerterms, publicreadingguide.

138local/138independentrelatedchecks11modules12new,0finalfail/error/skip,
notfullsuite. First127pass11fail: rawnewlinephrase, oldorder/count and8
negativefixtures usingfirstscalar-matchingrow now selectedA10notoldF/J.
Fixfullsourcerowidentity/approvedorder, do notrelaxphysicalvalues. Rawfailure
XML2025eef5retained. FinalCLIa6230114byteexact, testedtreee48fddab/archive
dd93dbc9. No productive physics/runtime/model coefficient changes.

RealLuaLaTeX36pagePDFallrendered/contactreviewed, detailed11/16/30,
captionQA foundA13/A15crosspage32/33and34/35. Notesmovedtoprose, short
captionsnow33/35sameastables. Final1--31pixelidentical, changed32--36
allreinspected, zeroover/underfull/missingchars. MDb9e29fca/PDF6f48c872.
InitialPDF/metadataretained, notsubmission-ready. Allbuild/snapshotsC.

RO11:36:01curvature28868827RUNNINGnode9/32CPU/54of272/noSCFfailure,
mean119s scenario18:49notguarantee. No newDFTinitiation/scheduler/physics/
domain/budget/privatearrays/+.005/independentlabels/material/CI/install/
release/PyPI/monitor/otherhost/submission. Fiveusereditsprotected.
Independentstability/forecasts/fullB0--B5/physicalholdout/scientificgain
stillmissing, fullgoalactive. Evidence:manuscript_atomic_treatment_E195_20261015.

### E205 2026-10-15 — PROGRESS: full offline regression and bounded manuscript-check repair

PreviousE204VERIFIED_WAIT; E196--E204 only rechecked the original live job.
Clean4d03befe export:3502passed/2failed/2skipped. OldE094section helper
still split3.4/3.5; absentending caused read-through to laterappendices,
potentialfalsepasses and139numbers instead16Cmmabounds. RestrictE.1--E.2,
uniqueanchors/nointerveningheadings; retainall8signedrows/4prosepartitions/
4Cmmaframes/16bounds/original5.0001e-7atol/zero-relative tolerance.
Sixnewchecks refusemissing/duplicateboundaries/rowborrow/all16boundmutations
despitecorrectvalueselsewhere. Target21passed; no scientificdatachanged.

Independenttestedtree1777a194/archivebb9dca41 fullsuite3510passed/0failed/
0errors/2skipped; all3506oldtestIDsretained,6new. Same skips:optional
ase.io.abacus unavailable/licensedVASP POTCARnotredistributed.314warnings
remain, notwarning-free. First/finalfullXML32a590b0/97083f53andtarget
e40f202fretained. Onlytest_hfo2_E094_manuscript_work.py executablechanges;
manuscriptb9e29fca/fiveuserhashes/runtime/physics/domain/budgetunchanged.
No DFTinitiation/schedulerwrites/independentlabels/install/CI/release/
PyPI/submission/monitor. Cscratchonly. RO12:28:26curvature28868827RUNNING
node9/32CPU/82of272/noSCFfailure; mean117.746s scenario18:41notguarantee.
Releasedstability/finalforecasts/fullB0--B5/physicalholdout/scientificgain
stillmissing. Fullgoalactive. Evidence:full_regression_E205_20261015.
