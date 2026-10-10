# Registered common-substrate endpoint matrix: actual progress

This index covers only the already registered five phases at two substrate
strains0/+1%, P0/E0. The +0.5% independent holdout remains unseen/uncreated.
These are endpoint screens, **not ten completed stable wells or G2 barriers**.
Original PBE/100Ry/full10auDZP/Gamma2x2x2 physical bytes are unchanged.
All atoms and third-row length/two tilts relax; substrate rows0/1 stay fixed.
No symmetry/mode is forced to retain a nominal phase.

| Condition | Endpoint | Evidence/status | Force eV/A | Open traction kbar |
|---|---|---|---:|---:|
|0|T|E047/28454923 completed; three labels P4_2/nmc|0.000499915|0.009594677|
|0|PO+|[E046/28453655](../clamped_PO_continuation_E046_20261009/README.md) completed; three labels Pca2_1|0.022284634|1.928434990|
|0|M|E048/28456312 completed; three labels P2_1/c|0.018206569|0.843275452|
|0|PO- T-pattern-preserving|E048/28456313 completed; three labels Pca2_1, variant gate separate|0.014888417|1.692852091|
|0|PO- T-pattern-reversing|E049/28464144 completed; three labels Pca2_1, ordered-gauge review below|0.014888417|1.692852095|
|+1%|T-derived seed|E050/28469832 completed; Ccce already present in strained seed, see geometry review|0.003364262|1.803737454|
|+1%|PO+|E050/28469833 completed; Pca2_1 at three tolerances|0.013274962|0.995806640|
|+1%|M|E051/28472493 completed; P2_1/c at three tolerances|0.007755026|1.452557036|
|+1%|two PO- candidates|registered seeds, not submitted|pending|pending|

T required one fresh SCF and no ionic step. PO+ reused the exactly identical
last clamped canary SCF, then required two fresh SCFs/steps with a fresh BFGS
Hessian. This is not an acceleration benchmark or free-cell cache reuse.
Endpoint criterion is atomic max<.03eV/A and open-traction norm<2kbar;
clamped reaction stress is informational. Phase labels at.001/.01/.05A and
1degree are not Hessian stability or electronic-polarization/variant proofs.

At this screened condition E(T)-E(PO+)=12.359975meV/HfO2. Free-cell separation
is81.321193meV/HfO2. The boundary switch is not a same-ensemble strain derivative;
neither gap proves a switching barrier or selective retention improvement.
Paired endpoints alone do not close the G2 candidate-network gate.

The preserving PO- seed subsequently completed0:0 at23:33:56CST after
7BFGS steps/8fresh SCFs. All8original HF raw E/F/stress and six physical
bytes are reviewed, phase labels Pca2_1at all3tolerances. The force/open
traction screens pass; this does not on its own establish electronic polarity,
the structural pattern registration, a Hessian minimum or G2 barriers.
Its actual energy is-9782.973971976650eV/Hf4O8; no enforced energy equality
with PO+ or symmetry repair is applied. Actual cost17.427534SCF core-hours
vs18.008889allocation core-hours. The72observable export and recorded
zero-DFT audit invocation are preserved in the phase directory. The earlier
23:29running-status receipt remains historical, not the current job state.

M completed0:0 at23:42:52 after17steps/18fresh SCFs. All18original raw
SCFs/physical bytes pass; P2_1/c at three tolerances, with no symmetry
restoration. Its endpoint energy is-9783.343992919006eV/Hf4O8, or
-92.540494meV/HfO2 relative to the same clamped PO+. This is a well-energy
gap, not its escape barrier. The152observable export and reproducible
zero-DFT audit invocation are preserved under M. Phase labels and physical
screens still do not certify all-variable stability or electronic polarization.

The remaining zero-strain reversing PO- seed was manually reviewed and
submitted as E04928464144 at23:50:34CST, same immutable tested source, all
six physical bytes,32MPI,20steps/21SCF maximum and2hour cap. M and preserving
jobs had actually ended before submission; only one study calculation was
then running. This is the existing fifth endpoint, not an additional strain,
orientation, independent prediction or automatically submitted G2 matrix.

E049 completed0:0 at2026-10-10T00:11:23CST after7BFGS steps/8fresh SCFs.
All8actual raw SCFs, six physical bytes,32MPI and three Pca2_1 labels pass.
Energy is-9782.973971976642eV/Hf4O8; actual SCF cost11.061641core-hours,
allocation11.102222core-hours. Its72observable export and zero-DFT replay
are preserved; no endpoint is resubmitted and no physical input is changed.

## Ordered endpoint correspondence is not a new phase or a new path

The [read-only variant audit](zero_strain_variant_audit.json) compares all five
screened endpoint **representations**, not five independent physical phases.
Composing the already registered parent inversions gives translation(.5,.5,0)
in the research axes, with source-to-target atom indices
`[2,3,0,1,9,8,11,10,5,4,7,6]`. This operation is declared from the parent,
not selected by a best fit to the relaxed structure or barrier labels.

The two minus seeds match under that operation to1.7e-15A; their independently
relaxed endpoints match to1.85e-11A. Yet the same ordered indices differ by
1.36150A under PBC, so the strict ordered endpoint cache correctly rejects
them as identical. Raw terminal energies differ by7.28e-12eV/cell, forces
covary after the declared permutation to3.0e-10eV/A, and stresses to
2.38e-12eV/A3. These are observed pair residuals, not a general numerical
error bound or permission to reuse a +1% endpoint without a separate gate.

Both minus endpoints retain opposite signs of the dominant **geometric**
T-pattern coordinate(-.924790/+.924790A), and their oxygen-minus-Hf mean
displacement along research x is negative while PO+ is positive. These are
ordered-gauge pattern observations, not electronic polarization, irrep labels
or mode-energy decomposition. Non-Gamma pattern signs depend on the parent
origin; a single endpoint sign is not an intrinsic phase identifier.
Translation equivalence of endpoints does not
prove whole paths equivalent, nor distinct topology: retain both original
path correspondences/lifts and do not permute a production chain to merge them.
No source geometry, atom order, inputs or runtime is rewritten by the audit.

## Registered +1% training condition: E050

After every earlier study job was confirmed terminal in sacct, the already
prepared T and PO+ seeds were reviewed against their original hashes,
substrate/P0/E0 and all six physical bytes, then submitted once as28469832/
28469833. Both started00:25:50CST on node6/node28;00:27:06 snapshotRUNNING.
Each has1node32MPI/threads1,20BFGS steps/≤21SCF and2h maximum, same.03/2/.02
and100Ry/full10auDZP. The original1330-pass immutable production source is
reused read-only, not replaced by the new offline analysis code. Exact handles,
preflight and invocation are in [E050 receipt](submission_E050.json).
Expired completed IDs may be absent from squeue; their sacctCOMPLETED0:0
records were verified before submission. Unrelated jobs remain untouched.
The +0.5%holdout is ungenerated/unread; no new phase, G2 band or G3 is submitted.

Both E050 handles subsequently COMPLETED0:0. PO+ ended00:38:02 after6steps/
7SCFs, same Pca2_1. T-derived ended00:40:45 after6steps/7SCFs. All14SCFs
have raw full-EFS/32MPI/original-six-byte reviews and128observable files.
SCF costs are6.358583/7.935382coreh for PO+/T-derived respectively; both are
training endpoints, not independent forecasts. Their well gap is14.145883meV/
HfO2, not a measured switching/escape barrier or TS response.

**Do not require P4_2/nmc for the +1% nominal T seed.** The registered substrate
contains the original long axis and one of the two formerly equivalent short
axes; stretching its two rows while opening the third makes those short axes
unequal. The seed already has Ccce at all three tolerances *before DFT*, and
all seven evaluated geometries retain Ccce. The[geometry review](strain_p0100/T/strained_T_geometry_review.json)
shows only the original geometric T pattern changes(.831575→.885753A), with
off-pattern residual≤1.33e-14A in the fixed original metric. This supports a
metric-lowered strained-T descendant, not a numerically restored tetragonal
phase or a newly certified bulk polymorph/Hessian minimum. The first new
replay test failed because it incorrectly expected the free-T group; corrected
expectation plus a before-DFT metric/pattern regression preserve the observed
structures without changing any input. Single orientation, no new phase task.

E051 submits only the existing +1% M seed after PO+'s raw/phase review. Job
28472493 started00:42:08CST with the same32MPI/20steps/21SCF/2h recipe and
physical bytes. T had actually ended before it was submitted, so this was one
study allocation, not three. Receipt[in E051](submission_E051_M.json);
neither running source nor any unrelated job is touched. No auto restart.

E051 subsequently COMPLETED0:0 at01:18:38CST,16steps/17fresh SCFs.
All17raw full EFS, original six physical bytes, native32MPI and unchanged
substrate boundary pass the zero-DFT replay. P2_1/c is retained at all three
tolerances; force.007755026eV/A and open traction1.452557036kbar pass.
Energy is-9783.405592825406eV/Hf4O8, or-114.746973meV/HfO2 relative to
the same +1% PO+; this is a well gap, not an escape barrier. Actual SCF
cost19.369124core-hours and allocated cost19.466667core-hours are separate.
The144observable bundle, invocation, raw audit and completion receipt are
preserved under strain_p0100/M. At10:16:42CST squeue was empty: there is
no surviving study job, no failure or unknown-live resubmission. Eight of
the ten registered endpoint representations pass the physical screens;
both +1% minus seeds remain unsubmitted. G2 barriers, G3 and the uncreated
+0.5% holdout remain missing. No endpoint count is a JCTC acceptance gate.

E048 is exactly two manually reviewed **existing** endpoints, each1HFnode/
32MPI,20BFGS steps/at most21fresh SCFs,2hour cap, original .03/2/.02 gates.
Maximum two study allocations are retained; unrelated user jobs are not stopped
or included as this study's evidence. The legacy generic CLI returns1 for a
completed nonconverged step-limit segment; a Slurm failure alone must therefore
be reconciled with its actual summary/raw logs, not treated as a reason to
repeat an unknown/live task or relabel an SCF failure. No automatic restart,
matrix, CI, new material/strain/orientation or holdout submission occurs.

Per-endpoint receipts and exact observable bundles remain under condition/
phase directories. PO+ continuation retains its existing stable case link;
do not duplicate its raw data here. Licensed/proprietary basis/pseudo/charge
files stay on HF and are checked there, not claimed from these local exports.

The [E049-E051 delivery receipt](validation_E049_E051_20261010.json) identifies
the exact clean staged source tree:1354passed/2skipped/0failed, plus660actual
observable files byte-matched. The current13page working manuscript was
compiled and changed pages11-13 reviewed; it remains scientifically incomplete.
The earlier validation_delivery.json is retained as a historical milestone.
