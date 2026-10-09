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
|0|PO- T-pattern-reversing|E049/28464144 RUNNING at23:51:45CST; max20steps|pending|pending|
|+1%|T/PO+/M/two PO- candidates|five registered seeds, not submitted|pending|pending|

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
