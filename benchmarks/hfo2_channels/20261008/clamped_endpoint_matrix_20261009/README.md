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
|0|M|E048/28456312 submitted23:00:09CST; max20steps|pending|pending|
|0|PO- T-pattern-preserving|E048/28456313 submitted23:00:09CST; max20steps|pending|pending|
|0|PO- T-pattern-reversing|registered seed, not submitted|pending|pending|
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
