# E140: audited +1% preserving child40 cap, not G2 convergence

28737465 ends natively COMPLETED/0:0 at08:50:08CST on2026-10-12, after40
updates. The same-boundary final residual is0.11414663193850914eV/Angstrom.
This is a successful finite step cap, not an ordinary0.10 pass, a failed SCF,
a stationary saddle, a continuous MEP/barrier precision certificate, or a
finished training label. G2 remains6/8; reversing28764704 remains running.

The unchanged E103 driver executes its existing `audit-child --parent
28709790` logic, only after actual native completion. It verifies280fresh
SCFs,338original E055 runtime files, six physical input/asset hashes, native
DSIZE32, electronic convergence and complete energy/forces/stress. Nine
latest-frame cache evaluations require zero new SCFs. Endpoints are not
recomputed. PBE/100Ry/full10auDZP/Gamma2^3/P0/E0/ordinary.10/noCI, clamped
rows0/1/open third-row tilt and the original optimizer parameters are unchanged.

Terminal audit SHA256:
df34573384f76650e123691dbb508a52d7d8ee529e6e50dca0fce5e4e6835c02.
Latest-seed manifest SHA256:
1a222cc63b89344816f0829bbf079793fa06306d3d4d13f5cc3c2b401f5fbf3f.
The driver SHA remains6c64c0155800082e28cbbb280511ac761929bdb6ccc5757fa8d68d156b32eac7;
the wrapper SHA remainsab914119aaa7918bd0cb466ed8dd86daafecb20609a077a5d249d89a3228eed0.
No existing receipt/journal, running source, physical asset or old figure is
overwritten. Audit outputs are in the new `child_terminal` namespace, not
the original starter's `terminal` namespace. The audit initiates zero DFT
and zero scheduler writes; the independent reversing job continues DFT.

## Cost and reviewed next action

All280fresh calls total29,003.678637seconds and257.81047732-core-hour
transport proxies. These are call-time measurements, not the allocation
accounting itself. The allocation runs00:45:12→08:50:08,29,096seconds,
258.631111allocated CPU-hours. These two different cost measures must not
be silently interchanged. Preparation/previous allocations are additional.

The residual falls1.325704→0.114147 in this segment, with the last ten
completed updates still decreasing0.137561→0.114147. The E139 diagnostic
rules out spring dominance at step37 but does not prove a different optimizer
would help. No single rebound or step cap justifies terminating the physical
candidate. A single further40-update/12h geometry/cache continuation is
justified by this healthy decreasing trend, with the same physical/optimizer
settings and no new independent channel. Nominal upper budget:280new SCFs,
about8.1h/258core-hours at the measured per-call rate, wall ceiling12h/384
allocated CPU-hours; actual convergence time is not predictable from this.

This audit **does not submit that continuation**. Existing E103 dispatch
accepts only the original ten-step starters; its accepted journal must not
be reset/replayed or its parent identity silently changed. A tested explicit
child40 ancestry registration and native spare-slot check are required first.
The study ceiling remains two committed one-node/32MPI allocations on
hfacnormal01. No automatic rolling dispatch, G3/reference/curvature selection,
physical+0.005holdout, CI or new monitoring automation.

## Portable replay

`completed_HF` is a byte-preserved copy of the completed native audit bundle,
including nine retained raw terminal SCFs and the exact latest seed. Local
tests reparse those nine E/F/stress records and replay the original residual;
they do not re-audit all280SCFs, the338HF runtime files or the undistributed
pseudopotentials/orbitals. Native all-call evidence is the pinned receipt.

```sh
python -m pytest -q tests/test_hfo2_E140_preserving_cap.py
```

Actual local and clean-archive regressions each pass127tests, zero failures,
errors or skips. Tested tree44fc496b9430b3b24a555aefecf7842e225639d0 has all
6097tracked files present before independent tests;78raw-evidence/JUnit files
are byte-equal and4source/metadata texts newline-normalized equal. No native
evidence bytes are normalized. Archive SHA256:
6ec76e75abe4f5b2a47a555a0b01eda6f39b0d4c8aea1385422c334514315e45.
Local JUnit619c1a939c32789a4800c198263b741ad62d80b8ce869e847215d5b2335b1dc3,
10.038s; independent JUnit70b97e66a3930c528bdc2cb55271070bc624319ba4e76b8074e5f090f4bff5b6,
7.942s. These are related tests, not a new full-suite or HFpytest run.
Only delivery metadata, independent JUnit and the ledger are added after the
tested snapshot; scientific data and test source remain unchanged. Native
portable bundle:77files/1,250,797bytes. Original five dirty tracked files and
all older untracked user work remain untouched and unstaged.

At09:06:33CST only reversing28764704 remains RUNNING32/hfacnormal01, latest
complete step17/.144601, stderr0. The preserving cap is terminal and is not
already resubmitted. Full JCTC objective remains active and incomplete;
manuscript labels are not promoted from this unconverged cap.
