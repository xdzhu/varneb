# E054: first matched clamped G2 pair, reviewed geometry continuation

P=0, common T-plane substrate at training strain ε=0; fixed cell rows0/1,
third vector's three components and all12atoms released. Every SCF retains
the original ABACUS100Ry/full10-auDZP/six physical files/32trueMPI/thread1.
No CI, physical adjustment, +0.5%holdout or extra independent channel.

## Actual completed first segments

| Channel | Prior job | Completion CST | Initial→step10 ordinary fmax (eV/Å) | Step10 sampled forward/reverse (meV/f.u.) |
|---|---|---|---|---|
| T-pattern-preserving PO flip | 28574708 | 2026-10-10 14:12:02 | 3.541469→1.466699 | 287.182317 /287.217576 |
| PO→M | 28574709 | 2026-10-10 13:45:24 | 0.994612→0.615276 | 114.427812 /206.968306 |

Both SlurmCOMPLETED0:0, material`max_steps_reached`, not ordinary-converged
MEPs. Do not rank mechanisms using these unfinished discrete peaks or use
them as H1 labels/TS energies. Last four original log rows decline for both;
neither was stopped for a rebound. Each segment has77fresh interior SCFs
(11complete frames×7interiors), no new endpoint SCF. Native log and all
six physical bytes/STRU/pinned EFS checked for all154calls. SCF transport
wall×32 sums84.231436+70.219580=154.451016coreh; this is not CPU profiling.
Allocation wall×32 totals155.822222coreh from actual scheduler elapsed.

`observations/` freezes complete step0/10 for both chains: exact POSCAR
lift, numeric E/F/stress trajectory, raw provenance and **same-boundary**
force replay. Export/analysis adds zero DFT calls, never modifies old runs.
Raw licensed PP/orbital files are not redistributed. Full physical-byte
checks ran on HF; portable data do not pretend those files exist locally.

`seeds/` prepares the identical audited step10 as a fresh geometry
continuation. Nine exact hash-pinned current-frame caches verified on HF;
after interior movement its cache invalidates normally, endpoints stay
fixed/cached. New FIRE state, not restored velocities/time step. The new
segment cap is20steps/12h, oneHFnode/32MPI per channel, at most2study
allocations; it remains the same two of the ≤8registered G2 chains.
Actual submission and live state are recorded separately after all gates pass.

At20:09:08CST, continuation28661019(node6)/28661020(node25)areRUNNING.
Both step0log rows reproduce old step10 exactly (1.466699/.615276),
first new interior SCF logs confirmDSIZE32. Fixed archive full regression
1401passed/2skipped/0failed/errors,666.41s; see`validation.json`and
`clean_regression_fixed.junit.xml`. Mean prior SCF wall123.0654/102.5935s
gives ~4.79/3.99h for20updates, roughlyOct11 00:00–02:00 at similar speed;
12h caps at08:07–08:08. This forecasts segment completion, **not**convergence.

## Validation and two resolved delivery-harness issues

1. Analysis here-document: PowerShell added a final CRLF to the `PY`
   delimiter. The Python body completed all audits, wrote complete receipts,
   then raised`NameError: PY`. No DFT or submission occurred. Preserve that
   nonzero exit; `verify_existing_HF.sh` independently rechecks18exact
   caches and all artifacts, writes a new passing verification receipt,
   without re-exporting or overwriting caches. A final shell comment fixes
   future delimiter transport.
2. First clean archive regression:1392passed/7failed/2skipped. All7failures
   correctly reject substrate SHA drift, not force/DFT changes. Native
   substrate is984bytes/LF/SHA`1b718f99…`; Windows Git archive of previously
   unprotected E053data changed it to1004bytes/CRLF/SHA`46dacfdc…`.
   Add `-text -eol` protection to both G2artifact namespaces; new tests check
   every E053manifest byte and substrate SHA in the clean archive. Do not
   normalize source files, change the physical contract or loosen the gate.

Initial analysis source tree`2c3aca9`, archiveSHA
`48c7c0e8550e3812d53d810c5e9d720f03693b0a7bc085856a663afc2e7712c2`.
Fixed archive tree`2b03239`, SHA
`170170ed6c3f8749c4930ca8af6a67fc17689dcdf57481e916577e40aea70d9c`.
Executable/analysis/Slurm files are identical between them; only packaging
attributes and archive tests differ. Fresh`source-fixed`never overwrites
the earlier analysis or E053production source. Full regression receipt and
submission gates must pass before any continuation is launched.

Observable archive SHA
`a8f9d9bd7909210084c4aa4a22c366ca7914b23fe53edcdc7d06e9f665e3c188`;
`audit_prepare.json` SHA
`d0415b847d038c68cf20d9ce30191e67db32988dccc807b7eaa5680919472abc`.
Existing +1% training wells are not rerun; the +0.5%holdout remains unread.
No G3sampling, independent prediction, manuscript-ready claim or release.
