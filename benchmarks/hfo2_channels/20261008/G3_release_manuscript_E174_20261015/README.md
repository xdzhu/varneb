# E174: completed ablation in the manuscript and a negative bounded optimizer test

PROGRESS, not scientific completion. The JCTC goal remains active. No new
DFT submission, restart, reservation, input change or held-out label access.

## Connected manuscript

The completed E173 eight-static frozen/released ablation now enters Methods
2.7, Results 3.7/Table7 and the terminology ledger. Both original centers,
both probe steps and all16model-point records remain. Maximum frozen/released
shape discrepancies, meV/12atomcell, are:

| Center | Step/A | Frozen | Released | DFT release variation |
|---|---:|---:|---:|---:|
| Preserving | 0.01 | 0.075880 | 4.671404 | 4.606959 |
| Preserving | 0.02 | 0.116535 | 4.720758 | 4.606959 |
| T | 0.01 | 0.416264 | 1.327177 | 1.168761 |
| T | 0.02 | 0.403160 | 1.346071 | 1.168761 |

The signed pointwise identity, not a subtraction of unmatched maxima,
separates frozen-model discrepancy from missed release-shape change. One T
point has cancellation, not evidence of global accuracy. These observations
do not uniquely identify anharmonicity or branch change, certify stable
conditional minima, or replace the strong-control/independent prediction
matrix. Eight static call-proxy/allocation costs8.136814/9.075556core-hours;
296curvature SCFs and all prior path/conditional preparation remain extra.

All113old table lines,9figure files,13external citations and the TeX wrapper
are preserved. Six new table lines are added. The old preparation-only
sentence and its prose test are superseded by actual completed-native
evidence; the sheet-preservation test retains the original113lines as a
multiset subset plus an exact119line total. No scientific threshold or
uncertainty test is weakened. The first PDF build exposes a43.94397pt long
inline-coordinate overflow. Separate math spans fix the line break without
changing any coordinate. Final32page LuaLaTeX build has no overfull/underfull
or missing-character warnings; all pages are contact-reviewed, Methods and
Table7 inspected at higher resolution. Build products stay outside Git;
`render_metadata.json` binds the source and actual PDF.

## Actual terminal material test

The already reserved arm28866132 finished COMPLETED0:0,1417s, one node6,
32CPUs,hfacnormal01. Exactly12new SCFs/attempts completed, with no backend
failure. It produces **no gradient candidate** before its original cap.
Seven exact old cache hits precede the12new calls. Best new trial released
norm0.0458552401eV/A remains above0.03; the old12-trial minimum is0.0337907153.
Trial extrema do not identify accepted iterates; this is not an equal-cold-
cost comparison, and two capped nonconvergences provide no speedup ratio.
New completed-call proxy12.0498025core-hours; allocation12.5955556. Old12
SCFs/296curvature SCFs and other preparation are additional, not free.
No default change, repeated arm, larger radius or extended budget is made.

`terminal_audit.json` reparses all12new plus12old native E/F/stress,
six-asset/DSIZE32/SCF/reservation/contract records and strictly reproduces both
optimizer request sequences using cached E/g only. It rechecks original380
runtime files and all9immutable execution-bundle files. `new_cache.zip`
retains new result JSONs/contract/reservations and the terminal arm summary;
native logs remain on HF with SHA-bound validation. This is a failed bounded
acceleration hypothesis, not an ABACUS calculation failure or convergence.

## Audit failures and exact-coordinate resolution

The first audit setup expects the old helper at a wrong external path; it
never reaches a calculator. Copying the existing helper fixes only setup.
The initial full audit then rejects new trial9: its full norm is0.5A but
`basis.T@full` gives0.5000000000000001A. Native E/F/stress/validation and fixed
substrate all agree at all12records; this one round-trip guard failure is
retained in `native_roundtrip_diagnostic.json`, not rounded away. An initial
diagnostic print also fails on NumPy bool JSON serialization; explicit bool
conversion fixes the diagnostic only.

The new read-only audit reconstructs the **original internal optimizer
proposal** by strict, original-order cache-only replay. Its unchanged strict
guard must pass and `basis@original_internal` must equal the recorded full
coordinate bitwise. No epsilon, inward move, raw-result mutation or physical
domain change is permitted. The round-trip failure remains explicit. A
regression rejects even a one-ULP altered full coordinate or inadmissible
original internal proposal. Recorded path aliases are accepted only after
exact real-directory identity inside the permitted arm is established.

An audit replay launched with default-login thread variables unset does not
reproduce the old exact coordinate sequence; enforce the original
OMP/MKL/OpenBLAS=1 environment before replay. Optional threadpoolctl is absent;
nothing is installed. A subsequent standalone solver import lacks its parent
package; loading the pinned solver as `vcneb.feasible_relaxation` repairs only
the auditor's import, with a dedicated regression. Final raw-source and
both replay checks pass. No running worker, native result or physical
parameter is repaired; every failed audit is a tooling failure, not a new
material failure or a justification for rerunning DFT.

Next finish and audit the original preserving face, retain failed nodes as
holes, and settle a training-only valid domain/representation before any
prospective strong-control predictions. Original totals remain T29/81,
preserving30/81, with8independent slots/face and physical+.005sealed. No
unregistered third Hessian center, new material or new full path is added.

Final validation:112local and112independent related tests across8modules,
27new,0failures/errors/skips; not a full-suite claim. The tested source tree
is633650930f108bd6ce58a85456a90a4d1a1e9553. Raw terminal proof/cache/diagnostic
bytes survive independent Git extraction. The05:37:54CST snapshot has the
original preserving job still RUNNING:22/25processed,17gradientcandidates,
holes13/14/19/20/21 and139complete SCFs. This is not a final surface audit.
