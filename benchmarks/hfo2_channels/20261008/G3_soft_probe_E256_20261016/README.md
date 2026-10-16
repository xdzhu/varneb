# E256: four finite-amplitude statics along the measured omitted direction

Planned2026-10-16; submitted2026-10-17 as **28912687** on hfacnormal01.
This is the explicitly bounded follow-up to E254/E255, not a rerun of the
NEB, a replacement independent test, a new Hessian center or a stable3D PES.

## Fixed design and accounting

Four literal full-coordinate structures use the original preserving Q0
anchor plus t times the frozen E255 soft direction, in the order
+0.05,-0.05,+0.10,-0.10A. The original Q1/Q2 and substrate control remain0.
No B B^T reconstruction, structural repair or new relaxed center is used.
All six physical input assets and the binary are unchanged:
ABACUS/PBE/100Ry/full10auDZP/Gamma2x2x2/P0/E0. Native inert serialization
passed with maximum radius0.326290072A and minimum separation2.066735547A,
inside the original0.5A/1.6A guards. This is not another electronic setting.

Charge all four as new training-structure slots, although Q is unchanged:
preserving35 ->39/81, remaining42; T34/81 unchanged. At most four new
attempts, one per structure; failures remain charged, no retry/rollover.
The original native center is reused, never recalculated. Original296/E181
272 derivative calls and all earlier preparation remain additional.
The original16 local test slots and physical+0.005 condition stay sealed.
The original2D necessary-sign refusal and E185/E188 forecasts are retained.

Both original quadratic predictions were stored before these labels.
The new energies/gradients will check finite-amplitude curvature and work
consistency, not be fitted and relabeled as original independent success.
All signs/amplitudes remain regardless of agreement. No global branch,
full stationary index/TS, error-bound or predictive-advantage claim follows.

## Actual execution gate

Native preflight checked original E161 runtime380 files/binary/six assets,
reparsed the cached native anchor and wrote four inert initial decks. It
did not start DFT. Job28912687 was submitted held, resource-checked and
released on the same handle. One node/32 CPUs, genuine MPI32 per SCF,
OMP/MKL/OpenBLAS1, one-hour walltime, at most two committed study allocations.
The successful original ABACUS environment and mpirun command are reused.
No running source was overwritten or software installed.

The initial literal held check rejected Slurm's NumNodes=1-1 spelling.
Read-only inspection confirmed both bounds1, squeue nodes1, ReqTRES node1
and NumCPUs32. held_validation.json records that exact normalization;
no resource request was changed, broadened or resubmitted. Actual status
must be read from Slurm; COMPLETED alone would not certify the four SCFs.

- config.json:f7dc1053bff96bbd141d0a7e2e5c14d949f91f69832852d61bf7af9b932de7a7.
- bundle_files.json:ed53bc642ad6d7c1e3255a3d27a0648d8c0d93ae33b1e60687cf59b416dbc8b2.
- native_preflight.json:e05b34d2f9d90d65d6b701017f857c4a0009060c78daadfcc046b42f2105b005.
- accepted.json:bfed8d352ea4d3f4b50217d338cd4a9791c4a84d8f4084bafcbbb92a4039c3b6.
- held_validation.json:b133b031ab286d2c374d695ecd47c6e1f33e1fab622b5c8921c4ba0925a73bd2.
- released.json:183697b23e8356b9df16e5021fc3d6193936cd3cd481e25c508df5a4a4ae6c8f.

HF namespace: /public/home/iai806/abacus/agent-runs/20261016-varneb-G3-soft4-E256-r1.
Its bundle is immutable; probes/ is fresh output. This case is site-specific
provenance, not a recommendation to copy hard-coded cluster paths for users.

## Tests and the preserved first failure

24 local/24 independent related checks, three modules/eight new, zero final
failures/errors/skips; not a full repository run. They check exact literal
points/predictions/budget, immutable policy, refusal of modified targets and
local/unapproved execution, plus existing actual-curvature diagnostics.
These are not four completed DFT calculations.

The first independent export had23passes/one metadata failure because Git
archive converted protocol.md LF to CRLF. Local and staged protocol bytes
were9bc8300c..., whereas the exported file was2e3b999c.... The new explicit
protocol.md -text -eol attribute preserves the **same original bytes**;
no worker algorithm, config, forecast, physical input or strict hash check
was changed. First failure XML and source archive remain available, not
discarded. The actual submitted bundle/config hashes are unchanged.

Final tested tree33df68aacd630c67c7b20cf7ddbc483882b452f5, archive
736eb8b812ec9015a2c92601224b683709064a78b99b0d5fbd2ad6142dda1ca0.
First tree7b85774ac795d959410720c56636c1f19a0d56e3, archive
fe57866d78b6f945538cfba65a0e75ae38da89875ed79e47b89c49955e22656e.
Scratch and exports are on C:, not the nearly full E: drive. No licensed
pseudopotential/orbital or executable bytes are redistributed.

```text
python worker.py prepare --bundle THIS_CASE --config FRESH_CONFIG.json
python -m pytest -q tests/test_hfo2_soft_probe4.py
```

prepare is inert; native preflight and execution require the exact fresh
HF namespace and original runtime, with Slurm/accepted-handle guards.
Source preparation/held admission is progress, not finished material
validation or JCTC readiness. The full research goal remains incomplete.

Initial metadata snapshot,2026-10-17 00:15:34CST: RUNNING/node9/32CPUs,
one completed result record, no failure records, producer completion absent.
This is not the final four-record native audit or completed material result.

## E257 completed result (2026-10-17)

Job28912687 completed normally (0:0,478s,node9,32CPUs). The read-only
postprocess.py reparsed all four new native SCFs and the cached center,
verified the original runtime/binary/six physical assets and genuine MPI32,
and matched the producer records. native_complete_audit.json is the final
audit; the initial RUNNING snapshot above remains an unmodified historical
record. Native audit SHA256:
41321300ae00dc46a526ea14069fb541f6d02811fcac8dcccdf5fda72ace9928.
analysis_E257.json SHA256:
761bb4bd2784fca366cfbb32abc4800647b8545e8281b5ccd03cf93028ceb318.

| Signed amplitude (A) | Energy change (meV/cell) | Directional gradient (eV/A) |
|---|---|---|
| +0.05 | -0.050532852 | -0.000841328 |
| -0.05 | -0.050260353 | +0.000840951 |
| +0.10 | -0.068631462 | +0.000135813 |
| -0.10 | -0.067975881 | -0.000136940 |

The four-formula-unit cell gives only0.0126--0.0172meV/HfO2 lowering.
This is not a significant whole-path barrier reduction. At amplitudes
0.05/0.10A, energy secants are -0.040317/-0.013661eV/A2, whereas gradient
secants are -0.016823/+0.001364eV/A2. These are distinct finite-amplitude
quantities outside a harmonic profile. The outer restoring gradient and
an energy below the origin can coexist in a shallow anharmonic well;
neither a stationary minimum nor an energy-force inconsistency is proved.

Opposite-endpoint work nearly cancels and cannot validate center-to-wing
lowering by itself. Retrospective adjacent-interval trapezoids leave
0.029245/-0.029491meV/cell residuals next to the origin and less than
0.000461meV/cell on the two outer intervals. These are explicitly post-label
diagnostics, not a new blind acceptance rule or a total physical-error
bound. The origin, coefficients, frozen predictions and signs were not
changed to close these residuals.

The first read-only native audit refused producer paths spelled with
bundle/../ even though the strict resolved files, hashes and native records
matched. first-native-audit-failure.txt preserves that refusal. The second
auditor accepts only a strict same-file resolution inside the declared
namespace and preserves the producer's lexical source spelling; no field,
hash or numerical check was removed and no DFT was repeated. The original
and corrected remote audit sources are separately retained.

New cost:3.756877 completed-call proxy core-hours and4.248889 allocation
core-hours. All prior costs remain additional. Training stays T34/81 and
preserving39/81. The original16 independent points are uncomputed, the
physical+0.005 holdout is sealed, and the original2D sign-admission refusal
is retained. A third-coordinate validation would require a separately
approved pre-label protocol, not a retrospective success label.

Manuscript methods/results/discussion and AppendixE.7/TableA17 now report
these actual measurements. See results_validation_E257.json for separate
result-test and PDF receipts; setup validation.json above is not rewritten.

Result checks:157 local/157 independent-export related checks,14 modules,
10 new, zero final failures/errors/skips; not a full repository regression.
The initial local153pass/4fail XML is retained: manuscript wording/count
guards still expected the old pending text and148 rather than154 table
lines. All148 old table lines remain; only new measured TableA17 adds six.
No physical tolerance was changed. The independent CLI reproduced the
analysis JSON byte-for-byte. Tested tree8e067efa7e764af15c8c9da1710c53a5581f062c,
archive1de64b5fed6beaf014d9c03de2e39b62e5564d2d3ff6e69eb2fca6edda2b6d41.
The 38-page working PDF was compiled and all pages visually reviewed in
contact sheets, with individual inspection of new pages9,15,38; the receipt
distinguishes this from a final submission-quality or scientific sign-off.

```text
python postprocess.py analyze --audit native_complete_audit.json --audit-sha256 41321300ae00dc46a526ea14069fb541f6d02811fcac8dcccdf5fda72ace9928 --bundle . --output FRESH.json
```
