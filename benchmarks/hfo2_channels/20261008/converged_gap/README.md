# First ordinary G1 T--PO residual pass

Slurm28300425 ended COMPLETED/0:0 at2026-10-08 23:13:32 CST, after1h26m26s.
The production summary terminates by force threshold, not a step limit.
The chain has10 total images(8 moving interior), P=0, ordinary NEB with
fmax0.10 eV/Angstrom and no climbing image. Electronic inputs remain the
historical ABACUS/PBE/100Ry/full10auDZP/Gamma2x2x2 contract.

The [frozen final observation](gap_converged_step06/observation.json) matches
all ten images to exact completed, input/geometry-audited SCFs. Neither fixed
endpoint was rerun. Its maximum force replays as0.05988161569025288 eV/Angstrom
on HF and0.059881615690252875 locally. The unchanged T/PO endpoints have equal
ordered periodic geometry to the registered variants.

At step6 the discrete maximum is image3. T→PO and PO→T maxima are respectively
33.888971 and115.210164 meV/f.u.; their difference is the endpoint energy
change-81.321193 meV/f.u. The largest residual belongs to atomic motion at
image2. See [descriptive analysis](analysis.json) and the original
[production summary](production_summary.json).

This is ordinary residual convergence, not a full-variable first-order saddle
certificate or a sampling/error-converged absolute barrier. Reference-mode
projections remain descriptive, not energy partitions. No extra static/CI
task was submitted on the strength of this pass; the rest of G1 is incomplete.

The observation was exported read-only from the R14 production namespace by
the unchanged checked exporter in the separate R19 analysis archive. Original
raw logs and run/source directories remain on HF. Hash-bound snapshots,
full numeric E/F/stress and a SinglePointCalculator trajectory are retained
here to reproduce the ordinary residual with no DFT calls.

Combined delivery validation uses a clean Git-tree archive containing this
trajectory, not the dirty working tree: 905 passed, 2 skipped in89.01s.
The [final verification receipt](final_delivery_verification.json) identifies
the archive, raw data hashes and unchanged analysis/API bytes also tested
against actual HF sources. Post-test additions are documentation/receipts only.
