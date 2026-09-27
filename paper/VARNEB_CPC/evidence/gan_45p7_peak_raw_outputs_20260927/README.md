# GaN 45.7 GPa image-15 original output check

The three original, byte-preserved electronic output files here come from
the exact hf result roots recorded in
`../gan_45p7_remote_source_index_20260927.md`. Their SHA-256 values match
that source index. The paired **actual** QE/ABINIT/CP2K input geometries are
archived in `../gan_45p7_final_peak_inputs_20260927/` and match the final
trajectory peak geometry.

`scripts/audit_gan_image15_raw_outputs.py` reads the outputs and the three
archived final trajectories. QE and ABINIT both have SCF convergence and
normal-completion markers; ASE output parsing reproduces final-chain
image-15 energy, all four atomic forces, and six stress components exactly.
This verifies raw-file-to-trajectory consistency but is not an independent
electronic parser or a complete 29-image audit.

The CP2K `cp2k.out` is a **cumulative** text file with 53 `FORCE_EVAL`
energy entries. Its last visible energy is `−4600.267100295448 eV`, whereas
the final trajectory's image 15 is `−4600.2667338125 eV`, a difference of
`−0.000366482947 eV`. None of the 53 energies matches the final-chain energy
within `1e−8 eV`; the closest is record 28 at `−4600.266642524207 eV`,
still `+0.000091288293 eV` away. The last run has an SCF-converged marker,
but its `PROGRAM ENDED` marker is absent; the final visible completion marker
precedes that run. Thus the text tail is not a completed provenance record for
the final peak. This does **not** establish the production chain is wrong: the actual
last input geometry matches the chain, but the cumulative text file alone
cannot be joined to the exact final calculator return. We therefore also
froze the exact-geometry worker cache file
`cp2k_image15_worker_cache.npz` (SHA-256
`f9999926bb269008fe16815e032f1d628f1c337ae3ca3b40b4606b44cb20f1fd`).
Its filename is the SHA-256 cache key recomputed from the final image-15
species, positions, cell and PBC by `ThreadedCalculatorExecutor`; its energy,
forces and stress match the final chain exactly. This establishes the
serialized worker return, **not** independent raw CP2K electronic-output
agreement. The cumulative text-output discrepancy and its cause remain open;
do not relabel the last visible output as the final peak value.

The adapter formerly copied CP2K's output before closing the shell, which
could archive an unflushed tail. It now closes first and then copies the
output, and preserves an output already located in the image directory.
That prospective fix does **not** reconstruct the missing historical text or
explain the energy difference by itself. The historical CP2K raw-output gate
therefore remains open pending same-protocol source reconciliation; no cutoff
or other physical setting should be changed to force agreement.

The machine-readable result is `audit.json`. Endpoint static outputs, the
other 28 images, and executable/pseudopotential freezing remain open.
