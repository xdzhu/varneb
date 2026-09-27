# GaN image-15 saved geometries versus final chain

These five small VASP-format structures were copied byte-for-byte from the
image-15 directories in the remote result roots listed in
`../gan_45p7_remote_source_index_20260927.md`. **They have different
semantics.** ABACUS `POSCAR.final` and VASP `POSCAR` happen to match the final
chain image. For QE/ABINIT/CP2K, `structure.start.vasp` is deliberately
written once when `vcneb.backends.attach_image_calculators` first attaches
calculators. It is an **initial-chain snapshot**, not the input of the last
calculator call. Its difference from the final chain is expected and is no
evidence that the current electronic output was overwritten.

| Backend | Saved file | SHA-256 | `V(saved)-V(final)` (Å³) |
| --- | --- | --- | ---: |
| ABACUS | `15/POSCAR.final` | `f1050345c9750c8adb67bf3b82c4e75746e5d30dd95f75c7f00e295b27ed6e07` | 0 |
| VASP | `15/POSCAR` | `3e323ad457fb458106c1e67de5ff5f747e8229013be6b0ed5ef09f069f425865` | 0 |
| QE | `image_0015/structure.start.vasp` | `99ced52beab323983966e67bef15a76f29afda436d7841ccafe6d0334e389570` | +0.2822261453 |
| ABINIT | `image_0015/structure.start.vasp` | `eee69c323aac34a66eda681442766d44dbf62a229a4661766cfa4fe7b43250d5` | +0.2106990786 |
| CP2K | `image_0015/structure.start.vasp` | `4c23aa4c6de89b1466467136d6242fa1741a89f6e15293a2ac7a811515953798` | +0.2786819613 |

That correct input-geometry comparison is now in
`../gan_45p7_final_peak_inputs_20260927/`: all three actual input files
match the final chain image 15. The copied `structure.start.vasp` files
could not decide that question. Paired output energy, forces, stress and SCF
completion remain to be checked before the raw-output provenance gate closes.
