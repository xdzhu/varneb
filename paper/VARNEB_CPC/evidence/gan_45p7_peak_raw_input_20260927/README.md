# GaN peak-image raw directory versus archived final chain

These five small VASP-format structures are copied byte-for-byte from the
**last visible** calculator directory for image 15 in the remote result roots
listed in `../gan_45p7_remote_source_index_20260927.md`. They are not
interchangeable with the final chain's image 15. Read them with ASE and
compare against `../gan_45p7_final_chains_20260927/`.

| Backend | Remote source | SHA-256 | `V(raw)-V(final)` (Å³) |
| --- | --- | --- | ---: |
| ABACUS | `15/POSCAR.final` | `f1050345c9750c8adb67bf3b82c4e75746e5d30dd95f75c7f00e295b27ed6e07` | 0 |
| VASP | `15/POSCAR` | `3e323ad457fb458106c1e67de5ff5f747e8229013be6b0ed5ef09f069f425865` | 0 |
| QE | `image_0015/structure.start.vasp` | `99ced52beab323983966e67bef15a76f29afda436d7841ccafe6d0334e389570` | +0.2822261453 |
| ABINIT | `image_0015/structure.start.vasp` | `eee69c323aac34a66eda681442766d44dbf62a229a4661766cfa4fe7b43250d5` | +0.2106990786 |
| CP2K | `image_0015/structure.start.vasp` | `4c23aa4c6de89b1466467136d6242fa1741a89f6e15293a2ac7a811515953798` | +0.2786819613 |

The QE/ABINIT/CP2K last visible work directories were overwritten by a
geometry different from the converged chain's peak. Their current `.pwo`,
`.abo`, and `cp2k.out` files therefore **cannot** be used as direct raw-log
proof for the published final image-15 energy/force/stress. This observation
does not invalidate the saved trajectory, result cache, or barrier, but it
does keep the raw-output provenance gate open. Resolve it by locating an
immutable calculator log matching each final chain geometry, or by a
hash-pinned static re-evaluation of the final chain under the same backend
contract. Future runs should retain immutable per-evaluation input/output
directories rather than only a mutable `image_0015` directory.
