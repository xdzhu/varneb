# GaN image-15 actual calculator input geometry

These are the byte-preserved **actual input files** from the visible image-15
calculator directories, unlike the one-time `structure.start.vasp` snapshots
in `../gan_45p7_image15_geometry_sanity_20260927/`. The source hashes
were obtained on hf; ASE's `espresso-in`, `abinit-in`, and `cp2k-restart`
readers were used to compare each input with the corresponding archived
final-chain image 15.

| Backend | Original file | SHA-256 | Max cell difference (Å) | Max Cartesian position difference (Å) |
| --- | --- | --- | ---: | ---: |
| QE | `image_0015/espresso.pwi` | `b558fff1f3b16379a009f905514fa92c141062f24bb2859d1b323e57b190c767` | 4.9e-15 | 4.9e-11 |
| ABINIT | `image_0015/abinit.in` | `b4d4388996655c4d13aae3c786dd3bd15f01c81ab13a23198a2d74e13b2069f3` | 3.6e-15 | 2.5e-15 |
| CP2K | `image_0015/cp2k.inp` | `935170d769aa368ddc8f1bc3089223953b275ea5b52bc70168b84285ab93cff1` | 0 | 0 |

The actual input geometry **does** match the published chain peak for these
three backends. This reverses the earlier, incorrect inference from the
initial `structure.start.vasp` snapshots. It is still necessary to pair each
input with its actual output and independently compare final SCF energy,
forces and stress to the chain before claiming complete raw-output audit.
