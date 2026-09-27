# GaN 45.7 GPa five-backend remote source index (working)

Read-only index captured on 2026-09-27 from `ssh hf` and Slurm `sacct` for
the **completed** production jobs. It identifies original result trees and
pins representative input/output plus whole-chain files. It is **not** a
portable raw-data deposit or a per-image SCF/force/stress audit. In
particular, the image-15 files below demonstrate where the dominant peak's
raw evidence lives; they do not prove the other 28 images were individually
checked. The compact plotted values are in
`gan_45p7_multibackend_vcneb_20260924.json`.

Let `R=/public/home/iai806/abacus/agent-runs/20260921-varneb-material`,
`M=$R/cases/gan`, and
`V=/public/home/iai806/abacus/agent-runs/20260918-varneb-vasp-hf/gan_b4_b1_vasp_pbe_paw_qian`.
Each Slurm record was queried with `sacct -X --format=JobID,JobName,State,WorkDir,SubmitLine`;
the long `SubmitLine` fixes the endpoint, profile, worker-count, and source
paths for QE/ABINIT/CP2K. Jobs must not be reconstructed from a neighboring
failed or superseded run.

| Backend | Slurm job | Exact result tree on hf | Job state |
| --- | ---: | --- | --- |
| ABACUS | 27760826 | `$M/abacus_vcneb_45p7_mpi4_20260922_retry1` | COMPLETED |
| VASP | 27719610 | `$V/vcneb_hf_n29_w9_mpi32` | COMPLETED |
| QE | 27770529 | `$M/qe_vcneb_45p7_pd04_29i_20260924_retry2` | COMPLETED |
| ABINIT | 27770687 | `$M/abinit_vcneb_45p7_pd04_29i_20260924` | COMPLETED |
| CP2K | 27770714 | `$M/cp2k_vcneb_45p7_k4_cut800_29i_20260924_retry1` | COMPLETED |

## Core file SHA-256

All paths in this table are relative to the backend result tree above. The
`vcneb.traj` hash is for the remote trajectory, **not** the independently
archived CP2K `gan_cp2k_45p7_chain_step_0052.traj` hash
`790599c1e8d4a700067fd5a5c5cdd2cf5ca4ff2e5ad572d37b798fb6a4eb9fc1`.

| Backend | `vcneb_summary.json` | `vcneb_preflight.json` | `image_worker_manifest.jsonl` | `vcneb.traj` |
| --- | --- | --- | --- | --- |
| ABACUS | `547bd20822ee418a1c4c96063f414b575ef533ce1a9a463e816764868684cfee` | `8f0363330046597617d6f431c40a32af6cb895f4a180a1607eaf91ced0bef962` | `d31d2fcd3e55cd04de90b7530083548992e44ade4b746e3975a0a62dfa392772` | `a6edbbef515500f770f28af9ef102c4eecc0dd76f5cfe944d977d169c3e85e3d` |
| VASP | `bb01f0e6c10d3c4f4583fcc7859d5d93e11c3080062797189b61c899d508c16d` | `a9367ed806ca7ed64922a85e2e830b22981a7f1e78c185a623be351499d72a8c` | `dcf67dcb9c6d5ba90c03a47db45016ddea9877ac1238eb46374e9dbb03f519cc` | `106b9e7d2b95d91d027c2387af1f478092fc31030e5653305869e7661d4b6439` |
| QE | `c3b97ca9630d29f54da3395efa1e1a3ae9484a57980cb89a44c34a77887c6123` | `6445bf0b707f5509708719f39e4c0cd9aab3827c4a3edc07010d92af389592aa` | `4040a46531587fca0995f1b26ded1c1183c89d7e3019b73b8d22a648e1eccffc` | `d4a359a9327e5c87a84d54e6f378dab572f637667a7cb01619c806ed4915709a` |
| ABINIT | `8a2d292df08048dc4c0607b62a0435537a743fe4feda24e199d6d2e11a589c1a` | `d0307eb031b1d8b0f8285a5e5c0ac8efb73b597b18ab964c03696eb87695d013` | `6277cc71b231c33746892059e5a529779f615a9036f1f8c01a1dfc06289b2fce` | `cb9f86260e874733e78dfca9a31a1ee72f27cec897c1a8e5caa923b312a30416` |
| CP2K | `7a10b32340e36c32fac1a5f31921b0793eb115717a83ace8e66d7f6f9f85d0c1` | `fdb1fc2fe5961ffdd0e130ff4648f3f6823c73182927cf5532f78646a2ef0dc1` | `d46b20ab7d6c0bce7b8f1caeaec70dcd606ff70a3d61253c6242fbbc84d23669` | `0bd431e4dacfe0a9cf94d50cd2fe3a256bb377d7fba83cbbc23a20535393ffc1` |

## Peak image 15: original calculator files

The following hashes refer to the **last visible files** under each image-15
directory. The final-chain trajectory/snapshot and worker manifest still need
to be joined to the exact calculator call before these can certify the final
image-15 point. `sha256sum` was run on hf; no VASP POTCAR is redistributed.
The geometry sanity check in `gan_45p7_image15_geometry_sanity_20260927/`
shows ABACUS `POSCAR.final` and VASP `POSCAR` matching final image 15.
The QE/ABINIT/CP2K `structure.start.vasp` files differ by +0.282226,
+0.210699, and +0.278682 Å³ in volume, **but these are deliberately
initial-chain snapshots written once by `attach_image_calculators`, not the
last calculator inputs**. These differences cannot establish whether their
`.pwo`, `.abo`, or `cp2k.out` match the final chain. That requires parsing
the actual `.pwi`, `abinit.in`, and `cp2k.inp` geometries and paired outputs.
The former geometry check has now been done in
`gan_45p7_final_peak_inputs_20260927/`: all three actual inputs match the
final image-15 geometry within `1e-9 Å`. The paired output-value audit is
still open.

| Backend | File in result tree | SHA-256 |
| --- | --- | --- |
| ABACUS | `15/INPUT` | `dc6684ffa709bf3d953c647272589b05d844c1293477bd144226d495207cb9df` |
| ABACUS | `15/KPT` | `bda588d9f8664a3465bb3ec87b2d648db1cc2c2998199122d3b2ad71f4a4bc3d` |
| ABACUS | `15/STRU` | `3abddc183471053e595e3d8343d637a7cacec1a983b7022a379780b2ac599de8` |
| ABACUS | `15/OUT.ABACUS/running_scf.log` | `4a3bcf6b174ce2417100691c5e374f3564c839d6dbe165f0ac8802d466ed9e12` |
| ABACUS | `15/Ga.upf` / `15/N.upf` | `e1aad2af1f8e69635634fd601d945943e033a36a1bb404cd4a1134e67180de7d` / `153a989f135fffe2e0f091bad0bbd724c0df5c12042ad3365f24feb0a426c529` |
| ABACUS | `15/Ga_gga_10au_100Ry_2s2p2d1f.orb` / `15/N_gga_10au_100Ry_2s2p1d.orb` | `8d6e55bcac77e743bf95e4992543214fcdb74252f9b7eab2c356767191527c98` / `e58669063b2c9e941a58a965292b64fa0cc9da428f5d880bb98c6571400996dc` |
| VASP | `15/INCAR` / `15/KPOINTS` / `15/POSCAR` | `83ba34d4b14ff4ea641bd74dec26e462ea1cc24b575db79a07138ccdfd552772` / `b5215f3e7608c27f49d45e8129d3edca2cfaa3ba88b8c389d4b52604715712ee` / `3e323ad457fb458106c1e67de5ff5f747e8229013be6b0ed5ef09f069f425865` |
| VASP | `15/OUTCAR` / `15/POTCAR` | `e97fc5321d373436270087c6a365015dcd58e0ab44ae40d1e4829b560f761c04` / `f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168` |
| QE | `image_0015/espresso.pwi` / `image_0015/espresso.pwo` | `b558fff1f3b16379a009f905514fa92c141062f24bb2859d1b323e57b190c767` / `8c151175c2e85e087da3181a2c43bd403eec05e6855a5a48f5bd4003a7659cb8` |
| ABINIT | `image_0015/abinit.in` / `image_0015/abinit.abo` | `b4d4388996655c4d13aae3c786dd3bd15f01c81ab13a23198a2d74e13b2069f3` / `eb34bbdb7435490057eb51cfe67959449de96ccd06caa299bae1b4ba5be30100` |
| CP2K | `image_0015/cp2k.inp` / `image_0015/cp2k.out` | `935170d769aa368ddc8f1bc3089223953b275ea5b52bc70168b84285ab93cff1` / `439ceac169eb5186e09afa3f51123f15a5147e29070e337c9fc2d12d4044607e` |

## Endpoint and profile identity

The QE/ABINIT/CP2K parameter/factory JSONs are under
`$R/varneb_guarded_source/examples/material_profiles/`. The ABACUS source
and Slurm script are under `$R/varneb_abacus_45p7_source/`; the VASP Slurm
script is under the parent of `$V` in `cluster/`. Structure hashes refer to
the exact `B4/CONTCAR` and `B1/CONTCAR` paths in the `sacct` submit lines,
except VASP, whose paths are `$V/vcneb_input/initial/CONTCAR` and
`$V/vcneb_input/final/CONTCAR`.

| Backend | Initial B4 / final B1 structure SHA-256 | Explicit profile or source SHA-256 |
| --- | --- | --- |
| ABACUS | `61fb33013240bbbc33cd3b87913b1221e53b087ea098c17e381d2c6af3639a32` / `594c0fa9ed6a8e334d3523a10365072377134f81e3241f924432350f75a35d68` | `varneb_abacus_45p7_source/examples/run_vcneb_abacus.py`: `0292a969b7a36b76539758a47a965169a1f2a43d3cf698886cad31bc77c3995c`; Slurm script `c69a3e014b3cfac799c534ec4feab78abf18a6fd9744f6d5a79f80836d255402` |
| VASP | `6417133891af3799dc309d6c7c8c2942cc9db97885cfc48a3918aeb62d680c24` / `0ddfa9f0bda4888b31a9004b1491386a1d5ad88f968ef41d5a55a629c9b3ee7d` | `vcneb_input/initial/INCAR`: `7495a37adf2971fddbf02e6242bef5d5b7009e8a58e0d329e895e5f706780ae7`; KPOINTS `9c462028dce835e55efda9471dc4b03f15773413a90fcb68cc912736a640fd1f`; POTCAR hash above; Slurm script `2fef391a95549d64703589850bd2cc52f7229019c2a0735886a6ccb088b28951` |
| QE | `65ecaaf8540decb8c8ca10fb6ddecdde6c9c955363cf43ee51e59c56ad424bcd` / `8de54c2b6fd7b894a0b67078bcca123ef3ee508dde2bba7152e930707f70ed0e` | parameter JSON `456e7f6a0ea3cab3d25a5d2afd64122445ffd64f5e89a30b4040e14bfe712ad7`; factory JSON `01604befb28a231e708fbfb00db45e4999ff5c3345ab048bb11770f30a03447e` |
| ABINIT | `27274061bd12e6b1db693170b3a60103731c30dad7d072b277d3118cadbb973b` / `3c996227b5a2814b96026216f49e21129ae63bdc9018ca996d018a8fd9f86d22` | parameter JSON `8c92a8b30ecac9a6ccfee81aac3dbd893c1f20c836c2d22597c14c40c76fd990`; factory JSON `8967fec1f5e76a81e989fad51cdac09d7c762a0ddb17432ec2167b08ae5d287c` |
| CP2K | `7411c5ff4d1e29132f97c1e13818f032378b32f73d4e27898b1e64ef67e5197c` / `187402f1eca823878a33732369520c0ab263dce6d6c7fa3a9102ab92bc59de9f` | parameter JSON `e274b93390b400289fba6c84fd7cf56d58fe8c6fa640d9b152f5bb2795af5a17` |

The complete remote trajectories have since been independently re-read and
their last 29 images archived under `gan_45p7_final_chains_20260927/`;
`gan_45p7_final_chain_audit_20260927.json` proves that the plotted curves
are reproduced from those evaluated chains. Next raw-output gate: locate
an actual calculator input/output pair matching **each exact final image
geometry** (including endpoint static jobs). Peak actual-input geometry is
already matched for QE/ABINIT/CP2K, but paired raw energy/force/stress and
the other images are not yet independently audited. Only if records do not exist
should a hash-pinned same-protocol static re-evaluation be considered.
Then compare raw energy, force, stress, SCF completion and atom/cell identity
to the chain. Do not treat the final trajectory audit as this raw-log audit.
