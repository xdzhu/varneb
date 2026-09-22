# VARNEB 多后端生产修复状态

更新时间：2026-09-22（HF `hfacnormal01`）

本文只记录已由日志、作业状态或回放实验支持的结论。旧失败目录保留，修复使用独立目录，不覆盖原始轨迹。

## 统一判定边界

VARNEB 控制器只依赖每个 image 的 `energy`、`forces` 和完整 `stress`。后端适配器负责输入生成、外部程序启动、输出解析和 SCF 诊断；FIRE、BlockFIRE、SplitFIRE 等路径优化器不负责修复后端输出，也不把 SCF 重试伪装成路径迭代。

失败分类现在进一步按后端契约分层：ABINIT 的 `chkorthsy` 归为
`abinit_symmetry_failure`，缺失或未固定伪势归为 `pseudopotential_contract`，
能量/力/应力区块缺失或解析失败归为 `output_contract_violation`。这些分类在
`MPI_ABORT` 或 launcher 文本之前判定，避免把输入/输出契约错误误报为 MPI 故障，
也避免错误地修改 FIRE、NEB 图像数或路径步长。

生产路径仍使用普通 VC-NEB `fmax = 0.10 eV/A`、固定端点和无 CI。最终是否收敛以完整链的 `final_max_generalized_force_eV_per_A <= 0.10` 判定，不能以 Slurm exit code 或某一轮回弹判定。

HF 的 ABINIT 8.6.1 启动契约也已固定：该二进制链接 Intel MPI/Fortran
2017.4/2017.5，必须加载匹配的 `compiler/intel/2017.5.239` 与
`mpi/intelmpi/2017.4.239`，再由 `mpiexec.hydra -bootstrap slurm -n N abinit`
启动。HF 默认 Intel 2021 环境的 `srun ... abinit` 会向该二进制传入
`pmi_args`，而仅切换到 `mpiexec` 仍会因缺少 `libifport.so.5` 或 PMI ABI
不匹配而失败。独立 32-rank canary `27749997` 在匹配的 compiler/MPI 与
Hydra 组合下返回 `abinit 8.6.1`；对照 `27749996` 的 `srun --mpi=pmi2`
仍以 `PMI_KVS_Get returned -1` 失败。此前 `27749797/98`、`27749871/72`
的目录和日志均保留，未将启动层失败误判为 SCF 或 NEB 发散。

ABINIT 的伪势也是输入契约的一部分：包含 `pps` 的参数必须同时提供存在的
`pp_paths` 目录。缺少该目录会在 ASE 写入输入前给出明确的
`ABINIT parameters with 'pps' require explicit pp_paths`/目录不存在错误，
不再等到 32 个 rank 启动后才产生含糊的伪势解析失败。GaN HGH-LDA 试跑
`27750045/46` 正是捕获了这一缺口；补充 `abinit_hgh_factory_kwargs.json`
后的 `27750056/57` 进入真实端点计算；其中 B1 暴露了自动对称性中止，B4
随后为保持端点契约一致而取消，二者目录均保留。

另外，ABINIT 后端默认写入 `nsym=1`（恒等对称操作），避免由变胞路径或
POSCAR 末位舍入触发自动对称性识别的 `chkorthsy` 中止。需要利用更高对称性
时必须显式提供并审核 `nsym/symrel`，不能让生产路径隐式依赖自动分类。

## ABACUS GaN：结果解析契约修复

### 根因

GaN ABACUS 作业 `27741518` 在 step 0 的 image 5 报出 NumPy 不规则数组 `(30,)`。回放确认，ASE-ABACUS 3.23.1b1 在构造完整结果字典时会提前解析 eigenvalues；该 image 的 30-k-point eigenvalue 区块形状不规则，而 VCNEB 实际只需要能量、力和应力。

修复版第一次重跑 `27744183` 绕过 eigenvalue 后，image 5 的 header 又出现并行输出格式损坏，ASE 私有 header 正则在 `SELF-CONSISTENT` 区块抛出 `NoneType.group`。这仍是输出解析问题，不是 FIRE 或路径几何问题。

### 处理

- `vcneb/abacus.py` 新增只解析 VCNEB 三项的结果入口，避免可选 eigenvalue 解析。
- 完整 ASE 解析失败时，使用严格的最小日志解析器读取最后一组 `final etot`、`TOTAL-FORCE` 和 `TOTAL-STRESS`；结果仍由统一的 `ImageEvaluation` 检查数量、形状和有限值。
- 对旧损坏日志的回放结果为：energy `-4726.0740513 eV`、forces `(4,3)`、stress `(6,)`。

Guarded continuation `27749598` 在完整写出 step 0--9 后，于 image 10 因
ABACUS 日志的人类可读 `final etot is` 行被长行截断而被误判为无能量；同一
日志实际包含完整的 `!FINAL_ETOT_IS -4726.5466304132096411 eV`、力和应力
区块。此前错误分类为 `mpi_failure`，根因是最小解析器只识别前一种能量标记，
不是 MPI 或晶胞信赖域问题。解析器现在同时接受两种官方输出标记，并增加了
该真实截断格式的回归测试；原 `27749598` 目录保留，续算使用独立目录。

新的独立续算 `27755927` 已从 `chain_step_0009.traj` 启动；两次仅参数错误
的提交 `27755906/24` 在 DFT 启动前退出，不计入计算失败。续算仍沿用
`maximum_cell_step=0.05`、9 个并行 image worker 和每像 32 MPI。
- 修复已加入回归测试并推送；ABACUS 标记解析的实现提交为 `ec10b2a`，后续
  契约分类修复在 `6fc573d`。

该续算已写出 `chain_step_0000.traj` 和 `chain_step_0001.traj`，控制器最大广义力
由 `1.694169` 降至 `1.525129 eV/A`。两步快照中的端点最大应力都约为
`484.97/468.50 kbar`，说明这组
输入端点没有通过 ABACUS 自身的变胞静态门禁；它只能作为并行执行与输出解析
诊断，不能进入物理生产矩阵。快照中间像的首步能量最高点约为
`1.0701 eV`，第二步约为 `1.0609 eV`（相对初端），但在端点未通过门禁前不作为能垒结论。

为防止同类问题再次发生，`hf_gan_abacus_vcneb.slurm` 现在默认要求
`ENDPOINT_STATIC_SUMMARY`，并在启动 VC-NEB 前调用统一的 `validate_ase_static_gate.py`。
只有显式设置 `REQUIRE_ENDPOINT_STATIC_GATE=0` 才能绕过；绕过的结果必须标记为
诊断而不是生产结果。`relax_abacus_endpoint.py` 的摘要也补齐了通用端点身份、
状态和能量字段，可直接接入统一静态摘要构建器。

旧的 `27744907` 是单像晶胞信赖域失控的已取消诊断目录；它不再是待等待的生产
任务。当前唯一的 ABACUS 生产续算是 `27755927`，从安全快照在独立目录运行，
在其完整链收敛前不把 ABACUS 宣布为生产级路径已通过。

## CP2K BTO：SCF 与资源启动分层修复

### 根因

原生产 `27741430` 在 BTO image 4 的 CP2K 内层 SCF 达到 500 次仍未收敛。输出显示 Broyden mixing `ALPHA=0.20` 下残差长期振荡；这是电子自洽问题，不是 NEB 候选步问题。

两个探针失败也已区分：

- `27744198` 申请共享节点 1 CPU，却在内部要求 `srun --exclusive`，一直等待作业 step，随后取消；属于 Slurm 资源请求错误。
- `27744265` 让 ASE 的 `cp2k_shell.psmp` 直接启动 32 MPI，shell 握手失败；CP2K shell 适配器必须用单 rank step。

### 处理

新增 `examples/material_profiles/cp2k_bto_pbe_dzvp_stable.json`：

- `EPS_SCF = 1e-5`，`MAX_SCF = 1000`；
- `ADDED_MOS = 30`，300 K Fermi smearing；
- `DIRECT_P_MIXING`，`ALPHA = 0.05`；
- 生产资源采用独占节点，多个 one-rank `cp2k_shell.psmp` image worker，而不是把 32 MPI 直接交给 shell。

本轮端点重跑进一步发现并修正了一个输入单位错误：ASE 的 CP2K calculator
把 `cutoff` 解释为 eV，而 CP2K 文献 profile 通常以 Ry 给出。旧的稳定 profile
裸写 `400`，实际只有约 29.4 Ry，足以产生异常大的 Pulay 应力并把 BFGS
推向错误的晶胞。现在 `make_ase_cp2k_factory` 接受显式的 `cutoff_ry`，在
后端边界转换为 eV；BTO/GaN 稳定 profile 均改为 `cutoff_ry: 400`，并拒绝
同时提供 `cutoff` 与 `cutoff_ry`。已取消的 `27749624/25/27/28` 只作为这
一错误输入的诊断证据保留，修正版将在独立目录重新准备端点。

单结构探针 `27744287` 已成功；两个相同的旧 image-4 结构均在 101 次 SCF 内收敛，退出码为 0。BTO 生产 `27744298` 随后因端点静态审计不通过、且路径力连续多步发散而取消；独立目录和日志均保留。

GaN CP2K `27741431` 也因端点静态审计显示固定端点基线无效而取消；原目录保留，只有重新生成并通过端点门禁后才允许新的 profile 续跑。

通用 ASE Slurm 模板现已在源头加入两项防护：

- CP2K 默认对每个持久 `cp2k_shell.psmp` 进程使用单 rank，不能把分配给 image 的 32 MPI 直接交给 shell 握手；
- 生产 VC-NEB 必须提供 `ENDPOINT_STATIC_SUMMARY`，并在启动路径前验证两个端点的力 `< 0.10 eV/A`、应力 `< 0.10 kbar`。先生成摘要时显式使用 `STATIC_ONLY=1`；绕过门禁必须显式设置 `REQUIRE_ENDPOINT_STATIC_GATE=0`。

这些是输入契约层的预防检查，不会修改或“修复”已有轨迹。

端点准备现在也有统一入口 `scripts/relax_ase_endpoint.py`：它通过与生产路径相同的 ASE factory，对单个端点执行可变胞 BFGS（或 `--fixed-cell` 原子 BFGS），保存 `CONTCAR`、断点、力/应力/焓摘要。端点准备、静态门禁和 VC-NEB 因而成为三个可审计阶段，而不是在路径迭代中临时改变后端参数。

### 运行中 ABACUS 链的几何诊断

`27744907` 的 `chain_step_0000` 到 `chain_step_0003` 显示 image 19 的体积从 `35.188` 增至 `39.020`、`43.093`、`47.353 Å³`，而相邻 image 保持在约 `35 Å³`。第 3 步日志的 `fmax` 仍为 `5.131 eV/A`。这不是一次普通回弹，而是单像晶胞信赖域失控；历史轨迹和作业继续保留用于审计。

控制器新增 `maximum_cell_step`：它按当前晶胞到候选晶胞的相对变形范数检查每一步，并在超过阈值时于 DFT 调用前拒绝候选。通用 HF 模板默认 `0.05`，配合 FIRE 的 8 次几何回溯；用第 2→3 步实际快照回放时，image 19 的相对变形为 `0.05548`、体积比 `1.09887`，已被门禁准确拒绝。这一检查与初始路径的 `maximum_deformation` 不同，专门覆盖运行时单像发散。

`examples/run_vcneb_ase.py` 和 HF 模板现在支持 `--resume-snapshot`/`RESUME_SNAPSHOT`。续算会先逐字节验证快照的初末端点身份，再在新工作目录中重新绑定 calculator；因此可以从最后一个安全完整快照续算，而不覆盖失控链。

## 其他当前结果

| 案例/后端 | 作业 | 当前判定 |
|---|---:|---|
| GaN / QE | 27741412 | 已完成；最终广义力 `0.091018 eV/A`，按默认阈值可接受；能垒 `1.638143 eV`。旧 summary 未写显式 `converged`，审计以最终力为准。 |
| GaN / ABINIT（旧路径） | 27741516 | `max_steps_reached`，最终力 `9.537243 eV/A`，且 HGH-LDA 端点应力约 `-2791 GPa`；未收敛，不作为生产证据。 |
| GaN / LAMMPS | 27741574 | 最终力 `0.095441 eV/A`；仅是经典 Tersoff 对照，不能与 DFT/PBE 能垒等价比较。 |
| BTO / QE | 27741411 | 最终力 `0.052831 eV/A`，无内部势垒；按阈值可接受。 |
| BTO / VASP | 27741623 | 最终力 `0.095145 eV/A`，能垒 `0.042644 eV`；已收敛。 |
| GaN / ABACUS（新解析器重跑） | 27744907 | 已取消；运行至第 15 步，`fmax` 在 `3.17–3.27 eV/A` 平台并出现 `54190`、`26643 eV/A` 尖峰；image 19 单像体积发散，原目录保留为失败证据。 |
| GaN / ABACUS（guarded continuation） | 27749598 | 已写出 step 0--9 后因 image 10 的 `!FINAL_ETOT_IS` 解析缺口失败；原目录和完整快照保留，不作为物理失败。 |
| GaN / ABACUS（marker-parser continuation） | 27755927 | 已从完整 step 9 快照写出 `chain_step_0000`；端点应力 `484.97/468.50 kbar`，未通过 ABACUS 静态门禁，仅作解析/并行诊断，不进入生产矩阵。 |
| BTO / CP2K endpoint relaxation（旧单位错误） | 27749624, 27749625 | 已取消并保留；旧 profile 把 400 Ry 错写成 400 eV，诊断应力 1459–3788 GPa，不进入结果矩阵。 |
| GaN / CP2K endpoint relaxation（旧单位错误） | 27749627, 27749628 | 已取消并保留；同一 cutoff 单位错误，不进入结果矩阵。 |
| BTO / CP2K endpoint relaxation（400 Ry 修正版） | 27749725, 27749726 | 已完成；初端点力很小但应力约 `5.97 kbar`，终端点应力约 `4.35 kbar`，均未通过严格静态门禁。 |
| BTO / CP2K endpoint tightening | 27749774 | 已完成；初端点 `max_generalized_force=4.24e-4 eV/A`、最大应力约 `0.0094 kbar`，通过严格门禁。 |
| GaN / CP2K endpoint relaxation（400 Ry 修正版） | 27749727, 27749728 | 已完成；两端最大应力约 `6.84/4.82 kbar`，未通过严格静态门禁。 |
| BTO/GaN / CP2K endpoint tightening | 27755902, 27755903, 27755904 | 独立加严目录运行中，目标 `fmax=5e-4 eV/A`、应力 `<0.1 kbar`，不覆盖粗弛豫结果。 |
| GaN / ABINIT-HGH-LDA endpoint relaxation（错误 launcher） | 27749797, 27749798 | 已取消并保留；`srun` 注入 `pmi_args` 导致 MPI socket 等待。 |
| GaN / ABINIT-HGH-LDA endpoint relaxation（ABI 不匹配复现） | 27749871, 27749872 | 已取消并保留；Intel 2017/2021 混用，`mpiexec` 仍复现 `pmi_args`/socket 等待。 |
| ABINIT 32-rank launcher canary | 27749992, 27749995, 27749996, 27749997 | 仅诊断：缺 compiler runtime、错误 PMI 组合均失败；匹配 Intel 2017 + Hydra 的 `27749997` 成功返回 8.6.1。 |
| GaN / ABINIT-HGH-LDA endpoint relaxation（缺伪势路径） | 27750045, 27750046 | 预检后失败并保留；`pps=hgh` 未配 `pp_paths`，输入契约已补强。 |
| GaN / ABINIT-HGH-LDA endpoint relaxation（自动对称性试跑） | 27750056, 27750057 | B1 因 `chkorthsy` 末位晶格噪声失败，B4 为一致性取消；独立目录和日志保留。 |
| GaN / ABINIT-HGH-LDA endpoint relaxation（`nsym=1` 修正版） | 27750109, 27750110 | 已完成但静态门禁失败；端点体积塌缩至 `7.22/7.08 Å³`，应力 `59.3/60.8 kbar`。HGH-LDA 测试结果不进入 PBE/VASP 生产矩阵。 |

## 下一步

1. 审计 `27755927` 的 ABACUS 诊断续算；若再次失败，只看最小解析器的具体契约错误，不再改 FIRE 参数。生产 ABACUS 链必须先用同一套 100 Ry/DZP/PBE 设置完成端点门禁。
2. ABACUS `!FINAL_ETOT_IS` 解析修复已部署到独立源码，并从 `27749598` 的
   `chain_step_0009.traj` 续算；只在新目录中验证后续路径，不覆盖原始链。
3. 对 CP2K BTO/GaN 修正版端点先用 `STATIC_ONLY=1` 生成独立摘要，只有通过门禁后才允许新的生产路径，不覆盖已取消目录。
4. 对 ABINIT GaN 先完成 HGH-LDA 端点重建并通过同一门禁，再决定是否重跑；当前 `max_steps_reached` 结果不得进入生产矩阵。新的端点任务必须使用匹配 Intel 2017 + Hydra 启动契约。若要与 PBE 结果比较，必须另行获得并固定 PBE 赝势，不能把 LDA/HGH 结果标成 PBE。
5. 生成统一后端状态表、路径图和文献比较数据；未达到阈值或非同一物理模型的结果不得进入“已验证生产矩阵”。
