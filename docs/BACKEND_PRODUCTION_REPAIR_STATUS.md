# VARNEB 多后端生产修复状态

更新时间：2026-09-21（HF `hfacnormal01`）

本文只记录已由日志、作业状态或回放实验支持的结论。旧失败目录保留，修复使用独立目录，不覆盖原始轨迹。

## 统一判定边界

VARNEB 控制器只依赖每个 image 的 `energy`、`forces` 和完整 `stress`。后端适配器负责输入生成、外部程序启动、输出解析和 SCF 诊断；FIRE、BlockFIRE、SplitFIRE 等路径优化器不负责修复后端输出，也不把 SCF 重试伪装成路径迭代。

生产路径仍使用普通 VC-NEB `fmax = 0.10 eV/A`、固定端点和无 CI。最终是否收敛以完整链的 `final_max_generalized_force_eV_per_A <= 0.10` 判定，不能以 Slurm exit code 或某一轮回弹判定。

## ABACUS GaN：结果解析契约修复

### 根因

GaN ABACUS 作业 `27741518` 在 step 0 的 image 5 报出 NumPy 不规则数组 `(30,)`。回放确认，ASE-ABACUS 3.23.1b1 在构造完整结果字典时会提前解析 eigenvalues；该 image 的 30-k-point eigenvalue 区块形状不规则，而 VCNEB 实际只需要能量、力和应力。

修复版第一次重跑 `27744183` 绕过 eigenvalue 后，image 5 的 header 又出现并行输出格式损坏，ASE 私有 header 正则在 `SELF-CONSISTENT` 区块抛出 `NoneType.group`。这仍是输出解析问题，不是 FIRE 或路径几何问题。

### 处理

- `vcneb/abacus.py` 新增只解析 VCNEB 三项的结果入口，避免可选 eigenvalue 解析。
- 完整 ASE 解析失败时，使用严格的最小日志解析器读取最后一组 `final etot`、`TOTAL-FORCE` 和 `TOTAL-STRESS`；结果仍由统一的 `ImageEvaluation` 检查数量、形状和有限值。
- 对旧损坏日志的回放结果为：energy `-4726.0740513 eV`、forces `(4,3)`、stress `(6,)`。
- 修复已加入回归测试并推送，当前代码提交为 `1812352`。

新鲜生产重跑：`27744907`，目录 `cases/gan/abacus_vcneb_production_minparser`。在它完成前，不把 ABACUS 宣布为生产级路径已通过。

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
| GaN / ABINIT | 27741516 | `max_steps_reached`，最终力 `9.537243 eV/A`，未收敛；暂不作为生产证据。 |
| GaN / LAMMPS | 27741574 | 最终力 `0.095441 eV/A`；仅是经典 Tersoff 对照，不能与 DFT/PBE 能垒等价比较。 |
| BTO / QE | 27741411 | 最终力 `0.052831 eV/A`，无内部势垒；按阈值可接受。 |
| BTO / VASP | 27741623 | 最终力 `0.095145 eV/A`，能垒 `0.042644 eV`；已收敛。 |
| GaN / ABACUS（新解析器重跑） | 27744907 | 已取消；运行至第 15 步，`fmax` 在 `3.17–3.27 eV/A` 平台并出现 `54190`、`26643 eV/A` 尖峰；image 19 单像体积发散，原目录保留为失败证据。 |
| GaN / ABACUS（guarded continuation） | 27749598 | 已从 `snapshots/chain_step_0000.traj` 在独立源码/工作目录启动，启用 `maximum_cell_step=0.05` 和 8 次回溯；当前运行中，尚未产生第一步摘要。 |
| BTO / CP2K endpoint relaxation | 27749624, 27749625 | 初末端点独立 BFGS 准备任务，使用稳定 PBE/DZVP profile；运行中。 |
| GaN / CP2K endpoint relaxation | 27749627, 27749628 | 初末端点独立 BFGS 准备任务，使用稳定 PBE/DZVP profile；运行中。 |

## 下一步

1. 等待并审计 `27744907` 的 ABACUS 生产结果；若再次失败，只看最小解析器的具体契约错误，不再改 FIRE 参数。
2. 用 `snapshots/step_0000` 的完整快照在独立目录启动 guarded continuation，启用 `maximum_cell_step=0.05`；当前 `27744907` 已取消，原目录只作为失控证据保留。
3. 对 CP2K BTO/GaN 先用 `STATIC_ONLY=1` 生成独立端点摘要，只有通过门禁后才允许新的生产路径，不覆盖已取消目录。
4. 对 ABINIT GaN 先修正端点物理设置并通过同一门禁，再决定是否重跑；当前 `max_steps_reached` 结果不得进入生产矩阵。
5. 生成统一后端状态表、路径图和文献比较数据；未达到阈值或非同一物理模型的结果不得进入“已验证生产矩阵”。
