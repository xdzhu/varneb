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

单结构探针 `27744287` 已成功；两个相同的旧 image-4 结构均在 101 次 SCF 内收敛，退出码为 0。新 BTO 生产 `27744298` 正在运行，独立目录为 `cases/bto/cp2k_vcneb_production_stable`；截至最近记录尚处于高力早期迭代，不能提前判断收敛。

GaN CP2K `27741431` 仍使用原 profile 运行；若它完成或失败，将先审计 SCF，再决定是否用稳定 profile 独立重跑，不覆盖运行目录。

## 其他当前结果

| 案例/后端 | 作业 | 当前判定 |
|---|---:|---|
| GaN / QE | 27741412 | 已完成；最终广义力 `0.091018 eV/A`，按默认阈值可接受；能垒 `1.638143 eV`。旧 summary 未写显式 `converged`，审计以最终力为准。 |
| GaN / ABINIT | 27741516 | `max_steps_reached`，最终力 `9.537243 eV/A`，未收敛；暂不作为生产证据。 |
| GaN / LAMMPS | 27741574 | 最终力 `0.095441 eV/A`；仅是经典 Tersoff 对照，不能与 DFT/PBE 能垒等价比较。 |
| BTO / QE | 27741411 | 最终力 `0.052831 eV/A`，无内部势垒；按阈值可接受。 |
| BTO / VASP | 27741623 | 最终力 `0.095145 eV/A`，能垒 `0.042644 eV`；已收敛。 |

## 下一步

1. 等待并审计 `27744907` 的 ABACUS 生产结果；若再次失败，只看最小解析器的具体契约错误，不再改 FIRE 参数。
2. 等待 `27744298`，检查每个 image 的 CP2K SCF 次数、完整 force/stress 和最终广义力；若仍失败，保留稳定 profile 的探针证据，继续区分单 image SCF 与路径几何放大效应。
3. 对 `27741431` 只做状态和日志核验，不中断、不覆盖；终止后再决定 GaN CP2K profile 是否迁移。
4. 生成统一后端状态表、路径图和文献比较数据；未达到阈值或非同一物理模型的结果不得进入“已验证生产矩阵”。
