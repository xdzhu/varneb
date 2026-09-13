# Codex 线程交接说明（2026-09-13）

## 项目与目标

- 项目目录：`D:\Work\Code\vasp_neb`
- 项目正式名称：`VARNEB`，上游仓库为 `https://github.com/xdzhu/varneb`；算法仍称 VC-NEB，Python 导入包暂保留为 `vcneb` 以维持兼容。
- 总目标：按 `VCNEB_PROJECT_PLAN.md` 推进至完整可用，包括理论推导、核心算法与模式/方向约束升级、VASP/ABACUS 计算器适配、集群典型晶体相变验证、与文献/参考实现对比、CPC 论文草稿及可复现发布材料。
- 长时间计算约束：BTO 主线统一使用合肥 `hfacnormal01`，每个 image worker 使用 32 MPI；每次提交前检查 `sinfo`/`squeue`，避免抢占用户任务。旧 `cu*` 记录仅作历史基线。
- 发布状态：版本 `0.0.1` 已通过 Trusted Publishing 发布；提交 `35c0dbb` 与标签 `v0.0.1` 已推送，重跑后的 workflow `34763435900` 为 success，PyPI 可下载 `varneb==0.0.1`。`.github/workflows/publish-pypi.yml` 仅由版本标签、已发布 Release 或手动确认触发。

## 原线程中断原因

原线程 ID：`019edb54-7732-77f3-8820-42123f31b2a3`。

最后一次请求在发送阶段失败：`input[1102].arguments` 长度为 `2,492,367`，超过接口允许的 `1,048,576`。这是请求字段长度校验失败，不代表项目代码或集群作业本身失败。新线程不要粘贴完整历史；以本文件、当前工作树、作业状态和日志为准。

## 已完成/已确认

- 端点弛豫摘要判据已修复：`converged` 使用 `FrechetCellFilter` 的完整广义力，同时保留原子力与应力字段。
- 端点路径预检已通过：最短距离约 `1.994 Å`（最终路径预检约 `1.983 Å`），最大变胞形变约 `0.1238`；HfO₂ 模板门槛已调整为可配置，默认约 `0.25`。
- 通用端点脚本不再暗含固定 `mpirun -np 40`，launcher 由环境/调用方提供。
- HfO₂ 7-image 普通 VCNEB 预收敛 40 步已完成：内部像 3 有局部峰，初步正向能垒约 `1.052235 eV`，反应焓约 `-0.743338 eV`，`has_interior_barrier=true`；但残差约 `0.847 eV/A`，不能作为最终能垒。
- 从完整第 40 步链快照成功续算，确认是状态恢复而不是重新插值。

## 当前计算状态

- HfO₂ 作业 `27674272` 已完成并降级为 60-Ry 诊断，不作为生产能垒。
- BTO 端点 `27675330`/`27675331` 及其 5/7/9-image 路径（`27676251`/`27676310`/`27676513`）均为 `100 Ry + DZP-10au` 的 ABACUS 静态 SCF，由 ASE cell-filter 负责外层更新；完整结果已归档并继续作为 BTO 生产证据。原生 `cell-relax` 方法迁移不再回头重算 BTO。
- BTO 已有完整的 5/7/9-image、反向路径和精度对照（均为 `100 Ry + DZP-10au`）；这些结果继续作为 BTO 生产证据。此前误发起的 ABACUS-native 端点方法重算 `27676929` (C, failed) / `27676930` (T, cancelled) 不纳入证据，也不再触发 BTO 路径重跑。后续原生 `cell-relax` 迁移只在 HfO₂ 主线执行，除非另有明确需求。
- BTO C→T 反向 5-image 已收敛（`27676299` + `27676485`）：`0.0195577 eV/A`、反应焓 `-0.0872206 eV`、无内部势垒。
- CI 闸门已审计并 withheld：当前没有内部 saddle，不启动物理上无意义的 CI；记录见 `outputs/batio3_t_to_c_pbe100_dzp10au/bto_ci_gate.json`。
- BTO 电子精度对照已完成：6×6×6/SCF 1e-9 的端点与静态路径反应能约 `0.0724 eV`，较 4×4×4 生产设置低约 15 meV；详见 `bto_precision_sensitivity.json`。
- BTO 端点求解器等价性已单独落档：`outputs/batio3_t_to_c_pbe100_dzp10au/endpoint_optimizer_equivalence.md`；FIRE 收敛端点不因迁移到 BFGS/native `cell-relax` 而失效，BTO 不重复计算。
- HfO₂ 端点兼容性已核验：当前 T/PO 输入均为 `12` 原子、`Hf4O8`，元素顺序一致；T 是为匹配 PO 常规胞而构造的 `√2×√2×1` 12 原子四方共格超胞，PO 是 12 原子极性正交常规胞。证据见 `outputs/hfo2_t_to_po_pbe100_dzp10au/endpoint_structure_audit.json` 和 `validation/hfo2_t_to_po/mapping_report.json`。
- 发布前干净安装烟测已通过：在临时目标目录构建/安装 `vcneb` wheel，并从安装目录运行 toy VCNEB（`barrier_eV=0.250004`）及 `vcneb --help`。
- PyPI 项目简介已同步到 `pyproject.toml`：`VARiable-cell Nudged Elastic Band code with universal first-frinciples calculators`；README 顶部保留规范拼写，Homepage/Repository/Issues 元数据链接已更新。新增 `tests/check_release_metadata.py` 并接入 `.github/workflows/publish-pypi.yml`，发布前会阻止简介或触发策略漂移。此次仅推送代码，未打 tag，因此不会触发 PyPI 发布；已发布的 `v0.0.1` 不覆盖，下一次按需发布的新版本（需新版本号）会携带该简介。
- HfO₂ 100 Ry/Orb-DZP-10au 端点重弛豫曾改用 ABACUS 原生 `cell-relax`/`relax_method bfgs`：ASE 外部 BFGS `27676827`/`27676828` 已保留首步诊断后取消；原生初次作业 `27676868`/`27676869` 已分别完成约 19/23 个离子步，但因驱动解析 `STRU_ION_D` 字段 bug 失败，断点已保留。修复解析器（含标签、磁矩、晶格常数单位换算）后从断点续算 `27676976` (T) 与 `27676977` (PO)，均 32 MPI；两端均因 ABACUS BFGS trust-radius breakdown 保留 `STRU_ION_D`，已切换到同属 ABACUS 的 `bfgs_trad` 续算 `27677116` (PO) 与 `27677117` (T)。两者均在观察多个完整 cell cycle 后因数值不稳定安全取消，完整输出已归档：PO 末态 `fmax=7.55946 eV/A`、`max_stress=2856.90 kbar`、`pressure=1644.58 kbar`；T 末态 `fmax=0.132648 eV/A`、`max_stress=27.40 kbar`、`pressure=-26.68 kbar`。两端仍保持 12 原子 `Hf4O8`，这只是优化器/步长稳定性问题，不代表物理端点不同。下一步可在相同能量面和自由度下采用更保守的外部 ASE BFGS 或阻尼 native 重启；两端通过受力+应力门槛后才启动高精度普通 VCNEB。
- 当前有效的保守外部 BFGS 端点作业为 `27677331` (PO) 与 `27677332` (T)，均 32 MPI、100 Ry、Orb-DZP-10au、`maxstep=0.0005`；前两轮仅因远端脚本同步/轨道文件名错误在预检阶段退出，不计入 DFT 证据。外部与 native 的验收标准统一为 `fmax<0.02 eV/A` 且 `max|stress|<0.1 kbar`。
- 上述 `maxstep=0.0005` 作业已在多步稳定窗口后取消并完整归档（T 第 19 步 `0.020087 eV/A`，PO 第 15 步 `0.483720 eV/A`），没有丢失恢复状态。当前续算为 `27677466` (T) 与 `27677467` (PO)，从最后轨迹帧恢复，改用 `maxstep=0.005`、`fmax=0.01`、120 步预算；观察 3--5 步稳定性后再继续。

## 最新轮询（2026-09-13）

- T 续算 `27677491` 已完成并通过双门槛：最大广义力 `0.00049882 eV/A`、最大应力 `0.01841 kbar`，12 原子 `Hf4O8`。
- PO 续算 `27677467` 仍在 `node148` 上以 32 MPI 运行；最新完整步为 9，`fmax=0.294759 eV/A`，连续下降。它是较低对称性的正交 PO 常规胞，不是已经收敛的 T 端。
- T 的完整输出已归档于 `outputs/hfo2_t_to_po_pbe100_dzp10au/endpoint_relax_T_ase_bfgs_maxstep005_job27677491/`。

## 端点晋级与普通 VCNEB（最新）

- PO 续算 `27677874`/`27677918` 已完成，最终 `fmax=0.0004973 eV/A`、最大应力 `0.07076 kbar`；T/PO 均为 12 原子 `Hf4O8` conventional cell，已由 `endpoint_promotion_gate.json` 原子晋级。
- 7-image calculator-free preflight 已通过：最小距离 `2.02496 A`、最大 deformation `0.04841`，无折返；所有 image 的 ABACUS 计算器能力检查通过。
- 原普通 HfO₂ 作业 `27677945` 因并行解析竞态失败；修复后的旧布局续算作业 `27678004` 已安全 CANCELLED（保留其完整 step 15 轨迹作为分布式迁移的恢复源）：128-task controller allocation，4 个并行 image worker，每 worker 32 MPI；100 Ry、Orb-DZP-10au、2x2x2、FIRE、无 CI。该轨迹只属于预收敛过程。

## 并行恢复修复（最新）

- `27677945` 后续因 ASE-ABACUS 的进程级 `ase_sort.dat` 竞态在 image 3 解析失败；已保留 step 2 链快照、manifest 和错误日志，未把该作业当作物理失败。
- `vcneb/abacus.py` 已加入 identity species-order 的安全补丁，并新增 `tests/check_abacus_parallel_sort.py`；远程 ICU 测试通过。旧全局 `ase_sort.dat` 已移到 `.stale_job27677945`。
- 修复后的旧布局续算作业 `27678004` 曾从 step 2 恢复，128 task / 4×32 MPI / no-CI；在迁移到按 image 分布式布局后已安全取消，作为恢复源的 step 15 完整轨迹和日志均保留。

## 分布式 image 语义与资源优化（最新）

- `n_images` 的定义是总帧数，始终包含两个固定端点；因此“7-image”=2 个端点+5 个 interior image，不是 7 个中间帧。
- VCNEB manager 现在首次读取端点的 energy/forces/stress 后缓存，后续迭代只通过 image executor 评估 `image_indices=[1,2,...,n_images-2]`；worker manifest 会记录 `image_count=n_images-2`，用于审计端点未重复计算。
- 分布式模板 `cluster/hf_hfo2_vcneb_distributed.slurm` 已改为默认 `IMAGE_WORKERS=N_IMAGES-2`、每 worker 32 MPI；7-image 默认 5 workers/160 tasks，仍需 2 个节点（单节点 128 核不足），不再为两个端点启动 worker。
- 新验证作业 `27678218` 已在 `hfacnormal01` 的 `node[381-382]` 启动，`NumTasks=160`、`NumCPUs=192`，首个 manifest 已确认只包含 image 1--5；从 `27678004` 的 step 15 完整轨迹恢复。此前错误布局的 `27678137`/`27678176` 已取消，不纳入物理结果。
- `27678218` 已完成（Slurm `COMPLETED`、`01:08:38`、exit `0`）：7 总帧（5 个 interior image）普通 VCNEB 在 FIRE `fmax_target=0.05 eV/A` 下达到 `final_max_generalized_force=0.0455156 eV/A`。最终正向焓垒 `0.1567509 eV`、反应焓 `-0.3252848 eV`，内部峰为 image 2（`has_interior_barrier=true`，未启用 CI）；manifest 共 96 条记录且每条只含 image 1--5，5 个 worker、每 worker 32 MPI，7 个 ABACUS image 目录无错误退出。
- 远端 `scripts/audit_vcneb_result.py --max-min-distance 2.0` 审计返回 `status=ok`、`issues=[]`：最小路径距离 `2.026335 A`、最大形变 `0.04894`、几何有效；summary 记录 `endpoint_evaluation_policy=fixed_cached_once`，确认端点未在 VCNEB 迭代中重复派发。结果文件为 `vcneb_summary.json/.txt`、`vcneb.traj`、`vcneb_barrier.png`、snapshots 和 worker manifest。
- 为补齐 HfO₂ 的 image-count 对照，已提交普通无 CI 作业 `27678406`（5 总帧=3 interior，1 节点×96 ranks）与 `27678407`（9 总帧=7 interior，2 节点×224 ranks）；两者均使用 100 Ry/10 au、32-MPI worker、固定端点一次缓存，独立工作目录分别为 `vcneb_n5_fire_distributed_cmp` 与 `vcneb_n9_fire_distributed_cmp`。
- `27678406` 已连续运行到 step 50：step 36 的最佳 `fmax=0.098970 eV/A` 后连续 14 个完整 optimizer step 回弹至 `0.188661 eV/A`。这已超过约定的 3--5 步观察窗口，因此在保留原目录、51 个完整快照和 manifest 后于 2026-09-14 02:35 左右安全取消。受控续算 `27678689` 已从旧 `vcneb.traj` 的完整 step 36 在新目录 `vcneb_n5_fire_distributed_resume_best36` 启动，仍为 5 总帧/3 interior、1 节点×96 ranks、3×32 MPI、普通无 CI，只将 `MAXSTEP` 从 0.02 降为 0.005；预检已通过且未覆盖旧证据。该续算完成 step 10（`fmax=0.112725 eV/A`）后也呈连续缓慢上升，于 2026-09-14 03:01 左右安全取消；新目录保留 11 个完整快照，step 0/旧 step36 的 `0.098970 eV/A` 是该 5-image 分支最佳证据。结论是 5 总帧在当前 100 Ry/10 au、`k=0.2`、FIRE 路径下未收敛，不把它当作最终物理能垒。
- 截至 2026-09-14 03:50 左右，`27678407` 已完成 step 55（`fmax=0.065946 eV/A`，从 step 0 持续下降），仍为 RUNNING；manifest 只记录 `[1,2,3,4,5,6,7]`，当前未见 ABACUS 错误 task。
- 由于 7-image 普通路径已达到真实内部峰且 `fmax<0.05 eV/A`，CI 精修作业 `27678507` 已提交：7 总帧/5 interior、2 节点×160 tasks、每 worker 32 MPI、`MAXSTEP=0.01`、`FMAX=0.03`、`CLIMB_AFTER=0`，从 `27678218` 的完整 `vcneb.traj` 恢复，工作目录为 `vcneb_n7_ci_refine`。该作业是首次真实 HfO₂ CI，不用于 BTO。
- CI `27678507` 已于 2026-09-14 04:14 左右正常完成（Slurm `COMPLETED`、`03:12:03`、exit `0`）。它在经历 step 15--80 的回弹/平台窗口后继续推进，于 step 91 达到 `final_max_generalized_force=0.0295475 eV/A`（目标 `0.03`），没有因首次回弹提前停止。CI 结果的正向焓垒为 `0.1291722 eV`，反应焓 `-0.3252848 eV`，最高/攀爬 image 为 2，`has_interior_barrier=true`；`endpoint_evaluation_policy=fixed_cached_once`、`climb_after=0`。独立审计返回 `status=ok`、`issues=[]`，最小路径距离 `2.035829 A`、最大形变 `0.0489263`、最大应力 `1.76792 kbar`。image worker manifest 共 279 条、状态全为 `ok`，只含 interior `[1,2,3,4,5]`，5 个 worker、每 worker 32 MPI；结果目录为 `vcneb_n7_ci_refine`。
- 9-image 普通无 CI 作业 `27678407` 已于 2026-09-14 04:30 左右正常完成（Slurm `COMPLETED`、`03:47:07`、exit `0`）：9 总帧/7 interior，最终 `fmax=0.0497191 eV/A`，正向焓垒 `0.1562092 eV`，反应焓 `-0.3252848 eV`，最高内部峰为 image 3。独立审计 `status=ok`、`issues=[]`，最小路径距离 `2.027974 A`、最大形变 `0.0496014`、最大应力 `5.38872 kbar`；manifest 207 条记录全为 `ok`，只含 image 1--7，7 个 worker、每 worker 32 MPI。
- HfO₂ image-count 对照已归档于 `outputs/hfo2_t_to_po_pbe100_dzp10au/image_count_comparison.md`；收敛的普通 7/9-image 焓垒只差 `0.0005417 eV`（约 0.35%），但 5-image 分支未收敛，因此不纳入物理能垒比较。已从远端复制 7-image CI 与 9-image 普通路径的 summary、audit、preflight 和 worker manifest；完整 ABACUS scratch 仍保留在 hfacnormal01。
- `scripts/compare_vcneb_images.py` 已用半个局部反应坐标段作为离散峰位默认容差，7/9-image 普通对照报告 `status=ok`、能垒 spread `0.0005417 eV`；机器可读结果为 `outputs/hfo2_t_to_po_pbe100_dzp10au/ordinary_image_comparison.json`。
- 新增 `outputs/hfo2_t_to_po_pbe100_dzp10au/hfo2_validation_provenance.json`，集中记录 12 原子端点、100 Ry/10 au、32-MPI worker 资源、5/7/9-image 与 CI 作业 ID、结果和本地/远端 artifact 位置。
- 代码工程化增量已在本地回归全通过并推送为 `31318f5`（功能主体 `0735990`）：新增可选 exact-state image cache、calculator namespace 锁定、原子 cache 写入、损坏条目回退以及 manifest 命中/未命中 provenance；远端源文件已同步，当前运行作业不受影响。
- 后续工程化修复已提交：`ThreadedCalculatorExecutor` 在并发 batch 的某个 image 失败时仍收集并原子保存已完成 sibling image，随后再返回 batch 错误；新增恢复回归确认成功 image 不丢失、重启不重复计算，manifest 保留失败与后续 cache hit provenance。
- `run_vcneb()` 现在可通过 `failure_report=` 原子写出 optimizer/calculator 中断报告（异常类型、已完成步数、trajectory/snapshot 位置和恢复提示），同时保持原异常继续抛出；对应回归已通过。line-search/ABACUS SCF 专用分类仍是未完成项。
- `examples/run_vcneb_abacus.py` 与 `examples/run_vcneb_vasp.py` 已默认把该 failure report 写入各自 workdir 的 `vcneb_failure.json`，因此集群模板发生中断时可直接定位恢复入口。
- 已按当前 `docs/theory.md`、README、代码和回归证据同步 `VCNEB_PROJECT_PLAN.md`：广义坐标/单位/切线/收敛判据、模式投影顺序、约束释放、VASP 静态 image 语义、calculator 切换和 dry-run/validate-only 现标为完成；ABACUS 专用 SCF 分类、USPEX/旧实现对比、真实材料模式对照、论文图表等仍保留为未完成。

## 建议的新线程第一步

1. 新计算先检查 `sinfo`/`squeue`，并使用 `hfacnormal01` 的 32-MPI image worker 资源模型。
2. 普通路径若发生回弹，至少观察 3--5 个完整 step，保存窗口内最佳快照，再决定是否续算或审计。
3. CI 只在普通路径达到门槛且存在真实内部能量峰时启动；当前 BTO 结论是 withheld。
4. 每个长作业都归档完整 trajectory、snapshots、summary、preflight、worker manifest、节点/task 数和审计结果。
5. 继续按 `VCNEB_PROJECT_PLAN.md` 推进，并把每个结论绑定到当前文件、日志、测试或作业输出。

## 交接原则

- 以当前工作树和实时作业状态为事实来源；本说明只提供压缩上下文。
- 新线程提示词保持短小，引用本文件，不要复制整段历史日志。
- 任何新结论都要区分“诊断/预收敛”与“最终收敛结果”。
