# VCNEB 项目标准化任务清单

> 项目目标：构建一个纯 Python、计算器无关、可复现的变胞 NEB（VC-NEB）工具包，支持 VASP、ABACUS 以及后续可插拔的 ASE Calculator；同时支持普通 VC-NEB、爬山图像（CI-VCNEB）、模式引导路径和严格的模式/方向约束，用于研究晶体相变过渡态与变胞能垒。
>
> 计划版本：v0.1，日期：2026-09-12。当前工作树是研究开发版，不代表已经达到论文发布或生产计算标准。
>
> 代码上游：<https://github.com/xdzhu/varneb>。正式项目名称为 `VARNEB`；算法名称仍写作 VC-NEB，Python 导入包暂保留为 `vcneb` 以维持兼容。

## 1. 执行原则

- 先定义广义坐标、内积、应力符号和单位，再扩展算法；所有公式都要有数值回归测试。
- 先用解析势和廉价计算器验证算法，再使用 DFT。不能把单次 DFT 跑通当作算法正确。
- 将“路径初始猜测的模式引导”和“优化过程中只允许沿模式运动”作为两个不同功能实现、测试和描述。
- 保持计算器适配层薄：VCNEB 核心只依赖统一的结构、能量、原子力和应力接口。
- 每个长任务都必须能恢复、能记录、能复核；每个科学结论都必须保留输入、代码版本、参数、节点和输出摘要。
- 集群规则：后续长时间计算统一从 `ssh hf` 进入合肥超算并通过 Slurm 提交；`cu17`、`cu22`–`cu26` 仅保留为历史试算记录，不再作为后续生产任务入口。每次提交前检查 `sinfo`/`squeue`，按体系规模、计算器和实测并行效率选择 `ntasks`，不机械固定 40 核，也不申请小体系的整节点独占。所有计算节点共享程序环境，代码与结果目录需保留版本、作业号、节点、实际 task 数和输入摘要。

## 2. 当前基线状态

### 已完成

- [x] 已找到并分析旧实现目录：`/home/zhuxd/abacus/8.dielec/vcneb`。
- [x] 已检查 USPEX 10.6 压缩包；本地包主要是安装脚本、许可证/说明和 MATLAB Runtime 自解压二进制，未获得可审计的 VCNEB 源码，因此不把它当作可直接复用的代码依赖。
- [x] 已建立纯 Python VCNEB 核心原型，使用原子坐标和 cell 的扩展空间进行路径优化。
- [x] 已加入 cell 力/应力通道、cell mask、端点固定、轨迹恢复、快照追加和优化器接口。
- [x] 已加入 `Mode` 数据结构及文本、JSON、NPZ 读取；支持原子模式、质量加权、去平移和可选 cell 变形分量。
- [x] 已加入模式引导初始路径和模态投影诊断；端点保持严格不变。
- [x] 已加入 `atom_mask`，可冻结部分原子自由度。
- [x] 已有 VASP、ABACUS 和 toy/model 案例入口。
- [x] 已通过当前回归检查：fractional force、cell force、cell mask、atom mask、mode-guided path、trajectory resume、snapshot append、optimizer API、cache invalidation、parallel endpoint ownership、ABACUS command profile 和 validator flags。
- [x] 已完成解析 CI 鞍点诊断和固定 cell ASE CINEB 对照；前者恢复 0.25 eV 联合势鞍点，后者逐 image 能量一致且能垒差为 9.33e-6 eV。
- [x] 已完成 VASP 与 ABACUS 单 image 真实 calculator smoke；两者均返回有限 energy/forces/stress，并保留独立运行目录与输入验证记录。
- [x] 上述回归曾在 `235` 和 `cu05` 环境验证通过；后续正式验证统一迁移到本项目指定的 `cu17`、`cu22`--`cu26` 节点。

### 尚未宣称完成

- [x] 已完成固定随机扰动路径和不同 cell 参数化的对照：4 个种子分别测试 `linear/log_strain`，两阶段 CI 的 8/8 案例通过；六个独立应变有限差分误差表、解析基准与固定 cell ASE 对照也已完成。
- [~] ABACUS 已完成 HfO₂ 100 Ry/10 au 生产级 VCNEB（普通 7/9-image 与 CI）；VASP 真实端到端生产级 VCNEB 仍待完成。
- [x] HfO₂ T 相到 PO 相的 12 原子结构来源、原子一一映射、端点独立弛豫和生产精度路径/能垒已保存并通过审计。
- [x] 已实现基础严格模式子空间和 projected-update 约束，并加入端点子空间验证、方向冲突诊断和解析势 release-and-refine；真实材料对照仍未完成。
- [~] 已形成初版误差预算与资源效率记录（`outputs/hfo2_t_to_po_pbe100_dzp10au/error_budget_and_efficiency.md`）及论文结果表；完整投稿数据包和正式图表仍待补。

## 3. 阶段总览与里程碑

| 阶段 | 名称 | 目标 | 主要出口标准 |
|---|---|---|---|
| P0 | 基线冻结 | 固化当前原型和回归基线 | 可复现安装、测试全绿、版本快照 |
| P1 | 理论规范 | 固化 VCNEB 和模式约束的数学定义 | 公式、符号、单位、边界条件完整 |
| P2 | 核心算法升级 | 提高插值、切线、cell 处理和约束的可靠性 | 解析势/有限差分/对照测试通过 |
| P3 | 计算器与工程化 | VASP、ABACUS、ASE 接口统一且可恢复 | 两种 DFT calculator 端到端 smoke test |
| P4 | 数值验证 | 证明结果对图像数、弹簧、cell 尺度和优化器稳定 | 系统收敛矩阵和误差报告 |
| P5 | 材料案例 | 完成 HfO2 和至少一个钙钛矿相变案例 | 可复核的能垒、路径、应力和结构数据 |
| P6 | 模式/方向扩展 | 完成模式子空间、原子方向限制和组合约束 | 约束语义、投影力、释放约束测试清楚 |
| P7 | 论文与发布 | 形成论文、文档、示例和可发布代码 | 新机器可按文档复现核心结果 |

依赖主线为：`P0 -> P1 -> P2 -> P3 -> P4 -> P5 -> P7`；`P6` 从 P2 后即可并行，但必须在 P7 前完成。P3 的接口开发可与 P2 的解析测试并行，不能用 DFT 结果替代 P2 的数值证明。

## 4. P0：基线冻结

### 任务

- [x] 为当前工作树建立一个明确的开发版本号和变更摘要；记录 `core.py`、`modes.py`、测试脚本和示例入口。
- [x] 已在空闲的 `cu17` 上重新运行现有回归集，并将 Python、ASE、NumPy、SciPy、VASP/ABACUS 环境信息写入 `outputs/vcneb_p0_baseline_manifest.json`。
- [x] 历史 cu 试算均在启动前检查负载、进程和用户任务；后续 HF Slurm 任务按 `sinfo`/`squeue` 和实际 scaling 选择资源，忙节点不挤占。
- [x] 利用共享目录完成一次代码同步即可；节点无需安装 Git，运行目录和 manifest 均记录对应的 Git commit。
- [x] README 已将 toy、model、VASP/ABACUS DFT-smoke 与 production-template 示例显式标注，避免把模板误读为已收敛结果。
- [x] 检查当前 `/home/zhuxd/abacus/8.dielec/vcneb` 旧代码与工作树的差异，保留可借鉴算法说明和输入格式，但不复制无法验证的逻辑。
- [x] 规定结果目录命名：材料、端点、图像数、calculator、参数集、代码版本必须可从目录名或 manifest 读出。

### 产出与验收

- 产出：`baseline_manifest.json`、回归日志、环境清单、基线版本说明。
- 验收：新节点上从干净工作目录执行一条命令即可完成所有非 DFT 回归；失败时日志能定位到测试名和输入。

## 5. P1：理论与接口规范

### 5.1 广义坐标与几何

- [x] `docs/theory.md` 明确周期结构对象：原子分数坐标 `s`、cell 矩阵 `h`、元素/质量、周期边界和原子映射。
- [x] `docs/theory.md` 明确路径变量采用 `(s,F)` 并给出与 `(r,h)` 的转换、最小镜像分支和变胞时的笛卡尔位置定义。
- [x] `docs/theory.md` 固化扩展空间内积、原子与 cell 的联合欧氏度量，并记录 `cell_scale` 的长度单位、默认值和敏感性要求。
- [x] `docs/theory.md` 规定 cell 插值的正体积/非奇异域、正定 log-strain 条件和整体旋转处理；初始化阶段有几何硬审计。
- [x] 已评估线性 cell 插值与对数/指数晶格插值；`interpolate_vcneb()` 提供 `linear`、正定 deformation 的 `log_strain` 和自定义回调，并保留线性方式作为兼容默认选项。

### 5.2 力、应力与 NEB 方程

- [x] `docs/theory.md` 统一能量、原子力、应力和 cell 广义力的单位/符号约定，并由有限差分回归覆盖。
- [x] 从能量对 `(s, h)` 的导数推导原子广义力和 cell 广义力，明确应力到 cell force 的转换、体积因子以及张量转置约定，并用非对角 deformation/外压有限差分验证。
- [x] `docs/theory.md` 定义能量加权切线、真实力垂直投影、弹簧平行分量、固定端点和 climbing-image 规则。
- [x] `docs/theory.md` 明确最大广义力、原子力、cell force/stress、能量和路径几何等收敛量；`path_diagnostics()` 输出对应字段。
- [x] `docs/theory.md` 区分固定 cell、固定原子/方向、cell 分量 mask、软模式和硬子空间约束，并有对应 API/测试。

### 5.3 模式与方向约束语义

- [x] 模式引导、`projected` 更新和 `subspace` 硬约束已在 API 与 `docs/theory.md` 中分层命名和说明。
- [x] 已定义模式在 Cartesian/fractional/扩展空间中的转换、归一化、多模式正交化与质量加权输入处理；严格动力学质量度量的进一步验证仍单列。
- [x] 已支持原子方向 mask/投影矩阵，并在文档中区分方向约束与冻结原子。
- [x] `outputs/uspex_vcneb_mode_analysis.md` 已记录 ABINIT `iatfix/iatfixx/iatfixy/iatfixz` 与线性组合约束语义，并明确不把 ABINIT 输入格式暴露给核心 API。
- [x] 已固定“先投影广义力，再做 NEB 切向/弹簧分解”的顺序，并由模式投影与有限差分回归验证。

### 出口标准

- [x] `docs/theory.md` 已完成广义坐标、力/应力、切线、CI、约束、收敛判据、实现伪代码和公共参数表。
- [~] 核心公共参数的默认值、合法域和单位已在理论说明中定义；calculator-specific 经验参数与磁性/电子占据敏感性仍需各自案例补充。

## 6. P2：核心算法升级

### 6.1 插值与结构路径

- [x] 实现端点原子自动映射检查：按元素分组、用周期最小镜像几何代价求 permutation，并允许用户显式提供 mapping；`validate_atom_mapping()` 返回可审计位移报告。
- [x] 实现分数坐标最小镜像插值，避免跨周期边界产生长路径，并加入公共周期平移与 mapping 联合预处理。
- [x] 实现 cell 插值策略选择器：线性、正定 deformation 的对数应变插值和用户自定义插值，并对不适用的 log 路径明确报错。
- [x] 对每个中间 image 检查 cell determinant、原子重叠、最小距离和异常应变；异常时在初始化阶段明确失败。
- [x] 已加入 calculator-free 的路径几何审计，报告体积、周期 MIC 最短距离、deformation 和扩展坐标折返 cosine；可通过插值参数在初始化阶段硬拒绝碰撞/异常 deformation，也可通过 `fold_cosine_threshold` 拒绝折返。
- [x] 已加入通用 `linear`/`log_strain` 初始路径比较脚本，固定 image 数和端点映射后输出体积、最短距离和 deformation 分布，并可复用于 HfO2/BaTiO3；BaTiO3 预检已发现并修正 identity 映射导致的 O 原子近重合。
- [x] 增加基于模式的初始路径组合：结构插值 + 模式位移 + cell 模式，并保持端点严格一致；已补充 cell-mode 回归。

### 6.2 NEB/VCNEB 力与优化

- [x] 复核切线在能量单调、能量拐点和平坦区的行为；已覆盖单调段、非对称峰值/谷值和 improved tangent 的加权方向。
- [x] 复核弹簧力在扩展空间的定义，保证原子与 cell 部分使用同一反应坐标度量；已用非均匀 image 间距逐项回归。
- [x] 实现 CI-VCNEB，并在联合原子-cell 解析势上验证 climbing image 不引入端点漂移或错误的 cell 方向；真实 DFT CI 仍待收敛案例。
- [x] 优先复用成熟优化器 API；核心只提供广义坐标、梯度和约束投影，不重复实现通用 LBFGS/FIRE 数值细节。
- [x] 将“普通 NEB 松弛后再 CI”固化为 `run_vcneb(climb_after=N)`，并加入折返线段 cosine 诊断及硬拒绝选项。
- [x] 增加局部峰/内部势垒诊断：单调路径或内部峰低于端点时明确给出
  `has_interior_barrier=false` 和 `ci_warning`，避免把无势垒路径误报为过渡态。
- [~] calculator 异常、NaN/Inf、非法 cell、SCF/超时/MPI 关键词现在会保留带 `failure_category` 和恢复提示的 failure report；并发 image 失败批次会保留已完成 sibling 结果。line-search 自动重试和基于 ABACUS 原始 stdout 的专用分类仍待补齐。
- [x] 设计 image 级缓存和原子写入；中断后可以从最近快照继续，且不会混用不同参数集的结果；并发批次发生单 image 失败时，已完成的 sibling image 结果仍会先写入缓存。

### 6.3 严格模式子空间

- [x] 实现单模式和多模式的正交投影器；支持仅原子模式、原子+cell 联合模式。
- [x] 实现约束优化：每一步更新后回投影到允许子空间；约化变量优化仍可作为后续性能优化。
- [x] 实现方向约束投影矩阵，与 atom mask/cell mask 组合时通过 `direction_basis_conflicts()` 给出部分裁剪、完全失活和秩损失诊断。
- [x] 对比三种结果：无约束 VCNEB、模式引导但最终无约束、严格模式约束；已在耦合解析势上以 `2.6e-5 eV` 以内的能垒差完成对照，文档明确三者不能互相替代。
- [x] `examples/run_release_and_refine.py` 已提供先模式子空间搜索、再解除约束做全空间精修的可复现实例；真实材料对照仍待完成。

### P2 出口标准

- [x] 解析势上能从固定随机扰动初始路径收敛到已知 saddle/MEP；直接 CI 的错误折返反例已发现，并由“普通 NEB 松弛后再 CI”的 protocol 修正。
- [x] 对同一问题改变 `cell_scale` 后，经过一致的广义坐标定义，能垒和结构结果在合理范围内稳定；已覆盖 `4/5/6 A`。
- [x] 固定 cell 时与 ASE 标准 NEB 的能量、力和路径趋势一致，差异有可解释的坐标/优化器原因；七 image 对照的能垒差为 `9.33e-6 eV`。

## 7. P3：计算器无关层与工程化

### 7.1 Calculator contract

- [x] 固化最小 calculator 协议：输入结构、能量、原子力、应力/virial、单位和计算状态；新增 `vcneb.calculator` 预检与机器可读报告。
- [x] 实现 calculator capability 检查：是否支持 stress、是否支持 variable cell、是否能返回每 image 独立目录；支持重复目录检测。
- [x] 核心只依赖 ASE energy/forces/stress contract；VASP/ABACUS 命令、环境变量和 MPI launcher 均由 adapter/profile 与 cluster wrapper 注入，核心不拼接或执行外部命令。
- [x] 明确 stress 缺失时的行为：`run_vcneb()` 预检直接拒绝，不能静默把 cell force 当成零。
- [x] 运行时 calculator 异常会保留 image 编号、目录和命令上下文，并由回归测试覆盖。
- [x] 新增物理 HfO2 单 image ABACUS smoke 入口；已在集群完成真实 SCF、原子力和 stress 验证。

### 7.2 VASP 适配

- [x] 完成单 image 能量/力/stress 解析、目录隔离、重启文件策略和错误分类；真实 GaN VASP smoke 已记录到 `outputs/vasp_gan_single_image_manifest.json`。
- [x] 已明确 VASP 的 cell-relax 输入只用于端点预弛豫；VC-NEB image calculator 使用静态 `IBRION=-1`、`NSW=0`、`ISIF=2`、`ISYM=0`，cell 更新由 VCNEB 控制。
- [x] 完成极小晶胞的 VASP 单 image smoke，并保留 INCAR/KPOINTS/POTCAR 来源记录；完整材料案例仍待开展。
- [x] `run_vasp_single_image_smoke.py` 已在 `cu26` 用 40 核 VASP 6.3.2 实跑，energy/forces/stress 与静态输入验证均通过。

### 7.3 ABACUS 适配

- [x] 完成 STRU/KPT/INPUT 生成和结果解析；核对 stress 输出、单位和晶格方向，并提供参数化多 image 入口；HfO2 单 image 与 7-image ABACUS smoke 已通过。
- [x] 已加入 native endpoint 只读闸门审计器 `scripts/audit_native_endpoint.py`，对退出码、收敛标志、原子数/化学计量、最大原子力和应力缺失或超阈值直接判定为证据不足。
- [x] 已加入 `relax_abacus_endpoint.py`，使用 ASE `FrechetCellFilter` 完成 HfO2 两端点的低精度、宽松阈值独立原子/cell 弛豫 trial 并记录。
- [~] ABACUS command profile、32-MPI worker 宽度、Slurm 退出码和真实生产作业已验证；超时与 SCF 不收敛的专用分类/重试策略仍待补齐。
- [x] `run_vcneb()` 只依赖 ASE energy/forces/stress contract，VASP 与 ABACUS 示例通过替换 calculator factory 切换，核心路径参数不变。

### 7.4 可恢复运行

- [x] driver 已支持保存 manifest、每 image 结构、energy、force、stress、cell 指标、真实/弹簧/NEB 力分解和收敛指标；每轮广义切线仍作为后续增强项。
- [x] 已加入可选 image-level 并发执行器；主控制器通过独立 Slurm job steps 并发 image calculator，真实 BTO `27675909` 完成 4×32 MPI smoke 并生成正常 summary/trajectory。
- [x] 并发 executor 支持显式 per-image retry；fail-once 回归确认只重算失败 image。
- [x] 已加入 image-worker JSONL manifest（每个控制器 evaluation batch 记录状态、耗时、重试次数和失败信息），并加入 `--resume-step`/`RESUME_STEP` 精确恢复完整 chain snapshot；`ThreadedCalculatorExecutor` 现在支持可选的跨作业 exact-state image cache、namespace 锁定、原子写入、命中/未命中 provenance 和损坏条目回退。
- [x] 已提供 `validate_vcneb_inputs.py` dry-run、ABACUS `--validate-only`、VASP/ABACUS single-image smoke，以及 `--n-images` 路径入口；发布前元数据回归也已加入。
- [x] 生成机器可读的 summary JSON、人类可读的文本摘要，以及由 `scripts/export_vcneb_metrics.py` 导出的逐 image CSV；论文级图表仍作为发布增强项。

### 出口标准

- [x] VASP 与 ABACUS 各至少完成一个单 image smoke test；结果分别记录在 `outputs/vasp_gan_single_image_manifest.json` 和 `outputs/abacus_hfo2_smoke_manifest.json`。
- [x] 断开后恢复不会丢失已完成 image，也不会重复覆盖有效结果；新增失败批次恢复回归，确认成功 image 缓存可复用且 manifest 保留命中/失败 provenance。
- [x] calculator 失败时错误信息包含 image 编号、输入目录和命令上下文；建议动作由上层 launcher 根据 calculator 类型补充。

## 8. P4：数值正确性与收敛验证

### 8.1 单元和微分验证

- [x] 对原子坐标、cell deformation 分量和联合焓扰动分别做中心有限差分；已补充三种对角应变和三种对称剪切应变的六方向扫描。
- [x] 验证能量梯度、应力符号、cell force、mask 后梯度和模式投影梯度。
- [x] 覆盖非正交 cell、周期边界跨越、原子排列变化和极小体积保护。
- [x] 增加容差随有限差分步长变化的表格，覆盖零压力与有限外压，避免只测一个“碰巧通过”的步长。

### 8.2 解析模型

- [x] 构造含原子坐标和 cell 变量的可控多井势，已知 minimum、saddle 和能垒；toy 与 HfO2 12 原子 synthetic model 均已端到端运行。
- [x] 构造含耦合项的 toy model，验证路径不能被错误地拆成独立原子和 cell 两条路径；`ToyPhaseTransition` 的双井方向同时含原子与 cell strain。
- [x] 在解析模型上对比不同弹簧常数、不同 image 数、FIRE/LBFGS 和 `cell_scale`；无弹簧与真实 DFT 组合仍待补充。
- [x] `examples/run_vcneb_convergence.py` 已在指定集群完成 54 组收敛矩阵并记录结果。
- [x] 收敛矩阵首次运行暴露并修正了测试端点漏设最终 cell strain 的问题；修正版 54/54 通过。
- [x] 解析双井模型验证 CI-VCNEB 将最高 image 推向已知 saddle，并通过 `saddle_diagnostics()` 给出负切线曲率和力残差；真实 DFT 鞍点诊断仍待收敛路径。

### 8.3 参考实现与理论对照

- [x] 固定 cell：与 ASE NEB/CINEB 在同一 calculator、同一端点和同一 image 数下比较；解析模型逐 image 能量一致，能垒差为 9.33e-6 eV。
- [x] 变胞：依据 Qian 等 VCNEB 论文的广义坐标、cell 力和弹簧思想完成首轮逐项对照；有限变形 stress measure 和大应变差异仍需补充。
- [x] `outputs/uspex_vcneb_mode_analysis.md` 已依据公开 USPEX VCNEB 手册对比输入语义/用户流程，并明确安装包无可审计源码、不声称复现内部实现。
- [~] 已依据历史迭代日志整理旧 `/home/zhuxd/abacus/8.dielec/vcneb` 的实现/失败模式对照（`outputs/legacy_vcneb_comparison.md`）；旧目录当前不可访问，定量同条件重跑仍待补。

### 8.4 收敛矩阵

- [x] image 数：已覆盖 5、7、9。
- [x] `cell_scale`：已覆盖 4、5、6 A 三组。
- [x] 弹簧常数：已覆盖 0.05、0.10、0.20 eV/A^2，记录路径收敛情况。
- [x] 电子结构精度：BTO 已完成 100 Ry/DZP-10au 的生产与 6³/1e-9 加严对照，并记录能垒/端点反应能差异。
- [x] 优化器和步长：已覆盖 FIRE/LBFGS，记录迭代数、收敛状态和残余力。

### P4 出口标准

- [x] 所有梯度/应力微分测试通过，并有零压与有限外压的误差表和 provenance manifest。
- [x] 解析模型的能垒、端点反应能、saddle 位置和负切线曲率在目标容差内恢复；完整 Hessian 诊断仍不在当前实现范围内。
- [x] BTO 已完成 5/7/9 image 收敛矩阵、反向路径和电子精度对照（100 Ry/DZP-10au）；不因端点优化器方法迁移重复计算。ABACUS-native cell-relax 作为 HfO₂ 主线的端点规范。

## 9. P5：材料案例与集群运行

### 9.1 案例选择与结构准备

- [x] HfO2 T 相到 PO 相：已固定为 12 原子共格端点对（T 为 √2×√2×1 四方超胞，PO 为 12 原子正交常规胞），元素顺序一致；100 Ry/10 au 高精度端点已分别收敛并通过端点晋级审计。
- [x] 备用/对照案例：BaTiO3 四方→立方相变已确定为当前主案例，端点原子映射与变胞路径已通过预检。
- [x] BTO 已保存结构来源、原子映射、端点能量、轨迹、summary、manifest 和审计；误发起的 native cell-relax 方法重算已终止，不覆盖既有生产证据。
- [x] 已加入并运行端点独立原子/cell 弛豫 trial；HfO2 两端点的 ABACUS-native cell-relax 生产结果已通过原子力、应力和化学计量审计。
- [x] 已对 12 原子 T/PO 端点完成原子映射、最短距离、cell 形变和路径段预检；高精度 7-image 路径复核显示几何有效、无折返，最短路径距离 `2.026335 A`。

### 9.2 分层运行策略

- [x] Level A 的 model/低精度 DFT smoke 分支已确认输入、路径、计算器、每 image 目录、输出解析、恢复和 CI 独立分支；HfO2 100 Ry/10 au 生产精度 7-image 路径已完成。
- [x] Level B：BTO 已有生产精度无约束 VCNEB 基线；CI 仍由内部 saddle 闸门决定，不因端点优化器切换重跑。
- [x] Level B 首轮已在 cu17 完整执行 7-image/40-step FIRE，并从完整轨迹尝试 LBFGS 恢复；两条分支均明确记录为未收敛，不能作为物理能垒。
- [x] 已补做输入预检失败诊断和 `FIRE(maxstep=0.05)` 保守步长重试；预检在 ABACUS 启动前正确拒绝错误 basis 路径，重试仍在前七步出现残余力单调增长，因此暂停继续消耗 DFT 资源，转入路径/应力诊断。
- [x] 已对最新完整重启轨迹做 7-image 零步静态力/应力诊断；image 2 的最大原子力/应力与 image 3 的最高焓错位，结果已进入 manifest，下一步优先改进机制路径。
- [x] Level C：BTO 已完成 5/7/9 image、反向路径及 6³/1e-9 对照；方法迁移记录为后续 HfO₂ 规范，不重复 BTO 生产矩阵。
- [x] Level A 的 model/低精度 DFT smoke、Level B 诊断 trial 与 HfO2 生产路径均已记录 Git 版本、节点、核数、输入、输出和结果摘要；HfO2 7-image CI 精修与 9-image image-count 扩展均已完成并独立审计。

### 9.3 集群执行规范

- [x] 后续入口已切换到 `ssh hf`；合肥超算使用 `hfacnormal01` 分区、共享目录 `/public/home/iai806` 和 Slurm，ABACUS/ASE 环境已完成健康检查。
- [x] 已加入 `cluster/hf_batio3_endpoint.slurm` 与 `cluster/hf_batio3_vcneb.slurm`：默认 16 task、无 `--exclusive`，`srun` 使用实际的 `${SLURM_NTASKS}`；可通过 `sbatch --ntasks=32` 等方式按实测放大。
- [ ] 每次提交前检查 `sinfo`/`squeue` 和账户任务；记录 partition、node、`ntasks`、耗时和资源效率，优先选择能满足 SCF/NEB 收敛的最小有效并行度。
- [ ] 不在网关节点上启动 DFT，不在本机 WSL 上启动长任务；合肥超算生产任务使用 Slurm，不使用 `qsub`/PBS。
- [ ] 每个案例独立目录运行，避免多个 image 共用临时文件或覆盖同名结果。
- [ ] 启动前做磁盘空间、环境、可执行文件、赝势/PAW、K 点和 MPI 健康检查。
- [ ] 运行中记录每轮完成 image 数和预计剩余时间；失败时先保留目录和 stdout/stderr，再决定是否重算。

### 9.4 科学结果判据

- [x] 给出端点相对能量、最高 image/saddle 能量、正向能垒、反应焓和 cell 演化；BTO 反向 5-image 路径也已独立运行并归档，HfO2 反向路径仍未纳入当前主线。
- [x] 代码和 ABACUS driver 已支持给出最高 image 的残余原子广义力、cell 广义力/应力、真实力与弹簧力分解；HfO2 生产路径 summary 已记录最终数值。
- [~] summary 已可通过 `scripts/export_vcneb_metrics.py` 导出逐 image 反应坐标、焓、晶格长度/角度、体积、应力和 NEB 力 CSV；`scripts/plot_vcneb_metrics.py` 已生成论文级焓垒/晶格/体积/力四联图，关键键长和模式投影仍待补齐。
- [x] 已用 7/9 总帧复算 HfO₂ 普通 VCNEB，并用可复现比较器报告能垒和离散峰位变化；两组结果均保留，不挑选更漂亮的一组。
- [ ] 对自旋、磁性、电子占据、对称性破缺和可能的中间亚稳相做敏感性检查。

### P5 出口标准

- [~] HfO2 T->PO 已有 100 Ry/10 au、7 总帧/5 内部帧的可恢复收敛普通 VCNEB（job 27678218）、CI 精修（job 27678507，焓垒 `0.1291722 eV`）和 9 总帧/7 内部帧普通对照（job 27678407，焓垒 `0.1562092 eV`），均通过独立审计；5-image 分支未收敛，独立重复路径和更广泛敏感性仍待完成。
- [x] BTO 已完成第二材料案例所需的最小验证矩阵、可恢复性、方向检查和结果归档。
- [ ] 结果足以支撑论文中的“方法可用性”图表，但暂不把单个案例称为普适性证明。

## 10. P6：模式引导与方向限制升级

### 功能分层

- [x] `Mode-guided initialization`：用声子/软模/用户模式改善初始路径，只影响初始 images。
- [x] `Mode-projected VCNEB`：每轮将允许的广义力和更新投影到模式子空间，适用于研究指定机制。
- [x] `Directional constraints`：已提供原子方向和 cell 方向的 basis 构造入口，并通过 `direction_basis_conflicts()` 诊断部分裁剪、完全失活和秩损失。
- [x] `Release-and-refine`：已加入解析耦合势示例，先做 projected 模式约束搜索，再解除约束做全空间 VCNEB 精修并报告能垒差异；真实材料对照仍待完成。

### 验收

- [x] 一个可解析验证的单模式模型，证明禁止方向的位移和力均为零。
- [x] 一个多模式非正交输入，证明正交化和投影结果与显式线性代数一致。
- [x] 已加入解析势的 `projected` 约束搜索到全空间 release-and-refine 示例；真实材料对照仍待完成。
- [ ] 一个真实材料案例，比较无约束、模式引导、严格模式约束和释放后精修。
- [x] `docs/theory.md` 和 README 明确约束路径得到的是受限路径上的鞍点，不能自动解释为全空间一阶鞍点。

## 11. P7：论文、软件包与发布

### 11.1 论文结构

- [~] Introduction：晶体相变能垒、固定 cell NEB 的局限、现有 VCNEB/USPEX/ABINIT/ASE 生态和纯 Python 的需求（论文草稿已写）。
- [x] Theory：广义坐标、cell 度量、应力到 cell force、切线、弹簧力、CI、约束和收敛判据。
- [~] Software：核心数据模型、calculator contract、VASP/ABACUS adapter、重启、并行 image 运行和模式 API（论文草稿已写，VASP 生产级路径仍待补）。
- [~] Examples：HfO2 T->PO、钙钛矿案例、模式引导与释放精修、能垒和路径结构；解析势、有限差分、ASE fixed-cell 对照、收敛矩阵和公开方法差异已纳入草稿，论文图表与旧实现定量对比仍待补。
- [~] Conclusions/Availability：应力精度、cell 参数化、原子映射、磁性/电子态、多路径问题、计算成本、版本、许可证、输入、结构、脚本、manifest 和复现命令已纳入草稿；完整发布包仍待补。

### 11.2 软件发布

- [x] 正式项目/发行名称已定为 `VARNEB`，上游仓库为 `xdzhu/varneb`；Python 导入包继续使用 `vcneb` 兼容名；许可证为 GPL-3.0-or-later，版本策略已落到 0.0.1。
- [x] 已整理 `pyproject.toml`、核心依赖和 plot/dev 可选依赖，并提供 `varneb --version` 入口（兼容保留 `vcneb` 别名）；VASP/ABACUS 继续作为 ASE calculator 运行时配置，不强制打包进核心依赖。
- [x] 已完成最小安装烟测（wheel 安装、toy VCNEB、CLI 帮助）；API/adapter/故障排查内容已并入 README 与 validation protocol。
- [x] GitHub About 与 `pyproject.toml` 的 Summary 已统一为 `VARiable-cell Nudged Elastic Band code with universal first-frinciples calculators`；已发布的 PyPI `0.0.1` 不可覆盖，下一次发布必须递增版本号，workflow 会把新简介写入 wheel/sdist。
- [ ] 将核心测试放入持续集成；真实 DFT 作为可选的集群复现实验，不要求 CI 内运行（尚未按要求开启 GitHub CI）。
- [ ] 提供最小可运行案例、HfO2 案例输入模板、模式文件模板和结果解析脚本。
- [~] 干净环境安装与从快照恢复测试已通过；干净节点运行仍以 HfO2 100 Ry/10au 端点和后续路径为最后验收项。

### 论文/发布出口标准

- [~] 论文主要 claim 已在 `paper/claim_evidence.md` 映射到代码、测试、输入摘要和结果文件；正式投稿图表与完整原始输出发布包仍待补。
- [ ] 新用户不阅读内部源码，仅按 README 就能完成 toy、model 和至少一个 calculator smoke test。
- [ ] 代码、数据和论文中使用的参数一致；没有手工修改但未记录的结果。

## 12. 立即执行顺序

下一轮按以下顺序推进，避免过早消耗大量 DFT 资源：

1. 在 `cu17` 上重新跑 P0 基线并生成 manifest；若节点不可用，切换到 `cu22`。
2. 完成 P1 的理论规范，重点核对 cell force/stress 的有限差分定义和 `cell_scale`。
3. 完成解析变胞多井势与联合原子-cell toy model，建立 P2/P4 的硬验收门槛。
4. 实现严格模式子空间和方向投影，补齐“模式引导”和“模式约束”的对照测试。
5. 在不改变核心 API 的前提下完成 VASP/ABACUS calculator contract 和单 image smoke test。
6. 准备 HfO2 T/PO 端点，做端点独立弛豫和原子映射审计。
7. 先用 Level A 跑 HfO2 粗路径，再按收敛矩阵推进 Level B/C；所有长任务只在指定计算节点运行。
8. 根据验证结果冻结 v0.1 API，再开始论文图表和软件发布材料。

## 13. 项目级 Definition of Done

只有同时满足以下条件，才把项目标记为“正确可用”而不是“原型可运行”：

- [ ] 理论定义、代码实现和测试中的坐标、单位、应力符号完全一致。
- [x] 原子力和 cell force 均通过有限差分；固定 cell 与 ASE 对照通过。
- [x] 解析模型能恢复已知 saddle/能垒，且 CI、弹簧、约束、缓存和恢复功能均有回归测试。
- [x] VASP 与 ABACUS 各有真实 calculator smoke test，且 ABACUS 已完成 HfO₂/BaTiO₃ 生产级 VCNEB。
- [x] HfO₂ T->PO 与 BaTiO₃ T->C 均完成端点审计、5/7/9-image（或等价）收敛矩阵和路径对照；HfO₂ 7/9 结果见 provenance manifest。
- [ ] 模式引导、严格模式约束、方向限制和释放后全空间精修均有清晰定义与实证案例。
- [~] HfO₂ 端点、5/7/9-image、CI 的代码/输入摘要、作业资源、结果、审计和 artifact 已集中在 `hfo2_validation_provenance.json`；完整论文 claim 映射与所有原始输出发布包仍待完成。
- [x] 用户可在不安装 MATLAB/USPEX 的情况下，从纯 Python 环境完成最小案例，并按文档切换 calculator；wheel 安装与 toy 恢复烟测已通过。

## 14. 每次迭代必须记录的字段

- 任务 ID、日期、代码版本、执行人/自动化会话。
- 节点、核数、Python/ASE/NumPy/SciPy/DFT 版本。
- 输入结构来源、原子 mapping、模式文件 hash、calculator 参数和收敛阈值。
- image 数、弹簧常数、cell 度量、优化器、最大迭代数和随机种子（如有）。
- 运行状态、失败 image、重试次数、wall time、最终能量/力/stress 和收敛判据。
- 结果目录、manifest、日志、图表、结论和下一步动作。

相关背景分析见：

- `README_VCNEB.md`
- `VCNEB_ITERATION_LOG.md`
- `outputs/uspex_vcneb_mode_analysis.md`
- `outputs/vcneb_literature_review.md`
- `outputs/vcneb_p0_baseline_manifest.json`
- `docs/theory.md`
