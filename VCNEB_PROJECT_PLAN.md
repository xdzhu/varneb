# VCNEB 项目标准化任务清单

> 项目目标：构建一个纯 Python、计算器无关、可复现的变胞 NEB（VC-NEB）工具包，支持 VASP、ABACUS 以及后续可插拔的 ASE Calculator；同时支持普通 VC-NEB、爬山图像（CI-VCNEB）、模式引导路径和严格的模式/方向约束，用于研究晶体相变过渡态与变胞能垒。
>
> 计划版本：v0.1，日期：2026-09-12。当前工作树是研究开发版，不代表已经达到论文发布或生产计算标准。
>
> 代码上游：<https://github.com/xdzhu/vcneb>。当前项目暂称 `vcneb`，正式软件名称留待发布前确定。

## 1. 执行原则

- 先定义广义坐标、内积、应力符号和单位，再扩展算法；所有公式都要有数值回归测试。
- 先用解析势和廉价计算器验证算法，再使用 DFT。不能把单次 DFT 跑通当作算法正确。
- 将“路径初始猜测的模式引导”和“优化过程中只允许沿模式运动”作为两个不同功能实现、测试和描述。
- 保持计算器适配层薄：VCNEB 核心只依赖统一的结构、能量、原子力和应力接口。
- 每个长任务都必须能恢复、能记录、能复核；每个科学结论都必须保留输入、代码版本、参数、节点和输出摘要。
- 集群规则：`235` 只作为登录网关；实际计算只在 `cu17`、`cu22`、`cu23`、`cu24`、`cu25`、`cu26` 上执行，每个节点 40 核，直接登录节点运行，不使用 `qsub` 或 PBS 调度，也不在本机 WSL 运行 DFT/NEB 长任务。所有计算节点共享目录和程序环境，代码与结果目录只需同步到共享路径一次。

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

- [ ] 尚未完成随机路径覆盖和不同 cell 参数化的完整对照；六个独立应变的系统有限差分误差表、解析基准与固定 cell ASE 对照已完成首轮。
- [ ] 尚未完成 VASP/ABACUS 两套真实 DFT 的端到端生产级 VCNEB 收敛案例。
- [x] HfO2 T 相到 PO 相的 12 原子结构来源、原子一一映射和端点独立弛豫 trial 已保存；生产精度路径和科学能垒仍未完成。
- [x] 已实现基础严格模式子空间和 projected-update 约束，并加入端点子空间验证、方向冲突诊断和解析势 release-and-refine；真实材料对照仍未完成。
- [ ] 尚未形成可投稿版本的误差预算、效率统计、软件发布包和论文结果表。

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
- [x] 以后每次启动计算前检查 `cu17`、`cu22`--`cu26` 的负载、进程和用户任务；选择空闲节点运行，忙节点不挤占。
- [x] 利用共享目录完成一次代码同步即可；节点无需安装 Git，运行目录和 manifest 均记录对应的 Git commit。
- [ ] 将 toy、model、VASP、ABACUS 四类示例分别标注为 `unit`、`model`、`DFT-smoke`、`production-template`，避免用户误把模板当成已收敛结果。
- [x] 检查当前 `/home/zhuxd/abacus/8.dielec/vcneb` 旧代码与工作树的差异，保留可借鉴算法说明和输入格式，但不复制无法验证的逻辑。
- [x] 规定结果目录命名：材料、端点、图像数、calculator、参数集、代码版本必须可从目录名或 manifest 读出。

### 产出与验收

- 产出：`baseline_manifest.json`、回归日志、环境清单、基线版本说明。
- 验收：新节点上从干净工作目录执行一条命令即可完成所有非 DFT 回归；失败时日志能定位到测试名和输入。

## 5. P1：理论与接口规范

### 5.1 广义坐标与几何

- [ ] 明确周期结构的基本对象：原子分数坐标 `s`、cell 矩阵 `h`、元素/质量、周期边界和原子映射。
- [ ] 明确路径变量是 `(s, h)` 还是 `(r, h)`；给出两者转换、最小镜像处理以及 cell 改变时原子笛卡尔位置的定义。
- [ ] 固化扩展空间内积，例如原子位移与 cell 变形的加权 Frobenius 度量；记录 `cell_scale` 的物理单位和默认选择依据。
- [ ] 规定 cell 插值的合法域：体积为正、避免奇异 cell、处理整数等价变换和不必要的整体旋转。
- [x] 已评估线性 cell 插值与对数/指数晶格插值；`interpolate_vcneb()` 提供 `linear`、正定 deformation 的 `log_strain` 和自定义回调，并保留线性方式作为兼容默认选项。

### 5.2 力、应力与 NEB 方程

- [ ] 统一能量单位、原子力单位、应力单位和 cell 广义力的符号约定。
- [x] 从能量对 `(s, h)` 的导数推导原子广义力和 cell 广义力，明确应力到 cell force 的转换、体积因子以及张量转置约定，并用非对角 deformation/外压有限差分验证。
- [ ] 定义切线选择、真实力的垂直投影、弹簧力的平行分量、端点处理和 climbing-image 规则。
- [ ] 明确 VC-NEB 的收敛量：最大真实广义力、最大原子力、最大 cell force/stress、能量变化和路径几何变化。
- [ ] 明确固定 cell、固定原子、固定原子方向、固定 cell 分量、软约束和硬约束之间的差别。

### 5.3 模式与方向约束语义

- [ ] 将模式功能拆为三种模式并在 API 中显式命名：`guided_initial_path`、`projected_dynamics`、`hard_subspace_constraint`。
- [ ] 定义模式在笛卡尔空间、分数坐标空间和质量加权空间的含义；定义多个模式正交化、归一化和线性组合规则。
- [ ] 支持原子方向 mask/投影矩阵；记录它与“冻结原子”的差别。
- [ ] 调研并记录 ABINIT GeoConstraints/directional constraints 的可复用语义，但不把 ABINIT 输入格式强行暴露给核心 API。
- [ ] 定义约束后的力是先投影再做 NEB 切向分解，还是先做 NEB 力分解再投影；选择一种并用有限差分验证。

### 出口标准

- [ ] 完成一份 `docs/theory.md` 或等价理论说明，包含所有符号、方程、单位和伪代码。
- [ ] 任何公共参数都能在理论说明中找到定义；任何实现中的经验参数都标注默认值、范围和敏感性。

## 6. P2：核心算法升级

### 6.1 插值与结构路径

- [x] 实现端点原子自动映射检查：按元素分组、用周期最小镜像几何代价求 permutation，并允许用户显式提供 mapping；`validate_atom_mapping()` 返回可审计位移报告。
- [ ] 实现分数坐标最小镜像插值，避免跨周期边界产生长路径。
- [x] 实现 cell 插值策略选择器：线性、正定 deformation 的对数应变插值和用户自定义插值，并对不适用的 log 路径明确报错。
- [ ] 对每个中间 image 检查 cell determinant、原子重叠、最小距离和异常应变；异常时在初始化阶段明确失败。
- [x] 已加入 calculator-free 的路径几何审计，报告体积、周期 MIC 最短距离和 deformation；可通过插值参数在初始化阶段硬拒绝碰撞/异常 deformation。
- [x] 已加入 HfO2 `linear`/`log_strain` 初始路径比较脚本，固定 image 数和端点映射后输出体积、最短距离和 deformation 分布。
- [ ] 增加基于模式的初始路径组合：结构插值 + 模式位移 + cell 模式，并保持端点严格一致。

### 6.2 NEB/VCNEB 力与优化

- [ ] 复核切线在能量单调、能量拐点和平坦区的行为；对比简单切线、能量加权切线和 improved tangent。
- [ ] 复核弹簧力在扩展空间的定义，保证原子与 cell 部分使用同一反应坐标度量。
- [x] 实现 CI-VCNEB，并在联合原子-cell 解析势上验证 climbing image 不引入端点漂移或错误的 cell 方向；真实 DFT CI 仍待收敛案例。
- [ ] 优先复用成熟优化器 API；核心只提供广义坐标、梯度和约束投影，不重复实现通用 LBFGS/FIRE 数值细节。
- [ ] 增加 line search 失败、calculator 异常、SCF 不收敛、NaN/Inf 和 cell 奇异的可恢复处理。
- [ ] 设计 image 级缓存和原子写入；中断后可以从最近快照继续，且不会混用不同参数集的结果。

### 6.3 严格模式子空间

- [x] 实现单模式和多模式的正交投影器；支持仅原子模式、原子+cell 联合模式。
- [x] 实现约束优化：每一步更新后回投影到允许子空间；约化变量优化仍可作为后续性能优化。
- [x] 实现方向约束投影矩阵，与 atom mask/cell mask 组合时通过 `direction_basis_conflicts()` 给出部分裁剪、完全失活和秩损失诊断。
- [ ] 对比三种结果：无约束 VCNEB、模式引导但最终无约束、严格模式约束；文档中明确三者不能互相替代。
- [ ] 加入约束释放功能，使用户可以先在模式子空间寻找路径，再切换到全空间做最终精修。

### P2 出口标准

- [ ] 解析势上能从随机初始路径收敛到已知 saddle/MEP。
- [ ] 对同一问题改变 `cell_scale` 后，经过一致的广义坐标定义，能垒和结构结果在合理范围内稳定。
- [ ] 固定 cell 时与 ASE 标准 NEB 的能量、力和路径趋势一致，差异有可解释的坐标/优化器原因。

## 7. P3：计算器无关层与工程化

### 7.1 Calculator contract

- [x] 固化最小 calculator 协议：输入结构、能量、原子力、应力/virial、单位和计算状态；新增 `vcneb.calculator` 预检与机器可读报告。
- [x] 实现 calculator capability 检查：是否支持 stress、是否支持 variable cell、是否能返回每 image 独立目录；支持重复目录检测。
- [ ] 将核心与 VASP/ABACUS 命令行、环境变量、MPI 命令完全隔离；命令 profile 只负责启动和结果解析。
- [x] 明确 stress 缺失时的行为：`run_vcneb()` 预检直接拒绝，不能静默把 cell force 当成零。
- [x] 运行时 calculator 异常会保留 image 编号、目录和命令上下文，并由回归测试覆盖。
- [x] 新增物理 HfO2 单 image ABACUS smoke 入口；已在集群完成真实 SCF、原子力和 stress 验证。

### 7.2 VASP 适配

- [x] 完成单 image 能量/力/stress 解析、目录隔离、重启文件策略和错误分类；真实 GaN VASP smoke 已记录到 `outputs/vasp_gan_single_image_manifest.json`。
- [ ] 明确 VASP 的 cell relaxation 输入只用于端点预弛豫还是也用于 image calculator；避免将普通变胞弛豫误用为 NEB image 更新。
- [x] 完成极小晶胞的 VASP 单 image smoke，并保留 INCAR/KPOINTS/POTCAR 来源记录；完整材料案例仍待开展。
- [x] `run_vasp_single_image_smoke.py` 已在 `cu26` 用 40 核 VASP 6.3.2 实跑，energy/forces/stress 与静态输入验证均通过。

### 7.3 ABACUS 适配

- [x] 完成 STRU/KPT/INPUT 生成和结果解析；核对 stress 输出、单位和晶格方向，并提供参数化多 image 入口；HfO2 单 image 与 7-image ABACUS smoke 已通过。
- [x] 已加入 `relax_abacus_endpoint.py`，使用 ASE `FrechetCellFilter` 完成 HfO2 两端点的低精度、宽松阈值独立原子/cell 弛豫 trial 并记录。
- [ ] 验证 ABACUS 命令 profile、MPI 进程数、退出码、超时和 SCF 不收敛处理。
- [ ] 让同一套 VCNEB 输入只更换 calculator 配置即可切换 VASP/ABACUS。

### 7.4 可恢复运行

- [x] driver 已支持保存 manifest、每 image 结构、energy、force、stress、cell 指标、真实/弹簧/NEB 力分解和收敛指标；每轮广义切线仍作为后续增强项。
- [ ] 支持从中断 image 级恢复、从最后完整 iteration 恢复和只重算失败 image。
- [ ] 增加 dry-run、validate-only、single-image 和 N-image smoke test 命令。
- [x] 生成机器可读的 summary JSON 和人类可读的文本摘要；CSV/Markdown 汇总仍作为发布增强项。

### 出口标准

- [x] VASP 与 ABACUS 各至少完成一个单 image smoke test；结果分别记录在 `outputs/vasp_gan_single_image_manifest.json` 和 `outputs/abacus_hfo2_smoke_manifest.json`。
- [ ] 断开后恢复不会丢失已完成 image，也不会重复覆盖有效结果。
- [x] calculator 失败时错误信息包含 image 编号、输入目录和命令上下文；建议动作由上层 launcher 根据 calculator 类型补充。

## 8. P4：数值正确性与收敛验证

### 8.1 单元和微分验证

- [x] 对原子坐标、cell deformation 分量和联合焓扰动分别做中心有限差分；已补充三种对角应变和三种对称剪切应变的六方向扫描。
- [x] 验证能量梯度、应力符号、cell force、mask 后梯度和模式投影梯度。
- [x] 覆盖非正交 cell、周期边界跨越、原子排列变化和极小体积保护。
- [x] 增加容差随有限差分步长变化的表格，覆盖零压力与有限外压，避免只测一个“碰巧通过”的步长。

### 8.2 解析模型

- [x] 构造含原子坐标和 cell 变量的可控多井势，已知 minimum、saddle 和能垒；toy 与 HfO2 12 原子 synthetic model 均已端到端运行。
- [ ] 构造含耦合项的 toy model，验证路径不能被错误地拆成独立原子和 cell 两条路径。
- [x] 在解析模型上对比不同弹簧常数、不同 image 数、FIRE/LBFGS 和 `cell_scale`；无弹簧与真实 DFT 组合仍待补充。
- [x] `examples/run_vcneb_convergence.py` 已在指定集群完成 54 组收敛矩阵并记录结果。
- [x] 收敛矩阵首次运行暴露并修正了测试端点漏设最终 cell strain 的问题；修正版 54/54 通过。
- [x] 解析双井模型验证 CI-VCNEB 将最高 image 推向已知 saddle，并通过 `saddle_diagnostics()` 给出负切线曲率和力残差；真实 DFT 鞍点诊断仍待收敛路径。

### 8.3 参考实现与理论对照

- [x] 固定 cell：与 ASE NEB/CINEB 在同一 calculator、同一端点和同一 image 数下比较；解析模型逐 image 能量一致，能垒差为 9.33e-6 eV。
- [x] 变胞：依据 Qian 等 VCNEB 论文的广义坐标、cell 力和弹簧思想完成首轮逐项对照；有限变形 stress measure 和大应变差异仍需补充。
- [ ] 依据 USPEX VCNEB 公开手册对比输入语义和用户流程；不声称复现其内部实现，因为当前获得的安装包不含可审计源码。
- [ ] 对比旧 `/home/zhuxd/abacus/8.dielec/vcneb` 实现的结果和失败模式，保留可复现实验而非凭印象判断。

### 8.4 收敛矩阵

- [ ] image 数：至少 5、9、13 或按路径长度调整。
- [ ] `cell_scale`：至少三组具有物理解释的数值。
- [ ] 弹簧常数：弱、中、强三组，并记录路径是否出现 image 堆积。
- [ ] 电子结构精度：粗略、生产、加严三组；检查能垒误差与 SCF 噪声。
- [ ] 优化器和步长：至少两种成熟优化器，记录迭代数、失败次数和 wall time。

### P4 出口标准

- [x] 所有梯度/应力微分测试通过，并有零压与有限外压的误差表和 provenance manifest。
- [x] 解析模型的能垒、端点反应能、saddle 位置和负切线曲率在目标容差内恢复；完整 Hessian 诊断仍不在当前实现范围内。
- [ ] 至少一套真实计算案例对 image 数、cell_scale、弹簧和电子精度表现出可解释的收敛趋势。

## 9. P5：材料案例与集群运行

### 9.1 案例选择与结构准备

- [ ] 主案例：HfO2 T 相到 PO 相，优先使用 conventional cell、12 原子、元素顺序一致的端点。
- [ ] 备用/对照案例：选择一个原子数少、cell 明显变化且原子可一一对应的钙钛矿相变，例如 BaTiO3 或 SrTiO3 的结构相变；最终材料和相名以可靠结构来源核实为准。
- [ ] 为每个案例保存原始结构来源、空间群/晶格参数、原子映射、端点能量和独立弛豫设置。
- [x] 已加入并运行端点独立原子/cell 弛豫 trial；HfO2 两端点在低精度、宽松 `0.1 eV/A` 阈值下通过，生产精度端点仍待加严。
- [ ] 检查相变路径是否包含必要的原子重排、剪切、体积变化和模式分量，避免端点只靠一个不物理的线性缩放连接。

### 9.2 分层运行策略

- [x] Level A 的 model/低精度 DFT smoke 分支已确认输入、路径、计算器、每 image 目录、输出解析、恢复和 CI 独立分支；生产精度仍待完成。
- [ ] Level B：生产精度 VCNEB，完成无约束路径和 CI-VCNEB。
- [x] Level B 首轮已在 cu17 完整执行 7-image/40-step FIRE，并从完整轨迹尝试 LBFGS 恢复；两条分支均明确记录为未收敛，不能作为物理能垒。
- [x] 已补做输入预检失败诊断和 `FIRE(maxstep=0.05)` 保守步长重试；预检在 ABACUS 启动前正确拒绝错误 basis 路径，重试仍在前七步出现残余力单调增长，因此暂停继续消耗 DFT 资源，转入路径/应力诊断。
- [x] 已对最新完整重启轨迹做 7-image 零步静态力/应力诊断；image 2 的最大原子力/应力与 image 3 的最高焓错位，结果已进入 manifest，下一步优先改进机制路径。
- [ ] Level C：加密 image、提高电子精度和改变初始路径，验证能垒与 saddle 的稳定性。
- [x] Level A 的 model/低精度 DFT smoke 与 Level B 诊断 trial 已记录 Git 版本、节点、核数、输入、输出和结果摘要；Level C 仍待开展。

### 9.3 集群执行规范

- [ ] 登录流程固定为先 `ssh 235`，再登录目标计算节点；所有节点共享 `/home/zhuxd` 下的项目目录和软件环境。
- [ ] 计算节点只使用 `cu17`、`cu22`、`cu23`、`cu24`、`cu25`、`cu26`；每节点按 40 核配置 MPI/并行参数。
- [ ] 每次使用节点前先检查负载和现有任务；若有任务在运行则换到其他空闲节点，不通过提高优先级或超额并行挤占资源。
- [ ] 不在 `235` 上启动 DFT，不在本机 WSL 上启动长任务，不使用 `qsub`/PBS。
- [ ] 每个案例独立目录运行，避免多个 image 共用临时文件或覆盖同名结果。
- [ ] 启动前做磁盘空间、环境、可执行文件、赝势/PAW、K 点和 MPI 健康检查。
- [ ] 运行中记录每轮完成 image 数和预计剩余时间；失败时先保留目录和 stdout/stderr，再决定是否重算。

### 9.4 科学结果判据

- [ ] 给出端点相对能量、最高 image/saddle 能量、能垒（正向和反向）、反应坐标和 cell 演化。
- [x] 代码和 ABACUS driver 已支持给出最高 image 的残余原子广义力、cell 广义力/应力、真实力与弹簧力分解；HfO2 生产路径仍需在收敛后生成最终数值记录。
- [ ] 给出路径结构图、晶格参数/体积/剪切角曲线、关键键长和模式投影。
- [ ] 至少用两种初始路径或两组 image 数复算；若得到不同 saddle，必须报告而不是挑选更漂亮的一组。
- [ ] 对自旋、磁性、电子占据、对称性破缺和可能的中间亚稳相做敏感性检查。

### P5 出口标准

- [ ] HfO2 T->PO 至少有一套可收敛、可恢复、可重复的 ABACUS 或 VASP VCNEB 运行。
- [ ] 至少一个第二材料案例完成相同的最小验证矩阵。
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
- [ ] 文档明确：约束路径得到的是受限路径上的鞍点，未必是全空间的一阶鞍点。

## 11. P7：论文、软件包与发布

### 11.1 论文结构

- [ ] Introduction：晶体相变能垒、固定 cell NEB 的局限、现有 VCNEB/USPEX/ABINIT/ASE 生态和纯 Python 的需求。
- [ ] Theory：广义坐标、cell 度量、应力到 cell force、切线、弹簧力、CI、约束和收敛判据。
- [ ] Software：核心数据模型、calculator contract、VASP/ABACUS adapter、重启、并行 image 运行和模式 API。
- [ ] Examples：HfO2 T->PO、钙钛矿案例、模式引导与释放精修、能垒和路径结构；解析势、有限差分、ASE fixed-cell 对照、收敛矩阵和旧实现/公开方法的差异作为本节的小节、表格或图展示。
- [ ] Conclusions/Availability：应力精度、cell 参数化、原子映射、磁性/电子态、多路径问题和计算成本，以及版本、许可证、输入、结构、脚本、manifest、原始输出和复现命令。

### 11.2 软件发布

- [ ] 选定正式名称、Python 包名、许可证、版本策略和引用方式；当前暂用项目代号 `pyVCNEB`，名称未最终确定。
- [x] 已整理 `pyproject.toml`、核心依赖和 plot/dev 可选依赖，并提供 `vcneb --version` 入口；VASP/ABACUS 继续作为 ASE calculator 运行时配置，不强制打包进核心依赖。
- [ ] 增加最小安装测试、API 文档、tutorial、calculator adapter 文档和故障排查页。
- [ ] 将核心测试放入持续集成；真实 DFT 作为可选的集群复现实验，不要求 CI 内运行。
- [ ] 提供最小可运行案例、HfO2 案例输入模板、模式文件模板和结果解析脚本。
- [ ] 发布前做一次干净环境安装、干净节点运行和从快照恢复测试。

### 论文/发布出口标准

- [ ] 论文中每个主要 claim 都有代码、输入和结果文件对应关系。
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
- [ ] 原子力和 cell force 均通过有限差分；固定 cell 与 ASE 对照通过。
- [ ] 解析模型能恢复已知 saddle/能垒，且 CI、弹簧、约束和恢复功能有回归测试。
- [ ] VASP 与 ABACUS 至少各有一条真实 calculator 路径通过 smoke test；至少一种完成生产级 VCNEB。
- [ ] HfO2 T->PO 和一个第二材料案例完成端点审计、收敛矩阵和重复路径检查。
- [ ] 模式引导、严格模式约束、方向限制和释放后全空间精修均有清晰定义与实证案例。
- [ ] 论文结果、代码版本、输入、原始输出、节点信息和结论可以一一对应复现。
- [ ] 用户可在不安装 MATLAB/USPEX 的情况下，从纯 Python 环境完成最小案例，并按文档切换 calculator。

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
