# VARNEB 收敛加速策略与创新边界

## 结论先行

目前已经有第一性原理证据支持的贡献，是一个面向可变胞 NEB 的
**feasibility-aware local/staged optimization framework**，而不是“我们重新发明了
FIRE”。它把四个环节统一起来：

1. 每个内部 image 保存独立的 FIRE 状态和局部步长；
2. 原子自由度和晶胞自由度可以分块控制；
3. 粗阶段使用按 image 数归一化的较大信赖步长，接近收敛后切换到保守步长；
4. 每次提交昂贵 DFT 之前，对候选晶胞做原生几何可行性检查；失败时只缩短并重置失败 image 的候选步，不污染整条 band。

这套设计在 BTO 与 HfO2 的消融中减少了 calculator 工作，并在 VASP 长链上展示了
“拒绝危险候选、接受半步并继续”的恢复能力。其新意在于针对 variable-cell
DFT-NEB 的整体组合和可审计语义；不能把其中任何一个已有的优化器名字单独写成
“首次提出”。

## 已实现且有 benchmark 证据的策略

### 1. ImageScaledFIRE：按 image 数归一化全局步长

普通 FIRE 的 `maxstep` 约束整条 band 的欧氏范数。内部 image 越多，平均分到每个
image 的位移越小。`ImageScaledFIRE` 用 `sqrt(N_interior)` 修正全局步长，使单个
image 的有效位移尺度不随 image 数人为缩小。

证据：

- BTO n=9：8 步、772 s、约 65 次 image evaluation，最终
  `fmax=0.0781 eV/A`；普通 global FIRE 为 28 步、2669 s、约 205 次；
- HfO2 7-image 转移测试：30 步、3818 s、约 157 次，普通 FIRE 为 42 步、5242 s、
  约 217 次。

边界：这是合理的尺度修正，不是全新的 NEB 原理。论文中应称为
“image-count-normalized step control”，并与 global-cap control 做同起点消融。

### 2. BlockFIRE：逐 image 独立的 FIRE 状态

每个内部 image 维护自己的速度、时间步和位移上限。某个 image 的晶胞候选失败时，
只缩短该 image 并重置其状态，其余 image 不被最困难的局部 image 拖慢。

证据：

- BTO n=9，保守 `per-image maxstep=0.01`：5 步、478 s、约 44 次 evaluation、
  `fmax=0.0899 eV/A`；相对普通 FIRE，时间减少 82.1%，估计 calculator evaluation
  减少 78.5%；
- 同一 BTO 起点的 `per-image=0.02` 版本也在 10 步、963 s 收敛。

边界：image-local FIRE 不是 VARNEB 首创，且在 HfO2 n=20 T→PO 上曾平台化在约
`0.188 eV/A`。因此不能声称“对所有路径普遍加速”；可写成在 variable-cell DFT
中有效的局部状态隔离机制，并同时报告负 transfer。

### 3. SplitFIRE：原子/晶胞分块控制

原子位移与晶胞应变的刚度和尺度不同，`SplitFIRE` 为两类广义坐标使用独立步长、
速度和回退。BTO n=9 的 atom `0.02` / cell `0.01` 版本为 9 步、863 s、约 72 次
evaluation、`fmax=0.0995 eV/A`。

边界：当前证据表明它优于 global baseline，但没有证明它稳定优于最保守的
BlockFIRE；它更适合写为分块控制消融，不宜单独作为主创新。

### 4. StagedFIRE：粗到细的信赖步长

先用 image-count-normalized 的粗步长快速消除大尺度几何误差；当 force 降到预设
切换值（当前实验默认约 `0.12 eV/A`）时，重置动量并切换到 `0.02` 级别的保守步长。
切换复用已有力，不额外调用 calculator。

证据：

- HfO2 7-image：约 23 个粗阶段更新加 1 个精修更新，约 134 次 evaluation、3513 s；
  global FIRE 为 217 次、5242 s；
- 成功链的 barrier、reaction energy 和几何审计均保留在可接受范围内。

边界：分阶段信赖域思想本身并非全新；我们的贡献是把它和 variable-cell 候选验证、
image-local 状态及后端失败语义结合起来。

### 5. Feasibility-aware local backtracking：提交 DFT 前的局部门控

候选晶胞先经过原生晶格探针和最短距离/体积检查。完整步若会触发 VASP Bravais
判定或晶格病态，则不提交该候选；只在同一 image 上尝试半步等已授权缩短，成功后
继续。原始候选、拒绝原因和接受的替代步全部写入轨迹。

证据：

- HfO2 T→PO 的完整候选被 native VASP lattice gate 拒绝，半步通过并在一次接受更新
  后达到 `0.099347 eV/A`；
- exact-contract、零 backtracking 的对照在同一 image 处停止，说明收益来自候选
  步管理而不是改变 INCAR、SYMPREC 或物理路径；
- HfO2 PO→M 及 GaN 失败案例验证了“先查根因、同输入重试”和“几何失败不改物理参数”的
  语义。

边界：它首先是可靠性/安全收敛机制，而不是独立的能垒算法。论文可称为
“pre-DFT feasibility gate with local trust-region recovery”，并必须报告拒绝率和
回退次数。

### 6. Manager--worker 与端点缓存

manager 持有完整 band 状态，只把内部 image 发给 worker；固定端点只计算一次并缓存，
每轮仅分发内部 image。32 MPI/image-worker 的并行化显著降低 wall time，并避免端点
重复计算。

边界：这是重要的软件架构和可扩展性贡献，不是新的收敛数学。应报告总 image、内部
image、worker 数、实际 SCF/cache miss 和 wall time，不把并行度造成的加速写成算法
收敛步数减少。

## 已验证的整体 claim

最稳妥的论文表述是：

> VARNEB introduces a feasibility-aware, image-local and staged optimization workflow
> for variable-cell first-principles NEB. By isolating difficult images, normalizing
> displacement scales, and rejecting infeasible cell proposals before expensive
> calculations, it reduces calculator work on reproducible BTO and HfO2 ablations and
> provides auditable recovery on long VASP chains.

中文可表述为：

> VARNEB 面向可变胞第一性原理 NEB，提出一种具有几何可行性门控的逐 image、分阶段
> 收敛工作流；它通过隔离困难 image、校正 image 数导致的步长缩放，并在 DFT 前拒绝
> 不可行晶胞候选，在 BTO/HfO2 消融中减少 calculator 工作，同时为长链 VASP 计算提供
> 可审计的局部回退机制。

这里的“提出”指整体工作流和其可审计实现，不指 BlockFIRE/FIRE/信赖域的基本思想。

## 尚未实现，不能当前写成已验证创新

### A. 原子--应变耦合模态预条件器（最有潜力的理论创新）

计划把端点 Γ 力常数、弹性曲率和原子--应变交叉曲率组成广义 Hessian，求解
`K e = lambda G e`，用正定谱截断或低秩耦合近似作为只改变更新速度的预条件器。
关键是严格区分路径度量 `G` 与预条件器 `P`，证明不改变驻定带。

当前状态：理论设计已写入路线图，尚无代码、解析势和 DFT 消融证据。因此现在只能
称为“planned method direction”。如果完成解析势正确性、BTO/HfO2 DFT 调用数下降和
相同 barrier/path fingerprint 三重验证，它将是比 FIRE 组合更强的算法创新。

### B. Metric-consistent variable-cell path metric

当前 `cell_scale` 仍有超胞表示依赖。候选的强度化 Hencky/原子--应变度量可解决
复制、旋转、置换和晶格 gauge 的一致性，但必须先完成有限差分、超胞复制和等价 gauge
测试。当前不能把它写成已经验证的改进。

### C. 自适应 image 插入/删除

按段长、切线角、局部能量插值误差和峰值包络决定加密，而不是固定“越多越好”。这
已有 AutoNEB/dyNEB 等先例；VARNEB 的潜在差异在于同时使用可变胞几何误差和能垒误差，
并在最终以冻结 image 集重认证。当前仍是 planned/heuristic，不是已证实主 claim。

### D. 多保真 continuation / surrogate

先低成本后端粗搜，再把几何、mapping 和 gauge 传给 VASP 等高层级后端，清空缓存后
重新认证全部 image。已有相关先例；只有完成误差门槛、不可混缓存和高层级重认证后，
才能进入论文。

## 论文中建议的层级

1. **主创新（现在可以写，但需谨慎）：** variable-cell、calculator-agnostic、可审计的
   feasibility-aware local/staged workflow；用 BTO/HfO2 两套同起点消融支撑。
2. **软件创新：** manager--worker 并行、端点缓存、统一 calculator contract、输入/几何
   门控、精确重启和失败语义；强调可复现与节省节点小时。
3. **未来方法学主创新：** metric-consistent path geometry + atom--strain coupled-mode
   preconditioner；完成后再提升论文理论层级。
4. **辅助科学贡献：** Γ 模/广义模态投影和低维 PES，用于解释为何某条路径快或慢；它是
   机制分析，不等同于收敛加速算法。

## 下一轮 benchmark 的最低标准

- 固定相同端点、mapping、初始 chain、calculator 参数和 `0.10 eV/A` 阈值；
- 同时报告 outer steps、实际 SCF/cache-miss 次数、节点小时、拒绝/回退次数；
- 比较 Identity、ImageScaledFIRE、BlockFIRE、SplitFIRE、StagedFIRE 及组合；
- 至少 BTO、HfO2 和一个 VASP 长链；
- 使用 5--10 步窗口判断平台或反弹，不能因一次回弹提前停止；
- 要求相同路径 fingerprint、barrier/reaction energy 在离散/电子误差内一致；
- 保留失败 transfer，不能只保留最快成功分支。

详细原始数据和已有数值见 [`convergence_acceleration_benchmark.md`](convergence_acceleration_benchmark.md)。
