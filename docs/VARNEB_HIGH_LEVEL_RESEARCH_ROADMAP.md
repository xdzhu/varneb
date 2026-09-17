# VARNEB 高水平方法与论文研究路线图

**版本：** 2026-09-16  
**状态：** 研究规划；其中标为“候选”的公式和算法尚不是已验证的软件能力  
**目标：** 把 VARNEB 从“可运行的多后端 VCNEB 实现”推进为一个具有明确理论命题、可量化加速、可解释相变机制和完整证据链的方法学工作。

## 0. 核心判断

VARNEB 不应把“重新实现 VCNEB”当作理论创新。Qian 等的 VCNEB 和 Sheppard 等的 G-SSNEB 已经建立了原子与晶胞自由度共同参与路径优化的基本框架；预条件 NEB、AutoNEB、dyNEB 和高斯过程加速也分别已有先例。我们的高水平论文必须回答一个更窄、但更扎实的问题：

> 能否在一个**表示一致的可变胞路径度量**下，把路径投影、原子--应变耦合模态、收敛预条件和自适应计算统一起来，并且在多个第一性原理后端上证明它既不改变目标路径，又显著减少昂贵的计算器调用？

建议把最终工作组织成四层贡献：

1. **统一计算器契约。** 同一 VARNEB 核心只消费能量、原子力、应力和可追溯的静态 image 目录；ABACUS、VASP、QE、OpenMX 等只是 adapter。
2. **表示一致的可变胞路径几何。** 给出超胞复制、整体平移/旋转、晶格基变换和同元素置换下的协变或不变性，并修正当前经验 `cell_scale` 的理论缺口。
3. **不改变驻定路径的模态预条件。** 用原子--应变耦合曲率改善优化条件数；路径几何的度量 `G` 与只影响更新速度的预条件器 `P` 必须严格分离。
4. **模式分辨的相变机制。** 不只报告能垒，还回答哪些原子软模、均匀应变和两者的耦合驱动了路径，以及低维模式势能面在多大程度上重构了 VCNEB 路径。

自适应 image 和多保真计算是重要的加速层，但只有在完成严格消融后才进入主 claim。它们不能代替第 2、3 层的理论主线。

## 1. 当前基线与必须冻结的事实

### 1.1 已有能力

- 普通 NEB 到 VCNEB 的扩展坐标、焓、应力到晶胞广义力、improved tangent、普通/CI 力、固定端点、重启、精确状态缓存和 manager--worker 执行已经实现。
- 公共默认收敛阈值为 `0.10 eV/Å`。更严格阈值只能作为显式研究协议，不能暗示其必然具有更高物理意义。
- `n_images` 的当前代码语义是**总 image 数**。因此现有 BTO `n_images=7` 是“2 个固定端点 + 5 个内部 image”，每轮只计算 5 个内部 image。后续 CLI 和论文应同时写出 `n_images_total` 与 `n_images_interior`，不再只写含混的“7-image”。
- ABACUS 已有 BTO T→C 和 HfO2 T→PO 的材料路径证据；BTO 已有 5/7/9-total-image 敏感性结果和 cubic 端点 Γ 点力常数/本征矢投影。
- VASP 已完成 BTO T→C 的 7-total-image 普通 VCNEB：5 个内部 image 在 cu17 上逐 image 使用 40 MPI，端点由已验证静态缓存固定；最终广义力 `0.0951479 eV/Å`，达到 `0.10 eV/Å`，路径无内部势垒，反应焓约 `0.0426441 eV`，审计无错误。该远端证据仍需同步为仓库内 manifest 和论文 source data。
- QE adapter、preflight 和静态收敛门槛已存在，但按当前研究决定暂停真实 QE 生产路径。
- OpenMX 尚未实现通用静态能量/力/应力 adapter，现阶段不得写成已支持后端。

### 1.2 当前理论债务

当前优化器坐标为

\[
x(Q)=\left(\operatorname{vec}(s h_0),\;c\operatorname{vec}(F-I)\right),
\qquad c=\Omega_0^{1/3}\;\text{(default)}.
\]

若把同一均匀路径复制成含 `m` 倍原子的超胞，原子块的欧氏范数按 `m^{1/2}` 增长，而当前晶胞块按 `m^{1/3}` 增长，二者相对权重改变 `m^{-1/6}`。因此现有 `cell_scale` 敏感性检查并不能证明路径对超胞表示不敏感。G-SSNEB 已明确讨论原子/晶胞自由度随体系尺度的平衡；修补这一点首先是**正确性工作**，不是可单独宣称的创新。

### 1.3 论文 claim 的成熟度

| Claim | 当前成熟度 | 升级条件 |
|---|---|---|
| 开源、可重启、证据可审计的 VCNEB | 可写 | 保持回归测试、发布工件和 source data 一致 |
| calculator-agnostic 核心 | 架构已成立 | 至少 ABACUS/VASP 两条相同 BTO 路径归档；QE/OpenMX 按实际状态表述 |
| 多后端材料验证 | 部分成立 | VASP 证据入库并重绘跨后端图；未跑的后端不列为验证完成 |
| 超胞/表示一致的路径度量 | 尚未成立 | 完成推导、解析测试、复制超胞测试和材料敏感性验证 |
| 更快收敛 | 尚未成立 | 预注册消融证明减少昂贵 calculator 调用且路径不变 |
| 原子--应变耦合模态机制 | 原子 Γ 模基础已成立 | 完成广义 Hessian、简并子空间跟踪和低维 PES 验证 |
| 自适应 image / 多保真加速 | 构想 | 独立消融、误差门槛和失败回退齐全 |

## 2. 与前人工作的边界

| 已有方向 | 前人已解决的核心 | VARNEB 不能宣称 | 仍可形成的缺口 |
|---|---|---|---|
| SSNEB / G-SSNEB / VCNEB | 晶胞与原子共同参与固--固相变路径 | “首次提出变胞 NEB” | 可审计的统一度量、现代 calculator contract、表示不变测试和模式接口 |
| 预条件 MEP 方法 | 利用近似曲率降低条件数 | “首次给 NEB 加预条件” | 面向原子--应变耦合空间、保持物理路径度量的预条件及其实证 |
| AutoNEB / dyNEB | 动态加 image、选择性优化未收敛 image | “首次自适应 image”或“首次跳过已收敛 image” | 以可变胞几何误差和能垒误差为共同驱动的安全策略 |
| GPR/代理模型 NEB | 用不确定性驱动昂贵真值调用 | “首次多保真/代理 NEB” | 面向异构第一性原理后端的不可混缓存协议和最终高层级重认证 |
| 声子/正规模态路径分解 | 用 Γ 模或低维模态解释结构变化 | “首次用声子解释相变路径” | 统一原子--应变广义模态，并让同一模态既用于解释又用于预条件 |

因此推荐的主创新不是五个功能的并列罗列，而是一个连贯命题：

> **VARNEB 在表示一致的扩展构型流形上定义物理路径，再用同一原子--应变耦合模态体系构造不改变驻定带的预条件，并输出可重构的模式分辨相变机制。**

多后端、自适应执行和可追溯缓存是使该命题可复现、可验证、可使用的工程体系。

## 3. 理论主线 A：表示一致的扩展路径度量

### 3.1 候选构型与应变变量

设第 `a` 个原子的分数坐标为 `s_a`，参考晶胞为 `h_0`，当前晶胞满足

\[
h=h_0F^T.
\]

候选新表示使用旋转不敏感的 Hencky 应变

\[
E=\frac{1}{2}\log(F^TF),
\]

并以

\[
Q=(s_1,\ldots,s_N,E)
\]

表示一个 image。这样整体刚体旋转不进入应变块；但晶格基变换下的协变性仍必须单独证明，不能由旋转不变性自动推出。

### 3.2 候选强度化度量

以每原子体积定义长度尺度

\[
L_0=(\Omega_0/N)^{1/3},
\]

候选路径线元为

\[
d\ell^2=
\frac{1}{N}\sum_{a=1}^{N}\|ds_a\,h_0\|^2+
\alpha^2L_0^2\|dE\|_F^2.
\]

其中 `α` 是无量纲、必须报告的原子/应变相对权重。对于同一均匀路径的整数超胞复制，`N` 与 `Ω` 同时乘以复制倍数，原子均方位移和 `L_0` 均保持不变，因此该候选度量具有明确的复制不变性目标。

需要验证的不是“不同尺寸超胞得到完全相同机制”。更大超胞可能允许小胞不存在的成核或畴壁机制。应验证的是：**把同一均匀路径仅作表示性复制时度量和离散 VCNEB 残差不变；若路径改变，应能归因于新增物理自由度，而不是归一化伪影。**

### 3.3 商空间与 gauge

路径距离应在以下等价操作下保持物理一致：

- 周期整体平移；
- 整体刚体旋转；
- 同元素原子置换；
- 等价晶格基的整数幺模变换；
- 均匀超胞复制。

实际实现不必一次求解完整的抽象商空间，但必须给出确定性的对齐顺序：

1. 枚举/求解允许的晶格对应；
2. 去除整体旋转；
3. 联合同元素 assignment 与周期平移选择最短 gauge；
4. 在统一 gauge 下计算 `E`、路径距离和切线；
5. 将所选 mapping、translation、lattice transform 和 hash 写入 manifest。

### 3.4 Riemannian NEB 表述

令候选度量矩阵为 `G(Q)`，定义

\[
\langle u,v\rangle_G=u^TGv,
\qquad
\|u\|_G^2=u^TGu.
\]

切线必须在 `G` 下归一化，`\langle\tau,\tau\rangle_G=1`。焓的 Riemannian 梯度为

\[
\operatorname{grad}_G H=G^{-1}\frac{\partial H}{\partial Q},
\]

真实力的法向投影写为

\[
f_\perp=f-\tau\langle\tau,f\rangle_G.
\]

弹簧长度、improved tangent 和 CI 反转也必须使用同一个 `G`。不能只改距离而继续用旧欧氏点积，否则“新度量”只是局部补丁而不是一致算法。

### 3.5 工作包 T1：度量正确性

**代码设计**

- 新增 `vcneb/metric.py`：`PathMetric` 协议及 `LegacyEuclideanMetric`、`IntensiveHenckyMetric`。
- `VCNEB` 内部只通过 `metric.displacement()`、`inner()`、`norm()`、`tangent()`、`force_from_covector()` 操作广义坐标。
- 旧 `cell_scale` 路径保留为兼容模式；新度量在完成验证前必须显式启用，不直接改变 `v0.0.1` 行为。
- summary 记录 `metric_name`、版本、参数、参考结构 hash 和 gauge 选择。

**必须测试**

1. 原子与晶胞梯度的中心有限差分，包括非正交晶胞、剪切和非零压力。
2. 固定晶胞极限与普通 ASE NEB 的切线和残差对照。
3. 整体平移、旋转、原子 permutation 和幺模晶格基变换的等价性。
4. `1×1×1 → 2×1×1 → 2×2×1` 均匀复制的路径长度、单位切线、投影残差和 barrier 不变测试。
5. 小胞与大胞允许不同非均匀机制时，测试不得错误地把物理差异判为失败。

**晋级门槛**

- 解析势有限差分相对误差 `≤10^-6`（绝对小量附近另用绝对误差门槛）。
- 等价 gauge 下的路径长度和残差只在浮点容差内变化。
- 复制路径的每 image 强度化残差一致；DFT 数值噪声必须单独估计。
- 在这些条件满足前，论文只能把该度量写为 proposal，不能写为已验证方法。

## 4. 理论主线 B：不改变驻定路径的原子--应变模态预条件

### 4.1 为什么要分开 `G` 与 `P`

`G` 定义什么是路径长度、切向和法向，因此会影响“要找的路径”。预条件器 `P` 只定义如何更快接近同一个驻定带。推荐更新形式为

\[
\Delta Q=P(Q)R_{\mathrm{VCNEB}}(Q),
\]

其中 `R_VCNEB` 已完全按 `G` 构造。若 `P` 在活动子空间中正定且可逆，则

\[
P R_{\mathrm{VCNEB}}=0
\quad\Longleftrightarrow\quad
R_{\mathrm{VCNEB}}=0,
\]

因此预条件器只改变收敛轨迹，不改变目标驻定带。任何同时改变 `G` 和 `P` 的实验都无法判断加速来自哪里，必须拆成独立消融。

### 4.2 广义原子--应变 Hessian

在平衡端点或经验证参考结构附近，定义

\[
K=
\begin{bmatrix}
K_{uu} & K_{uE}\\
K_{Eu} & K_{EE}
\end{bmatrix}.
\]

- `K_uu`：原子位移 Hessian，可由 Γ 点力常数构造；
- `K_EE`：应变--应变曲率，与弹性刚度及外压条件有关；
- `K_uE`、`K_Eu`：内部应变耦合，分别由原子力对应变、应力对原子位移的有限差分得到。

广义模态满足

\[
K e_m=\lambda_m G e_m.
\]

纯原子 Γ 声子只是 `K_uu` 的质量加权问题，不能自动替代上式。平移零模、整体旋转、固定 mask、虚模和数值近零模必须在明确的活动空间内处理。

### 4.3 候选预条件器

可从谱截断形式开始：

\[
P=\sum_m
\frac{e_me_m^TG}{\max(|\lambda_m|,\lambda_{\min})},
\]

并配合最大 `G`-范数步长、信赖域或自适应时间步。对不稳定参考端点，不能直接用带符号逆 Hessian 沿任意负曲率放大更新；首版应使用 `|λ_m|`、谱移位或正定近似，另将虚模用于路径初始化和解释。

需要比较三类低成本近似：

1. **对角块预条件：** `K_uu` 与 `K_EE` 分开，忽略耦合；
2. **低秩耦合预条件：** 只保留最软的若干原子/应变模及 `K_uE`；
3. **更新型预条件：** 端点曲率初始化，随后用 band 历史作受限 BFGS/L-BFGS 更新。

### 4.4 工作包 A1：预条件实现与消融

**代码设计**

- 新增 `vcneb/preconditioners.py`：`IdentityPreconditioner`、`BlockCurvaturePreconditioner`、`LowRankModePreconditioner`。
- `run_vcneb()` 接收 `metric=` 与 `preconditioner=`；优化器日志分别记录原始残差、预条件残差、接受步长和信赖域裁剪。
- 预条件器只能读取结构/曲率 artifact，不得直接改变 calculator 结果或路径能量。
- 重启文件保存预条件器版本和状态；版本不匹配时安全回退到 identity，而不是静默复用。

**消融矩阵**

| 编号 | 路径度量 | 预条件 | 自适应 image | 目的 |
|---|---|---|---|---|
| B0 | legacy | identity | off | 当前基线 |
| B1 | intensive | identity | off | 只测新度量的正确性/影响 |
| B2 | intensive | block | off | 原子/应变块曲率收益 |
| B3 | intensive | low-rank coupled | off | 耦合模态的额外收益 |
| B4 | intensive | coupled + history update | off | 最终加速候选 |

**主评价指标**

- 第一指标：昂贵 calculator 的 energy/force/stress 调用次数；
- 第二指标：总 SCF 迭代数、节点小时、峰值失败重试数；
- 正确性约束：同一收敛阈值、barrier/reaction enthalpy、路径距离、最高 image 身份、最小原子距离和最终残差；
- wall time 只作次级指标，因为排队、文件系统和 MPI 抖动会污染比较。

**可发表的加速门槛**

- 在预注册的至少两个材料路径和多个确定性初始链上，昂贵 calculator 调用数的中位数降低 `≥30%`；
- 所有成功分支达到相同 `0.10 eV/Å` 验收，且 barrier 差异不超过预先测得的电子结构/离散 image 误差；
- 不增加几何失败、SCF 失败或错误路径分支比例；
- 若只在解析势上加速、DFT 上无稳定收益，只能作为实现探索，不能进入 abstract 的主 claim。

## 5. 科学主线 C：从 Γ 声子到原子--应变耦合路径机制

### 5.1 三种必须区分的对象

1. **端点正常模：** 只在驻点 Hessian 上定义，可有稳定或虚频模。
2. **路径位移/切线投影：** 非驻定 image 上是集体坐标分解，不应直接称为该 image 的声子。
3. **广义原子--应变模：** `K e=λGe` 的本征方向，包含均匀应变和内部位移耦合，是本项目拟新增的分析对象。

### 5.2 路径投影量

对参考结构 `Q_0` 和 `G`-正交模式 `e_m`，定义

\[
a_m(s)=\langle e_m,Q(s)-Q_0\rangle_G,
\]

\[
p_m(s)=|\langle e_m,\tau(s)\rangle_G|^2.
\]

`a_m` 描述累计结构振幅，`p_m` 描述局部路径方向由哪些模式驱动。必须同时报告截断基的重构残差

\[
\epsilon_K(s)=
\frac{\|\Delta Q-\sum_{m=1}^{K}a_me_m\|_G}{\|\Delta Q\|_G}.
\]

对简并或近简并模，追踪的是子空间 projector 和主角，而不是任意相位/任意基下的单根本征矢。跨 image 跟踪可用相邻子空间 overlap 的最大匹配和 SVD/Procrustes 对齐。

### 5.3 模式势能面

首批只做能直接回答机制问题的 1D/2D PES：

\[
H(q_1,q_2)=H_0+
\frac12\kappa_1q_1^2+
\frac12\kappa_2q_2^2+
\gamma q_1q_2+
\beta_{12}q_1^2q_2^2+\cdots.
\]

建议同时保留：

- **frozen-mode PES：** 其余自由度固定，辨认裸耦合；
- **relaxed-mode PES：** 正交补空间放松，检验低维坐标是否足以描述路径；
- **VCNEB overlay：** 把每个 image 投影到 `(q_1,q_2)` 平面，比较路径是否沿低谷演化。

不以“最低频率排序”自动选择模式。候选模式应由端点位移、路径切线参与率、对称性和前向选择的重构改进共同筛选。

### 5.4 工作包 M1--M3

**M1：广义 Hessian 构造**

- 新增 `vcneb/generalized_modes.py`；读取标准化原子力常数、弹性矩阵和内部应变耦合。
- 为原子位移、对称应变、应力和单位约定建立单一 schema。
- 在解析耦合势上验证 `K_uE=K_Eu^T`、广义正交性、平移零模和有限差分步长收敛。

**M2：路径模式跟踪**

- 扩展现有 `vcneb/phonons.py`，输出 `a_m(s)`、`p_m(s)`、简并子空间 overlap、累计解释率和重构残差。
- BTO：从 cubic 不稳定 Γ 子空间出发，加入 tetragonal 端点模态，定量分离 Ti--O 极化位移与 tetragonal strain。
- HfO2：以 T/PO 端点和最高内部 image 为锚，判断 barrier 是否由单一软模、多个模式协同，或显著的原子--应变耦合控制。

**M3：低维 PES 与反事实路径**

- 对排名最前的 1--3 个广义模式做 frozen/relaxed 扫描。
- 比较仅用模式子空间生成的路径、释放后的 VCNEB 路径和无约束基线。
- 若低维模式不能重构路径，应如实报告多模/非谐机制；不得为了图形简洁强行归结为单软模。

**晋级门槛**

- 完整基重构误差达到数值精度；截断基必须报告 `ε_K(s)`。
- 简并子空间结论不依赖本征矢相位或基选择。
- 选定低维坐标能在明确误差范围内重现实路径能量轮廓，或给出其失败的物理原因。
- “某模式驱动相变”至少需要同时得到路径参与率、低维 PES 和结构可视化三类证据支持。

## 6. 加速层 D：自适应 image 与多保真计算

### 6.1 自适应 image 数量

固定 image 越多不一定越好：成本随内部 image 数增长，过密 image 还会引入高度相关的 SCF；太少则可能漏掉弯曲路径和狭窄势垒。image 的增加/删除应由误差指示器驱动，而不是由一个固定经验数决定。

候选每段指标：

\[
\eta_i=w_d\,\Delta\ell_i+
w_\theta(1-\cos\theta_i)+
w_H|H_{i+1}-2H_i+H_{i-1}|+
w_R\max(R_i,R_{i+1}).
\]

- `Δℓ_i`：`G` 下的段长；
- `θ_i`：相邻切线夹角；
- 二阶焓差：检测窄峰和强非线性；
- `R_i`：局部 VCNEB 残差。

策略：先用少量 image 得到粗路径，只在高 `η_i` 段插入；仅当相邻两段都平滑、能量插值误差低且删除后不改变峰值包络时才删除。每次拓扑改变必须迁移 optimizer 状态或安全重启，并保留旧完整链快照。

**安全门槛**

- 同时满足最大段长、最大切线角、局部焓插值误差和最终残差门槛；
- 至少一次使用更密固定 image 作盲审参考；
- 对 HfO2 的 barrier 差必须落在固定 7/9-total-image 离散误差内；
- 若动态图像改变了路径分支，只允许用于初始路径搜索，最终结果必须以冻结 image 集重认证。

### 6.2 多保真 continuation

推荐先实现“分阶段 continuation”，不立即实现在线 GPR：

1. 用较便宜但物理一致的 calculator/参数形成路径；
2. 只把几何、mapping 和 gauge 传给高层级 calculator；
3. 清空能量/力/应力缓存，在目标层级重新计算全部固定端点和内部 image；
4. 在目标层级继续松弛到同一收敛门槛；
5. source data 同时保存低/高层级轨迹，不混用 barrier。

不同后端、赝势、cutoff、k mesh 或泛函的缓存 namespace 必须完全隔离。若低层级路径落入错误拓扑分支，高层级 continuation 应有回退到多初始链或直接 DFT 基线的机制。

在线 GPR/机器学习代理只有在上述协议稳定后再立项，因为其不确定性校准、应力学习、晶胞表示和训练数据持久化会显著扩大论文范围。

## 7. 多后端工程路线

### 7.1 通用契约，而不是后端名单

每个 adapter 必须满足同一静态 image 契约：

- 输入：结构、独立工作目录、不可变 calculator 参数和 provenance；
- 输出：总能量、Cartesian force、完整 stress tensor、退出状态和输入/赝势 hash；
- 禁止后端自身改变 image 几何：VASP `IBRION=-1/NSW=0`，QE `calculation='scf'`，ABACUS 静态 SCF，OpenMX 使用等价静态设置；
- 能量、力、应力单位和符号必须由有限差分/静态审计验证；
- 端点只计算一次并缓存，迭代只调度内部 image；
- manager 拥有路径状态，worker 只拥有单 image calculator 调用。

### 7.2 后端优先级

1. **VASP：** 立即把已完成 BTO 路径的远端 summary、manifest、endpoint identity 和 source data 入库，更新过时 capability 表。
2. **ABACUS：** 继续作为 HfO2 加速消融和 BTO 模态分析的主生产后端；不重跑已完成端点。
3. **OpenMX：** 在理论/预条件主线稳定后实现。先完成单 image energy/force/stress 约定与有限差分，再做 BTO 7-total-image 正向路径。
4. **QE：** adapter 保留，真实路径按当前决定暂停；除非重新给出明确优先级，不占用近期计算预算。

“universal first-principles calculators”应解释为一个可扩展协议，而不是声称所有 DFT 程序已经材料级验证。论文必须把 `adapter implemented`、`static validated`、`material-path validated` 分成三个等级。

## 8. 预注册验证矩阵

### 8.1 层级 0：解析与合成势

| 算例 | 主要问题 | 必须输出 |
|---|---|---|
| 固定胞双势阱 | 退化到普通 NEB 是否正确 | ASE NEB 力/切线逐分量对照 |
| 原子--应变耦合双势阱 | `K_uE`、广义模态和预条件是否正确 | 解析 Hessian、本征值、路径和迭代数 |
| 复制超胞势 | 度量是否表示一致 | 1×/2×/4× 复制的长度、残差、barrier |
| 旋转/置换/gauge 变体 | 商空间对齐是否确定 | 相同 path fingerprint 与数值结果 |
| 窄峰/弯曲二维 PES | adaptive image 是否漏峰 | 与高密固定 band 的盲审误差 |

### 8.2 层级 1：BTO

定位：低成本、多后端、软模主导、T→C 正向路径无内部势垒。

- 冻结相同端点、mapping、7-total-image 初始链和 `fmax=0.10 eV/Å`。
- ABACUS 与 VASP 比较的是路径形状、相对焓、体积/应变演化和模式坐标；不同赝势下不要求绝对能量逐 meV 相等。
- 用 BTO 验证 calculator contract、metric 表示一致性和模式分析；由于该正向路径无内部峰，不把它作为 CI 或 barrier 加速的唯一证据。
- OpenMX 完成 adapter 后只需补同一条正向路径，不扩展为新的材料课题。

### 8.3 层级 2：HfO2

定位：真实变胞、存在内部势垒、用于预条件和 adaptive image 的主压力测试。

- 复用已接受端点和已完成普通 VCNEB 路径，不重新做无意义的端点优化。
- 从同一完整链快照、同一 calculator 参数和相同随机性运行 B0--B4 消融。
- 报告调用次数、SCF 迭代、节点小时、收敛曲线、barrier、最高 image、cell 轨迹和失败率。
- 模式分析聚焦于原子--应变耦合是否解释鞍点附近路径弯曲和收敛慢方向。

### 8.4 可选层级 3：普通固定胞对照

若正文需要一个标准 NEB 的清晰内部势垒，可使用 BTO `+P→-P` 固定胞极化翻转，展示 VARNEB 在冻结 cell 时退化为普通 NEB。它只服务理论闭环，不取代 HfO2 主案例。

## 9. 统一评价、误差预算与停止规则

### 9.1 误差预算

每条材料结果至少拆分：

1. 电子 SCF 误差；
2. cutoff/k mesh/赝势误差；
3. image 离散误差；
4. VCNEB 残差误差；
5. metric/`α` 敏感性；
6. 端点和原子 mapping 不确定性；
7. 不同后端产生的模型差异。

算法 A 与 B 的 barrier 差小于上述组合误差时，应写“在当前误差分辨率内一致”，不能写“完全相同”。

### 9.2 反回弹停止策略

单步残差回弹不是停止理由。普通 NEB、VCNEB 和 FIRE/BFGS 都可能短期非单调。统一停止/失败判断应使用：

- 已达到目标 `fmax`，正常收敛；或
- 在预先设定窗口内无统计显著改善，同时步长/能量/几何显示平台；或
- 连续多步恶化并触发几何、SCF、非有限值或信赖域失败；或
- 达到计算预算后保留最佳完整链，并明确标为未收敛。

建议保存 `best_state` 与 `last_state`，平台判断使用至少 5--10 个完整优化步的窗口，并把阈值和窗口写入 summary。任何新预条件器都必须在相同停止规则下比较。

### 9.3 image 数的界定

初始 image 数由路径复杂度而非“越多越好”决定：

- 端点位移/应变近线性、焓单调且切线平滑：5 个内部 image 通常足够作首轮；
- 有窄 barrier、路径弯曲或多个局部极值：先 5 个内部 image，再由段长、切线角和能量曲率插入；
- 最终用更密固定链验证 barrier 和路径拓扑，而不是把更密链当作默认生产配置。

报告必须写成“`N_total = N_interior + 2`”。计算资源按 `N_interior` 分配，端点不进入每轮 worker 池。

## 10. 实施顺序、依赖与止损

### P0：证据和术语冻结（立即）

- [ ] 把 VASP BTO 远端结果同步到 `outputs/`/manuscript source data，并更新 capability 表。
- [ ] 将所有“7-image”改成“7 total = 5 interior”或显式字段。
- [ ] 为 baseline 生成不可变 run manifest 和 path fingerprint。
- [ ] 冻结本路线图中的 claim 词典：implemented / statically validated / material-path validated。

**完成定义：** 任何论文数字都能回到 committed summary、输入 hash 和代码 commit。

### P1：度量抽象与表示一致性

- [ ] 实现 `PathMetric` 和 legacy 兼容层。
- [ ] 实现 intensive Hencky 候选度量。
- [ ] 完成 gauge、超胞复制和固定胞极限测试。
- [ ] 在解析势与现有 BTO 静态轨迹上做 `α` 敏感性，不先启动新 DFT。

**止损：** 若幺模基变换/一般晶格对应无法可靠自动化，首稿把适用范围限制为已预对齐的同构晶胞，并将自动晶格匹配移到后续版本。

### P2：曲率与预条件

- [ ] 在解析势上实现 block 与 low-rank coupled preconditioner。
- [ ] 用已有 Γ 力常数构造 BTO 原子块；补最少数量的 strain/cross 有限差分。
- [ ] 先用便宜 calculator 或已有数据做迭代消融，再提交 HfO2 DFT 对照。
- [ ] 通过后再研究 history update 和信赖域。

**止损：** 若构造 `K_uE` 的成本超过节省的 DFT 调用，改用端点一次计算、低秩更新或经验 block；若 DFT 调用减少不足 30%，不写“accelerated”，只把它定位为稳健的尺度平衡。

### P3：广义模式和机制图

- [ ] 完成 generalized-mode schema、简并子空间跟踪和重构误差。
- [ ] BTO 生成原子模/应变/能量联合图。
- [ ] HfO2 生成 barrier 邻域的广义模式参与率和 1D/2D PES。
- [ ] 制作模式动画与可编辑 source data。

**止损：** 若低维基无法解释 HfO2，论文转而报告定量多模混合和非谐耦合，不强造单模故事。

### P4：adaptive image

- [ ] 先在窄峰合成势实现插入/删除和状态迁移。
- [ ] BTO 只做无势垒平滑路径 sanity check；HfO2 做主要收益测试。
- [ ] 最终冻结 image 集重认证。

**止损：** 若动态图像频繁改变路径分支或优化器状态不稳定，只保留“粗路径初始化后一次性加密”，不进入在线主循环。

### P5：OpenMX 与多后端收尾

- [ ] 建立 OpenMX 静态 adapter 和能力审计。
- [ ] 做单 image stress 有限差分和目录隔离。
- [ ] 仅在 preflight 全绿后跑 BTO 7-total-image 正向路径。
- [ ] QE 继续暂停，除非后续重新排序。

### P6：多保真（stretch goal）

- [ ] 先实现离线 continuation 和不可混缓存。
- [ ] 对 HfO2 比较低层级 warm start 与直接目标层级。
- [ ] 在线 GPR 仅在前述主线完成且论文篇幅允许时启动。

### P7：论文和发布

- [ ] 自动生成 claim-to-evidence 表、参数表和全部 source data。
- [ ] 完成 blind rerun：从 manifest 重建至少一条解析路径和一条材料路径。
- [ ] 更新手册、案例、API 文档和 PyPI 元数据；按需触发发布，不因普通 push 自动发布。
- [ ] 投稿前做一次“前人工作重叠审计”和一次“证据不足措辞审计”。

## 11. 论文的推荐收敛形态（约 8 页）

### 11.1 最推荐的单篇主线

**题眼：** *Metric-consistent and mode-aware variable-cell nudged elastic bands with universal first-principles calculators*。

正文只保留四个核心问题：

1. 为什么普通 NEB 不能直接处理晶胞自由度；
2. 为什么可变胞路径需要表示一致的共同度量；
3. 为什么同一原子--应变模态可以既加速收敛又解释机制；
4. 这一设计如何跨 ABACUS/VASP（及完成后 OpenMX）保持相同核心算法。

自适应 image 和多保真若证据不够强，放补充材料或 outlook，不稀释主线。

### 11.2 页面与图表预算

| 内容 | 建议篇幅 | 主图 |
|---|---:|---|
| Introduction + prior-art boundary | 0.8 页 | — |
| 普通 NEB → metric VCNEB 理论 | 1.6 页 | Fig. 1：构型、度量和投影 |
| 模态预条件与软件架构 | 1.1 页 | Fig. 2：`G`/`P` 分离和 manager--worker |
| 数值验证 | 1.0 页 | 小表：不变性和解析误差 |
| BTO 多后端与模式机制 | 1.3 页 | Fig. 3 |
| HfO2 barrier 与加速消融 | 1.5 页 | Fig. 4 |
| Limitations + conclusion | 0.7 页 | — |

详细推导、所有 convergence 曲线、5/7/9 image、adapter 参数和动画放补充材料。八页限制下，不要同时把 adaptive image、GPR、多起点全都写成等权主创新。

### 11.3 投稿前最低组合

高水平版本至少需要以下组合之一：

- **推荐组合：** 表示一致度量 + 耦合模态预条件的 DFT 加速 + BTO/HfO2 模式机制；
- **保守组合：** 表示一致度量 + 多后端材料验证 + 完整模式机制，但不宣称加速；
- **不建议作为最终稿：** 只有多后端 adapter、作业调度和已有 VCNEB 复现。这可以形成软件发布，却不足以支撑更高层级的方法学论文。

如果广义度量/预条件的数学与实验内容膨胀，允许形成两篇工作：第一篇 CPC 聚焦可靠软件和多后端证据；第二篇方法/材料论文聚焦 metric-preconditioned VCNEB 与模式机制。但在证据不足前不提前拆稿。

## 12. 论文 claim 的最终验收表

| 拟写句型 | 必须具备的证据 | 不满足时的降级表述 |
|---|---|---|
| “calculator-agnostic” | 核心无后端分支；至少两个真实 DFT 后端材料路径；第三方 adapter 协议 | “calculator-interface based; demonstrated with ABACUS and VASP” |
| “representation-consistent” | gauge、旋转、置换、晶格基、超胞复制测试 | “uses an explicit configurable extended-space metric” |
| “accelerates convergence” | 预注册 DFT 消融，调用数中位数降低 ≥30%，结果等价 | “improves conditioning on analytic benchmarks” |
| “mode-resolved mechanism” | 广义 Hessian、子空间稳定性、投影残差、PES overlay | “projects the path onto endpoint Γ modes” |
| “adaptive image efficiency” | 密固定链盲审、barrier/拓扑误差与成本对比 | “supports heuristic image refinement” |
| “universal” | 清晰的协议与可扩展测试，不等于所有后端均验证 | 标题保留软件愿景，正文限定已验证后端 |

## 13. 下一步的最小闭环

近期不要同时铺开所有方向。最小闭环按以下顺序推进：

1. 归档已完成 VASP BTO 路径并修正文档现状；
2. 用合成势实现 `PathMetric`，首先证明超胞复制与 gauge 一致性；
3. 在同一合成势实现 `G`/`P` 分离的耦合模态预条件；
4. 只在解析消融通过后，用 HfO2 已有链做 ABACUS DFT 加速对照；
5. 同步构造 BTO 的原子--应变广义模式和二维 PES；
6. 根据第 4、5 步的数据决定论文以“加速”为主还是以“机制”为主；
7. adaptive image、OpenMX 和多保真按依赖顺序加入，QE 暂停。

这一路线的关键不是堆功能，而是让一个理论对象——扩展空间的度量与广义模态——贯穿算法、实现、验证和材料物理解释。只有这样，VARNEB 的优势才不是“又一个 VCNEB 脚本”，而是一套可验证、可扩展、能产生新物理分析的相变路径框架。

## 参考文献边界

1. Qian *et al.*, “Variable cell nudged elastic band method for studying solid–solid structural phase transitions,” *Computer Physics Communications* **184**, 2111–2118 (2013), [doi:10.1016/j.cpc.2013.04.004](https://doi.org/10.1016/j.cpc.2013.04.004).
2. Sheppard *et al.*, “A generalized solid-state nudged elastic band method,” *Journal of Chemical Physics* **136**, 074103 (2012), [doi:10.1063/1.3684549](https://doi.org/10.1063/1.3684549).
3. Caspersen and Carter, “Finding transition states for crystalline solid–solid phase transformations,” *PNAS* **102**, 6738–6743 (2005), [doi:10.1073/pnas.0408127102](https://doi.org/10.1073/pnas.0408127102).
4. Makri, Ortner, and Kermode, “A preconditioning scheme for minimum energy path finding methods,” *Journal of Chemical Physics* **150**, 094109 (2019), [doi:10.1063/1.5064465](https://doi.org/10.1063/1.5064465).
5. Kolsbjerg, Groves, and Hammer, “An automated nudged elastic band method,” *Journal of Chemical Physics* **145**, 094107 (2016), [doi:10.1063/1.4961868](https://doi.org/10.1063/1.4961868).
6. Lindgren, Kastlunger, and Peterson, “Scaled and Dynamic Optimizations of Nudged Elastic Bands,” *Journal of Chemical Theory and Computation* **15**, 5787–5793 (2019), [doi:10.1021/acs.jctc.9b00633](https://doi.org/10.1021/acs.jctc.9b00633).
7. Garrido Torres *et al.*, “Low-Scaling Algorithm for Nudged Elastic Band Calculations Using a Surrogate Machine Learning Model,” *Physical Review Letters* **122**, 156001 (2019), [doi:10.1103/PhysRevLett.122.156001](https://doi.org/10.1103/PhysRevLett.122.156001).
8. Okenyi, Ratcliff, and Walsh, “Multi-phonon proton transfer pathway in a molecular organic ferroelectric crystal,” *Physical Chemistry Chemical Physics* **23**, 2885–2890 (2021), [doi:10.1039/D0CP04236F](https://doi.org/10.1039/D0CP04236F).

详细来源、使用范围和证据边界见同目录的 `VARNEB_HIGH_LEVEL_RESEARCH_ROADMAP.provenance.md`。
