# VCNEB 理论与实现规范（研究草稿）

本文档描述当前 `vcneb/core.py` 的数学约定，目标是让实现、测试和论文使用同一套符号。它是研究版规范，不把尚未完成验证的参数化或约束算法包装成最终结论。

## 1. 结构空间

设体系有 `N` 个原子。ASE 使用行向量存储 cell：

\[
 h = \begin{bmatrix} \mathbf h_1^T\\
 \mathbf h_2^T\\
 \mathbf h_3^T \end{bmatrix},
 \qquad
 \mathbf r_a = \mathbf s_a h,
\]

其中 `s_a` 是原子分数坐标的行向量，`r_a` 是笛卡尔坐标。取第一个 image 的 cell 为参考 cell `h_0`，以 deformation gradient `F` 参数化当前 cell：

\[
 h = h_0 F^T, \qquad F = (h_0^{-1}h)^T.
\]

因此，一个 image 的广义构型为

\[
 Q=(s_1,\ldots,s_N,F).
\]

当前代码允许分数坐标不包回 `[0,1)`，以便在周期边界上保留连续路径；写文件或用户明确要求时才进行 wrap。

### 原子映射和周期分支

NEB 的端点必须有相同原子数、元素顺序和可审计的原子映射。若采用最小镜像插值，则端点差为

\[
 \Delta s_a = s_a^{(1)}-s_a^{(0)}-n_a,
 \qquad n_{a\alpha}=\operatorname{round}(s_{a\alpha}^{(1)}-s_{a\alpha}^{(0)})
\]

仅对周期方向执行上述分支选择。映射错误会产生看似收敛但化学意义错误的路径，因此不能仅由能量结果反推 mapping 正确。

## 2. 优化器坐标和扩展度量

优化器不直接使用 `(s,F)`，而是使用具有长度量纲的向量

\[
 x(Q)=\left(\operatorname{vec}(s h_0),\;
 c\operatorname{vec}(F-I)\right),
\]

其中 `c = cell_scale`，默认取参考 cell 体积的立方根。两点之间的当前路径距离为

\[
 d_{ij}=\|x_i-x_j\|_2.
\]

这一定义把原子部分和 cell 部分放进同一欧氏空间，但 `c` 是度量选择，不是材料常数。改变 `c` 会改变 NEB 的路径几何和弹簧分配，必须作为收敛/敏感性参数报告。

当前实现的 cell 插值是在 `F` 中线性进行：

\[
 F(\lambda)=(1-\lambda)F_0+\lambda F_1.
\]

由于 `F_0=I`，只要中间 determinant 保持正值，该插值具有明确的数值含义。对大剪切、晶格等价变换或近奇异端点，后续需要增加应变/对数晶格参数化，并通过对照测试决定默认策略。

## 3. 能量、焓和广义力

在外部静水压力 `P` 下，路径优化使用焓

\[
 H(Q)=E(Q)+P V(Q).
\]

代码要求 calculator 提供能量、原子 Cartesian forces 和应力张量。压力、能量和应力使用 ASE 内部兼容单位；对于实际 DFT calculator，适配层必须显式记录单位转换。

### 原子块

令 `f_a` 为 calculator 返回的 Cartesian 力，满足 `f_a=-\partial E/\partial r_a`。固定 cell 时，分数坐标共轭力为

\[
 g_{s,a}=f_a h^T.
\]

将其变换到优化器原子坐标 `x_a=s_a h_0` 后，得到

\[
 g_{x,a}=g_{s,a}h_0^{-T}.
\]

这对应 `fractional_force()` 和 `_force_to_x()` 的原子分支。

### cell 块

设 calculator 应力为 `\sigma`，体积为 `V`。当前 row-vector 实现先构造

\[
 W=-V(\sigma+P I),
\]

再将它变换到 `F` 的共轭广义力：

\[
 g_F = W F^{-T}.
\]

这里的符号和转置严格对应 `cell_force()` 中的
`solve(F, W.T).T`。这条公式必须用能量有限差分重新验证，尤其要覆盖非正交 cell、非零压力和应力张量非对角元；不同 DFT 程序的 stress/virial 符号不能默认相同。

由于优化器 cell 坐标是 `x_F=c(F-I)`，cell 块的优化器力为

\[
 g_{x_F}=g_F/c.
\]

因此 `_force_to_x()` 返回的完整广义力是

\[
 g_x=(g_{x,1},\ldots,g_{x,N},g_{x_F}).
\]

优化器梯度采用 `-g_x`，与 ASE 的 `get_gradient()` 约定一致。

## 4. VC-NEB 离散力

对内部 image `i`，记 `H_i=H(Q_i)`，`x_i=x(Q_i)`。当前实现使用能量加权的 improved tangent：

\[
 \tau_i = \operatorname{normalize}\left(
 |H_{i+1}-H_i|(x_{i+1}-x_i)
 +|H_{i-1}-H_i|(x_{i-1}-x_i)
 \right),
\]

在能量沿路径严格上升或严格下降时退化为相应的单侧切线。若切线范数太小，代码回退到正向差分；仍为零时返回零切线。

真实广义力在切线法向的分量为

\[
 g_i^\perp=g_i-(g_i\cdot\tau_i)\tau_i.
\]

相邻 image 距离为

\[
 d_i^+=\|x_{i+1}-x_i\|,\qquad
 d_i^-=\|x_i-x_{i-1}\|.
\]

弹簧力为

\[
 g_i^\mathrm{spring}=k_i(d_i^+-d_i^-)\tau_i.
\]

普通 VC-NEB 力为

\[
 g_i^\mathrm{NEB}=g_i^\perp+g_i^\mathrm{spring}.
\]

这意味着真实原子力和 cell 广义力先统一变换到同一个 `x` 空间，再做切线分解；不能分别对原子和 cell 使用两个互不相容的 reaction coordinate。

## 5. Climbing image

当前实现选取内部焓最高的 image 作为 climbing image。其力替换为

\[
 g_i^\mathrm{CI}=g_i-2(g_i\cdot\tau_i)\tau_i,
\]

即反转切线方向的真实力，不再叠加弹簧力。最终是否是一阶鞍点仍需检查全空间梯度、受限方向和 Hessian/负曲率诊断；“CI 收敛”本身不足以证明鞍点类型。

## 6. 掩码与约束

若原子 mask 为 `m_{a\alpha}\in\{0,1\}`，cell mask 为 `m^F_{\alpha\beta}`，当前实现只允许活动分量参与：

\[
 g_x\leftarrow m_x\odot g_x,
 \qquad
 \Delta x\leftarrow m_x\odot\Delta x.
\]

非活动分量在 `set_x()` 中保持 image 当前值，并不被重置到端点 cell。因而 atom/cell mask 是硬的分量冻结，不是一般线性模式约束。

## 7. 模式功能的三个层次

### 7.1 模式引导初始路径

给定归一化模式向量 `v` 和基准线性插值 `Q_i^0`，当前 `mode_guided_path()` 使用内部 image 的端点为零包络：

\[
 Q_i=Q_i^0+A\,b(\lambda_i)v,
 \qquad \lambda_i=i/(M-1),
\]

其中 `b` 可取 `sin(pi lambda)`、`4 lambda(1-lambda)` 或线性包络。端点严格不变；后续普通 VC-NEB 仍在全空间优化。因此它是初始路径生成器，不是约束 MEP。

### 7.2 严格模式子空间

严格模式约束定义活动子空间 `S=span{v_1,...,v_K}`，并将每轮广义更新投影到 `S`：

\[
 g_S=P_Sg_x,
 \qquad P_S=V(V^TV)^{-1}V^T,
\]

或在约化坐标中直接优化。当前代码通过 `build_mode_basis()` 将模式写成
VCNEB 的扩展 `x` 坐标列向量，再用 SVD 得到数值正交基；`constraint_mode="subspace"`
要求端点位移属于该子空间，并在初始化和每次 `set_x()` 时回投影。
`constraint_mode="projected"` 则以每个内部 image 的初始坐标为仿射参考，只约束优化更新。
对于质量加权模式，`V^TV` 必须替换为相应质量/扩展度量内积；当前实现的
`build_mode_basis()` 在归一化阶段显式处理质量加权输入，但严格的动力学质量度量仍需单独验证。
`mode_guided_path()` 仍然只是初始路径生成器，不能替代上述约束接口。

### 7.3 方向和原子约束

方向约束可以表示为每个原子的投影矩阵 `P_a`，或更一般的稀疏线性约束矩阵。
当前 `build_direction_basis()` 已支持每个原子一个方向、多个方向组合以及 cell deformation
方向，并交由 `VCNEB(mode_basis=...)` 完成正交化和投影。新增
`direction_basis_conflicts()` 会报告 mask 删除的分量、完全失活的 basis 列以及 mask 前后的秩损失；
它不替代大体系所需的稀疏 projector，也不改变先投影再进行 NEB 分解的当前顺序。

## 8. 收敛判据

当前 `fmax` 是所有内部 image 的优化器广义力分量范数最大值：

\[
 f_\mathrm{max}=\max_{i,a}\left|g_{x,i,a}\right|.
\]

生产级结果还必须同时报告：

- 最大真实原子力和最大 cell force/stress；
- 最高 image 的能量/焓、正向和反向能垒；
- image 间扩展距离和 reaction coordinate；
- 弹簧力与真实法向力的分解；
- `cell_scale`、image 数、弹簧常数、优化器和电子结构收敛设置。

`VCNEB.path_diagnostics()` 将上述量按 image 写成结构化表格，另外给出体积、晶格
长度和角度、最大原子力、最大应力、cell 广义力，以及内部 image 的真实切向/法向力、
弹簧力和最终 NEB 残余力。它用于区分“路径残余力变小”和“原始 DFT 应力/原子力已经
收敛”这两个不同判据。

## 9. 必须完成的理论验证

以下验证完成前，本文件中的方程只能作为实现草稿：

1. 原子 Cartesian force 到 fractional force 的中心有限差分。
2. 六个独立应变分量和非对角 cell 分量的 cell force 中心有限差分。
3. 非零压力下焓梯度和体积导数验证。
4. 固定 cell 时与 ASE NEB 的 force/tangent/spring 对照。
5. 不同 `cell_scale` 下扩展度量、路径长度和收敛行为的敏感性分析。
6. 严格模式投影与约化坐标实现的数值等价性。

## 10. 文献和实现定位

算法思想以变胞 NEB 文献中的扩展构型空间、cell 自由度和广义弹性带为理论背景；实现层参考 ASE 的 optimizer/calculator 契约、ASE UnitCellFilter 的应力到广义力处理，以及公开的 USPEX VCNEB 用户语义。最终论文必须逐项说明本项目与这些实现的坐标、cell 参数化、约束和 calculator 适配差异，不声称在没有源码证据时复现 USPEX 内部实现。
