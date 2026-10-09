# 固定基底平面的端点与 VCNEB

这是机械边界的实现和验证，不是新的外延理论。普通 NEB 默认阈值仍为
0.10 eV/Å；此接口不设置或改变计算器的截断、赝势、轨道、k 点或 SCF 参数。

## 明确固定什么

ASE 晶格矢量是矩阵 H 的行，VARNEB 使用 `H = H0 @ F.T`。
夹持 **H 的第 0、1 行**，而不是只将对应应力分量设为零。
设这两根矢量的单位法向为 n，允许的变形严格满足
`delta F = outer(v, n)`，因为两根面内矢量与 n 正交。

- `allow_tilt=True`：v 有三个独立分量，第三根矢量可改变长度和两个倾斜分量。
- `allow_tilt=False`：v 平行于 n，只释放法向长度，原有面外倾斜不变。
- 两种选择都是不同的机械边界，因此参数必须显式提供。

固定两个非共线矢量排除了整体晶胞旋转；开放的面外剪切并非任意刚体转动。
对斜胞/旋转基底，正确约束一般不能写成几个 Cartesian 元素的二值 mask。
接口返回全一 mask **加上**正交子空间；不得只用返回的 mask。
它与 `symmetric_inplane_vcneb_boundary` 的二维材料边界相反：后者释放面内应变并固定真空。

## 使用同一边界优化端点和路径

```python
from ase.optimize import BFGS
from vcneb import VCNEB, ClampedPlaneFilter, clamped_plane_vcneb_boundary

# substrate_cell 的前两行来自事先规定的同一基底。
# initial/final 的前两行必须已经相同；它们不能是不同面内晶格的自由胞端点。
boundary = clamped_plane_vcneb_boundary(
    len(initial), substrate_cell, allow_tilt=True,
)
target = ClampedPlaneFilter(initial, boundary, cell_scale_A=cell_scale)
BFGS(target, logfile="initial.opt.log").run(fmax=0.03, steps=endpoint_step_cap)
# final 同样优化并审核相身份，不能直接复用自由胞能量。

neb = VCNEB(images, cell_scale=cell_scale, **boundary.vcneb_kwargs(images))
```

`vcneb_kwargs(images)` 在 VCNEB 构造、任何内部投影及计算器调用**之前**拒绝
不兼容的原始链。全部 images 与端点共用基底；运行后的候选也检查相同约束。
图像准备、周期 lift、最短距离和相身份检查仍需独立完成；本接口不猜测原子映射。
同一基底下的端点若弛豫为另一个相，应记录原相消失，不通过额外模式冻结强造该相。

端点 filter 可用 ASE BFGS，优化所有原子及开放的晶格分量，拒绝未明确组合的额外
ASE 约束。负体积等非法步在修改原结构前拒绝。它不自动缩步或切换计算参数。
开放分量需接近零外加应力，**夹持反力应力允许非零**：不能以整个应力张量都小于
2 kbar 作夹持端点门禁。原始应力、开放子空间梯度、原子受力分别留档。
NEB 残差、真实切向梯度和原始应力也不能混为同一个收敛量。

`boundary.open_traction(stress, pressure=0)` 返回 eV/Å³ 的开放应力向量。
释放第三矢量的全部三个分量时，向量为 `(sigma+P I) @ n`；仅开放法向时，
取该向量沿 n 的投影。用向量范数给出与坐标旋转无关的门槛，而非任意选择
一个 Cartesian 元素。第三矢量移动 w 的应力功为
`substrate_area * dot(w, open_traction)`（w 在开放方向内）。
该关系已在一般斜胞中与完整应力功公式交叉验证。

`ClampedPlaneFilter(..., candidate_validator=guard)` 可组合逐步几何检查；
检查的是实际变形后的晶胞和原子位置，在修改原 Atoms 前拒绝危险步。
默认不自动缩步，不改变任何计算参数，也不以对称化“修复”路径。

HfO₂ 有界流程见
[共同基底端点种子](../benchmarks/hfo2_channels/20261008/clamped_endpoint_seeds/README.md)：
5类结构×2种基底，共10个几何种子，均无计算结果。
`scripts.relax_clamped_ase_endpoint` 用 BFGS 分别检查原子受力和开放应力；
兼容已验证的 ASE 旧/新收敛入口，检查点复用当前几何的缓存而不重复调用计算器。
达到数值门槛后仍须独立审计相身份、变体以及真实SCF输入/输出，不能仅凭文件名标相。

## 应力功和势垒应变响应

`cell_work_derivative(stress, H, dH_dlambda, pressure=0)` 给出固定分数坐标的
`d(E+PV)/dlambda = V (sigma+P I) : (dH_dlambda.T @ inv(H.T))`。
压力正号代表压缩；ASE 应力拉伸为正，单位 eV/Å³；返回单位是 eV/参数。
对等双轴应变 ε，`dH_depsilon[:2] = H[:2]/(1+epsilon)`，第三行偏导为零。
图表使用 meV/f.u. 时需再乘 1000 并除同一胞的 f.u. 数。

只有已在开放自由度上驻定、且跟踪同一连续分支的端点/鞍点，才能通过包络定理
将此偏导解释为已优化势阱/鞍点的应变导数。普通 NEB 最高像不自动满足这点。
预测 `dB/dε` 要相减鞍点与**共同初态**的导数；独立 ε=+0.5% 的检验不能参与拟合。
夹持 ε=0 与全释放 P=0 是不同系综，不作跨系综有限差分。

## 同一边界的联合曲率探针

```python
from vcneb import ActiveJointCurvatureCoordinates

chart = ActiveJointCurvatureCoordinates.for_clamped_plane(
    reference, cell_scale_A=cell_scale, boundary=boundary,
)
probe = chart.displaced(delta)  # len(delta) = 3*N + boundary.cell_dofs，单位Å
# 用真实同参数DFT的F/stress调用chart.enthalpy_gradient，再构造局部Hessian。
```

这里直接复用端点/VCNEB的精确开放子空间，倾斜开放时为3个形变方向，
仅法向开放时为1个。倾斜方向一般不对称，不能用六个对称应变方向替代；
后者在一般基底上会改变面内矢量并采样另一个边界。传入不兼容参考结构或
评估结构时直接拒绝，不投影修复，不调用计算器。原子序和周期lift不变。

新接口无约束调用默认6个对称方向；原`JointCurvatureCoordinates`及其历史
源码哈希保持不变。显式`deformation_basis=np.empty((0,3,3))`
给固定胞纯原子坐标。给定基矢需Frobenius正交归一，输入不被自动归一化；
基矢在内部复制并只读。额外ASE原子约束被拒绝，不能悄悄改动宣称的原子自由度。
该接口只生成探针/计算梯度，尚不代表已取得HfO₂联合Hessian或认证任何鞍点。

复现几何/实现检查：`python -m scripts.verify_clamped_joint_coordinates`。
检查10个已有未弛豫HfO₂几何种子及Cu/EMT斜胞、旋转和已形变点的有限差分功；
不会运行DFT，也不会提交作业。EMT的非零压力只是分析校验，不修改HfO₂的P=0。

后续方向采样使用[联合曲率探针与推断范围](JOINT_CURVATURE_PROBES.md)：
少量选定方向也保留完整梯度及未采样方向的耦合，不把k×k投影当全Hessian。
采样和解析独立于计算器，配对物理梯度不能替换为NEB投影力或弹簧力。

## 当前验证范围

`tests/test_epitaxial_boundary.py` 包括一般斜胞/旋转基底、全三个面外自由度、
不兼容内部像和实际 HfO₂ 自由胞端点的预拒绝、反力与开放应力分离，以及 EMT 下的
能量—应力有限差分、已变形点的 filter 梯度、每步固定基底的 BFGS。
这些是几何/实现验证，不是 HfO₂ 应变势垒或有利应变窗口的材料证据。
通用命令行现在也支持**显式选择**此边界；旧的自由胞/固定胞默认行为不变。

## 生产入口与配置

端点必须先在规定基底下优化并完成相/变体审计。路径入口不替端点弛豫，也不把
自由胞端点投影成夹持端点。引用文件规定晶胞前两行；普通路径坐标尺度也必须
显式给出，并与端点和后续曲率坐标一致。HfO₂登记值为5.12968067458423 Å。

```bash
python -m vcneb.material_runner \
  --initial path/to/audited-initial/POSCAR \
  --final path/to/audited-final/POSCAR \
  --workdir path/to/fresh/clamped-path-preflight \
  --clamped-plane-reference path/to/T-substrate/POSCAR.seed \
  --clamped-allow-tilt true --cell-scale 5.12968067458423 \
  --cell-mode full --no-align-cells --cell-interpolation linear \
  --mapping identity --mic --require-continuous-periodic-lift \
  --n-images 9 --no-climb --fmax 0.10 --validate-only
```

9个总像包含7个内部像和两个固定端点。`--validate-only`不加载计算器符号、
不构造计算器、不运行DFT；仅写几何与边界来源收据。真正计算须另在合规资源
分配内提供原计算器/工厂及**原参数**，不由此命令提交自己。
`--resume-snapshot`和`--initial-chain`均先审核未经改动的所有像。
不兼容内部像、倾斜被禁止却出现倾斜、额外ASE约束均拒绝，不事后投影修复。
每个优化候选同时检查基底、所选几何门槛和已有晶胞步长限制；FIRE缩步仍只在
显式配置的次数内进行，不重试电子故障或改变物理输入。

公共`varneb prepare`/`varneb run ... --execute`使用同一入口。配置中增加：

```json
{
  "cell_mode": "full",
  "cell_scale_A": 5.12968067458423,
  "clamped_plane_reference": "path/to/T-substrate/POSCAR.seed",
  "clamped_allow_tilt": true,
  "align_cells": false,
  "cell_interpolation": "linear"
}
```

这是完整配置中的边界片段，不是独立可执行配置；相对引用路径按配置文件所在
目录解析。`true`允许第三矢量的三个方向，`false`只开放法向长度，必须显式选择。
禁止自动晶胞对齐与log-strain插值，避免旋转/插值改动基底。夹持路径仍需要真实
stress，不能当作纯固定胞NEB。当前不支持与全局模式子空间artifact隐式组合；
此组合会明确报错，需另外定义与验证交集，不静默截去自由度。

`vcneb_preflight.json`和最终summary记录引用文件SHA256、固定行、法向、
倾斜策略、开放自由度及尺度；既有目录不能悄悄改成另一机械边界。
`tests/test_clamped_plane_cli.py`验证公共配置/直接入口、原始链预拒绝、斜胞
逐步基底不变、串行/双worker端点只算一次，以及计算器参数/命令原样传递。
10个未弛豫HfO₂种子仅用于几何检查，不能称作G2计算结果。

外加等双轴应变的响应探针另外使用
[受控参数坐标](BIAXIAL_CONTROL_CURVATURE.md)：ε参与偏导，但不进入端点/NEB的
可释放变量；其梯度保留基底反力。固定ε时原39维HfO₂可动空间不变，新增的
第40坐标只用于受控响应采样，不能把它计入固定ε的TS负方向或当自由胞优化。
