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

## 当前验证范围

`tests/test_epitaxial_boundary.py` 包括一般斜胞/旋转基底、全三个面外自由度、
不兼容内部像和实际 HfO₂ 自由胞端点的预拒绝、反力与开放应力分离，以及 EMT 下的
能量—应力有限差分、已变形点的 filter 梯度、每步固定基底的 BFGS。
这些是几何/实现验证，不是 HfO₂ 应变势垒或有利应变窗口的材料证据。
当前是 Python API；通用命令行尚未自动选择此机械边界。
