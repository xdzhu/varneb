# GaN 600 eV 局域能量—应力一致性审计

本审计只重新读取 GaN B4→B1、45.7 GPa 的现有 VASP 原始结果；没有运行
DFT，也没有改 `ENCUT`、k 网格或其它电子参数。36 个 Hessian 扰动和中心
静态的 `INCAR/KPOINTS/POTCAR` 哈希均与原 600 eV 生产输入一致。
可复现脚本是
[`scripts/audit_gan_600eV_energy_stress_consistency.py`](../scripts/audit_gan_600eV_energy_stress_consistency.py)，
逐轴、逐输出的数值和来源哈希在
[`gan_600eV_energy_stress_consistency_20260928.json`](../benchmarks/numerical_integrity/gan_600eV_energy_stress_consistency_20260928.json)。

## 核对模型与单位

沿广义坐标 `q`（Å）比较直接中心差分

`[H(q+h)−H(q−h)]/(2h)`，其中 `H=E+PV`、`h=0.02 Å`，

与中心及两侧 VASP 力/应力给出的广义焓梯度的 Simpson 积分均值

`[g(q−h)+4g(q)+g(q+h)]/6`。

两者均为 eV/Å。对三个正应变坐标又独立核对
`g=(σ_ii+P)V_0/L`；应力和压力均先换成 eV/Å³，`V_0/L` 为 Å²。
这一步与程序的变胞梯度逐点一致，排除了一个简单的压力单位或
正应变梯度变换错误。它不排除其它尚未发现的数值问题。

| 坐标组 | Simpson 后最大绝对残差 (eV/Å) |
| --- | ---: |
| 12 个原子坐标 | 0.000592 |
| 3 个改体积的正应变 | 0.021057 |
| 3 个一阶保体积剪切 | 0.0000196 |

三个正应变的残差分别为 `−0.020527`、`−0.020687`、
`−0.021057 eV/Å`；按各自 `dV/dq` 换算，能量差分暗示的压力
比应力积分高 `0.305`、`0.307`、`0.313 GPa`，即约
`3.05–3.13 kbar`。这是**数值一致性偏差**，不是统计误差条，
也超过此前接受的 `1 kbar` 压力尺度。

中心及正/负扰动的 FFT 网格均为粗网格 `24×24×40`、细网格
`48×48×80`；最大平面波数在正应变方向分别变化为
`1234→1242→1249`、`1237→1242→1251`、
`1242→1242→1258`（负扰动→中心→正扰动）。
[VASP 官方说明](https://vasp.at/wiki/Pulay_stress)指出，固定 `ENCUT`
下晶胞变化可改变有限平面波基组并造成能量—体积曲线不光滑。
因此基组切换/Pulay 效应是与这些数据相容的解释，**但仅凭平面波数变化
不能证明它就是本例残差的唯一根因**。相同 FFT 网格亦不能证明
完整的平面波基组相同。

## 原位弛豫与新静态计算不是同一机械输入

另一个已审计的 B4/B1 端点对比使用相同的最终结构字节、KPOINTS、
POTCAR 和所列电子 INCAR 参数，却有不同的机械设置：原位 `ISIF=3`
弛豫使用 `PSTRESS=457 kbar`，新静态使用 `PSTRESS=0`。两者的
`ENCUT` 都是 **600 eV**。VASP 的 `external pressure` 行等于打印的
`in kB` 原始应力对角均值减去 `PSTRESS`；原位 B4/B1 因而分别
打印约 `−0.14/−2.06 kbar`，新静态约 `458.72/453.89 kbar`。
相对于目标 **457 kbar** 的应力残差必须由原始张量计算，不能将
两类 `external pressure` 行直接相减或当作端点收敛判据。
[VASP 的 PSTRESS 说明](https://vasp.at/wiki/PSTRESS)明确指出其
压力项及应力对角修正；该输出恒等式也已由
[`audit_gan_600eV_basis_history.py`](../scripts/audit_gan_600eV_basis_history.py)
对原始 OUTCAR 逐一检查。

原位变胞弛豫保持原始平面波集合，而从同一末态重新启动的静态计算
重建集合；这是 [VASP 的 ISTART](https://vasp.at/wiki/ISTART) 与
[Pulay stress](https://vasp.at/wiki/Pulay_stress) 文档描述的区别。
本例两个末态的 384 个 K 点均发生平面波计数变化，静态减原位
`TOTEN` 为 B4 `−7.720`、B1 `+7.105 meV/四原子胞`。这些结果与
基组历史效应相容，但**不**单独证明能量差的唯一来源，也不构成
提高截断能或更换 K 网格的依据。完整数值、输入哈希和机械参数见
[`gan_600eV_basis_history_audit_20260928.json`](../benchmarks/numerical_integrity/gan_600eV_basis_history_audit_20260928.json)。

## 论文与下一步边界

600 eV 的原子—应变 Hessian、局域二维切片仍可如实作为固定协议的
局部诊断；不能据此把最高像认证为严格驻定的全变量胞 TS。
如需进一步判别，可以在**同一个 600 eV 协议**下检查更密的正应变
能量/应力采样和不同差分步长的稳定性，并预先规定接受门槛；
在现有数据未通过前，不得用 800/1000 eV 结果补入生产路径、
局域 TS 能量或论文同协议图，也不得把提高 `ENCUT` 当默认修复。
