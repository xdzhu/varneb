# BTO 条件模式面数值门槛：平移等价性与网格试验（2026-09-26）

本记录只针对已审计固定 `Q=(0.6,0) sqrt(amu)·Å` 的同一 5 原子 BTO
结构，沿实空间 `z` 整体平移 `−0.02,−0.01,0,+0.01,+0.02 Å`。
这五个晶体结构在周期边界下物理等价；能量差测量数值原点敏感性，
不是模式势能面的物理起伏。**本项目 BTO/ABACUS 计算契约固定
`ecutwfc=100 Ry`，并使用按 100 Ry、10 au 生成的 DZP 轨道；
不得通过提高 `ecutwfc` 为这套轨道另立生产协议。**这项边界由
用户于 2026-09-26 再次明确。下面的 120 Ry 试算是未经事先确认的
越界诊断；保留原始记录仅为审计，不作为参数选择、误差定量外推或
论文材料结果的依据。生产案例始终是 PBE/100 Ry/10 au DZP/
`4×4×4` k 点。

## 100 Ry 基线

此前作业 `27779827` 的五点独立原始日志审计
`bto_q1q2_translation_probe_independent_raw_audit_v4_job27779827.json`
给出 `50×50×50` 的电荷/势 FFT 网格与
`0.0214675133 meV/BTO` 的能量跨度。净平移力约为零，
故同一等价路径上仍有不应被解释为物理势垒的离散能量变化。
100 Ry 中心点的 `scf_thr=1e−8` 与 `1e−10` 同几何复算并未改变
这一诊断的结论。此处的平移跨度**不是任意模式方向的严格误差上界**，
但已经足以提醒我们：小于同一尺度的曲率能量差不能未经额外检查
就用于强结论。

## 只改 `ecutrho` 的失败试验

按 [ABACUS 官方 INPUT 说明](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html)，
`ecutrho` 控制电荷/势的截止能，默认 `4×ecutwfc`；因此曾保持
`ecutwfc=100 Ry`、赝势、轨道、k 点与 SCF 设置不变，仅显式设为
`ecutrho=800 Ry`。独立目录、预检与源码哈希守门后，Slurm 作业
`27781493` 在 `hfacnormal01/node11`、32 MPI 下运行 45 秒后
`FAILED 1:0`。实际 `INPUT` 确实含 `ecutrho 800`，MPI 日志也显示
32 个进程，但在首个 SCF 开始前 ABACUS 3.10.0 LTS 触发
`ModuleBase::matrix::operator(): Assertion ic<nc`；`addr2line` 将调用链
定位到 `Potential::update_from_charge` 和 LCAO `before_scf`。
该次日志报告的电荷网格仍为 `50×50×50`，没有最终总能、力、应力。
所以这**不是 800 Ry 的物理计算结果**，也不能凭这一例断言
`ecutrho` 参数普遍失效；在该二进制/输入组合上须先调查越界根因，
不能直接重投或把它用于插值。原始失败材料保存在 hf 的
`/public/home/iai806/abacus/agent-runs/20260926-varneb-bto-ecutrho800-translation-z`
及本地 `outputs/batio3_t_to_c_pbe100_dzp10au/`
`bto_q1q2_ecutrho800_failed_job27781493_raw/`；
其中预检 JSON SHA-256 为
`88e4f5713c864450734e0a424ddade43024fc76880b28978a7b584a4e9c4a635`，
原始 `running_scf.log` SHA-256 为
`f3f0579e202f3d0fd46b5d860bcd77e24251f92e09dc7adc1b500f995e9b8ea6`。

## 已撤回作为决策依据的 `ecutwfc=120 Ry` 越界试验

ABACUS [官方文档](https://abacus.deepmodeling.com/en/latest/advanced/input_files/input-main.html)
也指出 LCAO 仍用 `ecutwfc` 处理局域赝势等平面波表示；提高它
不只是改变显示网格，因此本试验是**截止能/网格联合敏感性**，
不是纯网格因果隔离。作业 `27781497` 在同一固定 Q 中心实算成功，
电荷网格增大至 `60×60×54`，原始 SCF/最终能量/力/应力均完整。
但提交时的 Python 校验器错误地复用了写死 `100 Ry` 的生产样例
校验函数，于 DFT **完成之后**拒绝该点，使 Slurm 显示
`FAILED 1:0`。独立审计 `bto_q1q2_ecutwfc120_center_raw_audit_v2_job27781497.json`
核对了提交时源码、真实 `INPUT/KPT/STRU`、原始日志、`DSIZE=32`、
无熵项及能量/受力/应力；随后本地诊断代码改为按“基线 INPUT
仅允许 `ecutwfc` 一项不同”验证，避免再把生产参数当作全局硬编码。
**没有重算这个中心。**

以该已审计中心为基准，作业 `27781503` 只算四个平移点，
在 `node11` 用时 `3:38`，四次 ABACUS 调用均正常完成。
独立原始日志审计
`bto_q1q2_ecutwfc120_translation_independent_raw_audit_job27781503.json`
核对所有点的结构、k 点、赝势/轨道、活跃 INPUT 的唯一差异、
完整 SCF、最终能量/力/应力、共同的 `60×60×54` 网格以及原始数据
与评估器缓存一致性。五点能量跨度为
`0.0205936285 meV/BTO`，相对 100 Ry 的
`0.0214675133 meV/BTO` 只降低约 `4.07%`；两个 Simpson
能量减受力功残差为 `−0.00553` 与 `−0.01100 meV/BTO`，
没有一致地改善。更重要的是：这与本项目 100 Ry 轨道/计算协议
不一致，因此**不能**用其 4.07% 差值选择 BTO 生产 cutoff，
不能把两套能量混入同一个条件势能面，也不能把它当作同协议的
数值误差评估。先前从它推断“提高 cutoff 能否解决曲率问题”的
参数探索方向在本项目中撤回。生产计算参数未改动，后续也不以
改变 `ecutwfc` 或 `ecutrho` 作为这条 BTO 研究线的默认修复路径。
五点两套原始能量的精简公共数据表为
`benchmarks/numerical_integrity/bto_q060_rigid_z_ecutwfc100_120.csv`；
120 Ry 审计 JSON 的 SHA-256 是
`8f2d1a342dea68313ac7c5a07d6bcbcafdb3665f4870a7b7e6c127374559ddf8`，
100 Ry 审计 JSON 的 SHA-256 是
`f5aa72c9073f5e147a180bb045577588455e482abbe8dee3efdc10268df0fa82`。
公共 CSV 是原始日志审计的压缩摘录，不能代替原始计算归档。

## 对下一步的约束

1. 固定 `ecutwfc=100 Ry` 与对应的 10 au DZP 轨道，不再提交
   120 Ry 或其他 cutoff 改动的 BTO 作业；不为追求亚
   `0.01 meV/BTO` 的漂亮数字重复整张二维 DFT 网格。现有 59 点
   冻结切片可作为探索性图，但不是已认证的条件最小值面。
2. 若要声称小正交曲率/局部最小值，须在实际采用的同一计算协议下
   给出数值误差与曲率信号的分离；现有最软正交 Hessian 本征值约
   `0.0022 eV/(amu·Å²)`，`h=0.1 sqrt(amu)·Å` 所对应的二次能量
   信号仅约 `0.011 meV`，小于上述平移跨度。这个比较是风险筛查，
   不是严格的误差传播或物理曲率反证。
3. 对十至百 meV 的稳健特征，可事先声明分辨率和独立留出误差，
   在结果达到足够尺度分离时推进粗粒度机制图；对小于数值分辨率
   的局部曲率，应保持“未认证”。优先做计算器无关的受约束路径/
   软模坐标定义和结果分级，不让一次失败的 ABACUS 调参阻滞整个
   VARNEB 研究目标。
