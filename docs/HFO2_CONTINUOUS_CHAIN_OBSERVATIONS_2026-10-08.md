# HfO₂ 完整链观察：连续模式坐标与受力归因

这是 G1 迭代中的证据，不是最终 MEP、TS 认证或应变选择性结论。
原始生产源仍为 `cc536a2`，R12 源与原始 SCF 未修改；分析独立运行。
全部沿用 ABACUS/PBE/100 Ry/完整 10-au DZP/Gamma 2³/P=0。

## 1. 本次实际取得的证据

`benchmarks/hfo2_channels/20261008/chain_observations/` 收录七组完整快照：
T→PO 漏峰修复链的 step 0/3/6，PO→M 的 step 0/4/8/10，共 66 个 image 评估记录。
不是 66 次新增 SCF；全部复用原计算，分析 DFT 调用数为零。
每像保留有序几何、原 SCF 路径、INPUT/KPT/PP/轨道/STRU/原日志 SHA、完整 E/F/stress。
逐一匹配真实完成的 32-MPI SCF 后重放普通 VCNEB，峰值能量和 fmax 与生产日志吻合。
本地可直接从 SinglePointCalculator 轨迹复核，不依赖 hf 或专有赝势文件。

| 链/完整步 | fmax (eV/Å) | 正向离散垒 (meV/f.u.) | 反向离散垒 (meV/f.u.) | 最大残差所在像 |
|---|---:|---:|---:|---:|
| T→PO / 0 | 0.829145 | 44.659 | 125.980 | 1 |
| T→PO / 3 | 0.215938 | 39.619 | 120.941 | 3 |
| T→PO / 6 | 0.344031 | 36.560 | 117.881 | 1 |
| PO→M / 0 | 0.881649 | 104.915 | 176.300 | 7 |
| PO→M / 4 | 0.362648 | 99.708 | 171.093 | 5 |
| PO→M / 8 | 0.319968 | 96.632 | 168.017 | 5 |
| PO→M / 10 | 0.293325 | 95.199 | 166.584 | 5 |

这些均未达到普通 0.10 阈值。正/反垒之差分别保持端点差
−81.321193 和 −71.384600 meV/f.u.，不能换用 M/T 能量零点去讨论同 PO 初态的竞争。
离散最高像均为 3；最高像的负切向曲率也不能代替全变量驻点与 Hessian index。

## 2. 回弹是什么，不是什么

T→PO 的 step3→6 回弹来自最大残差从像3转到像1。
step6 该像原子残差为 0.344031，晶胞块为 0.030582，真实垂直力为
0.346935，弹簧最大矢量为 0.023657 eV/Å。因而不能把此次回弹主要归咎于
弹簧太强或纯胞应力，也没有依据调整 cutoff/SCF 参数或立即停链。
step7 的生产日志已回落到 0.300872；继续看完整健康段。

PO→M step10 的最大残差在像5：原子块 0.293325、晶胞块 0.086426、
真实垂直力 0.295873、弹簧 0.004319 eV/Å，主要仍是尚未松弛的原子方向。
18:45 Slurm 28275259 为 COMPLETED，但材料摘要为 `max_steps_reached/step_limit`；
正常作业退出不是路径收敛。28274895 仍运行，两个翻转任务仍等待依赖。

## 3. 单一初始周期规范，不做逐像折叠

新通用接口 `continuous_reference_coordinates(images, reference)` 先拒绝不连续 lift，
仅在 image0 相对参考选择每原子的整数格矢代表元，并对所有后续 images 应用同一表：

`u_i = (q_i + n − q_ref) H_ref`。

不改变结构、原子映射、DFT 输入或运行路径，不逐像 MIC 到参考态；否则真实绕行可能
被图中的幅度跳变或伪平滑掩盖。完整 `F=solve(H_ref,H_i).T` 与
Green 应变 `(F.T F−I)/2` 独立记录，不把胞旋转强行变成六个对称应变。
这批几何均在声明的 quarter-site 局部范围内；这不是对更长路径的保证。

`project_reference_basis` 显式选 Cartesian 或质量度量，先去除相应刚性平移。
Γ 投影保留原 Phonopy O=15.9994 amu，不替换成 ASE 默认 15.999；未重算频率。
声学子空间以真实平移重叠识别。简并本征矢的单独坐标依赖选基，完整子空间权重才稳定。

## 4. 对后续技术路线的实质影响

T→PO step6 最高像在三个旋转 T 几何模式中的位移平方范数占比为 **93.91%**；
PO→M step10 最高像仅 **52.84%**。全路径最小占比约 31.43%，且端点固定，
继续优化不会让这些固定端点自动落进三模平面。投影占比不是能量贡献。

因此保留三模作为变体追踪，但否决“所有竞争通道都可由同一三模二维面解释”的做法。
G3 不马上启动粗暴三模铺网格；先等待通道优化，再检查瓶颈残差和联合原子—胞稳定性。
若采用局部路径切线/横向方向或增维，需要明确它们并非母相声子模，且对未参与选基的
独立采样点验证预测。T Γ 全基可以精确重构原子位移，只证明表示完整，不证明局部
谐近似、条件稳定性、能量预测或 TS 认证。冻结与释放的区别仍需真实 DFT 对照。

这构成一个可复现的“低维表示何时失效”结果，尚不是 JCTC 核心创新闭环。
仍须同 PO 初态的翻转/退相变竞争、同机械边界矩阵、独立预测与消融。

## 5. 复核与验证

```powershell
python -m scripts.analyze_hfo2_chain_observations `
  --root benchmarks/hfo2_channels/20261008/chain_observations `
  --variants benchmarks/hfo2_channels/20261008/reference_variants `
  --gamma benchmarks/hfo2_channels/20261008/gamma_analysis/T_d0.01.npz `
  --output observation_replay.json
python -m pytest -q tests/test_continuous_projection.py tests/test_hfo2_chain_observation.py tests/test_hfo2_observation_analysis.py tests/test_hfo2_lifted_restart.py
```

当前 29 项聚焦回归通过，含真实 hf 数据重放、改动日志/几何/数值证据拒绝、
绕行保持、单一整数规范不变性、斜胞有限应变、质量度量与简并子空间检验。
全回归及 HF 兼容结果另写进本目标进展文件；不以分析脚本测试冒称材料链收敛。

最终干净归档全回归：791 passed、2 skipped，69.10 s；HF ASE3.23.1b1
复核66个记录的力/模式/应变一致，最大差为浮点舍入量级。见同目录
`hf_replay_check.json`。记录并修复了遗漏忽略轨迹和跨平台目录排序两处交付问题，
没有放松审计门禁或改材料参数。19:10时gap最新step9为0.233694，
预计约19:25结束10步健康段；全链收敛仍待实际证明。

## 6. 夹持 G2 链的独立重放入口（2026-10-10）

本页第1–5节为历史**自由胞**G1观察，不能直接用其默认重放器审核夹持G2。
新入口明确读取生产 `vcneb_preflight.json`、固定基底文件及其SHA，使用同一个
`clamped_plane_vcneb_boundary`、3个开放晶胞方向与登记尺度5.12968067458423Å。
原子和晶胞的NEB力都在此子空间内重放，非零夹持反力不当作自由胞残差。
首版只接受已登记首段Slurm脚本的精确SHA，不能把任意运行套进固定弹簧/度量。

```bash
python -m scripts.export_hfo2_clamped_observation \
  --workdir /path/to/completed/G2/band --step 10 \
  --source-job-id ACTUAL_SLURM_ID --output /path/to/fresh/observation \
  --production-script /path/to/immutable/source/cluster/hf_hfo2_clamped_chain_pilot_20261010.slurm
python -m scripts.prepare_hfo2_clamped_resume \
  --observation /path/to/fresh/observation --output /path/to/fresh/resume-seed
```

每像须匹配同原子序的精确周期几何、原六份物理文件与STRU、原始SCF日志SHA；
端点缓存另与生产预检中登记的原端点记录相核对。内部像只接受完成的
`call_audit.json`，完整E/F/stress重新从32-MPI原日志解析。拒绝缺像、缺日志、
哈希变化、陈旧几何、断裂lift及重放与优化日志不一致；不修改生产目录或调用DFT。
离散垒、普通残差和完整物理梯度分别留档，未收敛链不称MEP或认证TS。

几何续算使用独立 `make_clamped_resume_cached_factory`：九个像均需精确缓存，
当前完整帧不重算，内部像移动后正常失效再计算，固定端点仍只读取原缓存。
它不是FIRE动力学状态恢复；新段从同一几何初始化新的FIRE速度/时间步，记录此限制。
新段用新目录和经测试的独立源码归档，绝不覆盖原运行源码、轨迹或SCF。
健康步数上限允许审核后续算，计为原独立链，不追加+0.5%留出条件或新通道。
