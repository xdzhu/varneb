# BTO 条件双模面：四角与预留中心点审计（2026-09-28）

本记录只合并**已经完成**的五个固定 `(Q_z,Q_x)` 条件点；没有新提交
DFT，也没有改变 ABACUS 3.10.0 LTS / PBE / 100 Ry / Ba–Ti–O 10 au
DZP / 4×4×4 电子 k 网格。模式来自五原子 `1×1×1` 立方胞的 Γ
力常数。每个点开放第三条 Γ 软模、其余原子方向和六个对称应变，
去除整体平移。能量统一减去同一个立方 C 参考，单位为 eV/BTO。
这里不是 T→C 的 VCNEB 焓路径或势垒。

| 固定 `(Q_z,Q_x)` (√amu·Å) | 所选 `E−E_C` (eV/BTO) | `|Q_y|` (√amu·Å) | 正交梯度 (eV/(√amu·Å)) | 最大应力 (kbar) | 最低正交本征值：步长 0.05 / 0.10 (eV/(amu·Å²)) |
| --- | ---: | ---: | ---: | ---: | ---: |
| `(0.60,0.00)` | −0.102198868588 | 1.106952603 | 0.00114220 | 0.53654 | 0.00534909 / 0.00537664 |
| `(0.60,0.30)` | −0.104965260420 | 1.070583787 | 0.00240128 | 0.49161 | 0.00538550 / 0.00541597 |
| `(0.90,0.00)` | −0.108553741792 | 0.966913679 | 0.00199022 | 0.98189 | 0.00569877 / 0.00572021 |
| `(0.90,0.30)` | −0.109929070434 | 0.920822449 | 0.00174174 | 0.97648 | 0.00581667 / 0.00584450 |
| **预先指定留点** `(0.75,0.15)` | −0.106887865039 | 1.032331062 | 0.00150764 | 0.77176 | 0.00543824 / 0.00546533 |

所有五点达到本次条件优化的正交梯度 `0.003 eV/(√amu·Å)` 与原始
应力 `1 kbar` 门槛；它们**不是**普通 NEB 的 `0.10 eV/Å`
收敛判据。每点三起点中的 `Q_y=0` 驻定分支比 `±Q_y` 低能
分支高；`±Q_y` 在各点近简并。每点两档完整 16 维正交 Hessian
未见负本征值，最软混合方向另有四个原始静态点作能量／梯度交叉
检查。但部分 Hessian 对角的能量–力差分不一致量与最小本征值
同量级；最大对角差并非严格本征值误差上界。因此只能称**局域
稳定候选**，不能给出严格条件极小值或全域下包络证书。

## 立方端点对“二维条件下包络”的限制

把固定坐标取为立方 Γ 软模的 `(Q_z,Q_x)`、再对遗漏的 `Q_y`、其余原子
坐标和应变全部最小化，定义的是一个**下包络**，而不是自动包含真实
T→C 两端点的能量面。已归档的立方 C 端 Γ 本征频率前三项均为
`−8.9624663 THz`（虚频，约 `−298.956 cm⁻¹`）。按 Ti−Ba 的 z、x、y
位移锚定，这三方向在不稳定子空间中的投影满秩；`Q_y` 在固定 z/x
平面外的分量比例为 `1.0`。七像路径的 C 端投影为 `(0,0)`，但在
该点保持其他坐标不动、沿自由 `Q_y` 作足够小的扰动就有负二阶
能量变化。因此 C **不可能**是释放 `Q_y` 后的局域条件极小值；
完全最小化的二维下包络在 `(0,0)` 的能量必低于 `E_C`。

四个离端点的实算条件点独立显示同一分支风险：`±Q_y` 低支比
分别弛豫的 `Q_y=0` 驻定支低 `28.308–57.467 meV/BTO`。这不是
在 C 点实测的条件能差，也不能外推为整片面的分支差。它说明
**继续加密当前下包络不能解决端点定义问题**。若要在二维图上
同时保留 T→C 的 C 端，应显式固定 `Q_y=0` 对称支，再释放其余
坐标/应变并报告其横向不稳定性；或者保留三个软模坐标作三维面，
或采用以真实路径为中心的局部 `(s,q_\perp)` 切片。三者回答的
物理问题不同，不得把其中一种改名为完整“条件双模 PES”。

该判断可从已跟踪的 Γ 本征对、力常数、四点分支 CSV 和七像投影 CSV
离线复算。在仓库根目录运行（输出文件须选一个尚不存在的路径）：

```text
python -m scripts.audit_bto_soft_triplet_conditional_endpoint --eigenpairs outputs/batio3_t_to_c_pbe100_dzp10au/bto_cubic_phonopy_gamma_eigenpairs.npz --force-constants outputs/batio3_t_to_c_pbe100_dzp10au/bto_cubic_gamma_force_constants.npz --provenance outputs/batio3_t_to_c_pbe100_dzp10au/bto_cubic_phonopy_gamma_eigenpairs_provenance.json --force-sets outputs/batio3_t_to_c_pbe100_dzp10au/bto_cubic_gamma_FORCE_SETS --branches paper/VARNEB_CPC/figures/bto_conditional_four_point_stage_2026-09-27_source_data.csv --path-projection paper/VARNEB_CPC/figures/bto_frozen_soft_mode_landscape_source_data.csv --output tmp/bto_soft_triplet_endpoint_check.json
```

已归档输出为
[`bto_soft_triplet_conditional_endpoint_20260928.json`](../benchmarks/numerical_integrity/bto_soft_triplet_conditional_endpoint_20260928.json)。
此结论使我们暂不提交无判别力的致密二维下包络作业；BTO 的
ABACUS/PBE/100 Ry/10 au DZP、Γ 声子 `1×1×1` 和电子 `4×4×4` 契约
均未更改。

## 独立留点与分支身份

四角双线性预测的中心能量为 `−0.106411735308 eV/BTO`；预留中心
点的独立 ABACUS 值比预测低 `0.476130 meV/BTO`。这只是在
**一个已测单元的一个留点**处的误差，不是整个二维面的误差上界。

`Q_z=0.9` 两角原始最低能选择了 `−Q_y`，其正分支也从独立起点
收敛到近乎相同的能量。对原子 `y` 位移和 Voigt `yz,xy` 剪应变
执行镜面变换后，正、负分支坐标最大差分别只有
`1.68×10⁻¹¹`、`1.99×10⁻¹¹`（前 15 项为 Å，后六项为无量纲
应变）；能量差分别为 `5.46×10⁻¹²`、`9.09×10⁻¹² eV/BTO`。
这是当前点位上镜面对称分支的证据，不应把原始选中符号的改变
误判成不同物理相。比较结构时统一取已实际计算、已缓存审计的
`+Q_y` 分支；并未修改 DFT 结果。

通用 `vcneb.conditional_evidence` 筛选器仍忠实检查调用方给出的
原始结构，不会自动把任意大跳跃宣布为对称等价。本案例的对齐是
**有证据的显式选择**：先核对两套实际计算的终态在镜面下坐标和
能量均一致，再在这个有限网格上统一使用正分支代表元。未来材料
若没有对应对称操作与独立分支证据，不能套用这一处理。

在同一质量加权原子＋已声明应变度规中，正分支四角预测中心
`Q_y=1.016318130 √amu·Å`，实测为 `1.032331062`，差
`0.016012932`。完整 21 维坐标的留点误差为
`0.1254203 √amu·Å`，其中原子项范数 `0.0204857`、应变项
范数 `0.1237360`；这两个分项平方相加后开方才是总误差，
**不是**原子位移各差 `0.125 Å`。用同一正分支代表，四条相邻边
扣除预设 Q 步长后的离面度规跃迁分别是
`1.44167, 1.47337, 0.37687, 0.39851 √amu·Å`；直接混用
两角的原始 `−Q_y` 选择，会把前两条夸大为 `2.67087,
2.62008`。大跃迁本身也可能是连续但强烈的应变响应，需独立
门槛和更多点验证，不能仅凭数值断言不连续。

`vcneb.conditional_evidence.screen_conditional_interpolation` 要求
每点曲率不确定度、分支连续性阈值、每单元独立留点的能量和
**完整结构**误差阈值都事先明确。当前只有一个单元及一个留点，
材料级容差未在看到这五点之前完整预注册；也没有严格的最小
本征值误差界。故**不宣布门已通过，不画获认证的平滑条件 PES**。
论文中的现有冻结双模切片必须继续标为冻结切片；这五点可以
作为方法与分支问题的有限区域实证，但不能替代整个条件面。

## 来源与复核

无需访问 hf 就能重算上述五点算术和镜面核验：
[`bto_conditional_five_point_patch_2026-09-28.json`](../benchmarks/numerical_integrity/bto_conditional_five_point_patch_2026-09-28.json)
保存五个点的坐标、能量、代表元及原始结果文件哈希；执行
`python scripts/audit_bto_conditional_five_point_patch.py benchmarks/numerical_integrity/bto_conditional_five_point_patch_2026-09-28.json`
会给出能量与结构留点误差、对齐前后邻点跃迁。该紧凑快照**不能**
替代下面列出的 hf 原始 SCF、力、应力逐点审计；脚本也不运行 DFT。

- `(0.60,0.00)`：本地忽略目录 `outputs/batio3_t_to_c_pbe100_dzp10au/`
  的 `bto_transverse_soft_conditional_q060_result_job27783331.json`
  （SHA-256 `88cd18884802278422f2d323a99659059f5eb9ae701595b55b00e0fb1b7fbe8a`）；
  原始点、两步长及分支审计见
  `docs/VARNEB_BTO_Q1Q2_CONDITIONAL_PILOT_2026-09-25.md`。
- `(0.60,0.30)`：hf
  `/public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-q060q030/run-q060-q030/conditional_q060_q030_result.json`
  （SHA-256 `176d8ec47e79c4e99972efb92402c827963c18df836771e82b7045f99817b2b8`）；
  小型原始审计摘要在
  `benchmarks/numerical_integrity/bto_q060_q030_three_branch_evidence_2026-09-27.json`。
- `(0.90,0.00)`、`(0.90,0.30)` 与留点：hf
  `/public/home/iai806/abacus/agent-runs/20260927-varneb-bto-soft-gridpilot-qx030/`。
  五个相关数组/作业 `27787086, 27787309, 27787447, 27787714,
  27787746` 以及留点 `27787471, 27787487, 27791225` 在
  `sacct -X` 均为 `COMPLETED 0:0`；这仅是运行状态，点位结论
  另据原始 SCF/能量/力/应力审计。该目录下三份累计原始点审计
  `audit-q090_q000-all-points-through-soft-probes-27787447.json`、
  `audit-q090_q030-all-points-through-soft-probes-27787447.json`、
  `audit-q075_q015-all-points-through-soft-probes-27791225.json`
  的 SHA-256 依次为 `0106e76fa2327a47b712f5c600cfc27ded71fc2787feebf9579d7ba930a4464e`、
  `0456db088fa165d41cd662515c002c2459cb5fb2d87f2ca6f64b9d1fa57248c0`、
  `8003e7d677e8c6658295547c5dbf6fde7b82a7aacf7b94c62226ce293d201264`。
  它们的静态 `INPUT`、`KPT` SHA-256 相同，分别为
  `298af81a1f3822905f23fcdadff3e8f3b81b2fdc62b66e8716a922765443736c`、
  `1ad281c0f88e2f4269b222db91fc59e3365188569a96e8629b57c8f65b0613e7`。
  `replay-q090_q000-final.json`、`replay-q090_q030-final.json`、
  `audit-q075_q015-branch-replay-v2-27787487.json` 是只读缓存分支
  重放，不是再次 DFT；两档曲率与定向四点报告仍保留在同目录。

本记录使用 `abacus-agent-skill` 的端点／路径参数一致性、原始输出
优先和静态点不冒充相变势垒的审计规则；未使用该 skill 的 235
运行模板，因为本项目明确限制在 hf `hfacnormal01`。
