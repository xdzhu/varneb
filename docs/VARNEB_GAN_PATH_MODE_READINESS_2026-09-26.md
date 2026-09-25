# GaN B4→B1：路径坐标已审计，声子模式和变胞 TS 尚待认证

本阶段仅分析已收敛的 VASP 四方路径，不启动新的 DFT。原路径在
45.7 GPa、29 总像、2 个 GaN/计算胞下，最高像 15 的正向焓垒为
0.338491 eV/GaN。源轨迹、端点身份、原始焓力与 TS 证据边界见
[`VARNEB_GAN_TS_CANDIDATE_AUDIT_2026-09-26.md`](VARNEB_GAN_TS_CANDIDATE_AUDIT_2026-09-26.md)。

离线分析器 [`audit_gan_path_adapted_atomic_axes.py`](../scripts/audit_gan_path_adapted_atomic_axes.py)
从已归档轨迹抽取最后一条完整链，对齐原子顺序与周期规范，以初始
B4 胞作为参考，逐像移除原子整体平移；把质量加权的原子位移做 SVD，
并**单独**报告相对于初始胞的六个对称应变分量。路径弧长仍按原
VCNEB 原子–胞度量计算。完整机器审计（含基矢、逐像坐标和输入
SHA-256）在
[`gan_tetragonal_empirical_atomic_axes_audit_20260926.json`](../outputs/neb_literature_benchmarks/gan_tetragonal_empirical_atomic_axes_audit_20260926.json)，
逐像公开数据在
[`gan_b4_b1_tetragonal_empirical_atomic_strain_20260926.csv`](../benchmarks/numerical_integrity/gan_b4_b1_tetragonal_empirical_atomic_strain_20260926.csv)。
报告里的 `empirical_Q1–Q3` 单位为 √amu·Å；其符号由主载荷为正
固定，只有与同一基矢文件搭配才可解释。

| 检验 | 当前结果 | 证据边界 |
|---|---:|---|
| 前 1 / 2 / 3 个经验轴所捕获的路径原子位移平方范数 | 96.635% / 99.246% / 99.999999954% | **同链拟合**，描述性压缩 |
| 3 轴最大逐像留一残差 | 0.000340 √amu·Å | 插值稳定性检查，不是外部路径验证 |
| 末态 B1 改作参考时的 3 轴最大残差 | 0.000183 √amu·Å | 参考规范的对照，不能证明声子身份 |
| 相对变形梯度最大反对称分量 | 8.87×10⁻⁷ | 数值级非对称；以对称应变报告误差最大 2.19×10⁻⁶ Å |
| 端点原子身份 | 通过 | 4 原子 Ga₂N₂，顺序、周期 gauge 与静态端点一致 |

三维经验子空间能几乎重构*这条已采样路径*，但这是从路径本身
训练的轴，**不是 Γ 声子本征矢、软模或过渡态不稳定方向**；不能从
三个奇异值推断存在三种物理驱动模式，也不能把 `Q` 的投影幅度
换算为“模式能量贡献”。胞应变是独立坐标，不能偷偷并入原子 SVD。
局部最高像仍只是 TS 候选：普通 VCNEB 残差收敛不等于全变量胞
焓梯度为零，更不等于 Hessian 恰有一个负特征值。

下一材料级门槛分两条独立验证：首先，在**相同 VASP 计算契约**下
得到端点 Γ 力常数，检验力差分线性、声学和规则、数值噪声，再把
路径投影到真正的端点声子基；固定胞 Γ 结果仅覆盖原子子空间。
其次，以最高像 15 为种子做局部全变量胞驻点精化，并在原子位移
与对称应变空间认证 Hessian 指数和两侧下坡去向。这两条证据齐全前，
论文只宜表述为“收敛的有垒路径和候选最高像”。

重生成命令（需要仓库外已归档的大轨迹和摘要；脚本拒绝覆盖输出）：

```powershell
python -m scripts.audit_gan_path_adapted_atomic_axes `
  --summary validation/gan_b4_b1/gan_b4_b1_hf_vcneb_summary.json `
  --trajectory validation/gan_b4_b1/gan_b4_b1_hf_vcneb.traj `
  --generic-audit validation/gan_b4_b1/vcneb_audit.json `
  --route-analysis validation/gan_qian_suite/tetragonal_analysis_20260918.json `
  --ts-audit outputs/neb_literature_benchmarks/gan_tetragonal_ts_candidate_audit_20260926.json `
  --output NEW_AUDIT.json --csv NEW_DATA.csv
```
