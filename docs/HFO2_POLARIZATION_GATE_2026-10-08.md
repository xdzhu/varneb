# HfO₂ 翻转端点的电子极化验证

## 有限预算与实际版本

两个按母相映射构造的 PO− 已通过同参数 E/F/stress 检查，但核位移反号不替代电子极化。
G1 翻转链前，独立验证 PO+ 与两个 PO−。仅 **3 个 SCF + 9 个固定电荷 NSCF**；
先提交 PO+ 一个小任务，实测通过后再提交另两个。不生成 Born 矩阵、不重优化端点。

合肥二进制为 ABACUS v3.10.0、commit `f7cb1d3`。已读取该提交的
[LCAO 官方案例](https://github.com/deepmodeling/abacus-develop/tree/f7cb1d3/examples/berryphase/lcao_PbTiO3)
和 [Berry 实现](https://github.com/deepmodeling/abacus-develop/blob/f7cb1d3/source/module_io/berryphase.cpp)。
[官方说明](https://abacus.deepmodeling.com/en/latest/advanced/elec_properties/Berry_phase.html)支持 LCAO。
案例的材料、50 Ry 和 k 点数仅是协议证据，**不照搬物理参数**。

## 固定契约与可见差异

| 环节 | 保持不变 | 声明的差异 |
|---|---|---|
| SCF | 原始 STRU、PBE、100 Ry、Hf/O 赝势、完整10au DZP、2×2×2、收敛与混合参数 | 追加 `out_chg 1`、`out_bandgap 1` 输出开关 |
| Berry NSCF | 同一几何、同一 SCF 电荷、同一泛函/截断/赝势/轨道 | NSCF/读电荷/禁对称/Berry/R3；关闭未使用的力/应力/结构输出；积分网格2×2×2、2×2×4、2×2×8 |

NSCF 是固定电荷上新观测量的积分检查，不是重新定义路径能量的 SCF 参数。
其任何能量都不进入能垒。NEB 固定物理哈希 guard 不放宽；原计算器文件不改。
SCF 输出追加前后须复现：E ≤1e−5 eV/cell、F ≤1e−4 eV/Å、stress ≤0.02 kbar。
这是复现阈值，不重新收紧端点的力/2 kbar 标准。失败先查原因，不改物理参数“修复”。

## 事前验收与分支边界

1. 真实 DSIZE32；SCF 电子收敛；NSCF 完整结束；无 singleton。
2. 从 UPF 校验 Hf/O 价电子12/6：Hf₄O₈共96电子、非磁占据48带。
   `istate.info` 占据列含 k 权重，不把0.25误读为分数占据；各采样间接带隙 >0.1 eV。
3. 同电荷的2×2×4→2×2×8极化模差 ≤0.01 C/m²。只验证纵向积分，不宣称全k收敛。
4. 两个 PO− 与 PO+ 反号关系在原生日志给出的模数下残差 ≤0.01 C/m²。
   若 PO+ 属于自身反号不变的0/半模数类，则此检验不充分，标记 review_required。
5. 该版本非磁、离子价数皆偶数时打印约 `2e|a₃|/V` 的模数。记录物理量子
   `e aᵢ/V` 与原生模数并核验倍率；**不能把输出 P 除以2**。
   源码所用较旧 SI 常量的微小单位差异也保留，不冒称 DFT 数值误差。

这是同胞 R3 极化类验证，不给绝对自发极化，不自动选翻转路径分支，不是完整三方向
极化或全布里渊区绝缘性证明。后续连续绝缘路径须另追踪晶胞、量子与分支；
“最小绝对P”不能替代物理分支选择。

## 实现与运行

- `vcneb/polarization.py`：量子晶格、原生输出解析、显式模差、采样带隙。
- `scripts/hfo2_endpoint_polarization.py`：固定源/输出差异、实际审计、极化类验收。
- `cluster/hf_hfo2_endpoint_polarization_20261008.slurm`：hfacnormal01，1节点32 MPI，45分钟/端点，最多并发2。
- prepare 后方可在 allocation 内显式 RUN_DFT=1；任何已有工作目录都拒绝覆盖。
- 本机上游 ASE 无 ABACUS IO，夹具测试不冒充实际接口；提交前用 hf 自带 ASE
  真正读取源 STRU 预检，不另装软件。

成功只解锁两个翻转端点的电子解释；失效保留原始输出，按证据修复协议或缩小结论，
不以无止境提高计算精度兜底。
