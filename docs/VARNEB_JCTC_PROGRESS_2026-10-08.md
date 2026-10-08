# JCTC 研究接力：2026-10-08

目标 active；总研究目标尚未完成。现行路线为
`docs/VARNEB_JCTC_HFO2_RESEARCH_PLAN.md`，不要按旧规划恢复无限材料/后端矩阵。

## 已完成

- 科学问题收束为 HfO₂ 同初态的翻转/非极性逃逸通道选择性；形成通道另行解释。
- 明确最近文献已覆盖应变、三模耦合和翻转变体枚举；新增方法/预测仍需实证。
- 后端无关的稳定正交自由度曲率约化已实现，失稳/未分辨方向会拒绝消去；解析验证通过。
- 四份历史链摘要归一化和两条最终链几何审计完成。
- 普通/引导低垒链均保持同一有序初末态与最终周期绕行整数，中间几何不同；仍不贴未经核验的 X2− 标签。
- 两个峰值同参数真实 SCF 均复现：28237630_0、28238926_1。
- 一次 node26 Hydra bootstrap 失败保留在 28237630_1；未进入 ABACUS，通过单节点 allocation 内 fork 解决，不改物理参数。
- 本机完整工作树测试 696 passed / 1 skipped；用户旧 BTO 文档/图及未跟踪工作保留，没有纳入提交。
- 8 个新原子/应变探针全部完成并通过有效输入和逐点 SCF 审计；导数符号/量纲相容，保留有限步长差异。
- 研究/开发里程碑已正常推送 `76bff18`、`08f2b7b`、`8923165`，无 tag、CI 或 PyPI 发布。

## 当前计算与检查入口

数组 `28240257`，hf / hfacnormal01，每点 1 节点 32 MPI、30 分钟墙时、最多并发 3。
8 个局部成对探针：原子方向与 xx+yy 应变方向，±0.01/±0.02 Å 联合坐标步长。
100 Ry、完整 10-au DZP、Gamma 2³ 和其余电子参数全部不变。

源/作业根目录：

`/public/home/iai806/abacus/agent-runs/20261008-varneb-hfo2-channels-r3`

12:31:38 北京时间：8/8 `COMPLETED 0:0`，汇总重读逐点文件完成审计；各元素墙时 103–154 s。
此前预估 12:30–12:34，这批已完成；不是整项 JCTC 研究完成，也不是路径/TS 已认证。
先 `sacct -n -X -j 28240257 --format=JobID,State,ExitCode,Elapsed`，不重复提交。
只读/生成汇总：

```bash
cd /public/home/iai806/abacus/agent-runs/20261008-varneb-hfo2-channels-r3/source
PYTHONPATH=$PWD /public/home/iai806/.conda/envs/icu/bin/python -m scripts.score_hfo2_channel_work_probes \
  --root /public/home/iai806/abacus/agent-runs/20261008-varneb-hfo2-channels-r3/probes
```

源 archive 不覆盖，结果审计保留真实失败。输入/代码哈希与句柄在
`benchmarks/hfo2_channels/20261008/run_registry.json` 和 `work_probe_manifest.json`。

## 下一步（按证据，有限推进）

1. 八点已收齐：原子方向导数差 0.000573/0.001418 eV/Å；应变方向导数差 −0.00000434/−0.000484。最大能量/梯度曲率差 0.10225 eV/Å²；两探测方向曲率均正。只接受局部方向诊断，未检查混合曲率/完整正交稳定性。
2. 若一致性不分辨或失配，检查解析、坐标 Jacobian、SCF/基函数数值来源；不得调 ECUT 去追结论。
3. 若 G0 通过，登记母相对称畸变与局部 Γ 基的区别，现有 12 原子胞 Γ 只需 1×1×1；准备 M 与两类 PO− 同胞变体，验证电子极化分支。
4. 第一批仅三条新无约束路径（两翻转、PO→M），T↔PO 复用；所有新链明确 9 total = 7 interior+2 endpoints，普通 0.10，不提前 CI。
5. 同边界同 PO+ 初态比较有效后，才开展计划中的夹持矩阵、局部条件面和独立预测。

低垒候选点目前投影梯度不是零（原子 secant 方向约 −0.1011 eV/Å）；它是普通链的最高像，不是已经认证的 TS。
这不否定该普通路径满足 0.10 的 NEB 条件，因为 NEB 残差并不要求每像真实切向力为零；小垒和鞍点认证需另做诊断。
只有局部/全可动驻定与负曲率证据充分后才谈 TS，不锁定二维图中心、不强迫过渡结构充当端点。

## 仓库及论文边界

CPC 现稿未重写为已完成的 JCTC 稿；JCTC 叙事/claim–evidence 设计在
`paper/VARNEB_JCTC/RESEARCH_DESIGN.md`。没有独立预测时不写完成式结论。
每个有效代码/数据/论文里程碑显式 stage、测试后 commit/push；不要 add 全工作树，不 force push。
下次接力无需读取旧会话全文；上述计划、记录和真实 hf 输出即可。

## 接续里程碑：母相映射与 Γ 批次准备

- 本机实际软件为 ASE3.28.0 / pymatgen2025.10.7 / spglib2.7.0 / phonopy2.21.0，
  没有为了技能示例版本重装环境；相身份分析保留 spglib 弃用警告。
- 已保存两个明确反演操作、母相置换、生成 POSCAR 和结构敏感性审计。
  T/PO/M种子三档容差均为 SG137/29/14；两候选均 SG29。
- Γ 几何畸变反向，但电子 Berry 极化尚未验证。PO 在原 T 模向量上的投影
  约 5.7e-7 Å，因此不使用这个量去宣布“保持/翻转 X2−”。
- 计划新增19次静态：T胞两步长 Γ 位移16点 + 两反演候选 + M种子；不是新 NEB，
  不重复优化 T/PO。Γ1×1×1，全部电子文件保持原始哈希。
- 完整本机工作树回归 **708 passed / 1 skipped**（45.86 s），含解析谐振子
  数据驱动完整 Γ 汇总测试；这些测试不充当 DFT 物理结果。
- 作业必须从新的提交 archive 提交；准备/提交状态将追加到 run_registry。
