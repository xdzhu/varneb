# VARNEB 全程研究交接：供另一会话独立审查

交接日期：2026-10-22（北京时间）。研究结果截至 2026-10-21；本次仅整理文档与核验索引，未启动新 DFT、声子、NEB 或集群作业。

## 0. 交接目的、版本与禁止误读

这不是新的研究计划，也不是“论文已达 JCTC 水平”的证明。作者要求对自研究开展以来的工作做完整、可独立复核的总结，明确成果、失败、缺口和停滞原因。接手人可以否定原解释，不必沿用原代理设定的研究方向或期刊判断。

- 仓库：<https://github.com/xdzhu/varneb>；本机根目录：`D:\Work\Code\vasp_neb`。
- 本文科学证据基线：`3e6286fb5256a48c857cec3b3b1636f3f46981a2`，分支 `main`。本文及旧交接入口链接是基线之后的文档交接提交，提交号以 Git 历史为准。
- 本文路径均相对仓库根目录；Markdown 链接相对本文所在 `docs/`，便于 GitHub 和新克隆阅读。
- 本次开始时，本地 HEAD 与本地缓存的 `origin/main` 相同；交付时另以实际推送及远端分支核验为准。本次没有实时检查 hf 上账户内所有作业。
- **当前授权仅限交接整理和推送。不要因为历史计划中写过 active、待提交、next、remaining 就自动续交作业。**
- 科学证据的优先级：原生输出/归档及其身份对应 → 完整已评估链 → 从它们导出的分析 → 图表 → 论文叙述。不是“最新 MD 一律覆盖原始数据”。
- 不把另一会话的 HZO、掺杂、氧空位等工作并入本会话成果；目录里存在相应计划或日志，不表示本会话完成了这些研究。

### 0.1 本次必须保留的既有用户改动

整理开始前，下面五项已修改，均不属于本次交接，不应混入本次提交：

```text
docs/VARNEB_BTO_Q1Q2_CONDITIONAL_PILOT_2026-09-25.md
docs/VARNEB_HIGH_LEVEL_RESEARCH_ROADMAP.md
docs/theory.md
paper/VARNEB_CPC/figures/varneb_architecture.pdf
paper/VARNEB_CPC/figures/varneb_architecture.svg
```

工作区另有其他会话的未跟踪文件。本次不使用 `git add .`，不清理或重置它们。干净克隆中的内容与上述本地修改不应混同。

### 0.2 建议先读的五个入口

1. 本文第 1、7、9、10 节：总判断、HfO₂ 结果、失败原因、待审问题。
2. [当前 HfO₂ 中文阅读指南](../paper/VARNEB_JCTC/AUTHOR_READING_GUIDE.md)。
3. [HfO₂ 阅读 PDF](../paper/VARNEB_JCTC/reading/VARNEB_HfO2_manuscript.pdf)：现有材料研究，不是新预测方法已经成功的论文。
4. [最接近文献及本工作增量](../benchmarks/hfo2_channels/20261008/mechanism_synthesis_20261021/literature_comparison.md)。
5. [明确的投稿前判断](HFO2_PRESUBMISSION_VERDICT_2026-10-21.md)：当前不视为已通过 JCTC 科学验收；这是内部判断，不是编辑决定。

## 1. 总判断与结果状态标签

**软件和验证底座已经形成；BTO/GaN 的四套主要二维面已有实算数据；HfO₂ 得到了“内部原子容纳改变候选路径应变响应”的定量结果。但新低维预测优势、普适加速、新路径拓扑和普适机械选择规律都没有建立。**

持续停滞的首要原因是原代理没有及时把科学问题、否证条件和停止条件固定下来，诊断与局部模型验证不断扩张。不能主要归咎于超算资源，也不能用代码测试通过代替科学创新。

本文使用以下标签，可靠性仅在声明范围内成立：

| 标签 | 含义 |
| --- | --- |
| 已完成／受限可靠 | 有可定位证据支持所列有限结论；不自动外推 |
| 已完成／负结果 | 实际执行后不支持原假设，必须保留 |
| 未建立 | 有相关分析，但不足以支持更强主张 |
| 未完成 | 没有同等级完整交付，不能用计划、目录或 canary 代替 |
| 历史失败／已修复 | 运行或方法错误已有诊断和后续修正，原失败仍不变成有效物理结果 |

三种“完成”必须分开：**Slurm 正常结束、数值达到项目阈值、论文科学问题得到充分回答**，三者互不等价。

## 2. 全程研究脉络与成果归属

| 阶段 | 实际工作 | 取得的进展 | 未解决或被否定的部分 |
| --- | --- | --- | --- |
| 基础 VCNEB | 联合原子/晶胞自由度、端点与插值、力变换、重启 | 建立可运行的软件基线与解析检查 | 实现已有 NEB/VCNEB，不是重新提出前人理论 |
| 材料与多后端 | BTO、HfO₂、GaN、CdSe，扩展 ASE 计算器 | 五 DFT 后端 GaN 闭环，若干文献对照 | 不是所有后端×所有材料全矩阵；B3/CdSe 有负结果 |
| 收敛工程 | 逐 image、分块、尺度修正、分阶段 FIRE 和几何门控 | BTO/HfO₂ 的外部 ABACUS 启动次数减少；长链续算救援 | 不能证明普适加速或相同鞍点盆地 |
| 模式与二维面 | Γ 模、投影、受限释放、联合曲率、真实 DFT 网格 | BTO/GaN 四套主面及图表、插值检验 | 冻结面不等于条件面，投影不等于能量分解 |
| HfO₂ 通道竞争 | 四通道、部分夹持、0/+1% 训练、+0.5% 独立检验 | 同初态比较及路径/模式/应变联系 | 简单共同模式空间及新模型优势未获支持 |
| HfO₂ 最终机制 | 匹配释放、更低完整链、三应变、两因素交叉 | 内部容纳使保持型应变增量 28.596→8.533 meV/HfO₂ | 原创性仍可能只是已知类路径的更充分松弛 |
| 稿件收束 | CPC 软件稿、HfO₂ 材料稿、主图与离线重放 | 有可评阅的完整材料稿，证据入口更清楚 | 不能据交付完成宣布 JCTC 科学目标实现 |

VASP/VCA 的 BST50/PTO/PZT50 实现与案例是作者从另一会话带入、随后梳理整合的成果，不应全部归功于本会话独立开发。此前提供的 ZStar 文稿仅用作格式参考，不是本工作的科学结果来源。

## 3. 软件：实现位置、成熟度与边界

### 3.1 核心代码索引

| 功能 | 关键代码 | 已实现的范围／限制 |
| --- | --- | --- |
| VCNEB 核心 | [vcneb/core.py](../vcneb/core.py) | 联合原子/晶胞路径、切线/弹簧/物理力；基于已有 NEB/VCNEB 思想 |
| 计算器契约与适配 | [calculator.py](../vcneb/calculator.py)、[backends.py](../vcneb/backends.py)、[abacus.py](../vcneb/abacus.py)、[vasp.py](../vcneb/vasp.py)、[qe.py](../vcneb/qe.py) | ASE energy/forces/stress；接口实现≠任意安装/赝势/输入均验证 |
| 配置与运行 | [cli.py](../vcneb/cli.py)、[config.py](../vcneb/config.py)、[material_runner.py](../vcneb/material_runner.py) | 显式参数、准备与执行分开；`prepare` 不跑 DFT |
| manager/worker、缓存与身份 | [executor.py](../vcneb/executor.py)、[provenance.py](../vcneb/provenance.py) | 独立 image 目录、固定端点缓存、完整链；并行降低墙时不等于减少数学迭代 |
| 优化器独立选择 | [optimizer_registry.py](../vcneb/optimizer_registry.py)、[block_fire.py](../vcneb/block_fire.py)、[scaled_fire.py](../vcneb/scaled_fire.py) | FIRE/逐像/分块/分阶段策略与后端独立；实验性策略保留负迁移 |
| 步与输入检查 | [step_control.py](../vcneb/step_control.py)、[vasp_contract.py](../vcneb/vasp_contract.py)、[vasp_lattice.py](../vcneb/vasp_lattice.py)、[periodic_path.py](../vcneb/periodic_path.py) | 候选几何和序列化检查，不通过时不能靠偷改物理设置掩盖 |
| 机械边界 | [epitaxial_boundary.py](../vcneb/epitaxial_boundary.py)、[slab_boundary.py](../vcneb/slab_boundary.py) | 部分夹持与 slab 边界区分；不是相同参数下自动选择的默认条件 |
| 模式与 Γ 声子 | [modes.py](../vcneb/modes.py)、[phonons.py](../vcneb/phonons.py)、[path_modes.py](../vcneb/path_modes.py)、[continuous_projection.py](../vcneb/continuous_projection.py) | Γ 本征矢/简并子空间/连续周期表示；模式幅度不等于电子极化或能量份额 |
| 条件释放与模式面 | [mode_surface.py](../vcneb/mode_surface.py)、[mode_evaluator.py](../vcneb/mode_evaluator.py)、[clamped_mode_coordinates.py](../vcneb/clamped_mode_coordinates.py)、[branch_relaxation.py](../vcneb/branch_relaxation.py) | 冻结/固定坐标释放、有限预算和分支；局域梯度通过不是全稳定性证明 |
| 曲率、响应与功 | [joint_curvature.py](../vcneb/joint_curvature.py)、[quadratic_reduction.py](../vcneb/quadratic_reduction.py)、[stationary_gap.py](../vcneb/stationary_gap.py)、[strain_work.py](../vcneb/strain_work.py)、[path_work.py](../vcneb/path_work.py) | 已有数学与研究诊断的实现；未形成独立验证的新预测理论 |

整体入口：[README](../README.md)、[架构](ARCHITECTURE_AND_USER_WORKFLOW.md)、[详细手册](USER_MANUAL.md)、[仓库布局](REPOSITORY_LAYOUT.md)。

### 3.2 软件语义与参数纪律

- 当前 `n_images=7` 是**总像数**，含两个固定端点，即五个内部像。历史对话里“7 个插帧”的口头用语不能替代实际配置。
- 默认普通 NEB 阈值为 `0.10 eV/Å`。这是项目停止准则，不是固定 meV 的能垒误差上界。
- 同胞、关闭晶胞自由度时退化为普通 NEB；不强制需要 stress。
- 默认先普通 NEB。无内部峰值的单调路径不应为做出峰而开启 CI。本次也不启动任何 CI 精修。
- ABACUS 核心研究保持 PBE/100 Ry/完整 10 au DZP 及原输入资产；GaN/VASP 模式面与主链使用 600 eV。历史 1000 eV 诊断不属于同参主证据。
- 参数应在同一计算器、同一物理协议内保持一致；不同后端的数值截断参数不能机械等同。
- 源仓库 Python 包名仍是 `vcneb`，公开项目/发行名是 `varneb`。

### 3.3 发行与使用边界

[pyproject.toml](../pyproject.toml) 的当前源码版本为 `0.0.3`。Git 有 `v0.0.1`、`v0.0.2` 标签；[0.0.2 发行记录](RELEASE_NOTES_v0.0.2.md)保存了 2026-09-29 对已发布 PyPI 文件的独立安装检查。该快照早于后续若干后端验证和研究功能，不能说等于最新主分支；没有把 `0.0.3` 的本地构建说成已发布。

[发布流程](RELEASE_PROCESS.md)及 [publish-pypi.yml](../.github/workflows/publish-pypi.yml)仅在显式标签/手动条件下发布，普通同步不发布。**本次不发布、不改变 CI、不创建 Release。**

## 4. 多后端材料验证与原始来源

### 4.1 最完整的共同案例：45.7 GPa GaN B4→B1

| 后端 | 历史作业 | 正向离散焓峰 eV/GaN | 最终广义力 eV/Å | 状态 |
| --- | ---: | ---: | ---: | --- |
| VASP | 27719610 | 0.338490866 | 0.096528 | 普通链完成；文献约 0.34，非同哈密顿量逐数复现 |
| ABACUS | 27760826 | 0.327365850 | 0.096844 | 普通链完成；不要误用早期 singleton/零压失败链 |
| QE | 27770529 | 0.329662882 | 0.098431 | 修复后的高压链完成；早期零压结果另归档 |
| ABINIT | 27770687 | 0.292422462 | 0.099105 | PBE/PSP8 生产链完成；不是旧 HGH-LDA 失败结果 |
| CP2K | 27770714 | 0.292755021 | 0.094092 | 链完成，但原生文本力/结束标记证据有缺口，见下文 |

可靠结论：五个计算器在各自固定输入下得到共同主峰拓扑，控制器具备真实多后端运行能力。**不是**完全同基组/同赝势、同能量面，也没有重新求每个后端自己的共存压力。

LAMMPS 有 ASE 接口与经典势测试，不能计入上述五套 DFT。OpenMX 未形成该类完整生产验证。“所有后端×BTO/GaN”的全矩阵未完成，不能写成所有材料均已验证。

### 4.2 从图表追到链和原始输出

| 证据层 | 位置 | 可审查内容与缺口 |
| --- | --- | --- |
| 五后端逐像绘图值 | [CSV](../paper/VARNEB_CPC/figures/gan_multibackend_validation_source_data.csv)、[汇总 JSON](../paper/VARNEB_CPC/evidence/gan_45p7_multibackend_vcneb_20260924.json) | 每后端29点、正反峰、端点差；不是原始 SCF 日志 |
| 已归档完整评估链 | [final_chains 目录](../paper/VARNEB_CPC/evidence/gan_45p7_final_chains_20260927/)、[CP2K final chain](../paper/VARNEB_CPC/evidence/gan_cp2k_45p7_chain_step_0052.traj) | SinglePoint E/F/stress 与有序几何；不是全部历史 calculator 调用 |
| 最终可见输入检查 | [全像输入检查](../paper/VARNEB_CPC/evidence/gan_45p7_full_image_input_contract_audit_20260929.md) | cutoff、网格、势/基组身份的一致性；不能外推为每次历史调用都审完 |
| hf 原始树与 SHA | [远端源索引](../paper/VARNEB_CPC/evidence/gan_45p7_remote_source_index_20260927.md) | 真实目录、job/summary/traj/input/output hashes；远端索引不等于完整数据已公开保存 |
| CP2K 原文重读 | [原能量/应力核对](../benchmarks/numerical_integrity/gan_cp2k_final_chain_raw_energy_stress_20260928.json) | 统一 Hartree→eV 因子后原能量/应力匹配；原文本未打印原子力且缺最后 PROGRAM ENDED，不能宣称完整原文力/结束证明 |
| 失败和修复记录 | [后端修复历史](BACKEND_PRODUCTION_REPAIR_STATUS.md) | 历史诊断，含已被后续结果取代的“下一步”；不应照此自动重跑 |

远端源索引定义：`R=/public/home/iai806/abacus/agent-runs/20260921-varneb-material`，`M=$R/cases/gan`。原结果树是：

```text
ABACUS: $M/abacus_vcneb_45p7_mpi4_20260922_retry1
QE:     $M/qe_vcneb_45p7_pd04_29i_20260924_retry2
ABINIT: $M/abinit_vcneb_45p7_pd04_29i_20260924
CP2K:   $M/cp2k_vcneb_45p7_k4_cut800_29i_20260924_retry1
VASP:   /public/home/iai806/abacus/agent-runs/20260918-varneb-vasp-hf/gan_b4_b1_vasp_pbe_paw_qian/vcneb_hf_n29_w9_mpi32
```

这些是历史定位，不是本次核实远端仍存在的声明。VASP 环境历史入口为 `source /public/home/iai806/Software/VASP/env.sh 6.3.2`；ABACUS 原安装在 `~/apprepo/abacus/v3.10.0LTS-intelmpi2025/`。无需新安装软件。本次不 SSH 执行计算。

## 5. 已有材料路径、文献对照和明确负结果

| 案例 | 实算结果 | 审查结论 |
| --- | --- | --- |
| BTO T→C / ABACUS | 单调，端点上升 0.087128902 eV/BTO，5/7/9总像控制 | 已完成；不是内部激活势垒，不是有限温度自由能 |
| GaN B4→B1 四方 / VASP | 0.338490866 eV/GaN；文献约0.34 | 有效机制/能量尺度对照，PBE/PAW 与原 PW91/USPP 不同 |
| GaN B4→B1 六方 / VASP | 0.385231877 eV/GaN；文献约0.39 | 已完成的独立分支对照 |
| GaN B3→B1 / VASP | 0.955290808 eV/GaN；文献约0.57 | 已完成／负结果：当前映射成单高峰，不是文献三事件机制 |
| CdSe RS→WZ cell mapping | 7.215583 meV/atom；文献2.4 | 已完成但非定量复现；PBE/PW91、映射和小峰精度有差异 |
| CdSe atomic mapping | 无同等级完整交付 | 未完成，不能以 cell 分支结果代替 |
| 早期 VASP/HfO₂ 20像 T→PO | 2.972970 meV/HfO₂；所用文献数据约2.5 | 早期特定协议对照；不可和后来 ABACUS 夹持结果混为一套 |
| 早期 VASP/HfO₂ 20像 PO→M | 83.623135 meV/HfO₂；文献约105.075 | 拓扑相似但不定量一致；是早期案例，不是最新同边界竞争矩阵 |
| BST50/PTO/PZT50 | T→C端点差31.569/221.640/27.146 meV/f.u.，无内部峰 | VCA/后端验证，来自另一会话后整合；未证明 VCA 鞍点定位 |

### 5.1 路径与文献图位置

- [逐路径比较图与 README](../outputs/neb_literature_benchmarks/)、[manifest](../outputs/neb_literature_benchmarks/manifest.json)：每路径 `*_source_data.csv`、PDF、SVG、PNG。
- [BTO 原案例](../outputs/batio3_t_to_c_pbe100_dzp10au/)：`bto_validation_provenance.json`、`bto_tetragonal_to_cubic_n7_static_audit_summary.json`、`bto_tetragonal_to_cubic_n7_static_audit_step3.traj`、端点身份和 Γ 数据已入 Git。
- [早期 ABACUS/HfO₂ 案例](../outputs/hfo2_t_to_po_pbe100_dzp10au/)：普通/CI/插值/像数对照及摘要。与上表早期 VASP20像、后文新夹持研究是三个独立证据集合。
- [VCA 说明](vasp_vca_validation.md)、[BST50](../results/bst50_vcneb_final/)、[PTO](../results/pto_vcneb_final/)、[PZT50](../results/pzt50_vca_vcneb_final/)：摘要/端点/检查；BST50另有完整绘图文件，后两者不要凭目录声称有同等图稿。
- [GaN 路线与映射方案](GAN_QIAN_2013_REPLICATION_PLAN.md)、[CdSe 方案](CDSE_SHEPPARD_2012_REPLICATION_PLAN.md)：含历史作业，不是当前状态仪表盘。

### 5.2 已知错误归因与原始数据缺口

1. [BTO 文献归因纠正](../paper/VARNEB_CPC/BTO_LITERATURE_ATTRIBUTION_AUDIT_2026-09-28.md)：文献2.1 kcal/mol是局部约束的 T→C 畸变成本，论文中的 NEB 讲氧空位迁移；不是 BTO T→C 的 NEB 文献基准。早期 `outputs/neb_literature_benchmarks/README.md` 中“restrained-NEB”等旧措辞不能继续引用。
2. 旧绘图脚本 [plot_literature_path_benchmarks.py](../scripts/plot_literature_path_benchmarks.py) 的默认摘要路径部分在 `tmp/plot_data/`。CSV/图/manifest已跟踪，**默认重放仍不能假定新克隆有全部 tmp 原数据**；原始摘要须按 manifest/远端索引另行定位。没有在本交接中虚构补全。
3. GaN 文献曲线为手工数字化，其横轴和本工作归一化弧长不是同一反应坐标；叠图不证明路径几何一一相同。
4. B3 三次事件要求的三段独立收敛路径与一致焓基准，没有完整交付；一条长链或一个高峰不能冒充三段拼接。

## 6. 加速策略与模式面：哪些成功、哪些不能主张

### 6.1 加速 benchmark 的有限成果

| 同起点、首次跨过共同0.10阈值 | ABACUS成功启动数 | 结论 |
| --- | ---: | --- |
| BTO global FIRE → BlockFIRE 0.01/image | 205→44，减少78.5% | 已有同起点消融，非所有体系的通用倍数 |
| HfO₂ global FIRE → staged | 217→134，减少38.2% | 已有转移收益，含粗/细两个启动段成本 |
| VASP GaN / CdSe 长链 | 分别1 / 4更新跨阈值 | 成功续算；CdSe基准未观测到最终跨阈值，不能给速度比 |
| HfO₂ 长链 BlockFIRE | 曾在约0.188 eV/Å平台 | 负迁移，不能删去 |
| 候选晶格门控+半步 | 完整步拒绝，半步通过并续算 | 可靠性机制；不等于已隔离出其普适收敛增益 |

统计的是**成功外部 ABACUS 进程启动**，不是电子步、并行倍数、完整每次SCF质量审计。历史基准运行的停止参数更严，但比较截断在首次共同0.10交叉处；不能拿更严基准的全程费用来夸大收益。

关键代码见第3节。完整依据：[证据审查](../paper/VARNEB_CPC/ACCELERATION_EVIDENCE_AUDIT_2026-09-27.md)、[原 Slurm 859行数据](../benchmarks/convergence/hf_slurm_abacus_steps_20260928.psv)、[启动数检查结果](../benchmarks/convergence/hf_slurm_abacus_launch_audit_20260928.json)、[首次阈值链](../benchmarks/convergence/first_crossing_chains_20260928/)、[原实验记录含失败](../benchmarks/convergence/autoresearch.md)。

重放代码：[启动数](../scripts/audit_acceleration_slurm_accounting.py)、[完整轮廓/几何比较](../scripts/audit_acceleration_path_equivalence.py)。后者显示链并非严格相同；路径等价容差也没有事先登记。可主张案例限定的启动数节省，不能主张相同MEP/同鞍点或首创FIRE思想。

### 6.2 四套主要二维面与数据、绘图器

`F`在本小节表示 `paper/VARNEB_CPC/figures/`。下列 stem 对应同名 `.pdf/.svg/.png`、`_source_data.csv`、`_qa.json`；不要把其他历史同名近似图替代它们。

| 物理对象 | 实算点数 | 主图stem | 审计输入／绘图器 |
| --- | ---: | --- | --- |
| BTO `Qy=0` 受限变胞面 | 27，9×3 | `bto_qy0_restricted_sheet_27_20260929` | [27节点](../benchmarks/numerical_integrity/bto_qy0_27_node_sheet_20260929.json)；[plot](../scripts/plot_bto_qy0_restricted_sheet27.py) |
| BTO 固定C胞双软模面 | 289，17×17 | `bto_frozen_soft_mode_289_20260929` | [289节点](../benchmarks/numerical_integrity/bto_transverse_soft_frozen289_20260929.json)；[plot](../scripts/plot_bto_transverse_soft_landscape.py) |
| GaN 峰附近原子—应变冻结焓切面 | 289，17×17 | `gan_600eV_local_joint_dft289_20260930_v2` | [289节点](../benchmarks/numerical_integrity/gan_600eV_ts_2d_17x17_refinement_20260930.json)；[plot](../scripts/plot_gan_600eV_local_joint_cut17.py) |
| GaN 中央路径—横向原子冻结焓切面 | 162＝18中心+144偏离 | `gan_600eV_atomic_transverse_162_20260929_v2` | [162节点](../benchmarks/numerical_integrity/gan_600eV_atomic_tube_q9_20260929.json)；[plot](../scripts/plot_gan_600eV_atomic_transverse_q9.py) |

图与CSV均在[figures目录](../paper/VARNEB_CPC/figures/)；[统一离线重建器](../scripts/rebuild_cpc_landscapes.py)固定了完整输入集合。加密不是只插值：BTO冻结面新增208个静态DFT，旧9×9预测新点最大误差0.488 meV/BTO；GaN局域面新增208点，误差0.00698 meV/GaN；GaN中央切面新增72点，误差0.0073 meV/GaN；BTO受限面15个新条件节点误差1.131 meV/BTO。

必须保留的物理边界：

- BTO Γ 力常数和本征矢采用5原子 `1×1×1` 胞，见[原 Γ 数据](../outputs/batio3_t_to_c_pbe100_dzp10au/bto_cubic_gamma_phonon_provenance.json)与[本征对](../outputs/batio3_t_to_c_pbe100_dzp10au/bto_cubic_phonopy_gamma_eigenpairs.npz)。没有把只算Γ不必要地扩胞。
- 立方 BTO 还有被两坐标省略的第三软方向；无约束全释放二维下包络不能同时把不稳定 C 保持为其原能量的极小点。`Qy=0` 是声明的限制，不是无约束面。[端点阻碍审查](../benchmarks/numerical_integrity/bto_soft_triplet_conditional_endpoint_20260928.json)。
- BTO 59点旧冻结面、79点路径适应冻结面、289点双软模面、27点条件面是不同对象；“漂亮旧图”不能靠同名替换物理定义。
- GaN 两面均使用45.7 GPa和原600 eV协议。中央切面不覆盖整个IS→FS，也不是完全条件释放的两声子面。
- GaN原子—应变联合Hessian有负方向，但[能量—应力一致性诊断](../benchmarks/numerical_integrity/gan_600eV_energy_stress_consistency_20260928.json)仍有法向体积应变差异，不能认证严格TS。既有比较不能靠擅自升ENCUT到1000 eV“修复”。
- 两段固定中心精修已完成；[协议与结果](GAN_B4_B1_FIXED_CENTER_TWO_SEGMENT_2026-10-04.md)、[链/manifest/raw audit](../paper/VARNEB_CPC/evidence/gan_45p7_split_20261004/)、[图CSV](../paper/VARNEB_CPC/figures/gan_45p7_split_local_mode_20261004_v4_source_data.csv)。拼接峰338.465 meV/GaN，切线余弦0.9636；小图域内只有中心像，虚线连邻像不是新增DFT点或已认证局域MEP。
- 所有面是0 K静态能量/焓，不是有限温自由能；投影、声子虚频、平滑等高线均不单独证明机制或精确鞍点。

## 7. HfO₂ 科学研究：最终事实、修正与阴性结果

### 7.1 核心问题、固定条件和坐标定义

问题：**同一部分夹持下，内部原子可以如何容纳重排，这怎样改变保持型翻转相对于离开极性相的候选竞争，以及对应变作用的判断？**

本节为ABACUS/PBE、100 Ry、完整10 au DZP、原六资产、Gamma2³、12原子Hf4O8＝4个HfO₂、零压零场。原T参考的x/z平面固定，y胞行与两项相关倾斜开放；普通NEB0.10/noCI。受限局域释放另有0.03 eV/Å目标，不能混作更严的NEB标准。

同一应变下各通道使用相同PO+初态能量。`Qx/Qy/Qz`是参考几何图样；“保持/反转”指某参考图样符号，两者都反转z结构代理。它们没有被认证为文献不可约声子标号。O–Hf物种平均位移代理不是电子极化；位移平方范数份额不是能量份额。

### 7.2 原八条训练链：可靠的有限比较

| 通道 | 0%峰 | +1%峰 | 变化（meV/HfO₂） |
| --- | ---: | ---: | ---: |
| PO→M | 98.458122 | 97.227719 | −1.230403 |
| PO→T | 84.416624 | 96.763650 | +12.347027 |
| 原保持型翻转 | 96.828183 | 125.424431 | +28.596248 |
| 原反转型翻转 | 338.968965 | 294.944634 | −44.024331 |

八条都通过项目普通残差标准。可靠认识：应变作用是通道相关的；只看高反转型下降，不能推断更低翻转路径总体受益。它们是**原候选集合**，不是搜索穷尽的最低势垒。

原生与汇总入口：[八链终态](../benchmarks/hfo2_channels/20261008/clamped_plus1_reversing_terminal_E148_20261014/)、[连续表示说明](HFO2_CONTINUOUS_CHAIN_OBSERVATIONS_2026-10-08.md)、[离线案例](../examples/cases/hfo2_competing_channels/README.md)。初始整链采用一次注册的周期代表，不逐像重新折叠来制造平滑模式轨迹。

### 7.3 新低候选改变了原解释

| 应变 | 原保持型峰 | 更低完整保持型峰 | PO→T峰 | T减新保持型 |
| --- | ---: | ---: | ---: | ---: |
| 0 | 96.828183 | 83.166281 | 84.416624 | 1.250343 |
| +0.5% | 98.530877 | 87.718545 | 90.418178 | 2.699634 |
| +1% | 125.424431 | 91.698937 | 96.763650 | 5.064714 |

单位meV/HfO₂。更低零应变局域构型已放回相同有序端点的完整自由链，9更新/63新SCF后通过0.10；不是只报孤立低点。以同一已固定源向两已有应变转移，+0.5%首轮即通过、7新SCF；+1%4更新、35新SCF。

**原保持型应变增量28.596248降为8.532656 meV/HfO₂。原“应变强烈硬化保持型、逃逸明显更容易”的候选解释不稳健。**但1–5 meV/HfO₂小差不足以认证全局翻转/逃逸排序、速率或器件寿命。

更低中心与原中心之差包含O_z交错约±0.200 Å、Hf_x约±0.068 Å；旧中心本来就有均值相消的横向O运动，不是“首次出现横移”。路径改变发生在相同物理设置下，没有将原链强制投影成新路径。原链也并非人为冻结图样，不能把它改写成有意受限的劣质基准。

### 7.4 匹配释放与最后两因素交叉

F＝原候选中心；A＝保留同四个观测量、固定整个胞的原子释放；J＝再开放允许的胞变量。固定胞A已降约13.599（0%）和17.172/17.175（+1%）meV/HfO₂；额外J约0.018和0.485/0.484。T通道三应变同定义总释放均小于0.25。

可靠认识：主要局域降能可在固定胞中发生，夹持没有消除内部原子容纳；不能全部归于开放胞松弛。但局域F/A/J不是完整路径势垒的唯一可加分解，局域梯度通过也不是全Hessian稳定证明。

最后四格实验保留向量q0/q1；**每个是四分量几何观测向量，不是两支声子或早期Q1/Q2**：

| 应变 | 向量 | 相对同应变初态局域能量 meV/HfO₂ | 释放梯度 eV/Å | 本次新SCF |
| --- | --- | ---: | ---: | ---: |
| 0 | q0 | 83.211101 | 0.011617 | 0，复用 |
| 0 | q1 | 114.553122 | 0.015790 | 10 |
| +1% | q0 | 92.398077 | 0.023149 | 4 |
| +1% | q1 | 107.767116 | 0.022660 | 0，复用 |

同应变图样差31.342021→15.369039；有限交互差−15.972982 meV/HfO₂消掉共同初态偏移。+1%持q0的局域点只比更低完整链峰高0.699140；持q1则高约16.068。这连接了保留图样、横向容纳、开放轴与候选响应，**没有**给出唯一Qx因果贡献或普适耦合常数。

固定q1的相对能量随应变下降主要涉及初态升高，不能直接叫该构型绝对软化。四分量同时变化，也不能只挑其中Qx解释全部作用。

### 7.5 独立检验：已完成的负结果，不是待继续调模型

- T局域8个独立候选在18新SCF内完成。张量二次插值最大原点对齐误差0.089077 meV/HfO₂，原曲率模型约0.129769/0.131911；原始未对齐误差的排名又不同。未证明新模型优于简单强对照。
- 修订保持型8候选在24新SCF内完成。非线性项略降RMSE，却增加最大误差；共同目标上0.01步长的谐/非线性最大误差0.033195/0.035625，RMSE0.027153/0.024916。未建立非线性优势。
- 真实+0.5%物理应变检验保留138条冻结记录，其中**6条数值预测、132条未能给出有效迁移预测**。不是138次验证成功。局域模型没有供应有效物理迁移，不能声称新预测机制通过。
- 中点从两侧初始化：保持型峰相差10.820844，T约−0.131935 meV/HfO₂。初始“中点非单调应变效应”的强解释因而失去稳健性；未偷偷替换原预测标签。
- 三模位移高覆盖不保证梯度覆盖；所释放补空间曾有负方向。原两坐标稳定条件面的拒绝不因后续更换目标而变成成功。
- +1%局部原0.5 Å窗口内两起点失败是预算/边界相关，不是“物理上不存在”。明确改为1.0 Å数值窗口后获得A/J；电子参数不变，旧失败仍保留。两应变数值窗口不同，不可隐瞒。

原预测与阴性结果入口：[T独立检验](../benchmarks/hfo2_channels/20261008/G3_T_independent_E258_20261017/)、[预测冻结](../benchmarks/hfo2_channels/20261008/G3_anharmonic_model_E265_20261018/)、[保持型独立执行](../benchmarks/hfo2_channels/20261008/G3_P8_execution_E266_20261018/)、[物理+0.5%检验](../benchmarks/hfo2_channels/20261008/physical_holdout_E268_20261018/)、[起点敏感性](../benchmarks/hfo2_channels/20261008/midpoint_reverse_E269_20261018/)、[冻结原子/胞交叉](../benchmarks/hfo2_channels/20261008/midpoint_factorial_E274_20261019/)。

### 7.6 原始归档、分析、图表一一对应

本小节B代表已跟踪的 `benchmarks/hfo2_channels/20261008/`；E编号只用于在历史中定位，不是研究创新的名称。

| 内容 | 关键证据位置 | 对应图／代码 |
| --- | --- | --- |
| 更低零应变完整连接 | [E273](../benchmarks/hfo2_channels/20261008/preserving_connection_E273_20261019/) 中 `native_terminal.tar.gz`、`native_terminal_receipt.json`；[原子运动分析](../benchmarks/hfo2_channels/20261008/preserving_connection_E273_20261019/atom_motion/analysis.json) | [analyze_hfo2_oxygen_motion.py](../scripts/analyze_hfo2_oxygen_motion.py)；新主图4来源 |
| 三应变完整链 | [E275](../benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019/) 的 `completed/strain_0050`、`completed/strain_0100`，各含 `native_terminal.tar.gz` 和 `collection_receipt.json`；[完整分析](../benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019/completed_analysis_full/analysis.json) | `analyze_terminal.py`、`plot_summary.py`；[主图9及来源](../benchmarks/hfo2_channels/20261008/preserving_strain_E275_20261019/figures/summary_three_strains/README.md) |
| 原匹配局域释放 | [E270](../benchmarks/hfo2_channels/20261008/matched_pattern_release_E270_20261018/)、[E271](../benchmarks/hfo2_channels/20261008/matched_pattern_seeded_E271_20261018/)、[E272](../benchmarks/hfo2_channels/20261008/preserving_pattern_branches_E272_20261018/) | T三应变、P零应变及原失败；不能用旧pending覆盖后来E276 |
| +1%局域缺口补齐 | [E276](../benchmarks/hfo2_channels/20261008/preserving_local_completion_E276_20261019/) 的 `completed/minus`、`completed/plus` 各有原生tar/receipt；`complete_analysis/`、`mechanism_analysis/` | `analyze_terminal.py`、`mechanism_data.py`、`plot_mechanism.py`；`figures/matched_controls/`主图10 |
| 最后两格交叉 | [E277](../benchmarks/hfo2_channels/20261008/strain_pattern_cross_E277_20261021/) 的 [原生tar](../benchmarks/hfo2_channels/20261008/strain_pattern_cross_E277_20261021/completed/native_terminal.tar.gz)、[receipt](../benchmarks/hfo2_channels/20261008/strain_pattern_cross_E277_20261021/completed/receipt.json)、[分析](../benchmarks/hfo2_channels/20261008/strain_pattern_cross_E277_20261021/analysis/analysis.json) | `analyze.py`、`plot.py`；`figures/`主图11、正文表5 |
| 最终原子机制综合 | [mechanism_synthesis](../benchmarks/hfo2_channels/20261008/mechanism_synthesis_20261021/) 的 `analysis.json`、`center_motion.csv`、`qa.json` | [plot_hfo2_mechanism_synthesis.py](../scripts/plot_hfo2_mechanism_synthesis.py)；[主图4 PDF](../benchmarks/hfo2_channels/20261008/mechanism_synthesis_20261021/hfo2_mechanism_synthesis.pdf)、同名SVG/PNG |

E275原作业28951226/28951227于10月19日完成，42新SCF；E276两作业28954771/28954772完成，27新SCF；E277两作业28995326/28995327于10月21日完成，14新SCF、15.52分配核时。它们是历史终态，不是待续算句柄。不能将这些局部成本相加当成全研究成本；早期路径/导数/失败预算另计，目前没有宣称完成全项目统一费用总账。

### 7.7 与已有文献的关系：当前最重要的原创性风险

[四篇最接近工作对照](../benchmarks/hfo2_channels/20261008/mechanism_synthesis_20261021/literature_comparison.md)和[来源/读取范围](../benchmarks/hfo2_channels/20261008/mechanism_synthesis_20261021/literature_comparison.provenance.md)分别对应：

1. Liu–Hanrahan 2019：自由/固定胞边界敏感性已经存在；本工作的部分夹持不等于其全固定T胞。
2. Behara–Van der Ven 2022：应变坐标、变体及shuffle关系已有研究；本工作几何Qx不是自动认证的文献irrep。
3. Qi–Singh–Rabe 2025：伴随非极性模式反转区分路径已有概念；记录明确部分判断只依据最终摘要，未冒称阅读全文。
4. Mu等2025：C/S是横向/直行、N/T是不过/穿过Hf层；本工作的“保持型/反转型”不是这些标签的同义词。

有序采样显示新旧保持型都有横向、不穿过已注册Hf层带的几何特征；完整发表轨迹/端点一一映射没有建立。**“更充分松弛了一个已知类别的路径”仍是足够解释，不能排除。**因此当前不主张新路径拓扑、新声子、普适选择规律或首创多模耦合。

## 8. 论文、主图与离线复核入口

### 8.1 两篇稿件必须区分

| 项目 | 当前来源与产物 | 状态 |
| --- | --- | --- |
| CPC软件稿 | [唯一TeX源](../paper/VARNEB_CPC/varneb_CPC.tex)、[证据表](../paper/VARNEB_CPC/MANUSCRIPT_EVIDENCE.md)、[投稿缺口](../paper/VARNEB_CPC/SUBMISSION_READINESS.md)；本机预览 `paper/VARNEB_CPC/varneb_CPC.pdf` **未被Git跟踪** | PDF不随本次push交付；从固定提交编译。作者元数据、CPC归档、最终干净版本/原始证据边界不是自动完成 |
| HfO₂材料稿 | [唯一科学正文源](../paper/VARNEB_JCTC/MANUSCRIPT_DRAFT.md)、[TeX薄封装](../paper/VARNEB_JCTC/varneb_JCTC_draft.tex)、[阅读PDF](../paper/VARNEB_JCTC/reading/VARNEB_HfO2_manuscript.pdf) | 当前题为 Atomic accommodation reshapes strain responses of competing hafnia pathways；可评阅，不能叫JCTC目标已实现 |
| 当前文稿构建记录 | [scientific_acceptance_20261021.json](../paper/VARNEB_JCTC/reading/scientific_acceptance_20261021.json) | 64页含文内附录；55项相关测试；明确 `JCTC_scientific_acceptance_passed=false`、非全仓测试 |
| 当前科学判断 | [10月21日意见](HFO2_PRESUBMISSION_VERDICT_2026-10-21.md) | 内部审视，不冒称真实期刊审稿；不建议仅凭当前证据默认首投JCTC |

HfO₂现有PDF：主图4在p11；表5在p22；图11在p23；讨论p24；结论p25；文内附录p27起；旧主图4完整移为附录F.3/图A9、p57。11幅主图；全部20个图目标含附录。64页不是科学成熟度指标。

本机CPC预览可能包含未提交的架构图版本，不将它补加Git来冒充冻结稿。新克隆能取得的是已跟踪TeX、图源和数据；以固定提交的独立编译为准。最新HfO₂阅读PDF则已经跟踪。文件在本机存在、已跟踪、已推送，是三种不同状态。

正文四部分仍为：路径竞争、机械约束、模式解释、独立检验及后续连接控制。详细阴性结果应保留，不必让读者先读全部诊断才理解主结果。旧草稿、旧编号、旧PDF收据只代表其历史版本。

### 8.2 不调用DFT的复核命令

以下从仓库根目录执行，依赖现有Python环境；图形命令需ASE/NumPy/matplotlib等。输出必须使用**全新目录**，不要覆盖研究证据。它们不运行SSH、Slurm或计算器。命令是给审查者的可选检查，不是本次全部重新执行的声明。

```powershell
# 仅标准库：全路径/最终交叉与原模型检验表，分别导出
python -m examples.hfo2_mechanism_results --output tmp/review_20261022_mechanism
python -m examples.hfo2_paper_results --output tmp/review_20261022_models

# 四套CPC二维面：验证源CSV和真实点数，不要求跨字体环境PNG字节相同
python -m scripts.rebuild_cpc_landscapes --output-dir tmp/review_20261022_landscapes

# 最终原子机制主图
python -m scripts.plot_hfo2_mechanism_synthesis --output tmp/review_20261022_figure4

# 最后交叉原生tar与两旧对角候选独立重放
python benchmarks/hfo2_channels/20261008/strain_pattern_cross_E277_20261021/analyze.py --output tmp/review_20261022_cross

# 已归档Slurm启动次数，非实时sacct，不跑ABACUS
python -m scripts.audit_acceleration_slurm_accounting audit --raw benchmarks/convergence/hf_slurm_abacus_steps_20260928.psv --output tmp/review_20261022_acceleration.json
```

聚焦最新科学结果的测试命令（最近已有55通过记录，不等于最新全仓都通过）：

```powershell
python -m pytest -q tests/test_hfo2_oxygen_motion.py tests/test_hfo2_strain_pattern_cross.py tests/test_hfo2_preserving_local_completion.py tests/test_hfo2_preserving_completion_mechanism.py tests/test_hfo2_terminal_manuscript.py tests/test_hfo2_mechanism.py tests/test_hfo2_mechanism_results_example.py
```

编译使用既有TeX Live，不安装新工具：

```powershell
# 进入各自目录；只生成build产物，不修改唯一源
Set-Location paper/VARNEB_JCTC
latexmk -norc -lualatex -interaction=nonstopmode -halt-on-error -outdir=build varneb_JCTC_draft.tex
# CPC需另回仓库根目录，再进入paper/VARNEB_CPC
# latexmk -pdf -interaction=nonstopmode -halt-on-error varneb_CPC.tex
```

编译成功仍须目视看图、字体、页码、单位；不自动把build产物覆盖已审阅阅读PDF。

### 8.3 新克隆与原始证据可用性

- 最稳妥独立审查从固定Git提交的干净克隆/`git archive`开始，不能只在共享脏工作区跑通就声称可复现。
- [源码快照工具与说明](../paper/VARNEB_CPC/CPC_SOURCE_SNAPSHOT.md)可固定提交导出；不是已经完成CPC/Mendeley归档。
- 历史CPC干净提交曾有605通过/1跳过记录；不能把它搬来代表10月21日全仓测试。最新55项也是限定范围。
- BTO/GaN主要图的CSV/审计JSON/部分完整链已跟踪；HfO₂最近E273/E275/E276/E277有跟踪原生压缩归档和receipt。**并非整个研究以来全部SCF原文都已入Git**。
- `.tar.gz`是未展开证据包，不是文本摘要。先核对receipt、成员清单与路径安全，再只在新目录展开；不覆盖其他快照或执行归档里的运行脚本。
- VASP POTCAR、DFT二进制、部分赝势/轨道身份只记录hash，内容留在授权环境；不能为了“审查完整”上传受限资产或凭hash伪造文件。
- 有些分析记录原hf绝对路径：能离线重放到什么层次，应以实际已跟踪归档为准；hash匹配不能证明未保存的数据当前可取。
- 本次不提供或声称全项目新增SCF/节点小时总计，没有将缓存命中、摘要数、预测条数冒充真实新计算量。

## 9. 执行失误、研究失败与停滞原因：应由原代理承担的部分

### 9.1 已记录的工程失误

| 问题 | 直接后果 | 后续边界 |
| --- | --- | --- |
| ABACUS名义32MPI实为32 singleton、错误绑核 | 重复电子步输出、同目录争写、速度异常 | 不是ABACUS物理上慢；错误轨迹不能续作生产证据 |
| 高压遗漏传压 | 部分QE/CP2K等实际零压 | 零压达残差也不能并入45.7GPa对照 |
| CP2K 400Ry误写400eV | 端点应力异常、无效尝试 | 单位是输入责任，不应靠更多NEB修复 |
| BTO标签/相身份不符 | T/C文件反置或退化同相 | 文件名不能代替结构和软模/极性识别 |
| 对齐后复用旧端点检查 | 有效几何与已审计hash不一致 | 对齐后的真正输入才是检查对象 |
| VASP Bravais冲突过度依赖调容差 | 没及时锁定序列化末位噪声和实际输入 | 17有效数字往返、原生探针比一味调SYMPREC更针对根因 |
| ABINIT MPI/编译运行环境及赝势协议错误 | PMI等待、错误模型/端点 | 旧HGH-LDA失败不能叫PBE验证 |
| 历史1000eV诊断及早期参数纪律不足 | 与原600eV对照不等价 | 当前同参主证据排除它；作者多次纠正是必要的 |
| 单次力回弹就停或过早称发散 | 截断正常优化过程 | 应看完整健康窗口；几何真发散与受力回弹分开 |
| 打包遗漏忽略数据、残留tmp依赖 | 共享树跑通而新克隆缺数据 | 已修部分CPC输入；其他历史源缺口仍必须披露 |

具体材料协议和日志索引见[后端修复记录](BACKEND_PRODUCTION_REPAIR_STATUS.md)、[VASP输入契约](vasp_input_contract.md)、[GaN分支计划](GAN_QIAN_2013_REPLICATION_PLAN.md)。这些根因不能全部归为集群随机故障。

### 9.2 科学与组织上的主要问题

1. **主线过多且停止条件不清。** 多后端、加速、TS、二维面、预测、应变选择性轮流成为“最后一个关键缺口”；原代理未及时分清主结论、辅助模块和可放弃假设。
2. **诊断膨胀替代科学推进。** 身份、预算、原生重放都必要，但修好它们不等于得到新科学问题的答案；曾用不断补检查回应原创性不足。
3. **最接近文献对照做得太晚。** 横向氧路径、非极性伴随反转、应变敏感性早有研究；不能在算了很多之后才判断“新”在哪里。
4. **对创新过早乐观。** 接口、多图、非加性四格、FIRE组合都不是自动方法创新；新局域模型的阴性结果应更早促成收缩。
5. **起点/分支敏感性改变结论。** 新低链大幅削弱原强应变硬化解释，说明原物理叙述曾先于充分的同问题对照；不应隐瞒这种修正。
6. **交付状态与论文达标混淆。** “有限任务完成”不能表达成“JCTC目标达成”。最新判断已明确否定这种等同，但此前汇报未一直保持清楚区分。
7. **文档组织没有跟上。** 一份旧日期handoff不断追加，长计划内混存pending/完成；大量E编号和历史诊断进入阅读路径，作者难以读懂物理问题。改名`papar→paper`并不等于完成仓库整洁化。

这些不是对每一项计算都无效的判断：软件、多后端、二维面、最后连接链的有限结果仍有价值。问题是不能用工作量、审计严密或篇幅替代科学增量。

## 10. 可靠、失败、没做：给审查者的最终清单

### 10.1 目前可保留的主张

- 统一控制器有多后端真实材料运行证据，后端与优化器选择分离。
- 已测试BTO/HfO₂同起点对照节省外部ABACUS启动，限定案例、不保证相同MEP。
- BTO/GaN四套二维面有真实采样，冻结/受限释放和图域有明确边界。
- HfO₂同部分夹持下混合原子容纳会改变保持型候选峰与表观应变硬化；有局域干预及完整链连接，不止投影叙述。
- 两因素四格局域能量非加性，且不由共同初态平移解释；只对这四候选成立。
- 原有局域模型/物理应变独立检验没有证明新方法优势，失败结果完整保留。

### 10.2 不应继续当成果的主张

- 首次提出VCNEB、逐像FIRE、信赖域、模式投影、多模耦合或横向氧路径。
- 所有ASE计算器无条件无缝可用、所有后端×材料全矩阵已完成。
- BTO单调端点上升就是内部鞍点活化能，或0K能量面等于有限温相变自由能。
- GaN/B3或CdSe已定量复现文献、三事件已完成三段独立计算。
- 二维插值像素是新DFT数据、投影路径必在面上、中心星号已认证TS。
- HfO₂当前通道穷尽、全局最低势垒已排序、可以推器件寿命/成核/速率。
- 新预测模型、普适加速理论、唯一模式能量贡献或可迁移机械选择定律已建立。
- 论文编译、55项测试或任务goal完成证明JCTC达标。

### 10.3 尚未做成，但本次不授权自动补算

- 已发表HfO₂路径与本工作完整有序端点/中间像的一一对应。
- 支撑细小翻转/逃逸差值的连续峰/全变量驻点和误差证据。
- 任意有限q声子模式作为坐标的材料级完整验证；已有Γ案例不能替代。
- 所有新低局域分支的全Hessian/全局条件下包络证明。
- 具有清楚成本优势且独立预测成功的新模型或预条件NEB方法。
- 全后端BTO/GaN矩阵、CdSe atomic完整交付、GaN B3三独立段。
- 本研究线的双层hBN/MPB等建议案例；不能混入其他会话结果冒充完成。
- 全研究原生输出公开全集、统一费用总账、CPC归档及最终作者/基金等元数据。

以上**不是新增任务单**。先判断它们是否对选定论文主张必要；不必要的应明确放弃，而不是全部补齐。

### 10.4 独立审查最值得回答的三个问题

1. 新旧保持型在文献的路径类别中到底是什么？现有几何证据是否足以支持任何“新机制”措辞？
2. 28.596→8.533这一变化能否给出超出“另一个起点找到更低路径”的、明确且可证伪的物理认识？如不能，应怎样诚实定位材料稿？
3. 软件稿与材料稿各自最有力的论文是什么？哪些研究功能应从核心论证移出，哪些失败足以终止原方法主张？

审查不必预设JCTC，也不必保住原代理拟定的标题或论证。当前最严肃的问题是科学增量，不是缺一张漂亮图，更不是缺少更多计划。

## 11. 本次交接的完成范围

本文和旧交接入口只作文档交付。完成标准：关键链接/路径可定位，列明Git内与远端/受限证据，数值与现有源表一致，已提交并推送；实际提交号以Git为准。

未启动任何新DFT、声子、NEB、SCF、集群任务、自动监控、CI、发布或投稿；未因整理文档更新计算状态为“已验证”；未修改物理输入、原结果、论文或用户五项既有修改。没有以本文替代原始输出，没有声称本次重新跑过所有历史复核命令。

本次实际核验：

- 本文129个本地Markdown链接均可解析；文件目标已跟踪，目录目标有已跟踪内容。这个检查不等于把链接目录内所有未跟踪文件都纳入交付。
- 发现CPC本机PDF未跟踪，已在第8节明确本地预览／新克隆编译的区别，没有强行加入可能包含用户修改的PDF。
- 两个标准库离线导出命令实际通过：最终四格及路径上下文；原八链与T/P的19/23模型臂。最终交互差重算为−15.972982391758705 meV/HfO₂；两命令均报告0新DFT。临时导出文件放系统临时目录，不加入Git。
- 未重新运行全仓测试、原生归档全集重放、四面重绘或论文编译；文中的历史通过数均注明原范围，不把文档检查升级为新的科学认证。

给接手会话的最终提醒：**先审查现有论证是否成立，再决定研究如何收束。不要沿着旧计划把所有未完成事项顺次算下去。**
