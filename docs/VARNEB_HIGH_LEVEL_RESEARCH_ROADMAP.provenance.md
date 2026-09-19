# Provenance: VARNEB high-level research roadmap

**Artifact:** `docs/VARNEB_HIGH_LEVEL_RESEARCH_ROADMAP.md`  
**Prepared:** 2026-09-16  
**Purpose:** 记录路线图所依据的仓库证据、原始论文和推理边界；本文件本身不把候选方法提升为已完成结果。

## 1. 仓库内来源

| 来源 | 用途 | 边界 |
|---|---|---|
| `docs/theory.md` | 当前 `Q=(s,F)`、`cell_scale`、广义力、NEB/CI、模式与收敛语义 | 描述当前实现，不证明候选 intensive/Hencky metric |
| `docs/next_priority_plan_2026-09-16.md` | 多后端、BTO Γ 模、论文篇幅与工作依赖 | 其中 VASP 状态已被 2026-09-16 后续生产结果超越 |
| `VCNEB_PROJECT_PLAN.md` | 项目范围、测试和历史决策 | 部分实时运行状态可能过时，不能代替 committed summary |
| `paper/VARNEB_CPC/MANUSCRIPT_EVIDENCE.md` | claim-to-evidence 格式 | “VASP 无权威材料路径”一行已过时，需在 P0 更新 |
| `vcneb/core.py` | 当前 legacy metric、路径力和运行 API | 仅作为现状基线 |
| `vcneb/calculator.py`, `vcneb/abacus.py`, `vcneb/vasp.py`, `vcneb/qe.py` | calculator contract 和后端能力层 | OpenMX 尚未实现 |
| `vcneb/modes.py`, `vcneb/phonons.py` | 现有模式引导、子空间、Γ 模对角化和投影 | 不等于原子--应变耦合广义 Hessian |
| commit `cdde89c` | cu17 VASP 静态 baseline runner | VASP 生产结果在远端运行目录，仍需同步入库 |

## 2. 当前材料证据

### ABACUS

- BTO T→C 5/7/9-total-image 证据由仓库 `outputs/batio3_t_to_c_pbe100_dzp10au/` 下的 provenance、summary 和图源数据承载。
- HfO2 T→PO 普通/CI/控制分支由 `outputs/hfo2_t_to_po_pbe100_dzp10au/` 和相应 provenance 承载。
- BTO cubic Γ 点 force constants/eigenpairs 和路径投影已经生成；其物理结论仍受单端点 Γ 模、简并子空间和非驻定路径解释边界约束。

### VASP

2026-09-16 在 `235 → cu17` 完成 BTO T→C 7-total-image 普通 VCNEB：

- 5 个内部 image，每个静态调用 40 MPI，串行遍历 image；
- 固定且缓存两个端点，不在每次迭代重复计算；
- source commit：`cdde89ce42facf76fc3de815a5df40897ac8ae9f`；
- 最终广义力：`0.09514793146149961 eV/Å`；
- barrier/reaction enthalpy：`0.04264410000000396 eV`，无内部 barrier；
- 最短距离：`1.801506 Å`；最大 deformation：`0.0516796`；
- audit：`status=ok`, `issues=[]`。

远端主目录：

`/home/zhuxd/abacus/agent-runs/20260916-varneb-v/batio3_vasp_pbe_paw_static/vcneb_tetragonal_to_cubic_n7_cu17_serial_cdde89c`

端点静态记录：

- `.../static_initial_ecut600_cdde89c/vasp_static_summary.json`
- `.../static_final_ecut600_cdde89c/vasp_static_summary.json`

这些数值在同步入仓库前属于“已完成、待归档”的运行证据；论文 source data 不应只引用本 provenance 文本。

### QE / OpenMX

- QE：adapter 和 preflight 已实现，真实材料路径按研究决定暂停。
- OpenMX：尚无静态 energy/force/stress adapter，不得宣称 material-path validated。

## 3. 原始论文来源及其使用方式

### 3.1 变胞/固态 NEB

1. Qian *et al.* (2013), DOI [10.1016/j.cpc.2013.04.004](https://doi.org/10.1016/j.cpc.2013.04.004).  
   本地 PDF：`C:/Users/zhu/Zotero/storage/LI6AZDJW/Qian 等 - 2013 - Variable cell nudged elastic band method for studying solid–solid structural phase transitions.pdf`。  
   用于确认：VCNEB 已将原子与晶胞自由度放在焓面上共同优化，讨论 finite strain、fractional coordinates、可变弹性常数、CI 和切线。  
   边界：路线图不声称复现未公开实现细节，也不把基本 VCNEB 作为 VARNEB 新理论。

2. Sheppard *et al.* (2012), DOI [10.1063/1.3684549](https://doi.org/10.1063/1.3684549).  
   本地 PDF：`C:/Users/zhu/Zotero/storage/TQI7HTTA/Sheppard 等 - 2012 - A generalized solid-state nudged elastic band method.pdf`。  
   用于确认：G-SSNEB 以统一结构向量处理原子/晶胞自由度，并通过长度/Jacobian 缩放讨论体系大小敏感性；不同超胞还可能出现协同与成核机制转变。  
   边界：VARNEB 的 intensive metric 必须与这一先例逐项比较；复制不变性首先是正确性而非独立创新。

3. Caspersen and Carter (2005), DOI [10.1073/pnas.0408127102](https://doi.org/10.1073/pnas.0408127102).  
   用于普通 SSNEB 的历史边界和固--固相变路径依赖性。

### 3.2 收敛加速

4. Makri, Ortner, and Kermode (2019), DOI [10.1063/1.5064465](https://doi.org/10.1063/1.5064465); author manuscript: [arXiv:1810.02705](https://arxiv.org/abs/1810.02705).  
   用于确认：利用近似曲率预条件 NEB/string 和自适应时间步已有明确先例，并已在经验势与 DFT 上展示。  
   边界：VARNEB 只能把“原子--应变耦合、物理 metric 与 optimizer preconditioner 分离、跨后端验证”的新增组合写作候选贡献。

5. Kolsbjerg, Groves, and Hammer (2016), DOI [10.1063/1.4961868](https://doi.org/10.1063/1.4961868).  
   用于确认：先粗路径、再围绕重要区域增加 image 的 AutoNEB 思想已有先例。

6. Lindgren, Kastlunger, and Peterson (2019), DOI [10.1021/acs.jctc.9b00633](https://doi.org/10.1021/acs.jctc.9b00633); preprint: [arXiv:1906.10257](https://arxiv.org/abs/1906.10257).  
   用于确认：按 image 收敛状态选择性计算和对鞍点区域缩放收敛判据已有先例。

7. Garrido Torres *et al.* (2019), DOI [10.1103/PhysRevLett.122.156001](https://doi.org/10.1103/PhysRevLett.122.156001); preprint: [arXiv:1811.08022](https://arxiv.org/abs/1811.08022).  
   用于确认：GPR surrogate 和不确定性驱动真值调用可显著减少 NEB 函数评估。  
   边界：路线图把在线 GPR 列为 stretch goal，不把多保真 continuation 冒充首次提出的 surrogate NEB。

### 3.3 模式分解

8. Okenyi, Ratcliff, and Walsh (2021), DOI [10.1039/D0CP04236F](https://doi.org/10.1039/D0CP04236F).  
   用于确认：以 Γ 点 phonon basis 逐步逼近相变 MEP、比较低维基与完整路径已有材料先例；模式选择不能简单按最低频率排序。

9. Lima, Oliveira, and Esteves (2026), DOI [10.1016/j.carbon.2026.121271](https://doi.org/10.1016/j.carbon.2026.121271).  
   用于近期边界检查：G-SSNEB 路径与 phonon polarization 的相关分析已经被用于碳同素异形体相变。  
   边界：这进一步说明“把路径投影到声子”本身不足以作为 VARNEB 独立创新；原子--应变广义模态、误差重构和预条件统一才是研究缺口。

## 4. 推理与建议的标注

以下内容是本项目提出的研究假说，不是上述论文的直接结论：

- 使用每原子均方原子位移和 `L0=(Ω/N)^(1/3)` 应变尺度构造 intensive metric；
- 采用 Hencky strain 并在 `G` 下统一切线、投影、弹簧和 CI；
- 严格分离物理路径 metric `G` 与 optimizer preconditioner `P`；
- 构造原子--应变 block Hessian 并用其广义模态同时做预条件和机制分析；
- 用 `≥30%` DFT calculator-call reduction 作为“accelerated”主张的项目内晋级门槛；
- 用多指标 `η_i` 驱动可变胞 adaptive image，并在冻结 image 集上重认证。

这些假说必须经过路线图所列测试和消融，失败时按止损规则降级表述。

## 5. 检索说明

- 检索日期：2026-09-16。
- 优先使用本地作者 PDF、DOI 页面、期刊/机构页面和 arXiv 作者稿。
- 路线图是定向方法设计，不是系统综述；投稿前仍需针对“Riemannian/metric NEB”“variable-cell preconditioning”“phonon-strain coupled path analysis”再做一次完整查新。
- 未把网页摘要中的宣传性加速倍数直接转写为 VARNEB 的预期结果；项目 claim 只由自身预注册基准决定。
