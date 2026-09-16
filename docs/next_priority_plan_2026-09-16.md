# VARNEB 下一阶段优先级清单（2026-09-16）

## 本轮结论

当前核心已经通过 ASE 的能量--原子力--应力契约实现了**架构上的**计算器无关性，
但实证仍不完整：ABACUS 有 BTO/HfO2 完整路径，VASP 只有单 image smoke，QE 只有
无 DFT preflight。因此论文不能把 ``calculator-agnostic`` 写成跨程序的材料验证结论。

现有 ``Mode``/路径投影功能也不是声子分析。它能读入用户给定的模式、生成模式引导
初始路径及投影路径；它不计算力常数、频率、虚频、Hessian 或沿路径的正规坐标分解。
下列任务把这两个缺口作为主线。此次文档更新只定义工作和验收标准，不启动 DFT 作业。

## 优先级与依赖

| 优先级 | 工作包 | 必须产物 | 后续依赖 |
|---|---|---|---|
| P0 | 论文理论重构与术语冻结 | 方程、符号表、claim 边界 | P1--P4 |
| P1 | 多后端真实能力 | QE adapter、VASP/QE BTO 正向 7-total-image 证据 | 论文“calculator-agnostic”主张 |
| P2 | 路径--声子模式分析 | 通用 force-constant 读取、投影、BTO 示例 | 机制图和声子结论 |
| P3 | 结果与论文整合 | 三后端 BTO 图、HfO2 主案例、模式图 | 投稿草稿 |
| P4 | 可选案例增强 | 一个有明确问题的补充案例 | 补充材料，不阻塞初稿 |

P0 应先完成；P1 的代码和 P2 的 calculator-free 分析实现可以并行开发。只有 P1
后端 preflight 全绿，才能提交相应的 BTO 生产路径。P3 只接收通过各自验收门槛的
数据，不把 smoke 或中途轨迹当成结果。

## P0：把论文的理论链条补完整

目标不是宣称新 VCNEB 理论，而是从普通 NEB 严格、自然地导出当前实现的变胞版本。

- [ ] 重排 CPC 稿为：1 Introduction；2 Conventional NEB；3 Calculator-agnostic
  VCNEB；4 Software workflow and reproducibility；5 Verification and material cases；
  6 Limitations and conclusions。
- [ ] 第 2 节先定义固定 cell 的构型向量、improved tangent、真力法向投影、弹簧
  切向分量及 CI 替换，并说明端点固定、image 数包含端点、普通 NEB 先于 CI 的原因。
- [ ] 第 3 节再将 ``R`` 扩展为 ``Q=(s,F)``，给出 ``r=s h``、
  ``h=h0 F^T``、共同度量、焓 ``E+PV`` 和由应力得到 cell 广义力的推导；明确
  ``cell_scale`` 是数值度量而非材料参数。
- [ ] 第 3 节单列 calculator contract：每个 image 只要求有限的 energy、Cartesian
  forces、stress 与独立工作目录；核心不含任何 ABACUS/VASP/QE 输入语义。说明
  单位、stress 符号和静态 image 计算由后端 adapter 负责验证。
- [ ] 在正文/补充材料中严格区分：(i) 普通 NEB 的默认 ``0.10 eV/A`` 验收，
  (ii) 显式较严研究协议，(iii) endpoint relaxation 阈值；不得把三者混为同一精度。
- [ ] 加入一张从固定 cell ``R`` 到扩展 ``(s,F)`` 的小示意图和一个符号/单位表，
  替换空泛的软件叙述。

**P0 验收：** `docs/theory.md`、代码和论文的坐标、应力符号、image 计数与默认
阈值逐项一致；所有“universal/calculator-agnostic”文字均有明确的证据等级。

## P1：把计算器无关性变成可运行、可验证的多后端能力

### P1.1 先完成 QE adapter 与回归测试

- [x] 已增加 `vcneb/qe.py`，以当前 ASE `Espresso`/`EspressoProfile` API 建立每 image
  独立目录的 calculator factory；输入固定为静态 `pw.x` 计算，不能使用 QE 自己的
  `vc-relax` 更新 image。
- [x] 已增加 `examples/run_vcneb_qe.py` 及零 DFT `--validate-only` preflight。公开参数包括 profile/
  command、pseudo directory、pseudopotential map、`ecutwfc/ecutrho`、k mesh、SCF
  阈值、smearing 和 image worker 设置；`ibrav=0`、`tstress=.true.`、
  `tprnfor=.true.` 必须由 validator 检查。
- [x] QE preflight 还要求 Ba/Ti/O UPF 均在指定目录内、元素 metadata 与 mapping 一致、
  有可见 PBE 标记，并把每个文件的 SHA256 写入 manifest；它不替代 QE 自身 cutoff/
  赝势收敛研究。
- [x] 已为 QE factory、目录隔离、stress capability 和 CLI 默认值增加
  mock/fake-calculator 回归，不依赖本地或 CI 中的 DFT 可执行文件。
- [x] 已为现有 VASP factory 增加静态-image validator/test：`IBRION=-1`、`NSW=0`、
  `ISIF=2`、`ISYM=0`、独立目录和 stress 读取都必须由测试覆盖。
- [x] 更新 capability 表、README 和手册：ABACUS/VASP/QE 都是已实现 adapter；只有
  完成对应真实路径后才标为 material-path validated。

### P1.2 统一的 BTO 三后端验证（只做正向 7 total images）

冻结同一份已验证的 BTO T-to-C 端点、原子 mapping、log-strain 初始路径、
`n_images=7`（两个固定端点加五个 interior images）、普通 NEB、`fmax=0.10 eV/A`。
不重跑现有 ABACUS BTO 参考路径；只新增 VASP 和 QE 两个后端。

- [ ] 每个新后端先做一次单-image static preflight：energy、forces、stress 有限，
  结构和原子顺序未变，输入确为静态 image 计算。
- [ ] VASP 与 QE 各运行一条正向 7-total-image BTO 普通 VCNEB；控制器缓存端点，
  每轮只计算 5 个 interior images。每个 interior worker 使用 32 MPI ranks，因此
  worker 池为 `5 x 32 = 160` ranks；提交前仍按 `sinfo/squeue` 和实际程序并行效率
  确认资源，不启动 CI 或本地 DFT。
- [ ] QE 采用经自身 cutoff/赝势收敛验证的 PBE 参数（100 Ry 仅在其赝势适用时保留）；
  VASP `ENCUT` 取 POTCAR 推荐值及独立收敛检查，不能把 ABACUS 的 100 Ry 机械换算。
  当前 hfacnormal01 已确认 QE 7.0 `pw.x` 可用，但其随包 pseudo 目录没有 Ba/Ti/O
  PBE UPF；VASP module 与此前已知 VASP 可执行路径亦不可见。因此实际 VASP/QE BTO
  作业尚未提交，等待已授权的 VASP 安装位置和经审查的同一套 Ba/Ti/O PBE UPF。
- [ ] 统一比较路径单调性、最高 image 身份、cell/volume 演化、reaction energy、
  最大广义力、最小距离和 wall time。不同赝势/实现的绝对能量不是逐 meV 对齐要求。

**P1 验收：** 三个后端均留下可审计 manifest；VASP/QE BTO 路径在相同物理设置下
无内部势垒、几何有效、达到 `0.10 eV/A`；论文才可称“在 ABACUS、VASP 和 QE 上
演示了同一 calculator contract 的端到端路径”。若任一后端失败，只报告 adapter/
smoke 状态，不提升该主张。

## P2：路径上的声子模式拆解（本阶段的科学重点）

### 先冻结物理语义

- [ ] 将“路径模式分解”和“声子”分开命名。非平衡 NEB image 不是驻点，不能把其
  位移投影直接叫作 phonon；它应称为相对参考结构的**正规坐标/集体模式投影**。
- [x] 对平衡端点（及经充分验证的 CI 鞍点）计算或导入 Gamma 点力常数，构造
  mass-weighted dynamical matrix、频率和本征矢。虚频只在该驻点 Hessian 的明确
  约束空间内解释。
- [~] VCNEB 的 homogeneous cell 变化单独以对称应变分量/strain amplitude 报告。
  原子 Gamma 声子与原子--cell 联合 Hessian 不应混称；后者需要外压下的完整
  广义 Hessian，是后续扩展而非首版交付。

### P2.1 通用、calculator-free 分析层

- [x] 增加 `vcneb/phonons.py`（或等价 analysis module）：读取标准化 `.npz` 力常数
  / eigenpairs，执行质量加权对角化、声学平移投影、频率单位/符号处理、模态正交性
  与 eigenvector phase 对齐。
- [x] 定义跨 cell 路径坐标：以一个固定参考 cell、已审计 mapping 与周期 translation
  gauge 表示原子位移，再计算
  `q_nu(lambda) = e_nu^T M^(1/2) u(lambda)`；同时输出投影残差、累积解释方差和路径
  切线与 mode 的重叠，不能只画单个漂亮的 mode。
- [~] 已有质量加权、平移投影、已知线性组合、映射与 gauge 合成测试；模式符号
  翻转不影响振幅、mapping/gauge 变换不改变投影、已知线性组合能被精确恢复。
- [~] 已提供 `examples/analyze_path_gamma_modes.py`，导出 JSON/NPZ、完整基重构残差和
  简并子空间切线重叠；CSV、模式动画与
  可编辑图仍待真实 BTO 数据后按图稿需要补充。把 `Mode` 的现有输入格式
  与真正的 phonon provenance（结构、calculator、supercell、displacement、hash）
  明确区分。

### P2.2 BTO 首个真实示例

- [x] 已提供 ABACUS BTO 端点的最小 Gamma 点有限位移/力常数 workflow。预检作业
  `27699613` 已验证只生成 6 个位移目录且未运行 ABACUS；`27704219` 已以 32 MPI
  对 cubic BTO 的 6 个位移顺序完成静态 SCF，并由 Phonopy 组装 `FORCE_SETS` 和
  `FORCE_CONSTANTS`。每个
  displacement 均通过静态 force preflight 并独立记录。可使用 ASE/Phonopy 的有限
  位移执行，但 VARNEB 分析层只读取标准化结果，不绑定某一声子程序。
- [~] cubic BTO 的三重简并不稳定 $Gamma$ 子空间（-217.474 cm$^{-1}$）已投影到
  已完成 ABACUS T-to-C 7-image 路径；原子 mapping、非整数 gauge translation、刚性
  平移去除与完整基残差均有记录。tetragonal 端点力常数、Ti--O/strain/energy 合图仍待补。
- [x] 只有在残差和 mode-overlap 支持时，才写“路径主要由某软模主导”；否则报告
  多模混合，而不是强行归因。VASP/QE 的 BTO 路径可复用该分析格式，但无需先做三套
  昂贵声子计算。

**P2 验收：** 能从一个可追溯的 BTO force-constant 输入重现频率、特征向量与全部
路径投影；论文把它作为机制分析，不把 T-to-C 的无势垒路径误称作 CI 鞍点证明。

## P3：把证据变成八页以内的论文

- [ ] 主图控制为四张以内：固定 NEB 到 VCNEB 的概念/方程图；manager--worker 与
  three-backend contract 图；BTO 三后端正向路径及软模/strain 投影；HfO2 barriered
  VCNEB 与文献限定比较。冗余 image/收敛表进补充材料。
- [ ] BTO 的定位：低成本、同一 7-image 正向路径的后端互操作性与无势垒 CI gate。
  HfO2 的定位：12 原子、真实变胞、有内部势垒的主难例。两者不能互相代替。
- [ ] 更新 abstract、program summary、capability table 和 conclusions：跨后端结果
  只陈述已经实际完成的范围；声子结论只陈述 BTO 数据真正支持的层级。
- [ ] 增加完整 claim-to-evidence 表、后端参数表和 source-data manifest；所有图均由
  提交的脚本从 summary/phonon artifact 重建。

## P4：可选案例增强（不阻塞初稿）

首选不是匆忙加入第三个昂贵相变，而是在 BTO 内补一个**固定 cell 的
`+P -> -P` 极化翻转**：它有真正内部鞍点，能展示第 2 节普通 NEB 与 VARNEB 固定-cell
极限/ASE-NEB 对照。开始前须确认相同端点定义、cell 约束和文献可比性。

若 P4 需要一个真正新的变胞体系，再做一个小型探索性筛选（例如 Si 的受压
diamond-to-beta-Sn），只在 mapping、压力边界、端点和计算成本均通过预检后立项。
不得为了“案例数量”重复 BTO/HfO2 或把未经收敛的路径写进正文。

## 推荐执行顺序

1. 完成 P0 的理论大纲、符号表和文字边界；同步修订文档中的阈值一致性。
2. 实现并测试 QE adapter；补 VASP BTO static validator；两者均先作单-image preflight。
3. 仅在 preflight 成功后，按 5 个 32-MPI interior workers 分别提交 VASP 与 QE 的
   BTO 正向 7-total-image 普通 VCNEB。
4. 并行实现 P2.1 的 calculator-free 模式分析与合成测试；随后只用 ABACUS BTO 做
   首个真实 Gamma-mode 数据集。
5. 用通过门槛的多后端路径和 BTO 模式结果重写 CPC 稿，再决定 P4 是否值得投入。
