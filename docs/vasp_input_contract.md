# VASP 输入契约：以 GaN 两次真实失败建立回归（2026-09-18）

目标是从输入生成源头预防可识别的错误，不是让每个 image 自动修改参数，
也不是通过强制对称化、去剪切或重新插值修复路径。当前补丁未触发 CI、PyPI
发布。早期契约补丁的cu17验证记录保留如下；当前运行位置已改为hf，新增的
候选全链门禁和六方续算记录见第5节，不重算已完成四方路径或端点。

## 1. 失败证据与边界

节点 cu17，VASP 6.3.2，40 MPI，Ga2N2，45.7 GPa，ENCUT=600 eV，
Gamma 8x8x6，生产 EDIFF=1e-7。29 **总像** = 2 个固定端点 + 27 个内部像。

原始失败结构保存在 `tests/fixtures/vasp_bravais/`，仅含 POSCAR，不含任何
受许可限制的 POTCAR 或 VASP 源码：

- `gan_image15.vasp`：来自 `recovery_isym_minus1/.../15/POSCAR`，ISYM=-1
  但默认 SYMPREC 下仍出现实空间 simple monoclinic / 倒空间 base centered
  orthorhombic 不一致。此前 SYMPREC=1e-8 的完整静态计算通过。
- `gan_image24.vasp`：来自 `recovery_isym_minus1_symprec1e8/.../24/POSCAR`，
  ISYM=-1、SYMPREC=1e-8，仍出现相同类型的拒绝。

两者具有正体积、正常 cell 条件数、合理最短距离；不能把它们归类为重叠原子
或发散结构。几何检查不应误拒绝这些物理有效的输入。

`LATTCHK` 的本机源码检查确认：VASP 分别分类实空间和倒空间并要求其对应关系
一致。因此关闭对称性使用不等于绕过全部晶格检查。

原样结构的启动探测结果（NELM=1，不是收敛 SCF，不提供生产能量）：

| 固定设置（ISYM=-1） | image 15 | image 24 |
|---|---|---|
| SYMPREC=1e-10 | 启动通过 | 失败：base-centered monoclinic / triclinic |
| SYMPREC=1e-12 | 启动通过 | 启动通过 |

这直接否决“只验证 image 15 即宣布整条路径修好”。更小容差不是普遍保证；
1e-12 仅是本机两个回归输入上通过的历史测试值，不再作为通用默认策略。
随后直接读取 `recovery_isym_minus1/.../vcneb.traj` 的最新完整 29-image 链，
固定 ISYM=-1、SYMPREC=1e-12，对端点与全部内部像启动探测：29/29 通过。
探测参数只在电子迭代预算上使用 NELM=1；不能作为完整 SCF 或未来路径保证。
两个原失败结构随后均完成完整 SCF，并通过解析结果后的第二次输入写入核对：
image 15 为 -21.91787247 eV，image 24 为 -21.65146756 eV；这些是单点能量，
不是能垒或路径收敛结果。image 15 与此前 1e-8 静态测试的能量在报告精度上相同。
上述历史策略的本地完整回归为 131 passed。证据保存在
`validation/gan_b4_b1/vasp_contract/{failed_policy_1e10,failed_inputs_1e12,complete_chain_1e12,full_scf_1e12}.json`。
完整 SCF 报告包含测试时实际适配器/策略/探测脚本 SHA256；调整默认策略后，
这些哈希只对应历史测试版本，不能冒充当前版本或新默认值的验证。
官方背景：[ISYM](https://vasp.at/wiki/ISYM)、[SYMPREC](https://vasp.at/wiki/SYMPREC)。

按用户指定，当前默认放宽为 SYMPREC=1e-4，ISYM=-1 保持不变。
该设置统一由后端、VASP 入口、静态验证器及启动探测脚本使用，GaN 续算脚本
也采用同一默认值。它是对称性识别容差，不改变 SCF 或 NEB 力收敛标准。
用户随后授权在 cu17 验证与续算：1e-4 下两个历史失败输入均完成完整 SCF，
并通过重复输入写入核对，能量与 1e-12 测试在报告精度上相同。
新证据为 `validation/gan_b4_b1/vasp_contract/full_scf_1e4.json`；
全链与续算进展见 `validation/gan_b4_b1/cu17_symprec1e4_20260918.md`。
不能将以上 1e-12 的历史测试结果当作新默认的通过证据，默认值调整本身不启动计算。

## 2. 已实现的源头契约

1. `prepare_vasp_static_parameters()` 统一静态设置，显式定义 ISYM=-1、
   SYMPREC=1e-4。允许显式覆盖为 ISYM=0 或其他正有限 SYMPREC，但必须在
   运行开始确定，端点静态结果与路径一致。整数标签不允许小数截断。
2. 每个 calculator 在挂接时冻结实际参数、INCAR/KPOINTS/POTCAR 指纹、
   原子种类和顺序。每次输入写入，以及 executor 缓存命中之前，都检查契约。
   动态修改参数或源文件会在启动外部进程前明确失败。
3. 每次 image 计算前检查原子/晶格有限值、三维周期性、正体积、晶格条件数、
   配置的最短距离（含自身周期副本）。默认仅排除重叠；材料距离门槛由调用方
   明确配置。对 VCA 检查物理位点，不误拒绝有意重合的组分表示。
4. 保留 ASE 的高精度 POSCAR 写入，并读回实际文件核对 cell、坐标和 species
   排序；核对实际 INCAR 静态/对称性标签及 POTCAR 指纹。原子化更新独立的
   `vasp_input_contract.json`，记录结构和实际输入 SHA256，不记录许可数据内容。
5. SCF 未收敛、缺失或非有限 energy/forces/stress 不允许进入优化链或 image
   结果缓存。Bravais 根错误优先于 MPI launcher 包装分类；不进行自动参数重试。

这些检查不等价于 VASP 的完整 Bravais 分类器；报告明确标记
`calculator_free_not_bravais_certification`，不能把预检成功说成 DFT 成功。
不配置候选门禁时，几何门槛越界仍明确停止并保留证据；新增的可选步拒绝接口
及使用条件见第5节，不能把普通几何检查误称为原生分类验证。

## 3. 结果与恢复一致性

- VASP 串行入口 `--image-workers 0` 对应一个执行 worker；每次仅启动一个
  image 的外部 DFT，仍可使用 40 MPI。不是多个 image 并行争用同一节点。
- 串行/并行统一使用持久化 exact-state 缓存和 JSONL manifest。缓存 namespace
  自动绑定有效参数、源输入内容指纹和 VCA 配置，用户标签不能绕过这些检查。
  串行模式逐像立即保存成功结果；一个像失败后不再启动后续像，manifest 明确
  记录未计算的索引。并行模式则回收已提交 worker 的有效结果以避免丢失已完成工作。
- 固定端点 summary 的复用现在检查有效 calculator 参数，而不只是源 INCAR
  哈希；相同源文件但不同命令行 SYMPREC 的旧 summary 不可误复用。
- 超过步数预算而未达到 fmax 的 summary 为 `step_limit_reached`，不是
  `completed`；仅通过默认 0.10 eV/A 或显式指定阈值才标记 completed/converged。
- 缓存与完整链快照是两类状态：失败迭代已完成的单点可以保留，但不能当作
  一次完整优化迭代，也不能拼成虚假的收敛轨迹。

## 4. 真实启动闸门及下一步验收

`scripts/probe_vasp_input_contract.py` 默认只生成检查过的输入，必须 `--run`
才执行真实 VASP。支持两个固定结构，以及 `--trajectory --n-images 29`
直接读取已有最新完整链，不重新插值。目录必须不存在，避免旧日志假阳性。
默认探测采用 NELM=1；通过必须同时满足进程返回 0 与实际开始电子迭代。
显式 `--full-scf` 使用完整源输入的电子迭代预算，要求有限且收敛的能量、力和
应力，并再次写入输入以验证解析结果之后没有破坏参数锁定。
输出永远标记 `production_cache_eligible=false`。

```bash
python scripts/probe_vasp_input_contract.py \
  --source-dir LICENSED_CASE/vcneb_input/initial \
  --structures tests/fixtures/vasp_bravais/gan_image15.vasp \
               tests/fixtures/vasp_bravais/gan_image24.vasp \
  --workdir NEW_PROBE_DIR --ncores 40 --vasp-bin VASP_BIN --run
```

补丁交付不等于 GaN 计算完成。验收顺序：两次历史失败输入 → 最新完整链全像
初始化 → 代表性输入完整 SCF 和契约重复写入检查 → 固定策略续算 → 最终收敛与
结果审计。不得再从单点启动成功推断整个未来路径已经无故障。

## 5. hf六方真实失败驱动的候选步门禁（2026-09-18）

**2026-09-19限制更新：** 27721448再次失败，旧probe对image7误判通过。
同源包装器、同许可对象、同输入及SYMPREC下，ifort `-O0 -march=core-avx2`
和通用`-O2`均返回4/4/4，而安装VASP所用组合`-O2 -march=core-avx2`返回4/11/4，
匹配实际失败。旧构建命令没有完整归档，不能假定是-O0，更不能归因仅缺-O2。
新增`scripts/build_hf_vasp_lattice_probe.sh`保存编译组合及源码/对象/binary哈希。
编译目标、ABI和对应VASP的真实回归均属诊断工具契约，不能只检查可执行文件存在。
新匹配组合仍无普适可靠性保证；六方0.919平台没有解决，没有重投生产。
详见[06:00接力证据](../validation/heartbeat_20260919_0600.md)。

### 5.1 manager正确、POSCAR末位舍入错误（09:29根因收敛）

继续逐字节比较发现，最后完整链的manager image7在匹配VASP构建的probe下为
`4/4/4`；ASE原固定小数POSCAR将晶格分量改变约$5\times10^{-16}$ Å后变为
`4/11/4`，实际VASP因此拒绝。manager基面角为120.006616°；其约0.007°剪切是
路径自由度，不能投影回120°。原先考虑的“六方规范化”据此否决，没有写入代码或路径。

`ExplicitPotcarVasp`现在只重写ASE已经生成并验证的POSCAR晶格三行，使用17位
有效十进制数字，要求每个binary64分量逐位往返相同；unit scale、旧写入与manager
相差超过1e-12 Å、非有限值或round-trip失败均硬退出。原子分数坐标、cell、物种顺序
和物理参数不变。`NativeVaspCandidateValidator`使用同一序列化函数预览，不再模拟
ASE旧固定小数格式。写入契约记录policy/digits和最终POSCAR哈希。

真实完整SCF canary 27723185从旧轨迹直接读取最后完整manager image7，不经中间
POSCAR；1节点32 MPI，VASP6.3.2/PBE/Ga_d+N/600eV/Γ8×8×6/ISYM=-1/
SYMPREC=1e-4。作业COMPLETED/0，SCF收敛，energy=-23.03083126 eV，最大原子力
0.155086 eV/Å，energy/forces/stress有限，重复写入通过。它不是NEB收敛结果。

同一完整链fresh FIRE无DFT重放：当前image7为4/4/4，首个完整候选的基面角
120.009445°、11/11/11，全链几何门禁通过。旧作业反复缩步使dt塌缩，后续大量
几何微步落到计算器容差内并复用内存结果，不能靠增加步数恢复。续算必须fresh FIRE，
缓存namespace绑定新写入源码/VASP binary，不能混用旧序列化缓存。

六方生产27721013在image10失败，实际ISYM=-1、SYMPREC=1e-4生效。
VASP6.3.2的POSCAR读取流程调用独立的晶格一致性检查，不受禁用对称性使用控制。
有效晶胞存在120.007104°的真实畸变；不把它投影回理想角度，不继续更改SYMPREC。

独立诊断包装器 `scripts/vasp_lattice_probe.f90` 链接用户已有且获许可的
VASP `lattlib.o`，不分发该对象或生成的二进制，也不改动VASP程序。
这是特定安装/Intel-classic ABI的可选接口，不是PyPI内置的VASP复制实现。
原生诊断须先对照对应VASP可执行文件验证：本例58个初始胞、3个失败候选
与实际初始化结果一致；半步image10随后完整SCF通过（27721349）。

`NativeVaspCandidateValidator` 在调用DFT前检查整条候选链的几何，以及
管理器精度和当前ASE POSCAR写回精度下的晶格分类。有效原子序、位置、晶胞
不能被验证器修改，原生工具失败/输出损坏/执行文件变化明确退出。

核心 `VCNEB.set_x()` 先在无calculator的副本上生成全链候选，全部通过后
一次性提交所有中间像；拒绝不能部分更新链、清除完整评估缓存或派发DFT。
端点保持不动。可选 `CheckedFIRE` 只处理 `CandidateStepRejected`，沿原位移
按固定因子有界缩短，接受后重置动量、减小dt；超过预算明确退出并保留记录。
SCF、程序、文件和其他错误不属于这个接口，不自动调整任何物理参数。
它是优化器的可计算性步控制，不是路径对称投影，也不保证更快收敛。

```bash
python examples/run_vcneb_vasp.py ... \
  --native-lattice-probe LICENSED_VALIDATED_LATTICE_PROBE \
  --optimizer FIRE --candidate-step-retries 8 --candidate-step-retry-factor 0.5
```

不提供原生probe时，默认候选重试为0；其他优化器不能启用该FIRE缩步策略。
`candidate_step_manifest.jsonl` 记录拒绝/接受、像索引、类别、位移hash及步长分数，
最终summary及失败summary均记录策略/历史。完整SCF缓存仍绑定原物理契约。

本例本机全套187项回归通过；hf无pytest，未安装或升级运行环境，改用真实失败
全链端到端回放，ASE3.23.1b1下原子性拒绝、半步接受、动量重置通过。
回放的输入力仅为人工delta/dt²，不是物理NEB力或路径结果；接受结构与已通过
SCF的半步canary差约1e-15Å。独立源码快照续算任务为27721448，旧B3任务
27721015运行快照不变。提交不等于收敛；最终仍需SCF/机制/能垒/轨迹审计。

证据目录：`validation/gan_qian_suite/failure_hex_image10/`，部署/缓存/提交记录：
`validation/gan_qian_suite/guarded_launch_20260918.json`。

## 6. 缓存命中的结果快照完整性（2026-09-18）

独立无DFT回归发现：executor从persistent cache返回有效结果时，新挂接的
live计算器可能仍没有results。直接交给ASE轨迹写入器会保存完整几何却漏写
已计算的energy/forces/stress；优化器使用的executor结果本身不因此改变。

核心现在以完整executor评估创建独立SinglePoint快照，保存实际energy、
笛卡尔forces、stress；不替换live计算器、不更改位置/晶胞、不触发额外计算。
评估不完整或几何在评估后变化，拒绝生成结果快照，不借用过期结果。
普通串行计算器保留原有ASE写入流程。公共几何目录导出仍允许未计算的几何，
不能将此类导出误标为完整电子计算。

新增3项回归检查：persistent命中保存全部字段且中间像不重算；缺失评估不启动
计算；外部几何修改不能保存旧结果。全套190项通过。当前修正只在本机，未
覆盖已提交的27721015/27721448源码；真实B3 step1全部29帧结果字段完整，
证据为 `validation/gan_qian_suite/snapshot_fields_20260918_2244.json`。

最终分析必须检查实际快照字段，不能仅因几何帧数正确就宣布可审计。
旧产物如确有缓存字段缺口，只能在匹配物理namespace与精确状态的前提下，
用原有效单点数据生成独立审计副本并保存出处；无匹配结果就明确失败，
不得修改路径、借邻近帧、插值能量/力，或把部分失败轮当完整收敛链。

### 6.1 实际续算产物的独立精确状态归档

2026-09-18六方27721448已运行，初始27像全部缓存命中，实际旧版step0快照
漏存results。本机快照修正没有热替换入这个活动运行包。独立归档工具
`scripts/restore_exact_cached_snapshot.py`要求原始完整已评估链作为交叉佐证：
两份preflight计算参数/输入哈希/用户namespace一致，缓存namespace一致；
每像cell/positions/numbers/pbc字节相同，精确cache key存在且字段完整有限；
缓存结果与原完整链能量/力/应力逐项相同。全链验证后才生成独立新traj/json，
原件保留，记录各数据来源和哈希，不调用DFT或改变几何，不证明路径收敛。

实际hf运行成功，独立副本29像结果完整。第一轮新SCF后的原始step1快照也已
实查全部结果完整，无须归档填补。6项新回归/全套202项通过；不可将这项
存储层操作描述为路径“自修复”或把字段存在当作精确鞍点证据。

## 7. 端点相身份和初猜周期绕行（2026-09-19）

CdSe非标准8原子超胞的FixSymmetry不自动施加完整母相晶格自由度约束。
其BFGS结束和小广义力不等于严格#225/#186：逐容差相检查失败必须保留，不删除
检查或扩大容差来获取期望标签。原始端点/summary保留并记录哈希。

`prepare_cdse_phase_candidates.py`仅对有来源、已完成的近相端点构造有界相候选；
主伸长改变<=1e-3、motif位移<=1e-3 Å，超出拒绝。结构重建发生在端点定义阶段，
从不作用于生产NEB像。候选不是新的收敛结果：检查作业必须重新完成静态SCF，
`verify_cdse_endpoint_statics.py`使用不带FixSymmetry的SinglePoint副本，检查原子力
和FrechetCellFilter virial广义力<=0.02。无需第二次DFT来执行该门禁；不通过则
不产生成功gate。生产同时要求afterok和成功gate；不重复原BFGS。
stationarity并不代替Hessian/局部稳定性分析。

参考motif还规定周期绕行分支。atomic初猜存在半胞位移，不能对带微小原点噪声的
端点再次np.rint选择最近周期代表：实际会从零绕行切换，造成1.674 Å碰撞。
显式零绕行使同一端点初猜最近距离2.536 Å，不降低1.8 Å门禁，不换原子映射。
新增真实结构/微小平移回归；计算参数和已运行源码不随修复改变。
