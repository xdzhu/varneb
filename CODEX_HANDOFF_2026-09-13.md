# Codex 线程交接说明（2026-09-13）

> **2026-09-19 09:42状态覆盖：** 六方根因已定位为ASE默认POSCAR末位舍入：
> manager image7原始binary64晶格是4/4/4，约5e-16 Å写盘差异才触发4/11/4；
> 未对真实120.006616度剪切做对称投影，也未改SYMPREC。通用VASP入口现以17有效
> 数字写晶格并要求逐位往返，preflight使用同一序列化。完整SCF canary 27723185
> COMPLETED/0。fresh-FIRE六方27723210在hfacnormal01 RUNNING；首轮27/27内部像
> 均新算、9 worker、每像一次成功、366.24 s，无Bravais错误，step0=0.919209；
> 29像Ga2N2的energy/forces/stress均有限。第一个非零位移step1也已真实重算27/27
> 新几何、每像一次成功、364.95 s，fmax降至0.757250；已跨过旧窄分类带，但尚未
> 达到0.10，不能声称路径收敛。B3 27721015 RUNNING，09:39 step37=0.196777，保持不动。
> 本机209项回归和diff-check通过，无CI/发布/推送。源包、作业和证据路径见
> docs/GAN_QIAN_2013_REPLICATION_PLAN.md及validation/gan_qian_suite/
> poscar17_recovery_20260919；目标未完成，不重复提交。
> 09:49四个生产均RUNNING、stderr空。按实测轮时和近期下降，仅作调度估计：CdSe两路
> 约10:30--11:30、B3约10:30--12:00、六方约11:00--13:00可能达到0.10；FIRE回弹
> 或后期变慢会延后，不能据此宣称完成。当前只剩活跃计算，按用户约定暂停目标等待接力。

> **2026-09-19 06:15 单次定时接力状态覆盖：** 已实际核验hf/sacct，不是提醒。
> B3 27721015仍RUNNING，FIRE25残差0.466259、单轮约18分钟；不因step24回弹停止。
> 六方27721448于05:23:39 FAILED：image7晶格分类错误；step26–102平台约0.919209，
> 不是0.10收敛。同源probe的-O0+AVX2/-O2通用均漏检，-O2+core-avx2匹配实际4/11/4；
> 不是简单“缺-O2”，旧构建flags没有完整记录。未重交未解决平台的六方生产。
> CdSe原端点27721577/78优化完成；79/81相身份检查FAILED，80/82依赖取消。
> 非标准超胞FixSymmetry不是完整母相晶格约束，另建有界相端点候选、保留原件；
> 必须完整静态SCF及未经对称投影的force/virial<=0.02才能进入生产，未重复BFGS。
> cell检查27722588 RUNNING(node124)→生产27722589 Dependency；atomic旧检查
> 27722590因半胞舍入绕行碰撞失败，91依赖取消。明确参考motif零绕行后，
> 独立atomic检查27722599 RUNNING(node363)→生产27722600 Dependency。
> 两路17像预检通过；cell已通过image0/1真实初始化，但完整SCF/路径仍待审核。
> v2/v3各独立目录，不覆盖B3/guard/cell运行源码；本机208项回归通过，无CI。
> 详情/哈希/耗时估计见validation/heartbeat_20260919_0600.md；没有新监控/发布/推送。
> 目标未完成，保持此前暂停状态；下次先核验上述四个新有效ID，不重复原六个提交。

> **2026-09-18 23:15左右实时覆盖：** 六方续算27721448已于23:04:54启动，
> node[39,478,491]/3节点291tasks，CheckedFIRE0→1残差4.505321→3.250342，
> 第一轮新SCF后chain_step_0001全29帧/4原子结果齐全有限。B3 27721015保持运行，
> FIRE2 fmax3.334191，step2全29帧/8原子结果齐全，未达到0.10。
> 续算step0真实出现27个缓存命中但原始快照缺results，未改生产源码/路径：
> 使用独立archive-audit目录工具，以精确cache key/匹配physics namespace，
> 对照原始hex step1逐项相同结果生成独立审计副本；29像几何字节一致、0DFT，
> 原始轨迹和core ca4d...保持不变。下载文件hex_resume_chain_step_0000_full_results
> 的traj/json包含全部来源哈希；不是收敛证明。202项本机回归通过。
> CdSe27721577/78仍RUNNING并持续电子迭代，四个后续任务Dependency，无失败证据。
> 研究目标仍active，不重复提交/不提前CI，不因为旧Priority条目去重启六方。

> **2026-09-18 CdSe 独立研究范围扩展：** 用户要求不要蹲守 GaN B3；保持27721015
> 正常运行、27721448排队，同时开展 Sheppard2012 CdSe RS→WZ 双映射。
> hf/hfacnormal01 已提交 RS/WZ BFGS 27721577/78 → 两路检查27721579/81 →
> 两路生产27721580/82。23:07核验两端RUNNING于node111/node122，真实VASP step各32MPI，
> 电子DAV迭代已开始；后续Dependency。提交/SCF启动不等于收敛。
> 独立不可变目录 `20260918-varneb-cdse-sheppard`，新目录部署最新core b65090...，
> 不覆盖旧GaN源码；196本机测试通过，实际hf两条种子17像原生/舍入晶格检查通过。
> Cd4Se4/8原子/4FU、0GPa、PBE/PAW455eV/MP10³；15中间像，无CI/fmax0.10；
> 每条生产2节点194tasks/6×32MPI，不派端点worker。原文PW91 Fig11小峰2.4meV/atom
> 与前文经验势区分，不能以0.10单独证明meV能垒精度或提前确认重建候选机制。
> 详见docs/CDSE_SHEPPARD_2012_REPLICATION_PLAN.md及validation/cdse_sheppard_2012。
> GaN旧目标与新增CdSe研究均未完成；下一轮分别检查记录中的真实ID，不重复提交。

> **2026-09-18 22:44状态覆盖：** B3 27721015仍RUNNING，FIRE1残差4.348028，
> 由5.530889下降；未达到0.10。实际chain_step_0001.traj全29帧/8原子都有有效
> energy/forces/stress，轨迹和字段审计已下载。六方27721448仍PENDING/Priority，
> squeue估计9月19日12:25启动（可变估计，不是保证）。保持作业，不重交。
> 本机新增缓存快照修正：从完整executor评估生成独立SinglePoint快照，不替换
> live calculator、不做DFT；过期几何拒绝，3项新回归，全套190项通过。
> **该IO修正尚未部署到已提交的运行/排队源码**：本机core SHA
> b65090fcd062f07adf3f3fa80f114000bfb4829592205771708c89702d76db31，
> guard运行包core仍ca4d3a2492dd933c290dc7038276ef7ceae6490b05195d6a79132eca10642547。
> 后续必要时只以精确状态/匹配物理namespace缓存生成独立审计副本，不能更改路径
> 或用邻近状态填补数据；目前B3该完整链没有缺口。目标仍active，无CI。

> **2026-09-18 22:31左右状态覆盖：** 已实施原子性全链候选门禁与可选CheckedFIRE，
> 仅在DFT前对明确候选拒绝有界缩步，接受短步后动量归零、dt降低；不修改物理契约。
> 本机全套187项通过；hf没有pytest，未安装/升级环境，真实失败链回放在ASE3.23.1b1
> 通过，接受半步与SCF canary几何差约1e-15Å。续算预检/端点缓存校验通过。
> 六方续算27721448已提交，最后查询PENDING/Priority；独立代码目录
> `/public/home/iai806/abacus/agent-runs/20260918-varneb-gan-lattice-guard`，不要覆盖运行源码。
> 3节点291核/9×32MPI，29总像，298剩余步，无CI；从27721013最后完整链继续，
> FIRE重新初始化，缓存复制到新目录，不导入诊断半步或不完整失败轮。
> 原代码包SHA695f2ea41e6812cb37a06f23c0bc8ac75ac25ce1c26aca7383ed8a517c0aaadd。
> 27721015仍RUNNING，完整初始链FIRE0 fmax=5.530889，后续SCF进行中。
> 提交和门禁通过不等于六方/B3收敛；目标仍active。作业记录guarded_launch_20260918.json
> 防重复提交，完整证据/适用边界见docs/vasp_input_contract.md第5节。

> **2026-09-18 22:13左右状态覆盖：** 27721349 COMPLETED/exit0/04:50，参考态与
> 半步image10均converged_static_passed、重复写入契约通过。参考重算energy/forces/stress
> 与原结果差均为0（输出精度内）；记录已下载。不准误读旧RUNNING/Priority条目。
> 半步只证实此次候选可计算，不能宣布整条六方路径收敛；控制器前置门禁/有界缩步
> 尚未实施，下一步实现并测试、用独立源码快照续算，不覆盖B3运行目录。
> 27721015最后查询仍RUNNING/06:26，目标仍active；保持既定物理契约和无CI。

> **2026-09-18 22:11左右状态覆盖：** 半步完整SCF验证27721349已RUNNING于node124，
> 参考态已 `converged_static_passed`，半步态尚在计算；此前Priority条目仅为历史状态。
> B3生产27721015仍RUNNING；本轮没有覆盖其运行源码、改参数或自动重算。

> **2026-09-18 22:07 CST 状态覆盖：** B3检查27721014 COMPLETED/exit0/28:37，
> 全29初始像与两端静态SCF通过；生产27721015已RUNNING于node[44-45,102]。
> scontrol逐项确认9个VASP step×32 MPI，各节点3个，不派发端点。
> 本轮得到原生晶格前置检查证据：调用已许可lattlib.o的独立诊断包装器，
> 不修改VASP/结构/参数；准确重现image10及两个等价换基失败，58个初始胞与
> 两条完整已评估链均一致。失败整步只有image10不一致，1/2、1/4、1/8步均全链一致。
> 这只是晶格分类验证，不是SCF成功；完整静态canary27721349已提交，最后为Priority等待。
> 不重复提交、不把该候选直接当生产断点。原生工具尚未接入生产控制器。
> 证据与下步要求见 `validation/gan_qian_suite/failure_hex_image10/native_lattice_preflight_20260918.md`。

> **2026-09-18 六方失败状态覆盖（晚于21:42条目）：** 27721013 FAILED/exit1，
> 用时17:43；image10在下一轮静态初始化发生实/倒空间Bravais分类不一致。
> SYMPREC=1e-4、ISYM=-1实际生效；保留两条完整已评估链，不是因受力回弹停止。
> 等价整数换基诊断27721271已完成报告，但三个输入全部初始化失败，绝非修复成功。
> 证据见 `validation/gan_qian_suite/failure_hex_image10/diagnosis.md`。
> 27721014最后查询仍RUNNING，已通过初始像0–23；27721015仍Dependency。
> 用户强调三次连续转变应分为三段：各段独立起点焓、相身份与鞍点审计；
> 原文单次长链优化与后续独立三段验证分开记录。目标仍active，不重复提交。

> **2026-09-18 21:42 查询路径更正：** 六方首条全链FIRE step0已于21:37完成，
> 初始fmax=5.969264 eV/Å，后续优化正常。实际优化日志为 `vcneb.opt.log`，
> 不是此前误查的 `vcneb.log`；不要由错误文件名推断无进展。目标仍active。

> **2026-09-18 21:39 状态覆盖：** B3/B1端点27721010/11均已COMPLETED/exit0，
> BFGS第2步收敛，8原子常规胞、45.0 GPa、#216/#225保持，数据已归档。
> 27721014已在node122启动，B3显式路径构造/preflight通过；27721015仍afterok等待。
> 六方生产27721013仍RUNNING，实际9×32 MPI/3节点，没有中断或参数修改。
> 端点证据见 `validation/gan_qian_suite/endpoints_45gpa_audit.md`，目标仍active。

> **2026-09-18 21:32 状态覆盖：** 六方检查27721012已 COMPLETED/exit0，生产
> 27721013已 RUNNING，node[364-366]，实际9个VASP step×32 MPI、每节点3个。
> B3/B1端点27721010/11仍RUNNING，B1 step1 fmax=0.037585>0.02；
> B3检查/生产27721014/15仍依赖等待。证据见
> `validation/gan_qian_suite/hexagonal_worker_launch_20260918.md`，目标仍active。

> **2026-09-18 21:21 GaN 新目标：** 保留四方 27719610 成功结果（0.33849 eV/GaN），
> 继续六方/B3 文献对照。已提交端点 27721010/27721011（RUNNING）、六方检查
> 27721012（RUNNING）；六方 VCNEB 27721013 afterok:27721012，B3 检查
> 27721014 afterok:27721010:27721011，B3 VCNEB 27721015 afterok:27721014。
> 新远端目录 `/public/home/iai806/abacus/agent-runs/20260918-varneb-gan-suite`，
> 仅 hf/hfacnormal01，9×32 MPI/3 节点/291 核生产，29 总像/27 中间像，无 CI。
> 端点优化尚未完成；目标 active。不要重复运行提交入口（launch.json 防重复）。
> 方案及原文三峰/三段区别见 `docs/GAN_QIAN_2013_REPLICATION_PLAN.md`；
> 作业、种子、指纹和当前状态见 `validation/gan_qian_suite/launch_20260918.md`。

> **2026-09-18 执行位置覆盖更新：** 用户要求停止 235/cu17 当前 GaN VASP 计算，
> 后续 VASP 测试、验证与生产使用合肥 `ssh hf`、Slurm `hfacnormal01`。
> 加载 `source /public/home/iai806/Software/VASP/env.sh 6.3.2`；PBE 数据集
> 位于 `$VASP_PSEUDO_ROOT/PBE/数据集/POTCAR`。不无声替换 Ga_d/Ga 等数据集。
> 当前默认 SYMPREC=1e-4，NEB fmax=0.10 eV/A；无 CI。迁移策略与状态见
> `docs/HF_VASP_WORKFLOW.md` 和 `validation/gan_b4_b1/hf_migration_20260918.md`。

## 项目与目标

- 项目目录：`D:\Work\Code\vasp_neb`
- 项目正式名称：`VARNEB`，上游仓库为 `https://github.com/xdzhu/varneb`；算法仍称 VC-NEB，Python 导入包暂保留为 `vcneb` 以维持兼容。
- 总目标：按 `VCNEB_PROJECT_PLAN.md` 推进至完整可用，包括理论推导、核心算法与模式/方向约束升级、VASP/ABACUS 计算器适配、集群典型晶体相变验证、与文献/参考实现对比、CPC 论文草稿及可复现发布材料。
- 长时间计算约束：BTO 主线统一使用合肥 `hfacnormal01`，每个 image worker 使用 32 MPI；每次提交前检查 `sinfo`/`squeue`，避免抢占用户任务。旧 `cu*` 记录仅作历史基线。
- 发布状态：版本 `0.0.1` 已通过 Trusted Publishing 发布；提交 `35c0dbb` 与标签 `v0.0.1` 已推送，重跑后的 workflow `34763435900` 为 success，PyPI 可下载 `varneb==0.0.1`。`.github/workflows/publish-pypi.yml` 仅由版本标签、已发布 Release 或手动确认触发。

## 原线程中断原因

原线程 ID：`019edb54-7732-77f3-8820-42123f31b2a3`。

最后一次请求在发送阶段失败：`input[1102].arguments` 长度为 `2,492,367`，超过接口允许的 `1,048,576`。这是请求字段长度校验失败，不代表项目代码或集群作业本身失败。新线程不要粘贴完整历史；以本文件、当前工作树、作业状态和日志为准。

## 已完成/已确认

- 端点弛豫摘要判据已修复：`converged` 使用 `FrechetCellFilter` 的完整广义力，同时保留原子力与应力字段。
- 端点路径预检已通过：最短距离约 `1.994 Å`（最终路径预检约 `1.983 Å`），最大变胞形变约 `0.1238`；HfO₂ 模板门槛已调整为可配置，默认约 `0.25`。
- 通用端点脚本不再暗含固定 `mpirun -np 40`，launcher 由环境/调用方提供。
- HfO₂ 7-image 普通 VCNEB 预收敛 40 步已完成：内部像 3 有局部峰，初步正向能垒约 `1.052235 eV`，反应焓约 `-0.743338 eV`，`has_interior_barrier=true`；但残差约 `0.847 eV/A`，不能作为最终能垒。
- 从完整第 40 步链快照成功续算，确认是状态恢复而不是重新插值。

## 当前计算状态

- HfO₂ 作业 `27674272` 已完成并降级为 60-Ry 诊断，不作为生产能垒。
- BTO 端点 `27675330`/`27675331` 及其 5/7/9-image 路径（`27676251`/`27676310`/`27676513`）均为 `100 Ry + DZP-10au` 的 ABACUS 静态 SCF，由 ASE cell-filter 负责外层更新；完整结果已归档并继续作为 BTO 生产证据。原生 `cell-relax` 方法迁移不再回头重算 BTO。
- BTO 已有完整的 5/7/9-image、反向路径和精度对照（均为 `100 Ry + DZP-10au`）；这些结果继续作为 BTO 生产证据。此前误发起的 ABACUS-native 端点方法重算 `27676929` (C, failed) / `27676930` (T, cancelled) 不纳入证据，也不再触发 BTO 路径重跑。后续原生 `cell-relax` 迁移只在 HfO₂ 主线执行，除非另有明确需求。
- BTO C→T 反向 5-image 已收敛（`27676299` + `27676485`）：`0.0195577 eV/A`、反应焓 `-0.0872206 eV`、无内部势垒。
- CI 闸门已审计并 withheld：当前没有内部 saddle，不启动物理上无意义的 CI；记录见 `outputs/batio3_t_to_c_pbe100_dzp10au/bto_ci_gate.json`。
- BTO 电子精度对照已完成：6×6×6/SCF 1e-9 的端点与静态路径反应能约 `0.0724 eV`，较 4×4×4 生产设置低约 15 meV；详见 `bto_precision_sensitivity.json`。
- BTO 端点求解器等价性已单独落档：`outputs/batio3_t_to_c_pbe100_dzp10au/endpoint_optimizer_equivalence.md`；FIRE 收敛端点不因迁移到 BFGS/native `cell-relax` 而失效，BTO 不重复计算。
- HfO₂ 端点兼容性已核验：当前 T/PO 输入均为 `12` 原子、`Hf4O8`，元素顺序一致；T 是为匹配 PO 常规胞而构造的 `√2×√2×1` 12 原子四方共格超胞，PO 是 12 原子极性正交常规胞。证据见 `outputs/hfo2_t_to_po_pbe100_dzp10au/endpoint_structure_audit.json` 和 `validation/hfo2_t_to_po/mapping_report.json`。
- 发布前干净安装烟测已通过：在临时目标目录构建/安装 `vcneb` wheel，并从安装目录运行 toy VCNEB（`barrier_eV=0.250004`）及 `vcneb --help`。
- PyPI 项目简介已同步到 `pyproject.toml`：`VARiable-cell Nudged Elastic Band code with universal first-frinciples calculators`；README 顶部保留规范拼写，Homepage/Repository/Issues 元数据链接已更新。新增 `tests/check_release_metadata.py` 并接入 `.github/workflows/publish-pypi.yml`，发布前会阻止简介或触发策略漂移。此次仅推送代码，未打 tag，因此不会触发 PyPI 发布；已发布的 `v0.0.1` 不覆盖，下一次按需发布的新版本（需新版本号）会携带该简介。
- HfO₂ 100 Ry/Orb-DZP-10au 端点重弛豫曾改用 ABACUS 原生 `cell-relax`/`relax_method bfgs`：ASE 外部 BFGS `27676827`/`27676828` 已保留首步诊断后取消；原生初次作业 `27676868`/`27676869` 已分别完成约 19/23 个离子步，但因驱动解析 `STRU_ION_D` 字段 bug 失败，断点已保留。修复解析器（含标签、磁矩、晶格常数单位换算）后从断点续算 `27676976` (T) 与 `27676977` (PO)，均 32 MPI；两端均因 ABACUS BFGS trust-radius breakdown 保留 `STRU_ION_D`，已切换到同属 ABACUS 的 `bfgs_trad` 续算 `27677116` (PO) 与 `27677117` (T)。两者均在观察多个完整 cell cycle 后因数值不稳定安全取消，完整输出已归档：PO 末态 `fmax=7.55946 eV/A`、`max_stress=2856.90 kbar`、`pressure=1644.58 kbar`；T 末态 `fmax=0.132648 eV/A`、`max_stress=27.40 kbar`、`pressure=-26.68 kbar`。两端仍保持 12 原子 `Hf4O8`，这只是优化器/步长稳定性问题，不代表物理端点不同。下一步可在相同能量面和自由度下采用更保守的外部 ASE BFGS 或阻尼 native 重启；两端通过受力+应力门槛后才启动高精度普通 VCNEB。
- 当前有效的保守外部 BFGS 端点作业为 `27677331` (PO) 与 `27677332` (T)，均 32 MPI、100 Ry、Orb-DZP-10au、`maxstep=0.0005`；前两轮仅因远端脚本同步/轨道文件名错误在预检阶段退出，不计入 DFT 证据。外部与 native 的验收标准统一为 `fmax<0.02 eV/A` 且 `max|stress|<0.1 kbar`。
- 上述 `maxstep=0.0005` 作业已在多步稳定窗口后取消并完整归档（T 第 19 步 `0.020087 eV/A`，PO 第 15 步 `0.483720 eV/A`），没有丢失恢复状态。当前续算为 `27677466` (T) 与 `27677467` (PO)，从最后轨迹帧恢复，改用 `maxstep=0.005`、`fmax=0.01`、120 步预算；观察 3--5 步稳定性后再继续。

## 最新轮询（2026-09-13）

- T 续算 `27677491` 已完成并通过双门槛：最大广义力 `0.00049882 eV/A`、最大应力 `0.01841 kbar`，12 原子 `Hf4O8`。
- PO 续算 `27677467` 仍在 `node148` 上以 32 MPI 运行；最新完整步为 9，`fmax=0.294759 eV/A`，连续下降。它是较低对称性的正交 PO 常规胞，不是已经收敛的 T 端。
- T 的完整输出已归档于 `outputs/hfo2_t_to_po_pbe100_dzp10au/endpoint_relax_T_ase_bfgs_maxstep005_job27677491/`。

## 端点晋级与普通 VCNEB（最新）

- PO 续算 `27677874`/`27677918` 已完成，最终 `fmax=0.0004973 eV/A`、最大应力 `0.07076 kbar`；T/PO 均为 12 原子 `Hf4O8` conventional cell，已由 `endpoint_promotion_gate.json` 原子晋级。
- 7-image calculator-free preflight 已通过：最小距离 `2.02496 A`、最大 deformation `0.04841`，无折返；所有 image 的 ABACUS 计算器能力检查通过。
- 原普通 HfO₂ 作业 `27677945` 因并行解析竞态失败；修复后的旧布局续算作业 `27678004` 已安全 CANCELLED（保留其完整 step 15 轨迹作为分布式迁移的恢复源）：128-task controller allocation，4 个并行 image worker，每 worker 32 MPI；100 Ry、Orb-DZP-10au、2x2x2、FIRE、无 CI。该轨迹只属于预收敛过程。

## 并行恢复修复（最新）

- `27677945` 后续因 ASE-ABACUS 的进程级 `ase_sort.dat` 竞态在 image 3 解析失败；已保留 step 2 链快照、manifest 和错误日志，未把该作业当作物理失败。
- `vcneb/abacus.py` 已加入 identity species-order 的安全补丁，并新增 `tests/check_abacus_parallel_sort.py`；远程 ICU 测试通过。旧全局 `ase_sort.dat` 已移到 `.stale_job27677945`。
- 修复后的旧布局续算作业 `27678004` 曾从 step 2 恢复，128 task / 4×32 MPI / no-CI；在迁移到按 image 分布式布局后已安全取消，作为恢复源的 step 15 完整轨迹和日志均保留。

## 分布式 image 语义与资源优化（最新）

- `n_images` 的定义是总帧数，始终包含两个固定端点；因此“7-image”=2 个端点+5 个 interior image，不是 7 个中间帧。
- VCNEB manager 现在首次读取端点的 energy/forces/stress 后缓存，后续迭代只通过 image executor 评估 `image_indices=[1,2,...,n_images-2]`；worker manifest 会记录 `image_count=n_images-2`，用于审计端点未重复计算。
- 分布式模板 `cluster/hf_hfo2_vcneb_distributed.slurm` 已改为默认 `IMAGE_WORKERS=N_IMAGES-2`、每 worker 32 MPI；7-image 默认 5 workers/160 tasks，仍需 2 个节点（单节点 128 核不足），不再为两个端点启动 worker。
- 新验证作业 `27678218` 已在 `hfacnormal01` 的 `node[381-382]` 启动，`NumTasks=160`、`NumCPUs=192`，首个 manifest 已确认只包含 image 1--5；从 `27678004` 的 step 15 完整轨迹恢复。此前错误布局的 `27678137`/`27678176` 已取消，不纳入物理结果。
- `27678218` 已完成（Slurm `COMPLETED`、`01:08:38`、exit `0`）：7 总帧（5 个 interior image）普通 VCNEB 在 FIRE `fmax_target=0.05 eV/A` 下达到 `final_max_generalized_force=0.0455156 eV/A`。最终正向焓垒 `0.1567509 eV`、反应焓 `-0.3252848 eV`，内部峰为 image 2（`has_interior_barrier=true`，未启用 CI）；manifest 共 96 条记录且每条只含 image 1--5，5 个 worker、每 worker 32 MPI，7 个 ABACUS image 目录无错误退出。
- 远端 `scripts/audit_vcneb_result.py --max-min-distance 2.0` 审计返回 `status=ok`、`issues=[]`：最小路径距离 `2.026335 A`、最大形变 `0.04894`、几何有效；summary 记录 `endpoint_evaluation_policy=fixed_cached_once`，确认端点未在 VCNEB 迭代中重复派发。结果文件为 `vcneb_summary.json/.txt`、`vcneb.traj`、`vcneb_barrier.png`、snapshots 和 worker manifest。
- 为补齐 HfO₂ 的 image-count 对照，已提交普通无 CI 作业 `27678406`（5 总帧=3 interior，1 节点×96 ranks）与 `27678407`（9 总帧=7 interior，2 节点×224 ranks）；两者均使用 100 Ry/10 au、32-MPI worker、固定端点一次缓存，独立工作目录分别为 `vcneb_n5_fire_distributed_cmp` 与 `vcneb_n9_fire_distributed_cmp`。
- `27678406` 已连续运行到 step 50：step 36 的最佳 `fmax=0.098970 eV/A` 后连续 14 个完整 optimizer step 回弹至 `0.188661 eV/A`。这已超过约定的 3--5 步观察窗口，因此在保留原目录、51 个完整快照和 manifest 后于 2026-09-14 02:35 左右安全取消。受控续算 `27678689` 已从旧 `vcneb.traj` 的完整 step 36 在新目录 `vcneb_n5_fire_distributed_resume_best36` 启动，仍为 5 总帧/3 interior、1 节点×96 ranks、3×32 MPI、普通无 CI，只将 `MAXSTEP` 从 0.02 降为 0.005；预检已通过且未覆盖旧证据。该续算完成 step 10（`fmax=0.112725 eV/A`）后也呈连续缓慢上升，于 2026-09-14 03:01 左右安全取消；新目录保留 11 个完整快照，step 0/旧 step36 的 `0.098970 eV/A` 是该 5-image 分支最佳证据。结论是 5 总帧在当前 100 Ry/10 au、`k=0.2`、FIRE 路径下未收敛，不把它当作最终物理能垒。
- 截至 2026-09-14 03:50 左右，`27678407` 已完成 step 55（`fmax=0.065946 eV/A`，从 step 0 持续下降），仍为 RUNNING；manifest 只记录 `[1,2,3,4,5,6,7]`，当前未见 ABACUS 错误 task。
- 由于 7-image 普通路径已达到真实内部峰且 `fmax<0.05 eV/A`，CI 精修作业 `27678507` 已提交：7 总帧/5 interior、2 节点×160 tasks、每 worker 32 MPI、`MAXSTEP=0.01`、`FMAX=0.03`、`CLIMB_AFTER=0`，从 `27678218` 的完整 `vcneb.traj` 恢复，工作目录为 `vcneb_n7_ci_refine`。该作业是首次真实 HfO₂ CI，不用于 BTO。
- CI `27678507` 已于 2026-09-14 04:14 左右正常完成（Slurm `COMPLETED`、`03:12:03`、exit `0`）。它在经历 step 15--80 的回弹/平台窗口后继续推进，于 step 91 达到 `final_max_generalized_force=0.0295475 eV/A`（目标 `0.03`），没有因首次回弹提前停止。CI 结果的正向焓垒为 `0.1291722 eV`，反应焓 `-0.3252848 eV`，最高/攀爬 image 为 2，`has_interior_barrier=true`；`endpoint_evaluation_policy=fixed_cached_once`、`climb_after=0`。独立审计返回 `status=ok`、`issues=[]`，最小路径距离 `2.035829 A`、最大形变 `0.0489263`、最大应力 `1.76792 kbar`。image worker manifest 共 279 条、状态全为 `ok`，只含 interior `[1,2,3,4,5]`，5 个 worker、每 worker 32 MPI；结果目录为 `vcneb_n7_ci_refine`。
- 9-image 普通无 CI 作业 `27678407` 已于 2026-09-14 04:30 左右正常完成（Slurm `COMPLETED`、`03:47:07`、exit `0`）：9 总帧/7 interior，最终 `fmax=0.0497191 eV/A`，正向焓垒 `0.1562092 eV`，反应焓 `-0.3252848 eV`，最高内部峰为 image 3。独立审计 `status=ok`、`issues=[]`，最小路径距离 `2.027974 A`、最大形变 `0.0496014`、最大应力 `5.38872 kbar`；manifest 207 条记录全为 `ok`，只含 image 1--7，7 个 worker、每 worker 32 MPI。
- HfO₂ image-count 对照已归档于 `outputs/hfo2_t_to_po_pbe100_dzp10au/image_count_comparison.md`；收敛的普通 7/9-image 焓垒只差 `0.0005417 eV`（约 0.35%），但 5-image 分支未收敛，因此不纳入物理能垒比较。已从远端复制 7-image CI 与 9-image 普通路径的 summary、audit、preflight 和 worker manifest；完整 ABACUS scratch 仍保留在 hfacnormal01。
- `scripts/compare_vcneb_images.py` 已用半个局部反应坐标段作为离散峰位默认容差，7/9-image 普通对照报告 `status=ok`、能垒 spread `0.0005417 eV`；机器可读结果为 `outputs/hfo2_t_to_po_pbe100_dzp10au/ordinary_image_comparison.json`。
- 新增 `outputs/hfo2_t_to_po_pbe100_dzp10au/hfo2_validation_provenance.json`，集中记录 12 原子端点、100 Ry/10 au、32-MPI worker 资源、5/7/9-image 与 CI 作业 ID、结果和本地/远端 artifact 位置。
- BTO 当前生产矩阵也新增 `outputs/batio3_t_to_c_pbe100_dzp10au/bto_validation_provenance.json`；它明确 100 Ry、Ba/Ti/O 全套 10 au DZP、5/7/9 总帧、反向路径、CI withheld 与 6³/1e-9 精度对照，并标明旧 `outputs/batio3/` manifest 仅为诊断记录。
- 已从 7/9-image 完成 summary 导出逐 image `vcneb_metrics.csv`，包含反应坐标、焓、晶格长度/角度、体积、应力和 NEB 力分解，可直接用于结构—能量图和论文表格。
- 项目计划的 DoD 已按现有证据更新：有限差分/ASE 对照、解析模型、VASP/ABACUS smoke、ABACUS 生产路径、HfO₂/BTO 收敛矩阵和纯 Python 最小案例均已勾选；真实材料模式约束、论文 claim 映射、VASP 生产级路径和更完整敏感性仍未完成。
- USPEX/ABINIT 语义对照条目也已按既有 `outputs/uspex_vcneb_mode_analysis.md` 证据勾选；旧 `/home/zhuxd/abacus/8.dielec/vcneb` 的可复现实验结果对比仍未宣称完成。
- 在检查 hfacnormal01 负载（`10454` idle CPUs）后提交的 HfO₂ 独立初始路径对照 `27678924` 已完成：7 总帧/5 interior、100 Ry/10 au、5×32 MPI/2 nodes、普通无 CI，唯一改变为 `CELL_INTERPOLATION=linear`；预检、终态审计和 linear/log-strain 比较均已归档。
- 为修复共享目录无 `.git` 时 provenance 的 `git_revision=null`，ABACUS driver 现在优先读取 `VCNEB_GIT_REVISION`，分布式模板也将其导出（未提供时记录 `remote-sync-unknown`）；后续提交任务应在 `--export` 中显式带当前源码 commit。
- 代码工程化增量已在本地回归全通过并推送为 `31318f5`（功能主体 `0735990`）：新增可选 exact-state image cache、calculator namespace 锁定、原子 cache 写入、损坏条目回退以及 manifest 命中/未命中 provenance；远端源文件已同步，当前运行作业不受影响。
- 后续工程化修复已提交：`ThreadedCalculatorExecutor` 在并发 batch 的某个 image 失败时仍收集并原子保存已完成 sibling image，随后再返回 batch 错误；新增恢复回归确认成功 image 不丢失、重启不重复计算，manifest 保留失败与后续 cache hit provenance。
- `run_vcneb()` 现在可通过 `failure_report=` 原子写出 optimizer/calculator 中断报告（异常类型、已完成步数、trajectory/snapshot 位置、诊断日志路径和恢复提示），同时保持原异常继续抛出；对应回归已通过。显式 line-search 的有界重试见后文。
- failure report 现在额外包含保守的 `failure_category`（`scf_nonconvergence`、`timeout`、`mpi_failure`、`nonfinite_evaluation`、`invalid_cell` 或 `calculator_or_optimizer_error`），分类会读取显式 image calculator 目录下日志的有限尾部，并由回归覆盖 ABACUS SCF 日志与关键词边界；它仍不替代完整 ABACUS/VASP stdout。
- `run_vcneb()` 新增显式 `BFGSLineSearch` 与有界 `line_search_retries`：仅对精确的 `LineSearch failed!` 在最后完整链状态重建 optimizer、缩小 `maxstep`/`stpmax` 后重试；其他 calculator/optimizer 异常仍直接抛出，回归覆盖失败一次后续跑与全局步数保持。
- `examples/run_vcneb_abacus.py` 与 `examples/run_vcneb_vasp.py` 已默认把该 failure report 写入各自 workdir 的 `vcneb_failure.json`，因此集群模板发生中断时可直接定位恢复入口。
- 已按当前 `docs/theory.md`、README、代码和回归证据同步 `VCNEB_PROJECT_PLAN.md`：广义坐标/单位/切线/收敛判据、模式投影顺序、约束释放、VASP 静态 image 语义、calculator 切换和 dry-run/validate-only 现标为完成；ABACUS 专用 SCF 分类、USPEX/旧实现对比、真实材料模式对照、论文图表等仍保留为未完成。
- README 的 Files 清单已为 `unit`、`model`、`DFT-smoke`、`production-template` 示例加上显式标签，并同步勾销 P0 对应计划项。
- `docs/theory.md` 已补齐可执行的 manager/worker VCNEB 伪代码和核心公共参数表（默认值、单位、合法域与 CI 时序）；calculator-specific 经验参数和磁性/占据敏感性仍保持未完成标记。

## 本轮新增（2026-09-14）

- 论文草稿已同步到当前生产证据并提交为 `e7d315a`：HfO2 普通 7/9-image、CI（含延长回弹观察）和 BaTiO3 5/7/9-image barrierless 矩阵均已写入；Program Summary 固定为 GPL-3.0-or-later、包版本 0.0.1。LaTeX 编译 exit 0（6 页）；输出 PDF 为本地生成物，未作为源码结果提交。
- 隔离环境 `python -m build` 已成功生成 `varneb-0.0.1` sdist/wheel，包内 `METADATA` 的 Summary 与 GitHub About 完全一致；许可证已改为 SPDX `GPL-3.0-or-later` 并保留 `LICENSE` 文件，不再触发旧式 license table 弃用警告。
- 当前 PyPI 上已存在的 `varneb==0.0.1` 仍显示发布时的旧 Summary；PyPI 不允许覆盖同版本文件，因此新的 GitHub About 只会在下一次递增版本（例如 `0.0.2`）并按需触发 Trusted Publishing 后生效。本轮未打 tag、未触发发布；`tests/check_release_metadata.py` 与隔离构建均已确认新 Summary 会进入 wheel/sdist。
- 已新增 `outputs/hfo2_t_to_po_pbe100_dzp10au/error_budget_and_efficiency.md`，分开记录 image-count、CI 路径选择、端点残差、电子精度和资源效率，避免把不同物理/数值因素合并成单一误差条。
- 已新增 `outputs/legacy_vcneb_comparison.md`：基于历史迭代日志记录旧实现的可证实差异和 MIC 失败模式；旧目录当前不可访问，定量同条件 benchmark 明确保留为未完成。
- 已新增 `scripts/plot_vcneb_metrics.py` 与 `tests/check_vcneb_plot.py`：从已归档的逐 image CSV 生成焓垒、晶格、体积和广义力四联图；回归通过，未调用任何 DFT 计算器。
- 新增 `tests/test_regression_scripts.py` 作为 pytest 收集入口，复用全部 `tests/check_*.py` 脚本；`python -m pytest -q` 已通过 `10 passed`，未启动 DFT。
- 新增 `examples/mode_template.json` 及 `tests/check_mode_template.py`，提供可复制的原子/可选 cell 模式 JSON 模板；回归入口现为 `11 passed`，仍不启动 DFT。
- 27678924 的线性 cell 插值对照已完成：Slurm `COMPLETED`、耗时 `01:51:16`、exit `0`，step 51 达到 `final_max_generalized_force=0.0495087 eV/A`；正向焓垒 `0.1596772 eV`、反应焓 `-0.3252848 eV`、最高内部峰 image 2。独立审计 `status=ok`、`issues=[]`，最小距离 `2.028435 A`、最大形变 `0.0496212`、最大应力 `3.88299 kbar`；manifest 共 159 条、全部 `ok`，只含 interior `[1,2,3,4,5]`、5 workers/32 MPI，端点策略 `fixed_cached_once`。结果已复制到 `outputs/hfo2_t_to_po_pbe100_dzp10au/vcneb_n7_fire_distributed_linear_job27678924/`。
- 新增 `--allow-duplicate-image-counts` 及插值/优化器字段到 `scripts/compare_vcneb_images.py`，并以回归覆盖同一 image 数的 variant 比较。HfO₂ linear 7-image 与 log-strain 7/9-image 的报告为 `outputs/hfo2_t_to_po_pbe100_dzp10au/cell_interpolation_comparison.json`：barrier spread `0.0034680 eV`、最高峰反应坐标 spread `0.04639`、状态 `ok`；linear 相对 log-strain 7-image 高 `0.0029263 eV`。metrics CSV 和四联图也已生成。
- 新增 `scripts/export_vcneb_structural_metrics.py` 与回归 `tests/check_geometry_metrics.py`：从最新完整 trajectory 导出指定 MIC 关键键长和扩展空间模式投影；HfO₂ linear 7-image 已生成 structural metrics CSV/JSON。端点位移投影明确不是声子软模，真实材料模式约束对照仍未宣称完成。
- 清理生产文档中的参数歧义：README 的 HfO₂ 示例已改为 100 Ry、全套 10 au DZP、2×2×2、SCF 1e-8；论文将旧 60 Ry/低精度 HfO₂ 仅标为历史 adapter smoke provenance，生产结果只引用 100 Ry/10 au。论文 LaTeX 使用 TeX Live 编译 exit 0（6 页）。
- HfO₂ smoke/fallback 示例默认值已统一为 100 Ry、全套 10 au DZP、2×2×2、SCF 1e-8、`mixing_beta=0.3`；新增 `tests/check_production_parameters.py` 锁定 BTO/HfO₂ 集群模板与示例的生产参数，完整回归为 `14 passed`。历史低精度 manifest 不被覆盖。
- `506d924` 已将显式 `BFGSLineSearch` 的有界重试参数贯通到 ABACUS/VASP 驱动及四个 Hefei Slurm 模板；默认 `LINE_SEARCH_RETRIES=0`，只有用户显式选择该优化器并设置预算时才启用，FIRE/普通 BFGS 行为不变。静态检查与完整回归仍为 `14 passed`。
- 本轮复核作业 `27674272`：父作业已不在队列（Slurm 对已结束父作业返回 invalid job id），可见数组子任务 `.497`--`.506` 均 `COMPLETED`、exit `0`；`hfacnormal01` 当前无本人运行任务。工作树中的未跟踪项均为既有 DFT 输出/归档目录，未纳入代码提交。
- `run_vcneb_abacus.py` 与 `run_vcneb_vasp.py` 现已提供统一的 `--mode`、`--mode-guided`、`--constraint-mode` 和模式归一化参数；四个 Hefei Slurm 模板也支持 `MODE_FILE`/`MODE_GUIDED`/`CONSTRAINT_MODE` 等环境变量。新增 CLI/模板回归后完整 pytest 为 `15 passed`，仅做接口验证，没有重跑 BTO 或启动 CI。
- 新增端点位移模式导出器 `scripts/derive_endpoint_mode.py`；修正为 VCNEB 参考 cell 下的扩展原子坐标（不是当前 cell 的 Cartesian 差），并由 `tests/check_endpoint_mode.py` 验证 12 原子 HfO₂ 严格子空间端点可连接。该模式只作结构位移诊断，不能称为声子本征模。
- 远端 `27679815` 的 `VALIDATE_ONLY=1` 严格模式预检通过（7 总帧、5 interior、2 节点资源，不调用 ABACUS）。首次真实模式引导作业 `27679826` 因远端核心文件滞后、在 DFT 前报 `run_vcneb() got unexpected keyword line_search_retries`，不计物理证据；同步 `vcneb/core.py` 后重试 `27679830` 正在运行，配置为 HfO₂ 7-image、5×32 MPI、模式引导后无约束普通 VCNEB、无 CI。已完成 step 0--11，`fmax` 从 `0.974822` 降至 `0.416943 eV/A`，未见回弹或 SCF/MPI 失败；须等终态 summary/audit 后再归档。

## 建议的新线程第一步

1. 新计算先检查 `sinfo`/`squeue`，并使用 `hfacnormal01` 的 32-MPI image worker 资源模型。
2. 普通路径若发生回弹，至少观察 3--5 个完整 step，保存窗口内最佳快照，再决定是否续算或审计。
3. CI 只在普通路径达到门槛且存在真实内部能量峰时启动；当前 BTO 结论是 withheld。
4. 每个长作业都归档完整 trajectory、snapshots、summary、preflight、worker manifest、节点/task 数和审计结果。
5. 继续按 `VCNEB_PROJECT_PLAN.md` 推进，并把每个结论绑定到当前文件、日志、测试或作业输出。

## 交接原则

- 以当前工作树和实时作业状态为事实来源；本说明只提供压缩上下文。
- 新线程提示词保持短小，引用本文件，不要复制整段历史日志。
- 任何新结论都要区分“诊断/预收敛”与“最终收敛结果”。
