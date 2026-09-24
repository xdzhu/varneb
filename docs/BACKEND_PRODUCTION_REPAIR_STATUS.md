# VARNEB 多后端生产修复状态

更新时间：2026-09-24（HF `hfacnormal01`）

本文只记录已由日志、作业状态或回放实验支持的结论。旧失败目录保留，修复使用独立目录，不覆盖原始轨迹。

## 2026-09-24：GaN 45.7 GPa 的 QE、CP2K 与 ABINIT 重建

这轮不再续用三个后端的零压旧链，而是从已验证的 GaN B4/B1 四原子种子出发，
分别在各后端自己的 PBE 契约下重新做 45.7 GPa 可变胞 BFGS。普通 VCNEB 仍用
`fmax=0.10 eV/A`；端点门禁使用用户指定的残余应力 `<1.0 kbar`。为避免 ASE
晶胞过滤器的力式停止条件提前终止，端点 driver 现分别检查原子力和相对于目标
静水压力的应力残差，不再用一个与体积相关的广义 `fmax` 代替压力判据。

三套端点均已完成并通过门禁，且 B4/B1 均保持 4/6 配位和 4/4 原子：

| 后端 | 端点作业 | B4/B1 残余应力 (kbar) | B4/B1 最大广义力 (eV/A) | 输入契约 |
|---|---|---:|---:|---|
| QE | `27770124/25`，B4 加严 `27770163` | `0.492/0.364` | `0.00302/0.00190` | QE 7.0，PseudoDojo NC-SR-PBE v0.4，100/600 Ry，4×4×3，单 rank `srun` |
| ABINIT | `27770099/27770101`，B4 加严 `27770164` | `0.763/0.225` | `0.00468/0.00117` | ABINIT 8.6.1，同源 PseudoDojo PSP8，1400 eV，4×4×3，Intel-2017 Hydra 8 rank |
| CP2K | `27770175/76` | `0.416/0.249` | `0.00283/0.00240` | CP2K 2024.1，GTH-PBE/DZVP，800 Ry、REL_CUTOFF 80，4×4×3，16 rank `mpirun` |

QE 的 32-rank `srun` 端点试投 `27770097/98` 在 `MPI_Init_thread` 失败；同机历史
成功记录和重试均证明该站点构建应使用单-rank `srun --exclusive ... pw.x`，并把
并行放在独立 image 上。QE 29-total-image 生产链 `27770529` 已通过端点哈希、
路径几何和计算器契约，9 个 interior-image worker 正常推进；初始 `1.5046 eV/A`
到 step 8 已降至 `0.4224 eV/A`，后续回弹按完整轨迹观察，不因单步上升停止。

CP2K 的结论需要区分启动器和 shell 生命周期。相同 GaN 静态输入的受控 benchmark
中，直接 `mpirun -np 8` (`27770135`) 用时 3:44，`mpirun -np 16`
(`27770136`) 用时 2:23，能量一致；原单-rank 端点在 12 分钟仍未完成首步，故取消
并保留。先前失败的是把多个 shell 交给错误的 `srun`/协议组合，而不是
`cp2k_shell.psmp` 永远只能单 rank。生产首投 `27770686` 又揭示 ASE 在 calculator
构造期立即启动 shell：29 个像会瞬间形成 464 rank，绕过 `IMAGE_WORKERS=5` 的并发
上限。CP2K factory 已改为惰性、一次 image 评价期间启动一个 16-rank shell，并在
复制输出后关闭；`--validate-only` 也不再实例化任何会启动外部进程的 calculator。
修正版生产为 `27770714`，首轮应只出现一个端点 MPI world，随后最多五个并发
interior worlds；须以实际 Slurm steps 和完整 step 0 为准继续审核。

ABINIT 已弃用先前不可比较的 HGH-LDA 端点，改用带审批清单、逐文件 SHA-256 的
PseudoDojo NC-SR-PBE v0.4 PSP8。生产 `27770687` 已完成 step 0 的全部像并进入下一轮；
9 个并发 image worker、每像 8 Hydra ranks 均返回完整能量、力和应力。作业 stderr
中的站点 ROCm modulefile 提示不影响 ABINIT 计算，但保留在审计记录中。

两次 QE 生产预提交 `27770190/27770301` 均在首个 DFT 前退出，根因是远端 guarded
source 的入口脚本与 `vcneb` 包 API 版本不一致；其目录保留但不计作物理失败。
远端现按整包同步并执行 `compileall` 与 calculator-free `--validate-only`，防止模板、
入口和包 API 的部分部署。三条生产链都使用 identity mapping、无整体平移、无晶胞
旋转、linear cell interpolation，以保持有效端点哈希与静态门禁完全一致。

BTO 的生产定义维持 **7 total images（5 interior）**。9 total images 不是统一硬要求；
只有出现未解析的路径曲率或分辨率依赖时，才做独立的 image-count 收敛检查。

## 统一判定边界

VARNEB 控制器只依赖每个 image 的 `energy`、`forces` 和完整 `stress`。后端适配器负责输入生成、外部程序启动、输出解析和 SCF 诊断；FIRE、BlockFIRE、SplitFIRE 等路径优化器不负责修复后端输出，也不把 SCF 重试伪装成路径迭代。

失败分类现在进一步按后端契约分层：ABINIT 的 `chkorthsy` 归为
`abinit_symmetry_failure`，缺失或未固定伪势归为 `pseudopotential_contract`，
能量/力/应力区块缺失或解析失败归为 `output_contract_violation`。这些分类在
`MPI_ABORT` 或 launcher 文本之前判定，避免把输入/输出契约错误误报为 MPI 故障，
也避免错误地修改 FIRE、NEB 图像数或路径步长。

生产路径仍使用普通 VC-NEB `fmax = 0.10 eV/A`、固定端点和无 CI。端点静态门禁的默认残余应力阈值为 `< 1.0 kbar`；高压案例比较的是相对于目标静水压力的残差，而不是绝对应力。最终是否收敛以完整链的 `final_max_generalized_force_eV_per_A <= 0.10` 判定，不能以 Slurm exit code 或某一轮回弹判定。

HF 的 ABINIT 8.6.1 启动契约也已固定：该二进制链接 Intel MPI/Fortran
2017.4/2017.5，必须加载匹配的 `compiler/intel/2017.5.239` 与
`mpi/intelmpi/2017.4.239`，再由 `mpiexec.hydra -bootstrap slurm -n N abinit`
启动。HF 默认 Intel 2021 环境的 `srun ... abinit` 会向该二进制传入
`pmi_args`，而仅切换到 `mpiexec` 仍会因缺少 `libifport.so.5` 或 PMI ABI
不匹配而失败。独立 32-rank canary `27749997` 在匹配的 compiler/MPI 与
Hydra 组合下返回 `abinit 8.6.1`；对照 `27749996` 的 `srun --mpi=pmi2`
仍以 `PMI_KVS_Get returned -1` 失败。此前 `27749797/98`、`27749871/72`
的目录和日志均保留，未将启动层失败误判为 SCF 或 NEB 发散。

ABINIT 的伪势也是输入契约的一部分：包含 `pps` 的参数必须同时提供存在的
`pp_paths` 目录。缺少该目录会在 ASE 写入输入前给出明确的
`ABINIT parameters with 'pps' require explicit pp_paths`/目录不存在错误，
不再等到 32 个 rank 启动后才产生含糊的伪势解析失败。GaN HGH-LDA 试跑
`27750045/46` 正是捕获了这一缺口；补充 `abinit_hgh_factory_kwargs.json`
后的 `27750056/57` 进入真实端点计算；其中 B1 暴露了自动对称性中止，B4
随后为保持端点契约一致而取消，二者目录均保留。

后续 `nsym=1` 的 ABINIT GaN HGH-LDA B4 端点虽满足宽松的原子力阈值，
其 4 原子晶胞却从 `39.244 Å³` 被 BFGS 压到 `7.222 Å³`，能量同时从
约 `-3666.866` 升至 `-3350.615 eV`；原始点 ASE/ABINIT 应力对角约
`+15.1 eV/Å³`（ABINIT 原输出约 `+2418 GPa` 的拉伸应力），与正常
GaN 量级严重不符。旧输入 `ecut=600 eV` 经 ABINIT 回显仅 `22.05 Ha`，
这是需要优先排除的截断能/Pulay 应力来源，而不是 NEB 优化器问题。
已提交单变量高截断静态诊断 `27761479`：保持 HGH-LDA/2×2×2 和同一
B4 起始结构，只将 `ecut` 提至 `2200 eV`（约 `80.85 Ha`）。此任务
不进行变胞、不生成 NEB；其结果只能判断旧应力异常是否对截断能敏感，
不能将 HGH-LDA 当成 PBE 生产验证。来源和判据见
`validation/backend_smoke/abinit_gan_ecut2200_stress_probe_20260923.json`。
首次 32-rank 尝试 `27761479` 在首个 SCF 前因 Slurm step 内存不足退出，
没有获得应力，故不能当作高截断物理失败。HF 拒绝 8 CPU/64 GB 的
高内存申请（站点限制单 CPU 内存）；重试 `27761497` 在独立目录保留
32 CPU 内存配额、实际只启动 8 个 Hydra MPI rank，输入完全不变。
重试已成功完成 SCF：同一 B4 几何在 `2200 eV` 下最大原子力
`8.607 eV/A`，应力对角约 `0.735/0.735/0.806 eV/Å³`
（约 `118–129 GPa`）。应力较 `600 eV` 的约 `15 eV/Å³`
降低约 20 倍，但仍远未可用；能量和力也对截断能强烈敏感。
故旧 HGH-LDA 端点不具备可比性，`2200 eV` 也尚未证明收敛。
下一步不能原样续算旧链；必须先有经验证的 PBE 伪势与截断能/应力
收敛契约，再重新端点优化及静态门禁。

另外，ABINIT 后端默认写入 `nsym=1`（恒等对称操作），避免由变胞路径或
POSCAR 末位舍入触发自动对称性识别的 `chkorthsy` 中止。需要利用更高对称性
时必须显式提供并审核 `nsym/symrel`，不能让生产路径隐式依赖自动分类。

## ABACUS GaN：结果解析契约修复

### 根因

GaN ABACUS 作业 `27741518` 在 step 0 的 image 5 报出 NumPy 不规则数组 `(30,)`。回放确认，ASE-ABACUS 3.23.1b1 在构造完整结果字典时会提前解析 eigenvalues；该 image 的 30-k-point eigenvalue 区块形状不规则，而 VCNEB 实际只需要能量、力和应力。

修复版第一次重跑 `27744183` 绕过 eigenvalue 后，image 5 的 header 又出现并行输出格式损坏，ASE 私有 header 正则在 `SELF-CONSISTENT` 区块抛出 `NoneType.group`。这仍是输出解析问题，不是 FIRE 或路径几何问题。

### 处理

- `vcneb/abacus.py` 新增只解析 VCNEB 三项的结果入口，避免可选 eigenvalue 解析。
- 完整 ASE 解析失败时，使用严格的最小日志解析器读取最后一组 `final etot`、`TOTAL-FORCE` 和 `TOTAL-STRESS`；结果仍由统一的 `ImageEvaluation` 检查数量、形状和有限值。
- 对旧损坏日志的回放结果为：energy `-4726.0740513 eV`、forces `(4,3)`、stress `(6,)`。

Guarded continuation `27749598` 在完整写出 step 0--9 后，于 image 10 因
ABACUS 日志的人类可读 `final etot is` 行被长行截断而被误判为无能量；同一
日志实际包含完整的 `!FINAL_ETOT_IS -4726.5466304132096411 eV`、力和应力
区块。此前错误分类为 `mpi_failure`，根因是最小解析器只识别前一种能量标记，
不是 MPI 或晶胞信赖域问题。解析器现在同时接受两种官方输出标记，并增加了
该真实截断格式的回归测试；原 `27749598` 目录保留，续算使用独立目录。

独立续算 `27755927` 曾从 `chain_step_0009.traj` 启动；两次仅参数错误
的提交 `27755906/24` 在 DFT 启动前退出，不计入计算失败。该续算采用
9 个 image worker、每像名义 32 MPI，但 `srun --mpi=pmix_v3` 与本机
Intel MPI/ABACUS 组合不兼容：每个 image 出现 32 次 `PMI server not found`，
ABACUS 进程表却只有 1 个进程。实际是 32 个 singleton 竞争同一个 image
目录，不能当作有效 DFT 结果。`27755927` 已于 2026-09-22 16:18 CST
取消，原始目录保留。
- 修复已加入回归测试并推送；ABACUS 标记解析的实现提交为 `ec10b2a`，后续
  契约分类修复在 `6fc573d`。

该无效启动曾写出 `chain_step_0000.traj` 和 `chain_step_0001.traj`，控制器记录的
最大广义力由 `1.694169` 降至 `1.525129 eV/A`，但多进程同目录写入使这些
数字不能用于收敛判断。更重要的是该旧链错误使用了零外压，不能作为
45.7 GPa 的 GaN B4→B1 文献对照。旧端点绝对应力约为
`484.97/468.50 kbar`；不能把它们同零应力比较来判定高压端点是否合格，
必须比较与目标 `457 kbar` 静水应力的残差。中间像的约 `1.07 eV` 峰值
不作为能垒结论。

安装包的 `case/abacus.slurm` 使用 `source .../scripts/env.sh` 后
`mpirun -np $SLURM_NTASKS abacus`。隔离 canary `27757793` 证明直接
`mpirun` 可形成真正的 4-rank MPI，但 `--cpus-per-task=1` 使四个 rank
全绑在同一 CPU，已及时取消。改为一个 Slurm task 拥有 8 CPUs 后，
canary `27757811` 完成同一 GaN 静态像的 4/8-rank 计算：分别用时
`170.44/124.79 s`，能量为 `-4727.3758080966627/-4727.3758080965881 eV`，
力与应力一致且完整。并发探针 `27757818` 证明两个 4-rank
`mpirun` worker 可在同一分配中同时完成。

正式 GaN 脚本因此改为单节点、10 个 Slurm tasks × 4 CPUs/task：
1 个 manager 加 9 个各自 `mpirun -np 4 abacus` 的 image worker，
27 个中间像分三批；端点不重复计算。解析器在读取结果前会拒绝
`PMI server not found` 输出，防止 singleton 回退的结果混入路径。
新生产作业 `27760826` 已观测到 9 个并行 image worker、36 个同时运行的
ABACUS MPI rank，且 27 个内部像每轮均返回完整结果。作业以 Slurm
`COMPLETED/0:0` 结束；step 44 的完整链最大广义力为
`0.096844 eV/A`，达到默认 `0.10 eV/A` 阈值。29 像焓均为有限值，
路径几何检查通过，日志未发现 `PMI server not found` 或 SCF 未收敛标记。
能垒为 `0.654732 eV/4-atom cell = 0.327366 eV/GaN`，反应焓为
`-0.033113 eV/GaN`，唯一内部峰位于 image 15。Qian 等在 45.7 GPa
报告约 `0.39 eV/GaN`；本结果低 `0.06263 eV/GaN`（约 16.1%），比较时
仍须注明赝势、基组和实现差异。机器可读审计见
`validation/backend_smoke/abacus_gan_45p7_vcneb_20260923.json`。

为防止同类问题再次发生，`hf_gan_abacus_vcneb.slurm` 现在默认要求
`ENDPOINT_STATIC_SUMMARY`，并在启动 VC-NEB 前调用统一的 `validate_ase_static_gate.py`。
只有显式设置 `REQUIRE_ENDPOINT_STATIC_GATE=0` 才能绕过；绕过的结果必须标记为
诊断而不是生产结果。生产 driver 还会比较静态摘要与经过原子映射、晶胞对齐和
周期平移后的实际两端结构哈希；任何不一致均在首个 DFT 路径迭代前终止。
`relax_abacus_endpoint.py` 的摘要也补齐了通用端点身份、状态和能量字段，可直接
接入统一静态摘要构建器。

旧的 `27744907` 是单像晶胞信赖域失控的已取消诊断目录；它不再是待等待的生产
任务。高压端点已在同一 ABACUS 3.10.0LTS/PBE/100 Ry/DZP 10 au
契约下重新以 BFGS 优化：B4 `27760738`、B1 `27760745` 均完成，
各为 4 原子；B4 保持四配位，B1 保持六配位。独立静态作业 `27760794`
完成，力为 `0.00622/0 eV/A`、相对于 45.7 GPa 的最大应力残差为
`0.357/1.017 kbar`，通过预先声明的力 `<0.10 eV/A`、残差 `<2 kbar`
门禁。29 total images（27 interior）的无 DFT 预检匹配两端静态哈希。
首次生产提交 `27760808` 因 Bash 空数组与 `set -u` 的启动问题在 DFT 前退出，
目录保留；修复并经启动 dry run 后，独立目录中的正式生产 `27760826`
于 2026-09-22 23:24 CST 启动，使用 45.7 GPa，现已按上述阈值收敛并完成审计。

## CP2K BTO：SCF 与资源启动分层修复

### 根因

原生产 `27741430` 在 BTO image 4 的 CP2K 内层 SCF 达到 500 次仍未收敛。输出显示 Broyden mixing `ALPHA=0.20` 下残差长期振荡；这是电子自洽问题，不是 NEB 候选步问题。

两个探针失败也已区分：

- `27744198` 申请共享节点 1 CPU，却在内部要求 `srun --exclusive`，一直等待作业 step，随后取消；属于 Slurm 资源请求错误。
- `27744265` 通过当时的 `srun`/MPI 组合启动 32-rank `cp2k_shell.psmp`，shell 握手失败；该结果只能否定这一启动组合，不能推出 CP2K shell 必须单 rank。2026-09-24 的直接 `mpirun` 对照已验证 8/16 rank 均可用。

### 处理

新增 `examples/material_profiles/cp2k_bto_pbe_dzvp_stable.json`：

- `EPS_SCF = 1e-5`，`MAX_SCF = 1000`；
- `ADDED_MOS = 30`，300 K Fermi smearing；
- `DIRECT_P_MIXING`，`ALPHA = 0.05`；
- 当时的诊断阶段采用多个 one-rank image worker 绕开错误的 32-rank `srun`；当前生产契约已由实测更新为每个活跃 image 一个 16-rank direct-`mpirun` world，并以惰性生命周期限制并发数。

本轮端点重跑进一步发现并修正了一个输入单位错误：ASE 的 CP2K calculator
把 `cutoff` 解释为 eV，而 CP2K 文献 profile 通常以 Ry 给出。旧的稳定 profile
裸写 `400`，实际只有约 29.4 Ry，足以产生异常大的 Pulay 应力并把 BFGS
推向错误的晶胞。现在 `make_ase_cp2k_factory` 接受显式的 `cutoff_ry`，在
后端边界转换为 eV；BTO/GaN 稳定 profile 均改为 `cutoff_ry: 400`，并拒绝
同时提供 `cutoff` 与 `cutoff_ry`。已取消的 `27749624/25/27/28` 只作为这
一错误输入的诊断证据保留，修正版将在独立目录重新准备端点。

单结构探针 `27744287` 已成功；两个相同的旧 image-4 结构均在 101 次 SCF 内收敛，退出码为 0。BTO 生产 `27744298` 随后因端点静态审计不通过、且路径力连续多步发散而取消；独立目录和日志均保留。

GaN CP2K `27741431` 也因端点静态审计显示固定端点基线无效而取消；原目录保留，只有重新生成并通过端点门禁后才允许新的 profile 续跑。

通用 ASE Slurm 模板现已在源头加入两项防护：

- CP2K 不再在 calculator 构造期启动所有持久 shell；每个活跃 image 惰性启动一个经 benchmark 验证的 16-rank direct-`mpirun` world，评价结束即关闭，实际并发由 `IMAGE_WORKERS` 限制；
- 生产 VC-NEB 必须提供 `ENDPOINT_STATIC_SUMMARY`，并在启动路径前验证两个端点的力 `< 0.10 eV/A`、残余应力 `< 1.0 kbar`。先生成摘要时显式使用 `STATIC_ONLY=1`；绕过门禁必须显式设置 `REQUIRE_ENDPOINT_STATIC_GATE=0`。
- 静态摘要中的两端 SHA-256 必须与映射/对齐后的实际 VCNEB 两端完全一致。该检查覆盖“源文件静态合格、但路径预处理改变了计算器实际输入”的缺口。对 CP2K 等有限实空间网格计算，生产模板允许显式设置 `ALIGN_TRANSLATION=0`，避免把静态门禁从一个网格原点带到另一个网格原点；这是可审计的路径定义，不是运行后的自修复。

这些是输入契约层的预防检查，不会修改或“修复”已有轨迹。

2026-09-22 的加严端点结果已通过统一门禁：BTO 初/终端最大广义力为
`4.24e-4/3.64e-4 eV/A`、最大应力为 `0.0094/0.0348 kbar`；GaN 初/终端为
`3.72e-4/4.82e-4 eV/A`、最大应力为 `0.0241/0.0496 kbar`。端点原子数分别
保持为 5/5 和 4/4；可复核记录见
`validation/backend_smoke/cp2k_bto_endpoint_gate_20260922.json` 与
`validation/backend_smoke/cp2k_gan_endpoint_gate_20260922.json`。

首次提交 `27756548/27756549` 在写门禁报告前退出，根因是新 `WORKDIR` 尚未
创建，而不是 CP2K SCF 或静态门禁失败。模板已在提交 `3190ccc` 中于门禁前
创建工作目录，并增加回归断言；修正版在独立目录重提为 BTO `27756554`
与 GaN `27756555`。两者当时仅通过了“源文件静态值”门禁并进入普通 VCNEB，
未覆盖首次失败目录；后续有效端点哈希审计分别否定了退化 BTO 路径和 GaN
平移后终态，因此这一历史门禁不能再作为生产通过证据。

端点准备现在也有统一入口 `scripts/relax_ase_endpoint.py`：它通过与生产路径相同的 ASE factory，对单个端点执行可变胞 BFGS（或 `--fixed-cell` 原子 BFGS），保存 `CONTCAR`、断点、力/应力/焓摘要。端点准备、静态门禁和 VC-NEB 因而成为三个可审计阶段，而不是在路径迭代中临时改变后端参数。

### CP2K BTO 端点相身份审计（2026-09-23）

`27756554` 以 `final_max_generalized_force=0.00963 eV/A` 正常结束，
但这**不是有效的 T→C VCNEB**。端点文件命名首先反置：
`cases/bto/endpoints/T_CONTCAR` 实际是
`4.028225×4.028225×4.028225 Å` 的立方结构，
`C_CONTCAR` 实际是 `3.990321×3.990321×4.239804 Å` 的
四方极性结构。旧 CP2K 初端从前者（立方）起步；旧终端从后者
（四方）起步，但在该 CP2K 契约下经 27 步 BFGS 弛豫至近立方
（`fmax=3.64e-4 eV/A`）。因此输入命名错误与四方端点坍缩是
两个独立事实，不能仅归因于其中之一。独立弛豫后的两端均为近立方：
初/终端 `c/a=0.9999976/1.0000231`，两端插值路径总长度仅
`1.80e-4 Å`。因此程序报告的 `2.61e-7 eV` 只是同相近重合端点间的
数值差，不能作为 BTO 相变能垒或多后端生产通过证据。静态力/应力门禁
只能证明端点驻定，不能证明它们保留了指定晶相；完整审计记录见
`validation/backend_smoke/cp2k_bto_endpoint_gate_20260922.json`。

为在 DFT 启动前挡住这种退化路径，通用初始链预检新增可选的
`--minimum-endpoint-separation`，按原子映射与平移对齐后的原子+晶胞扩展坐标
测量端点距离。HF 通用模板可通过 `MINIMUM_ENDPOINT_SEPARATION` 显式设置
案例阈值；对 BTO T→C 下一轮至少要求 `0.05 Å`，其数值低于已完成的
VASP BTO T→C 初始路径约 `0.345 Å` 端点距离，却远高于这次退化路径。
此通用几何检查不是晶相鉴定的替代品。新增计算器无关的
`scripts/validate_bto_phase_pair.py`，在 DFT 前检查五原子成分、四方
`c/a`、Ti 相对赤道 O 面的极性位移及立方端的近零畸变。实际旧输入
审计同时报出 T 标签缺少四方/极性、C 标签保留四方/极性。
`vasp_vcneb_production3/snapshots/step_0000/POSCAR_00` 与旧
`C_CONTCAR` 逐字节相同；因此从该 VASP T 文件提交的诊断
`27761267` 在识别为旧已完成弛豫的重复计算后，于运行 9 分钟时
及时取消，独立目录保留。下一步不是重跑同一参数，而是先对
CP2K 2×2×2 k 点、DZVP 基组和 300 K 电子展宽做受控精度/相稳定性
测试，再决定是否存在可用的 T 端点。上述 VASP T 文件及 CP2K
profile 的 SHA-256 分别为
`dc3f811c0d5b0c2988926d702eea534414edcf6ff91252934d3f0bc6d96f7ac1`、
`9260b8196f44955a24e6823616556da87b6fc04acb6199930c9da762909661ea`。
为避免再次依赖反置文件名，未改动原文件，另在 HF
`cases/bto/phase_verified_seeds_20260923/` 保存只读意义的正确命名副本：
`T_P4mm_seed.vasp`（spglib #99，`c/a=1.06252`，Ti 极性位移 `0.215 Å`）
和 `C_Pm3m_seed.vasp`（#221，`c/a=1`，Ti 极性位移 0）。两者是输入
种子而非 CP2K 已验证端点；哈希与来源见
`validation/backend_smoke/bto_phase_seed_audit_20260923.json`。
已启动受控的 CP2K 静态对照 `27761420`：只把原 2×2×2 k 网格改成
4×4×4，其余 PBE/DZVP/400 Ry/300 K/SCF 参数不变，在两个正确标记的
T/C 输入种子上各做一次独立静态计算。该测试用于判定 k 网格对能量、力
和应力的影响，不被当成 T 相已稳定或整条 VCNEB 已验证；配置和输出
边界见 `validation/backend_smoke/cp2k_bto_k4_phase_probe_20260923.json`。
结果已完成：旧 2×2×2 对同一输入种子给出
`E_T-E_C=+0.159293 eV`，新 4×4×4 给出 `-0.063233 eV`，仅 k 网格
变化就翻转两相静态能量排序。新 T/C 静态原子最大力仍为
`0.1252/0.1221 eV/A`，不能据此宣布端点收敛；已在独立目录提交
同一 4×4×4 契约下的 T/C 可变胞 BFGS `27761495/27761496` 均已完成。
相身份门禁确认 T 端 `c/a=1.06183`、Ti 相对赤道 O 面偏移 `0.22056 Å`，
C 端 `c/a≈1`、对应偏移约 `5.8×10^-10 Å`。随后 `27763724` 对映射和
平移对齐后的实际端点做独立静态计算：T/C 最大原子力为
`0.00446/0.05874 eV/A`，残余应力为 `0.914/0.0695 kbar`，按统一
`1.0 kbar` 阈值通过。曾按旧 `0.10 kbar` 标准提交的 T 端加严作业
`27763854` 在阈值修正后及时取消，未用其结果替换已验证端点。完整记录见
`validation/backend_smoke/cp2k_bto_k4_endpoint_gate_20260923.json`。

HF 通用模板同时显式传递 `PRESSURE_GPA` 到端点静态门禁和 VCNEB 控制器。
CP2K GaN `27756555` 使用旧模板的零外压与零应力端点，已在
`fmax=0.09612 eV/A` 时结束，但不能接收为生产结果。旧静态门禁验证的是
平移前终态（力 `4.82×10^-4 eV/A`、应力 `0.0496 kbar`）；路径自动平移后
实际终态哈希改变，CP2K 实空间网格下重新计算得到力 `0.18296 eV/A`、
应力约 `6.04 kbar`。这既不满足新 `1.0 kbar` 门禁，也证明旧门禁与有效
输入不一致。代码已加入有效端点哈希强制匹配，旧链只保留为零压诊断；
`validation/backend_smoke/cp2k_gan_effective_endpoint_mismatch_20260923.json`
记录了两套端点哈希和数值。当前单点探针只改变 CP2K `REL_CUTOFF`，用于判断
整体平移敏感性是否来自实空间网格精度；未通过前不重跑整条链。
同样，已完成的 QE GaN `27741412` 的原始 `vcneb_summary.json` 明确记录
`pressure_eV_per_A3=0`，其 `1.638143 eV` 势垒只属于该零压计算契约；
即使最大广义力 `0.091018 eV/A` 达标，也不可与 45.7 GPa 文献值
直接并列为同条件偏差。

### 运行中 ABACUS 链的几何诊断

`27744907` 的 `chain_step_0000` 到 `chain_step_0003` 显示 image 19 的体积从 `35.188` 增至 `39.020`、`43.093`、`47.353 Å³`，而相邻 image 保持在约 `35 Å³`。第 3 步日志的 `fmax` 仍为 `5.131 eV/A`。这不是一次普通回弹，而是单像晶胞信赖域失控；历史轨迹和作业继续保留用于审计。

控制器新增 `maximum_cell_step`：它按当前晶胞到候选晶胞的相对变形范数检查每一步，并在超过阈值时于 DFT 调用前拒绝候选。通用 HF 模板默认 `0.05`，配合 FIRE 的 8 次几何回溯；用第 2→3 步实际快照回放时，image 19 的相对变形为 `0.05548`、体积比 `1.09887`，已被门禁准确拒绝。这一检查与初始路径的 `maximum_deformation` 不同，专门覆盖运行时单像发散。

`examples/run_vcneb_ase.py` 和 HF 模板现在支持 `--resume-snapshot`/`RESUME_SNAPSHOT`。续算会先逐字节验证快照的初末端点身份，再在新工作目录中重新绑定 calculator；因此可以从最后一个安全完整快照续算，而不覆盖失控链。

## 其他当前结果

| 案例/后端 | 作业 | 当前判定 |
|---|---:|---|
| GaN / QE | 27741412 | 零压链已完成；最终广义力 `0.091018 eV/A`，按默认阈值可接受，能垒 `1.638143 eV` 仅属零压契约。旧 summary 未写显式 `converged`；不能与 45.7 GPa 文献直接比较。 |
| GaN / ABINIT（旧路径） | 27741516 | `max_steps_reached`，最终力 `9.537243 eV/A`，且 HGH-LDA 端点应力约 `-2791 GPa`；未收敛，不作为生产证据。 |
| GaN / LAMMPS | 27741574 | 最终力 `0.095441 eV/A`；仅是经典 Tersoff 对照，不能与 DFT/PBE 能垒等价比较。 |
| BTO / QE | 27741411 | 最终力 `0.052831 eV/A`，无内部势垒；按阈值可接受。 |
| BTO / VASP | 27741623 | 最终力 `0.095145 eV/A`，能垒 `0.042644 eV`；已收敛。 |
| GaN / ABACUS（新解析器重跑） | 27744907 | 已取消；运行至第 15 步，`fmax` 在 `3.17–3.27 eV/A` 平台并出现 `54190`、`26643 eV/A` 尖峰；image 19 单像体积发散，原目录保留为失败证据。 |
| GaN / ABACUS（guarded continuation） | 27749598 | 已写出 step 0--9 后因 image 10 的 `!FINAL_ETOT_IS` 解析缺口失败；原目录和完整快照保留，不作为物理失败。 |
| GaN / ABACUS（错误 launcher/外压） | 27755927 | 已取消；名义 32 MPI 实际为 32 个 singleton 写同一 image 目录，且错误使用零外压。轨迹不得作物理结果；高压端点应比较目标压力应力残差。 |
| GaN / ABACUS MPI canary | 27757793, 27757811, 27757818 | 前者暴露所有 rank 绑同一 CPU；后两者分别验证正确的 4/8-rank 静态计算和两个并行 worker。均非 VCNEB 生产结果。 |
| GaN / ABACUS 45.7 GPa 端点 | 27760738, 27760745, 27760794 | 两端 BFGS 优化及独立静态均已完成；4/4 原子，B4 四配位、B1 六配位，45.7 GPa 压力残差和原子力通过门禁。 |
| GaN / ABACUS 45.7 GPa 普通 VCNEB | 27760826 | 已完成并达到 `fmax=0.096844 eV/A`；能垒 `0.327366 eV/GaN`，单一内部峰 image 15，完整能量/力/应力及路径几何审计通过。 |
| BTO / CP2K endpoint relaxation（旧单位错误） | 27749624, 27749625 | 已取消并保留；旧 profile 把 400 Ry 错写成 400 eV，诊断应力 1459–3788 GPa，不进入结果矩阵。 |
| GaN / CP2K endpoint relaxation（旧单位错误） | 27749627, 27749628 | 已取消并保留；同一 cutoff 单位错误，不进入结果矩阵。 |
| BTO / CP2K endpoint relaxation（400 Ry 修正版） | 27749725, 27749726 | 已完成；初端点力很小但应力约 `5.97 kbar`，终端点应力约 `4.35 kbar`，均未通过严格静态门禁。 |
| BTO / CP2K endpoint tightening | 27749774 | 已完成；初端点 `max_generalized_force=4.24e-4 eV/A`、最大应力约 `0.0094 kbar`，通过严格门禁。 |
| GaN / CP2K endpoint relaxation（400 Ry 修正版） | 27749727, 27749728 | 已完成；两端最大应力约 `6.84/4.82 kbar`，未通过严格静态门禁。 |
| BTO/GaN / CP2K endpoint tightening | 27755902, 27755903, 27755904 | 已完成；旧源结构静态值通过当时的 `<0.1 kbar` 门禁，但 GaN 后续发现有效对齐终态哈希不同，故旧 GaN 门禁不得继续复用。 |
| BTO / CP2K 4×4×4 phase-correct endpoints | 27761495, 27761496, 27763724 | T/C 相身份、实际映射端点力与 `<1.0 kbar` 残余应力门禁均通过。 |
| BTO / CP2K corrected ordinary VCNEB | 27763906 | 9 像链以 `fmax=0.081828 eV/A` 收敛；无内部峰，`0.060805 eV/5-atom cell` 只是端点焓差，不是激活能垒。完整审计见 `validation/backend_smoke/cp2k_bto_k4_vcneb_20260923.json`。 |
| BTO / CP2K ordinary VCNEB（退化旧输入） | 27756554 | 数值收敛但两端均近立方，路径总长仅 `1.80e-4 Å`；保留为负例，不作为 T→C 生产结果。 |
| GaN / CP2K ordinary VCNEB（旧对齐） | 27756555 | 零外压链数值达到路径阈值，但有效对齐终态未通过静态门禁（`0.183 eV/A`、约 `6.04 kbar`），不接收为生产结果。 |
| GaN / CP2K endpoint identity retry | 27764113, 27764399, 27764401 | 前三次尝试均在 DFT 前被契约拒绝：先后暴露自动映射重排、`log_strain` 与关闭晶胞对齐不相容等输入组合问题；不产生物理结果。 |
| GaN / CP2K source-consistent production retry | 27767011 | `MAPPING=identity, ALIGN_TRANSLATION=0, ALIGN_CELLS=0, CELL_INTERPOLATION=linear` 已通过端点哈希预检并进入 29 像生产计算；结果待完整链审计。 |
| BTO/GaN / CP2K first VCNEB submission | 27756548, 27756549 | 仅运行 4 秒即退出；门禁报告父目录缺失，已作为启动契约失败保留，不覆盖 retry-1。 |
| GaN / ABINIT-HGH-LDA endpoint relaxation（错误 launcher） | 27749797, 27749798 | 已取消并保留；`srun` 注入 `pmi_args` 导致 MPI socket 等待。 |
| GaN / ABINIT-HGH-LDA endpoint relaxation（ABI 不匹配复现） | 27749871, 27749872 | 已取消并保留；Intel 2017/2021 混用，`mpiexec` 仍复现 `pmi_args`/socket 等待。 |
| ABINIT 32-rank launcher canary | 27749992, 27749995, 27749996, 27749997 | 仅诊断：缺 compiler runtime、错误 PMI 组合均失败；匹配 Intel 2017 + Hydra 的 `27749997` 成功返回 8.6.1。 |
| GaN / ABINIT-HGH-LDA endpoint relaxation（缺伪势路径） | 27750045, 27750046 | 预检后失败并保留；`pps=hgh` 未配 `pp_paths`，输入契约已补强。 |
| GaN / ABINIT-HGH-LDA endpoint relaxation（自动对称性试跑） | 27750056, 27750057 | B1 因 `chkorthsy` 末位晶格噪声失败，B4 为一致性取消；独立目录和日志保留。 |
| GaN / ABINIT-HGH-LDA endpoint relaxation（`nsym=1` 修正版） | 27750109, 27750110 | 已完成但静态门禁失败；端点体积塌缩至 `7.22/7.08 Å³`，应力 `59.3/60.8 kbar`。HGH-LDA 测试结果不进入 PBE/VASP 生产矩阵。 |

## 下一步

1. ABACUS GaN `27760826` 已完成审计；保留 45.7 GPa 的 `0.327366 eV/GaN` 结果和完整 provenance，不从 `27755927` 的 singleton 轨迹续算。
2. ABACUS `!FINAL_ETOT_IS` 解析修复已部署，但原 `27749598` 的
   `chain_step_0009.traj` 须先审计其 MPI 启动与端点契约，不能直接视作安全物理续算点。
3. CP2K BTO 的 4×4×4 T/C 端点已通过相身份与有效端点门禁，可在部署哈希契约后的独立目录启动 T→C VCNEB。GaN CP2K 先完成 `REL_CUTOFF` 平移不变性探针；有效端点仍不合格时不得重跑整链。若要与 45.7 GPa 文献比较，仍须重新构造同压力端点、静态门禁与独立生产目录。
4. 对 ABINIT GaN 先完成 HGH-LDA 端点重建并通过同一门禁，再决定是否重跑；当前 `max_steps_reached` 结果不得进入生产矩阵。新的端点任务必须使用匹配 Intel 2017 + Hydra 启动契约。若要与 PBE 结果比较，必须另行获得并固定 PBE 赝势，不能把 LDA/HGH 结果标成 PBE。
5. 生成统一后端状态表、路径图和文献比较数据；未达到阈值或非同一物理模型的结果不得进入“已验证生产矩阵”。
