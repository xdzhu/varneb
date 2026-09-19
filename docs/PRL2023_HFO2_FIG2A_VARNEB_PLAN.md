# PRL 2023 HfO2 图 2a：VARNEB + VASP 复现计算方案

## 1. 方案结论

本项目不使用 USPEX。路径优化由仓库自身的 VARNEB 完成：`vcneb/` 优化分数坐标与可变晶胞形变自由度，VASP 只对 VARNEB 给定的每个 image 计算静态总能、原子力和应力。

生产计算尽量复用论文的数值设置：PBE、600 eV、`4 x 4 x 4` Monkhorst-Pack、电子收敛 `1e-6 eV`、零外压；普通相变路径使用 20 个总 images，极化翻转路径使用 21 个总 images。VARNEB 的 `--n-images` 包含两个固定端点，因此 20 总像对应 18 个内部像，21 总像对应 19 个内部像。

赝势按本项目约定采用 VASP PAW-PBE 常规数据集：

- `PAW_PBE/Hf/POTCAR`；
- `PAW_PBE/O/POTCAR`；
- POTCAR 拼接顺序为 `Hf O`；
- 不使用 `Hf_pv`、`Hf_sv` 或其他半芯赝势。

论文没有报告具体 POTCAR，因此这是本复现方案的显式选择，必须在结果中记录 TITEL、发行版本和 SHA256，不能把它表述成论文已知参数。

## 2. 已核实的本仓库能力

当前实现满足本任务所需的关键条件：

- `vcneb/core.py` 同时优化原子和晶胞扩展坐标；
- VASP calculator 必须返回 energy、forces、stress，缺失应力会在预检阶段失败；
- VASP images 强制为静态调用：`IBRION=-1`、`NSW=0`、`ISIF=2`；
- 支持 FIRE、BFGS、LBFGS 和 BFGSLineSearch，本方案选择 FIRE；
- 支持零外压、普通NEB、climbing image、断点续算、逐image缓存和内部像并行；
- `--n-images` 是包含两个端点的总像数；
- `--initial-trajectory` 可保留作者路径的原子顺序、周期绕行和完整20/21-image初猜；
- 输入契约会冻结 INCAR/KPOINTS/POTCAR 指纹，并检查每个 image 的POSCAR往返精度。

## 3. 需要计算的路径

作者图 2a 数据表共有九条路径。相间路径一律20个总images；极化翻转为21个总images。

| 编号 | 路径 | 原子数 | f.u. | 总images | 内部images | 优先级 |
|---|---|---:|---:|---:|---:|---|
| P01 | T <-> PO | 12 | 4 | 20 | 18 | 第一批 |
| P02 | T <-> PO' | 12 | 4 | 20 | 18 | 第二批 |
| P03 | PO <-> M | 12 | 4 | 20 | 18 | 第一批 |
| P04 | T <-> M | 12 | 4 | 20 | 18 | 第二批 |
| P05 | PO' <-> M | 12 | 4 | 20 | 18 | 第二批 |
| P06 | PO' <-> PO | 12 | 4 | 20 | 18 | 第二批 |
| P07 | AO <-> T | 24 | 8 | 20 | 18 | 第三批 |
| P08 | AO <-> PO | 24 | 8 | 20 | 18 | 第三批 |
| P09 | +PO -> -PO，经T | 12 | 4 | 21 | 19 | 第三批 |

第一批先运行 P01 和 P03。P01 检验论文最关键的低势垒 T -> PO；P03 检验 M -> PO 的 169、189、196 meV/f.u. 数据分歧。第一批通过后再批量开展其余路径。

建议额外保留一个不并入图 2a 主结果的 P03b：构造保持 `X5y`、但反转/错配 `X2-` 的 M -> PO 映射，目标是复核 SM 图 S2 的约196 meV/f.u.分支。

## 4. VASP 统一输入

### 4.1 INCAR

所有端点和路径 images 使用同一电子结构设置。建议源 INCAR 为：

```text
SYSTEM = HfO2_PRL2023_VARNEB
PREC   = Accurate
ENCUT  = 600
EDIFF  = 1E-6
ALGO   = Normal
NELM   = 160

ISPIN  = 1
ISMEAR = 0
SIGMA  = 0.05

GGA    = PE
LASPH  = .TRUE.
LREAL  = .FALSE.
ADDGRID = .TRUE.
LMAXMIX = 4

IBRION = -1
NSW    = 0
ISIF   = 2
ISYM   = -1
SYMPREC = 1E-4

LWAVE  = .FALSE.
LCHARG = .FALSE.
```

说明：

- `IBRION/NSW/ISIF` 是 VARNEB 的计算器契约；晶胞和原子更新由外层 VARNEB/ASE 完成，不能再让VASP内部优化；
- 论文未报告 `ISMEAR/SIGMA/PREC/LREAL/LASPH/ALGO`，以上采用适合绝缘HfO2和可比应力的稳健默认值；
- 路径上可能出现较高能低带隙结构，因此主线使用 `ISMEAR=0, SIGMA=0.05`，不在不同images之间切换smearing；
- VASP中的 `EDIFFG` 不控制VARNEB收敛，故不作为路径参数；端点和路径阈值由外层优化器给出。

### 4.2 KPOINTS

主线计算严格采用论文报告的固定 Monkhorst-Pack 网格：

```text
HfO2 PRL2023
0
Monkhorst-Pack
4 4 4
0 0 0
```

12原子和24原子路径均使用这一网格，以保持对论文文字的最直接仿照。24原子Pbca路径另做一次较低密度网格的成本对照只能作为附加测试，不能替换主线 `4 x 4 x 4` 结果。

### 4.3 POTCAR与可追溯性

用许可目录中的常规 `Hf` 和 `O` PAW-PBE POTCAR生成源 `POTCAR`。提交计算前记录：

- 两段 `TITEL`；
- POTCAR总 SHA256；
- 每段赝势的 ENMAX；
- VASP版本和可执行文件路径；
- INCAR、KPOINTS、端点POSCAR和代码Git revision。

POTCAR不提交到Git，也不复制到公开产物目录。

## 5. 端点准备与优化

### 5.1 结构来源

优先使用作者 L23 仓库的结构和路径文件，而不是从数据库重新下载：

- T：P42/nmc；
- PO：Pca21；
- PO'：Pmn21；
- AO：Pbca；
- M：P21/c。

含AO的路径使用作者给出的24原子、8 f.u.相容晶胞；其余相间路径使用12原子、4 f.u.晶胞。不能在路径计算结束后才变更晶胞约定或归一化。

### 5.2 三阶段端点松弛

端点在 `pressure = 0 GPa` 下进行完全晶胞+原子优化，并保持目标空间群。使用仓库的 `scripts/relax_vasp_endpoint.py`，VASP仍为静态后端，ASE `BFGS(FrechetCellFilter)`负责更新结构。

建议按同一工作目录断点续算：

1. 预松弛：`fmax = 0.02 eV/Å`；
2. 中等收敛：`fmax = 0.005 eV/Å`；
3. 论文级端点：`fmax = 0.001 eV/Å`。

每阶段建议 `maxstep = 0.03 Å`，以300步为一个续算批次。最终端点必须同时检查：

- 外层最大广义力 `<= 0.001 eV/Å`；
- 原子最大力 `<= 0.001 eV/Å`；
- 残余应力与晶胞变化已稳定；
- 空间群仍为目标相；
- 五相相对能量顺序为 `M < AO < PO < PO' < T`。

命令骨架：

```bash
python scripts/relax_vasp_endpoint.py \
  --structure structures/T/POSCAR \
  --source-dir inputs/vasp_pbe_standard \
  --workdir runs/endpoints/T \
  --pressure-gpa 0 \
  --fmax 0.001 \
  --steps 300 \
  --maxstep 0.03 \
  --ncores 8 \
  --vasp-bin /path/to/vasp_std
```

后续批次增加 `--resume`。五个端点必须全部完成，不能用不同精度的端点混成一张能量图。

## 6. 初始路径构造

### 6.1 采用作者路径作为拓扑模板

首选方案不是重新猜原子映射，而是读取作者 `transition-poscars` 中的20/21个结构，保留：

- 原子排列；
- `X2-` 符号；
- 周期边界下的原子绕行；
- 作者给出的路径弯曲和晶胞演化趋势。

由于生产端点要在本项目的常规PAW-PBE势能面上重新优化，不能直接把作者原始首末结构固定为端点。准备脚本应把作者路径的内部偏离量移植到本项目的已松弛端点上，生成一个端点严格相等的 `.traj`：

1. 读取作者首末端点和全部中间结构；
2. 在作者端点之间建立相同原子顺序的线性分数坐标/晶胞基线；
3. 提取每个作者image相对基线的原子和晶胞偏离；
4. 在本项目松弛端点之间建立新基线，并按 `lambda` 加回偏离；
5. 强制第0和最后一个image精确等于本项目端点；
6. 检查最短原子距离、正体积、晶胞条件数、路径折叠和空间群端点。

生成的路径通过 `--initial-trajectory` 输入，因此运行时不再同时使用 `--mapping auto`、`--mic`、`--align-translation` 或 `--mode-guided`。

### 6.2 后备方案

若某条作者路径无法可靠移植，则按作者端点的原子顺序使用：

- `--mapping identity`；
- `--cell-interpolation linear`；
- 不启用 `--mic`；
- 显式检查 `X2-` 和 `X5y` 模式符号。

禁止仅凭自动最近邻映射决定 M <-> PO 路径，因为这正是论文势垒差异的来源。

## 7. VARNEB生产参数

主线普通VCNEB设置：

| 参数 | 值 | 理由 |
|---|---:|---|
| pressure | 0 GPa | 本征零压路径 |
| optimizer | FIRE | 当前实现中对长链和变量晶胞较稳健 |
| spring `k` | 0.10 | VASP驱动默认值；全项目固定 |
| total images | 20或21 | 与作者xlsx/POSCAR一致 |
| ordinary VCNEB fmax | 0.10 eV/Å | 第一阶段快速获得稳定、平滑的近似MEP |
| steps | 每批300 | 允许安全断点续算，不预设一次跑完 |
| climb | 第一阶段关闭 | 先获得平滑普通MEP，避免早期选错CI |
| CI refinement | 第二阶段，0.03 eV/Å | 对路径最高能区和鞍点进行论文量级精修 |
| cell interpolation | 由initial trajectory给出 | 保留作者路径拓扑 |
| image retries | 0 | SCF/MPI错误不静默掩盖 |

论文采用“各image RMS force `<0.03 eV/Å`”；本项目当前优化停止量是最大广义力，因此采用分阶段标准：普通VCNEB先以最大广义力 `<=0.10 eV/Å` 获得稳定路径，再开启CI并收敛到最大广义力 `<=0.03 eV/Å`。最终报告还应另算每个image的原子/广义NEB力RMS；只有CI精修结果满足相应RMS标准时，才作为与论文势垒的最终定量对照。

普通阶段命令骨架：

```bash
python examples/run_vcneb_vasp.py \
  --initial runs/endpoints/T \
  --final runs/endpoints/PO \
  --initial-trajectory paths/T_to_PO_20.traj \
  --workdir runs/paths/T_to_PO_n20 \
  --n-images 20 \
  --pressure-gpa 0 \
  --optimizer FIRE \
  --k 0.10 \
  --fmax 0.10 \
  --steps 300 \
  --no-climb \
  --image-workers 5 \
  --ncores 8 \
  --vasp-bin /path/to/vasp_std
```

未收敛时使用相同参数和 `--resume`，不要再次传入 `--initial-trajectory`。普通路径达到 `0.10 eV/Å` 并通过路径诊断后，可复制为独立CI工作目录，或在同一可追溯链上以 `--resume`、去掉 `--no-climb`、`--fmax 0.03` 进行CI精修。普通阶段结果用于筛查路径拓扑和异常，不作为最终精确势垒。

## 8. 计算顺序

### 阶段 A：输入与单点基线

1. 生成统一的 INCAR/KPOINTS/POTCAR；
2. 在 T、PO、M 三个12原子端点上做真实静态SCF；
3. 在 AO 24原子端点上做真实静态SCF；
4. 检查SCF收敛、有限energy/forces/stress、POTCAR标签和输入指纹；
5. 对一条20-image初始链执行 `--validate-only`；
6. 对最畸变的一个内部image做真实静态SCF canary。

### 阶段 B：端点优化

依次完成 T、PO、PO'、M、AO；核对晶格常数、空间群和相对能量。若端点能量排序错误，停止路径生产，先检查POTCAR、k点和结构映射。

### 阶段 C：第一批路径

1. P01 T <-> PO，20总像；
2. P03 M <-> PO，20总像；
3. 必要时P03b另一模态映射。

先普通VCNEB，后可选CI。第一批同时验证作者原始路径移植脚本和20-image并行配置。

### 阶段 D：其余12原子路径

运行P02、P04、P05、P06。各路径独立工作目录和cache namespace，不跨参数复用image缓存。

### 阶段 E：24原子和极化翻转路径

运行P07、P08、P09。P07/P08仍使用 `4 x 4 x 4`，势垒按8 f.u.归一化；P09用21总像，势垒按4 f.u.归一化。

## 9. 并行与资源建议

12原子、600 eV、64个k点的单image可先以8 MPI核计时。若一个节点有40核，建议：

- `image-workers = 5`；
- 每个VASP image使用8核；
- 20-image路径的18个内部像约分4批完成一次完整链评估；
- 21-image路径的19个内部像同样约4批。

24原子路径先比较每image 8核与16核的实际壁钟时间，再固定生产配置。不能让多个image共享同一VASP目录；Slurm并行时每个VASP命令必须获得排他的MPI资源。端点由控制器固定并缓存，不重复参与每一步内部像并行。

## 10. 收敛与结果验收

### 10.1 电子结构数值检查

在正式全路径前，对T、PO、M至少进行：

- k网格 `3 x 3 x 3`、`4 x 4 x 4`、`5 x 5 x 5` 单点比较；
- `ENCUT = 500, 600, 700 eV` 单点比较。

主线结果仍固定为论文的600 eV与`4 x 4 x 4`。若600/4与更高设置的相对能量差超过1 meV/f.u.，应把它记录为论文参数的数值敏感性，并另建高精度验证支线，不能偷偷改写主线参数。

### 10.2 端点验收

- 空间群正确；
- 晶格常数与Table S1接近，初始目标误差 `<0.5%`；
- 相对能量排序正确；
- M为零点时尽量接近 `AO 72.9, PO 84.3, PO' 142.9, T 166.3 meV/f.u.`。

### 10.3 路径验收

- 普通VCNEB最大广义力 `<=0.10 eV/Å`，用于确认路径拓扑和进入CI阶段；
- 最终CI路径最大广义力 `<=0.03 eV/Å`；
- 每image RMS力另行计算并与论文阈值比较；
- 路径没有折返、重复images、负体积或异常近距离；
- 最高内部image高于两个端点，或明确报告为barrierless/unresolved；
- 相邻image间距没有严重集中；
- 普通VCNEB和CI结果分别保存；
- 以总能差除以4或8 f.u.，报告正向和反向势垒。

目标对照值：

| 路径 | 主要论文目标 |
|---|---:|
| T -> PO | 2.5 meV/f.u. |
| T -> M | <5 meV/f.u. |
| T -> PO' | <5 meV/f.u. |
| PO' -> M | 28 meV/f.u. |
| PO' -> PO | 37 meV/f.u. |
| M -> PO，X2-保持 | 169 meV/f.u. |
| M -> PO，X5y保持 | 196 meV/f.u. |

作者公开xlsx可作为逐image形状参考，但其中 M -> PO 为189.425 meV/f.u.、T -> M为5.1、PO' -> M为25.5，且含Pbca的两个sheet端点文字标签反置。因此验收报告必须同时列“论文排版值”和“作者仓库派生值”，不能只选对我们最有利的一套。

## 11. 建议目录

```text
validation/prl2023_hfo2_fig2a_vasp_pbe/
  provenance/
  inputs/vasp_pbe_standard/
  structures/author_raw/
  structures/relaxed_endpoints/
    T/ PO/ PO_prime/ AO/ M/
  paths/initial/
  runs/endpoints/
  runs/P01_T_PO_n20/
  runs/P03_M_PO_n20/
  runs/P02_T_POprime_n20/
  runs/P04_T_M_n20/
  runs/P05_POprime_M_n20/
  runs/P06_POprime_PO_n20/
  runs/P07_AO_T_n20/
  runs/P08_AO_PO_n20/
  runs/P09_PO_switch_n21/
  reports/
```

## 12. 在真正提交计算前还需完成的实现工作

1. 增加一个作者多结构POSCAR到ASE trajectory的转换与端点校正脚本；
2. 为全部九条路径生成manifest，记录实际结构方向、原子数、f.u.和image数；
3. 在结果导出中增加逐image RMS原子力/广义NEB力，和论文收敛定义直接对照；
4. 生成统一的VASP输入模板，但POTCAR只在许可环境中组装；
5. 为第一批P01/P03生成 `--validate-only` 和真实SCF canary提交脚本；
6. canary通过后才提交长时间VARNEB生产任务。
