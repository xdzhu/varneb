# USPEX VCNEB 对照与模式引导路径设计

日期：2026-09-12

## 结论先行

`H:\ReSearch\VCNEB\USPEX\USPEX_v10.6.tar.gz` 已完成归档检查和解压。
该包不是可阅读的 USPEX 源码包，而是：

- `install.sh`；
- `README`；
- 两个许可证文件；
- 一个约 878 MB 的 MATLAB Runtime 自解压 ELF 文件 `USPEX_MATLABruntime.install`。

归档内没有 `.m` 源码、Python 源码或独立的 VCNEB 源代码目录。因此本项目不应从二进制中提取或改写 USPEX 实现。本文对 USPEX 的算法分析依据公开的 VC-NEB 论文和官方手册，代码设计保持干净实现和许可证边界。

## 1. USPEX 10.6 包的实际形态

解压目录：

```text
H:\ReSearch\VCNEB\USPEX\extracted_10.6\USPEX_v10.6\
  install.sh
  MCR_license.txt
  README
  USPEX_license.txt
  USPEX_MATLABruntime.install
```

`README` 明确要求 MATLAB Runtime R2016b，并列出 Python 依赖 `numpy scipy spglib pysqlite ase`。`install.sh` 会把安装结果放到 `application/archive/`，并设置 `MCRROOT`、`PATH` 和 `USPEXPATH`。这说明它是 MATLAB 编译后的应用发行包，而不是可供我们复用的 Python 库。

USPEX 许可证还明确禁止再分发、反汇编、反编译和逆向工程。因此我们只使用论文、公开手册和黑盒行为作为参考，不复制内部代码。

## 2. USPEX VCNEB 的算法结构

Qian 等人的 VC-NEB 把晶胞应变和原子分数坐标放到同一扩展构型空间：

```text
X = (epsilon_11, ..., epsilon_33, r_1, ..., r_N)
dim(X) = 9 + 3N
```

恒压时优化焓面：

```text
H = E + P Omega
```

晶胞部分由应力得到，原子部分由 Hellmann-Feynman force 得到；原子力通过晶胞度量张量映射到与应变变量兼容的扩展空间。NEB 力由真实焓梯度的路径垂直分量和弹簧力的路径平行分量组成。

论文和官方手册体现的 USPEX 流程为：

1. 端点晶胞和原子序列进行旋转规避与重排；
2. 端点或用户提供的中间结构生成初始 images；
3. 外部计算器计算各 image 的能量、原子力和应力；
4. 计算扩展空间切线、垂直真实力和弹簧力；
5. 更新晶胞应变和原子坐标；
6. 根据收敛标准重复。

USPEX 手册还提供：

- 固定或可变 image 数；
- 固定或可变弹簧常数；
- Steepest Descent 和 FIRE；
- 只放松原子、只放松晶胞或全放松；
- CI/DI image；
- 定期 restart 输出；
- 晶胞旋转规避；
- 端点原子序列匹配。

这些是我们必须达到的基线，而不是应当照搬的内部实现。

## 3. 与当前实现的逐项对照

| 能力 | USPEX VCNEB | 当前 `vcneb` | 判断与改进 |
|---|---|---|---|
| 扩展坐标 | 9 个应变 + 3N 分数坐标 | 变形梯度 `F` + 3N 分数坐标 | 数学上同一类；需正式固定 row-vector/column-vector 约定 |
| 恒压焓 | `E + P Omega` | 已支持 `pressure_gpa` 与焓 | 加入压力/应力符号的单元测试和结果元数据 |
| 原子力/应力 | 内部黑盒 | ASE calculator contract | 我们的 adapter 边界更清楚、可跨计算器 |
| 计算器 | VASP、GULP、QE 等有限集成 | VASP、ABACUS 已有 adapter | 继续加入 QE、LAMMPS/ML，统一错误与单位报告 |
| 初始路径 | 线性插值或用户中间结构 | fractional/cell 插值、MIC | 加强端点匹配、路径质量诊断和 mode-guided 初猜 |
| 晶胞旋转 | Euler + mirror 全局搜索 | 当前只做 polar rotation gauge | 加入可审计的 rotation/mirror 搜索和旋转量报告 |
| 原子匹配 | 内部自动重排 | 已有独立 fixture/mapping 脚本 | 做成稳定 API，不再依赖用户手工排序 |
| image 数 | 可变 image | 当前固定 image | 实现 split/merge，保留每次变更的 provenance |
| 弹簧 | 可变弹簧 | 当前 scalar 或 per-segment `k` | 增加能量自适应和 path-spacing 控制 |
| CI/DI | CI、DI、多 CI/DI | 基本 CI | 将 climbing policy 独立成策略对象 |
| 优化器 | SD/FIRE | ASE FIRE/BFGS/LBFGS | 保留 ASE，增加 step limit、预条件和状态恢复 |
| 重启 | USPEX restart 目录 | trajectory complete-snapshot resume | 增加 JSON manifest、calculator provenance 和 optimizer state |
| 并行 | image 级并行 | ASE `world` image ownership | 增加 Slurm array 和每 image 独立任务模板 |
| 可复现性 | 依赖版本/工作流较重 | 可以把配置与原始输入完整打包 | 这是 OpenVCNEB 的主要优势之一 |
| 许可证/分发 | 注册、禁止再分发和逆向 | 可采用宽松开源许可证 | 发布前仍需审计依赖许可证 |

## 4. 我们超过 USPEX 的具体方向

不把目标定义成“每一步都比 USPEX 快”，而是定义成更现代、更容易复现的研究工具：

### 4.1 计算器无关

核心只要求：

```python
energy, forces, stress = backend.evaluate(atoms)
```

后端负责输入文件、运行命令、单位转换、收敛状态、错误诊断和输出解析。VCNEB 核心不应该知道 VASP 的 INCAR，也不应该知道 ABACUS 的 `INPUT` 语法。

推荐后端顺序：

1. VASP；
2. ABACUS；
3. Quantum ESPRESSO；
4. LAMMPS/ML potential；
5. 用户自定义 Python callable。

### 4.2 方法显式化

不要把所有变胞 NEB 都叫作 VCNEB。配置中应明确：

- `method: qian_vcneb`：恒压焓面、Qian 型扩展梯度；
- `method: gssneb`：广义固体 NEB；
- `method: fdneb`：有限变形，显式选择 Cauchy、PK1、PK2 或 hydrostatic stress；
- `method: fixed_cell_neb`：固定晶胞基线。

这样才能在论文中公平比较不同度量、应力测度和外部边界条件。

### 4.3 结构准备成为一等公民

把以下步骤做成可报告、可复现的 API：

- 元素/原子一一对应检查；
- species-wise Hungarian mapping；
- fractional MIC 搜索；
- cell rotation/mirror 搜索；
- 端点体积、角度、应变变化报告；
- 每个 image 的最小原子距离和 cell determinant 检查；
- 对超胞和非原胞的映射记录。

### 4.4 集群优先

本机只负责配置、静态校验、可视化和日志汇总。每个 image 使用独立工作目录，计算器运行由 Slurm/PBS/srun/mpi 启动。任何 DFT smoke 和生产 VCNEB 都应在 235/cu05 或正式集群执行。

## 5. 指定振动模式的 NEB：正确的语义

“指定振动模式”不是单一标准算法，至少有三种不同含义，必须在软件中分开：

### A. 模式引导初始路径，推荐默认

先做普通端点插值，再在内部 image 上叠加端点为零的模式扰动：

```text
R(lambda) = R_linear(lambda) + A * s(lambda) * e_mode
```

其中 `s(0)=s(1)=0`，例如 `sin(pi*lambda)` 或 `4*lambda*(1-lambda)`。随后运行完全自由的 NEB。

优点：

- 不改变最终“全空间 MEP”的定义；
- 能生成不同 `+mode` / `-mode` 初猜，探索不同路径分支；
- 适合软模、八面体转动、极化模和离子迁移模式；
- 不要求计算器实现模式约束。

当前代码已经实现：

```python
from vcneb import Mode, mode_guided_path

mode = Mode.from_file("soft_mode.dat", n_atoms=len(initial))
images = mode_guided_path(
    initial,
    final,
    n_images=9,
    mode=mode,
    amplitude=0.20,
    envelope="sin",
    mic=True,
)
```

这里的 `amplitude` 单位为 Angstrom，模式在内部归一化；端点保持精确不变。

### B. 模式子空间约束，不应默认开启

如果只允许结构在一个或多个模式张成的子空间内移动，那么得到的是“约束 MEP”，不是完整势能面上的 MEP。它适合：

- 构造一维或低维能量曲线；
- 验证 Landau/软模模型；
- 分离一个已知机制的能垒贡献。

它必须明确报告 `constrained_mep: true`，并检查端点差向量是否落在所给子空间内。对于真实晶体相变，单个固定-cell Gamma 模式通常不足以描述晶胞应变、有限波矢模式或重构过程。

### C. 模式投影和后处理，必须提供

给定已收敛的 NEB/VCNEB 路径，计算每个 image 的位移在指定模式上的投影：

```python
from vcneb import project_path_onto_modes
q_mode = project_path_onto_modes(images, initial, [mode])
```

这对应 modal NEB 文献中用正常模投影分析迁移路径的思路。它回答的是“哪个模式参与了已经找到的路径”，不等于“强行沿该模式搜索路径”。

### D. 变胞模式

普通 phonon eigenvector 只有原子自由度；晶体相变还可能需要应变自由度。`Mode` 已预留：

```python
Mode(atomic=atomic_mode, cell=deformation_mode)
```

其中 `cell` 是无量纲变形梯度方向。后续需要定义正式的扩展模式文件格式，记录：

- atom order；
- q-point 和 commensurate supercell；
- frequency；
- mass weighting convention；
- atomic vector units；
- cell strain basis 和 stress measure。

## 6. ABINIT 约束的准确理解

ABINIT 当前文档提供两类几何约束：

### 6.1 分量固定

`iatfix` 固定原子的全部三个方向；`iatfixx`、`iatfixy`、`iatfixz` 分别固定某些原子的 x/y/z 分量。这与当前新增的：

```python
atom_mask = np.ones((n_atoms, 3))
atom_mask[atom_index, axis] = 0.0
```

语义一致。

需要注意，ABINIT 官方文档自己警告：方向固定在 `geoopt=viscous` 时按 Cartesian 坐标处理，而其他部分优化器/动力学情形按 reduced 坐标处理；使用方向固定时还会限制可用对称操作。我们的实现应始终显式声明坐标系，默认在 VCNEB 扩展 Cartesian metric 中投影。

### 6.2 线性组合约束

`nconeq`、`iatcon` 和 `wtatcon` 约束原子坐标的线性组合。投影后的力满足：

```text
sum_I,mu W_(mu,I) * F'_(mu,I) = 0
```

因此它可以保持两个原子的相对高度，也可以定义平面、质心等约束。但它不是“自动识别并沿任意 phonon eigenvector 做 NEB”。要实现任意模式子空间，需要我们在外部定义模式基、质量度量和正交投影，并明确这得到的是 constrained MEP。

## 7. 当前代码已落地的新增接口

### 原子分量约束

```python
chain = VCNEB(images, atom_mask=mask)
```

`mask.shape == (n_atoms, 3)`，值为 1 的分量参与优化，值为 0 的分量保持当前 image 的坐标，并从真实力、切线、弹簧和反应坐标中排除。

### 模式引导路径

```python
mode = Mode(
    atomic=np.asarray(...),
    cell=np.asarray(...),       # optional 3x3 deformation mode
    label="soft-mode-1",
    frequency=-2.4,
)
images = mode_guided_path(
    initial,
    final,
    n_images=9,
    mode=mode,
    amplitude=0.20,
    envelope="bell",
)
```

支持 `.dat/.txt`、`.json` 和 `.npz` 模式文件；支持质量加权输入转换和可选的平移分量去除。

### 模式投影

```python
modal_coordinate = project_path_onto_modes(images, initial, mode)
```

返回形状为 `(n_images, n_modes)` 的模态坐标，可用于比较 `+mode`、`-mode`、普通线性插值和最终收敛路径。

## 8. 还需要补强的功能

优先级从高到低：

1. 补一个可直接运行的模式文件示例；
2. 增加任意方向投影矩阵，而不仅是轴向 `atom_mask`；
3. 增加严格 `constrained_mep` 模式，并在端点不属于子空间时拒绝运行；
4. 加入多模式正交化、质量度量和 q-point/supercell 元数据；
5. 加入 mode ensemble：多个模式、正负振幅、不同 cell strain seed；
6. 完成可变 image split/merge；
7. 完成 Qian/GSSNEB/FD-NEB 三套显式 method class；
8. 在集群上做 VASP 和 ABACUS force+stress smoke；
9. 用 Ar/GULP、Si diamond-beta-tin 和 HfO2 T-PO 做跨方法 benchmark；
10. 发布 `OpenVCNEB` 的 `pyproject.toml`、CLI、文档、测试、输入样例和 DOI 归档。

## 参考资料

1. Qian et al., *Variable cell nudged elastic band method for studying solid-solid structural phase transitions*, CPC 184, 2111 (2013), DOI: 10.1016/j.cpc.2013.04.004.
2. USPEX VCNEB manual: https://uspex-team.org/online_utilities/uspex_manual_release/EnglishVersion/uspex_manual_english/vcneb.html
3. ABINIT GeoConstraints: https://docs.abinit.org/topics/GeoConstraints/
4. ABINIT Relaxation variables: https://docs.abinit.org/variables/rlx/
5. ABINIT TransPath: https://docs.abinit.org/topics/TransPath/
6. Gordiz et al., *Enhancement of ion diffusion by targeted phonon excitation*, Cell Reports Physical Science (2021).
7. Pham et al., *Correlated Terahertz phonon-ion interactions control ion conduction in a solid electrolyte*, arXiv:2305.01632.
8. `modalNEB` public analysis toolkit: https://github.com/kgordiz/modalNEB
