# VASP VCA 后端：实现、验证案例与适用边界

**状态：** VASP 后端已支持一个物理晶格位点由多个重合的 VCA 组分表示；已完成
Ba0.5Sr0.5TiO3（BST50）和 Pb(Ti0.5Zr0.5)O3（PZT50）的材料路径。

这是一份后端验证说明，不把 VCA 本身表述为 VARNEB 的新算法，也不把等价替位的
平均势近似外推为真实缺陷化学。

## 1. 实现契约

VARNEB 的物理 `Atoms` 仍保持一个原子对应一个物理晶格位点。对一个虚拟位点，例如
BST50 的 A 位 Ba/Sr，`VirtualCrystalCalculator` 在调用 VASP 前将它展开为相同分数
坐标上的 Ba 与 Sr 两个组分；其余原子不复制。VASP 读取与该展开顺序相同的 POTCAR 和
`VCA` 权重；返回时，重合组分上的力求和回映射到一个物理位点：

\[
\mathbf F_{\mathrm{virtual}}=\sum_{c\in\mathrm{components}}\mathbf F_c.
\]

能量和应力直接属于整个 VCA 构型。这样 VCNEB 核心始终处理不重叠的物理结构，而不是
把重合伪原子暴露给插值、mapping 或 image manager。

调用接口为：

```text
python examples/run_vcneb_vasp.py ... \
  --vca-virtual-symbol Ba --vca-components Ba Sr
```

或对 PZT50：

```text
--vca-virtual-symbol Ti --vca-components Ti Zr
```

运行前 adapter 现在会验证：

1. 一个物理 virtual symbol 恰好出现一次；
2. `INCAR` 有 `VCA` 向量，且长度与 POTCAR dataset 数一致；
3. POTCAR 中恰有一个、与 components 顺序一致的连续组分块；
4. 组分权重非负且和为 1，所有非虚拟 dataset 的权重均为 1；
5. 经过展开后的物理结构可被构造。

VASP 输入由 `ExplicitPotcarVasp` 在 ASE 写输入后重新恢复指定的源 POTCAR，避免 ASE
按默认赝势选择重写掉经审计的 dataset 顺序。每个 image 仍使用独立目录；端点静态
结果经过 fingerprint 和结构 hash 审计后缓存，后续迭代只计算内部 image。

## 2. 已验证的范围

当前验证的是**等价替位、单虚拟位点**的 VCA：

| 系统 | 虚拟位点 | PAW-PBE datasets | 路径 | 总 image / 内部 image |
|---|---|---|---|---:|
| BST50 | A 位 Ba/Sr = 0.5/0.5 | Ba_sv, Sr_sv, Ti_sv, O | tetragonal → cubic | 7 / 5 |
| PZT50 | B 位 Ti/Zr = 0.5/0.5 | Pb_d, Ti_sv, Zr_sv, O | tetragonal → cubic | 7 / 5 |
| PTO | 无 VCA，对照 | Pb_d, Ti_sv, O | tetragonal → cubic | 7 / 5 |

这三条路径均使用 VASP 6.3.2、PBE、`ENCUT=600 eV`、`4×4×4` k 网格、静态 image
调用（`IBRION=-1`, `NSW=0`, `ISIF=2`, `ISYM=0`）。端点使用
`FIRE(FrechetCellFilter)`，对称约束为 `FixSymmetry`；NEB 采用普通无 CI 的
FIRE、默认 `fmax=0.10 eV/A`。每一项中的 7 是总 image 数，不是 7 个中间插帧。

| 案例 | 端点最大原子力（T / C，eV/A） | 最终广义力（eV/A） | T→C 端点焓差（meV/f.u.） | 内部势垒？ |
|---|---:|---:|---:|---|
| BST50 VCA | 0.01750 / 0.00000 | 0.07780 | 31.569 | 否 |
| PTO | 0.00574 / 0.00000 | 0.09227 | 221.640 | 否 |
| PZT50 VCA | 0.02265 / 0.00000 | 0.09710 | 27.146 | 否 |

这里的数值来自 committed 的 summary/audit：

- `results/bst50_vcneb_final/`
- `results/pto_vcneb_final/`
- `results/pzt50_vca_vcneb_final/`

三个 audit 均为 `status=ok`，并明确给出 `has_interior_barrier=false`。因此 summary
中的 `barrier_enthalpy_eV` 字段在这些单调路径上等于终点上升，**不得**解读为活化能或
CI 鞍点。BST50 的完整 SVG/PDF/PNG 图和 source-data CSV 已归档；PTO/PZT50 目前归档
的是可复核的 endpoint summary、path summary 和 audit，图稿应由这些 summary 再生成。

另有 `results/bst50_vca_static_comparison.*`：它只是在既有纯 BTO 路径几何上做 VCA
静态单点比较，几何没有为 BST50 重新弛豫。它可用于说明平均 A-site 势改变能量轮廓，
不能替代上表的 BST50 VCA-VCNEB 结果。

## 3. 不应外推的范围

- 当前 adapter 只实现/验证一个虚拟位点。多位点、多组分混合需要扩展 manifest、
  force pullback 和测试，不能通过多次传参假装支持。
- 只验证了 Ba/Sr 与 Ti/Zr 的等价替位。Hf4+→Y3+ 属于异价替代：它需要明确的电荷补偿
  模型，通常是氧空位、显式缺陷超胞或另行验证的平均缺陷近似。不能把本 VCA 流程直接
  用作 Y:HfO2 的定量路径或氧空位迁移结论。
- VCA 平均化了化学无序，不能给出局域配位、短程有序、Y--VO 缔合、显式缺陷态或构型
  展宽。此类问题应以显式超胞/SQS 作基准。
- 三个当前路径均没有内部势垒，适合验证 calculator contract、端点缓存和路径几何，
  不足以单独验证 VCA 下 CI 的鞍点定位能力。

## 4. 可重复运行入口

- `scripts/setup_batio3_vca_case.py`：生成 BST50 的 Ba/Sr VCA 输入和 manifest；
- `scripts/setup_perovskite_vasp_case.py`：生成 PTO 或 PZT50 的五位点钙钛矿输入；
- `scripts/relax_vasp_vca_endpoint.py`：物理结构上的对称保持原子/晶胞端点弛豫；
- `cluster/cu17_bst50_vca_endpoint_relax.sh` 与
  `cluster/cu17_bst50_vca_vcneb.sh`：BST50 的端点和路径入口；
- `cluster/cu17_pto_pzt50_serial_workflow.sh`：在 cu17 上先完成 PTO，再完成 PZT50
  VCA；它仅调度 5 个内部 image 的串行 VASP 调用。

所有这些脚本都要求显式 `RUN_DFT=1`、受许可的 `VASP_BIN` 和确定的
`VCNEB_GIT_REVISION`；默认不启动计算。
