# VCNEB 完全体执行路线（2026-09-13）

本路线把当前项目从研究开发版推进到“可复核、可恢复、可发布候选版”。当前主材料切换为 BaTiO3 四方相（T）→立方相（C），ABACUS 统一使用 100 Ry 与 Ba/Ti/O 全套 10 au DZP 轨道；HfO2 暂作为后续高难度案例，不再阻塞当前主线。

## 总体验收目标

- 理论、坐标度量、应力符号和收敛判据在代码、测试和文档中一致。
- BTO T/C 端点在相同 calculator 设置下独立变胞弛豫并可复核。
- BTO 普通 VCNEB 在至少一组稳定的 image 数下达到目标广义力；CI 只在普通路径稳定后开启。
- BTO 的 image 数、弹簧、优化器和电子结构精度有最小收敛矩阵。
- 运行可从完整快照恢复，失败 image 不覆盖有效结果；所有作业有节点、核数、参数、代码版本和输出归档。
- HfO2 的旧 60 Ry 结果全部标记为诊断结果，不作为生产能垒。

## 阶段与闸门

### G0：配置纠正与基线（已完成/进行中）

- [x] 将 BTO 模板统一到 `100 Ry + Orb-DZP-10au`。
- [x] BTO 模板支持 `DIRECTION=tetragonal_to_cubic`，并隔离新端点/新路径目录。
- [x] 旧 HfO2 `60 Ry` 结果降级为诊断记录。
- [x] 已同步最新源代码和 BTO/HfO2 生产模板到 hf 共享运行目录；当前工作树尚未提交 Git，运行 manifest 以工作树状态和 Slurm 作业号为准。

### G1：BTO 端点重弛豫（当前任务）

- [x] BTO T/C endpoint 与 5/7/9-image 路径已完成（`100 Ry + DZP-10au`），现有结果继续作为生产证据；误发起的 native endpoint 方法重算 `27676929`/`27676930` 已终止，不纳入证据，不再重跑 BTO。
- [x] 两端均使用 100 Ry、4×4×4、SCF 1e-8、Ba/Ti/O 全套 10 au DZP。
- [x] 端点均收敛：cubic 最大广义力 `7.48e-4 eV/A`，tetragonal `9.05e-3 eV/A`，cell determinant 均为正；结构和摘要已归档到 `outputs/batio3_t_to_c_pbe100_dzp10au/`。
- [ ] 若某端点不收敛，保留完整日志，调整 optimizer/maxstep 或电子混合，不直接进入 NEB。

### G2：BTO 普通 VCNEB 路径

- [x] 已使用新端点、T→C 方向、7 images、普通 NEB、CI 关闭提交预收敛作业 `27675388`；其回弹分支已保留并停止。
- [ ] 检查每轮 fmax、image 能量、cell/应力、最短距离、segment cosine 和是否折返。
- [x] 已从完整外部轨迹恢复（`27675509`，32 MPI，`maxstep=0.002`），验证不重新插值；FIRE/LBFGS 两种续算均在约 `0.05 eV/A` 附近平台/回弹。
- [x] 已生成并归档 summary、trajectory、每-image 输入/输出、Slurm manifest/preflight 和人类可读摘要；静态审计 `27675643` 判定该直连路径无内部势垒，因此 CI 继续关闭。
- [x] 追加验证 `27675704`（32 MPI、`maxstep=0.003`）仅由 `0.029613` 降到 `0.029365 eV/A` 后回弹至 `0.031914 eV/A`；已安全停止并保留 `step0..2` 快照，未启动 CI。

### G3：image 数与数值收敛

- [x] 已完成同一 calculator/端点/路径参数下的 5→7→9 image 普通 VCNEB 收敛矩阵；三组均达到 `fmax≤0.02 eV/A`。
- [x] 已提交同设置的 5/9-image 普通分支（原始作业 `27675759`/`27675760`，32 MPI）；两条串行轨迹均保留并分别转入并行恢复。
- [x] 5-image 分支已按新规则从最佳 step 11 续算 12 步（`27676180`，4×32 MPI）；step 3 短暂回弹至 `0.036197 eV/A` 后继续下降，step 12 达到 `0.0231289 eV/A`，因此确认此前不能因单次回弹停止；仍未达到 `0.02`。
- [x] 5-image 随后从 step 12 完整轨迹续算（`27676251`，4×32 MPI），step 7 达到 `0.0191050 eV/A`，通过普通路径广义力门槛；审计确认几何有效、能垒 `0.0871289 eV` 且无内部势垒。
- [x] 7-image 从最佳完整 step 1 续算（`27676310`，4×32 MPI）；经历 step 1--5 的回弹窗口后重新下降，step 17 达到 `0.0196710 eV/A` 并通过普通路径门槛；审计确认几何有效、能垒 `0.0871289 eV` 且无内部势垒。
- [x] 9-image 串行分支在完整 step 9 后完成首段并行恢复（`27675981`，30 步，末步 `0.0739323 eV/A`），随后 `27676207` 再续 30 步至 `0.0289920 eV/A`；能垒仍为 `0.0871289 eV` 且无内部势垒，已从 step 60 继续提交 `27676513`。
- [x] 9-image 最终续算 `27676513` 达到 `0.0199074 eV/A`；5/7/9 汇总比较为 `barrierless-consistent`，能垒均 `0.0871289 eV`，最高 image 坐标差在协议容差内。
- [x] 以能垒差、最高 image 坐标、几何有效性和广义力共同完成 image 收敛判断；三组无内部势垒，未发现需 11/13 image 的高曲率欠采样证据。
- [x] 回弹处理规则已修正：单次或短暂回弹不再触发停止；每个分支至少继续观察 3--5 个完整 optimizer step，记录最佳完整快照。仅当连续回弹/平台且最佳值在观察窗口内不再改善时，才安全停止并转静态审计。
- [ ] 若高曲率区域局部欠采样，再做 11 或 13 images；不盲目全路径加密。

### G4：CI 精修与方向检查

- [x] 普通路径收敛后完成 CI 闸门审计；5/7/9 均 `has_interior_barrier=false`，因此 CI 精修有明确物理理由 withheld，记录于 `outputs/batio3_t_to_c_pbe100_dzp10au/bto_ci_gate.json`，未强行启动无意义的 CI。
- [x] 最高 image、路径坐标和几何审计已记录；由于不存在内部 saddle，不报告伪造的负切线曲率或 saddle 残差。
- [x] C→T 反向 5-image 普通 VCNEB 已完成（`27676299` + `27676485`，4×32 MPI）；在回弹观察后达到 `0.0195577 eV/A`，能垒 `0`、反应焓 `-0.0872206 eV`，无内部势垒，结果已归档。

### G5：电子精度与可恢复性

- [x] 已完成生产与加严电子精度对照：100 Ry/DZP-10au 下 4³/1e-8 与 6³/1e-9 的 5-image 静态路径及两端重弛豫；加严端点反应能 `0.0724417 eV` 与静态路径 `0.0723553 eV` 相差 `0.00009 eV`，但相对生产 4³ 结果约有 `0.01476 eV` 敏感性，记录于 `bto_precision_sensitivity.json`。
- [ ] 人为中断并从最后完整 iteration 恢复；验证不重新插值、不重复覆盖有效结果。
- [ ] 测试单 image 失败、SCF 异常、NaN/Inf 和 cell 奇异时的错误保留与重试边界。

### G6：核心工程化补缺

- [ ] calculator 命令/环境与核心算法彻底隔离。
- [x] 新增可选 image-level `ThreadedCalculatorExecutor`；主控制器统一推进
  VCNEB，Slurm 模板通过 `srun --exclusive` 并发启动独立 image worker；真实
  BTO smoke `27675909` 已验证 4×32 MPI job steps 与正常归档。
- [x] 已为 HfO₂ 生产线补齐同一控制器/worker 资源模板
  `cluster/hf_hfo2_vcneb_parallel.slurm`（默认 4×32 MPI，128 task）；仅在
  native 端点通过阈值后启用，不提前提交。
- [x] 并行 executor 提供显式 `image-retries`，并通过 fail-once 回归验证失败
  image 的局部重算；原子 snapshot、validate-only 和 calculator-free audit 已有
  独立覆盖。
- [~] 已完成跨作业从完整 chain snapshot 的精确恢复（`--resume-step`/`RESUME_STEP`）和每批 image-worker JSONL manifest；持久化 image 级结果缓存仍是后续优化项，不影响当前可恢复性验收。
- [x] 理论文档、BTO validation protocol、结果解析脚本已完成；发布前干净环境安装烟测通过（wheel 构建、toy VCNEB 与 CLI；当前发行名为 `varneb`，兼容保留 `vcneb` 入口）。

### G7：HfO2 回归与发布候选

- [x] BTO 外部 cell-filter 端点对应的普通 VCNEB、反向检查、image 收敛和电子精度对照已完成；不因方法迁移重复计算。HfO₂ 当前由 ABACUS 原生 `cell-relax`/BFGS 作业 `27676868`/`27676869`（T/PO）负责端点，尚未通过端点阈值。
- [ ] 对 HfO2 重新做端点审计、普通路径、image 收敛和 CI；旧 60 Ry 结果只作为历史诊断。
- [ ] 汇总 BTO、HfO2、解析模型、ASE 对照和 calculator smoke，形成论文表格与复现包。

## 当前停止条件

- 端点未达到阈值：停止 NEB，先修端点。
- 普通路径出现折返或最大广义力反弹：停止 CI，先修路径/优化器。
- image 数变化导致不同 saddle：不挑选结果，保留分支并报告多路径。
- 节点/资源未检查：不提交长任务。
