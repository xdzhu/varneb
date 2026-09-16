# BaTiO3 T→C VCNEB 验收协议

本协议固定当前主案例的计算设置和结果判据。它把“增加 image”定义为数值
收敛实验，而不是把 image 数量越大越好。

## 固定计算设置

- ABACUS Dojo-NC-FR，`ecutwfc = 100 Ry`。
- Ba/Ti/O 全部使用 `Orb-DZP-10au`：Ba `4s2p1d`、Ti `4s2p2d1f`、O
  `2s2p1d`。
- `4×4×4` k 点、`scf_thr = 1e-8`、`scf_nmax = 150`、PBE、静态 image
  calculator 提供 energy/force/stress。
- 端点先独立变胞弛豫，`fmax ≤ 0.01 eV/A`；NEB 使用相同端点，不在 image
  calculator 中再次做结构弛豫。
- 生产模板从后续作业起申请 32 MPI ranks、每 rank 1 CPU；已经运行的作业不
  中途改变资源。

## 跨计算器端点一致性闸门

ABACUS、QE 与 VASP 的每次 no-DFT preflight 都会在
`vcneb_preflight.json[endpoint_structures]` 写入有序物种、PBC、晶胞和包裹
分数坐标的 SHA256 结构指纹。QE/VASP 的 7-image 生产作业前，必须以
`scripts/compare_vcneb_endpoint_records.py` 将其初末端点分别与接受的 ABACUS
preflight 比较；任一端点不匹配即阻止生产提交。该规则有意保留原子顺序，避免
未经审计的重排被误报为同一条 NEB 路径。

QE/VASP 的 160-rank 模板还要求将这个通过的 JSON 路径作为
`ENDPOINT_IDENTITY_GATE` 传入；模板会在启动任一 image worker 前再次验证
`initial` 与 `final` 两项均为通过状态。

QE 还要求 `QE_PP_MANIFEST`。该 manifest 的 `approval_status` 必须为
`approved`，并为 Ba/Ti/O 逐项固定 UPF basename 与 MD5；驱动会同时检查
文件存在、元素 metadata、PBE 标记和 MD5。仓库中的 SSSP 候选登记不是
approved manifest，不能直接启动预检或生产作业。

VASP preflight 会记录初始目录 `INCAR`、`KPOINTS` 与许可 `POTCAR` 的完整
SHA256；同一份经过审核的输入必须复制到每一个静态 image 目录。

## Image 数量策略

第一轮只比较 5、7、9 个总 image（含两个端点），保持端点、calculator、
`log_strain` 插值、mapping、`cell_scale`、弹簧常数和 optimizer 完全一致。

1. 先跑 7-image 普通 NEB，确认路径没有折返、碰撞或异常 cell。
2. 以 5→7→9 顺序跑普通 NEB；普通路径稳定前不启动 CI。
3. 同时记录正向能垒、反应焓、最高 image 的 reaction coordinate、体积/晶格
   变化、关键 image 的原子力/应力和最大广义力。
4. 默认收敛门槛为：相邻 image 的能垒差不超过 `0.02 eV`，最高 image 的
   reaction-coordinate 位移不超过相邻 segment 长度的 `0.25`，且最高 image
   的结构/体积变化不出现新的局部极值。若 7→9 仍不满足，只在高曲率区增加
   11 或 13 image，不盲目把整条路径加密。

image 数偏少的证据是能量峰位置在 5→7→9 间明显移动、segment cosine 变差、
或能垒随加密单调漂移；image 数过多的证据是能垒和 saddle 已稳定但 wall time
近似按内部 image 数线性增加。最终报告保留全部分支，不挑选较低能垒的一支。

FIRE 的单步最大广义力回弹不视为失败。运行控制采用至少 3--5 个完整
optimizer step 的观察窗口，并保留窗口内最佳完整 chain snapshot；只有持续
回弹或平台且最佳值不再改善时才转入静态审计或更换优化器。

## CI 与恢复闸门

- 只有普通 NEB 满足广义力目标（当前预收敛目标 `0.02 eV/A`，最终可按路径
  噪声加严）且几何诊断通过，才从同一普通轨迹启动 staged CI。
- CI 报告必须包含最高 image、负切线曲率、垂直力、完整物理原子力/应力和
  `has_interior_barrier`；投影 NEB residual 不能替代完整物理力。
- 人为终止后只从最近完整 `chain_step_####.traj` 或完整 trajectory chain 恢复；
  恢复日志必须证明没有重新插值、没有覆盖原始 `initial-vcneb.traj`，且 image
  目录仍一一对应。

电子精度对照必须同时重评估端点和路径。当前 BTO 的 6×6×6、SCF 1e-9
对照与 4×4×4、SCF 1e-8 的拓扑结论一致，但绝对反应能相差约 15 meV，
因此最终能量报告不能只引用低精度端点。

## 结果归档

每个分支目录至少保存：端点结构、初始路径、每 image 输入/输出、optimizer
日志、完整 trajectory、chain snapshots、`vcneb_summary.json`、
`vcneb_preflight.json`、Slurm 作业号/节点/task 数、代码版本和
`scripts/audit_vcneb_result.py` 的审计结果。
