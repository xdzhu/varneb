# CdSe rock-salt → wurtzite：独立双映射案例

研究设置、结构来源和比较边界见
[复现计划](../../docs/CDSE_SHEPPARD_2012_REPLICATION_PLAN.md)。
本目录不含获许可的 POTCAR 或 VASP 二进制，也不把种子标为优化完成的端点。

入口：`scripts/prepare_cdse_sheppard_suite.py`、
`scripts/submit_cdse_sheppard_suite.py` 和 `cluster/hf_cdse_sheppard_*.slurm`。
实际源码采用独立不可变远端快照，不能覆盖 GaN 正在运行的目录。

RS/WZ 为相同 8 原子 Cd4Se4；两种显式映射各使用 15 中间像/17 总像。
普通 VCNEB fmax=0.10、无 CI；两端固定静态结果只缓存一次。
每条路径的计算器为 VASP/PBE/PAW，455 eV、MP 10³、0 GPa。
原文 PW91 的 2.4 meV/atom 小峰只是带设置差异的比较对象，
提交成功、初始化通过或达到 0.10 不能证明已准确复现这个小势垒。

六方/岩盐端点 BFGS 完成后才构造生产输入；每条路径分别完成全链晶格/输入检查
与两端完整静态 SCF，然后进入 6×32 MPI 的分布式 worker 计算。
后续结果须独立审核机制、相身份、峰和能量单位，不混用经验势结果。
