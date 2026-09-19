# GaN：Qian 2013 文献路径对照

此案例使用 VASP / PBE / Ga_d+N PAW，在文献压力下比较 B4→B1 四方/六方机制，以及 B3→B1。完整方案见 `docs/GAN_QIAN_2013_REPLICATION_PLAN.md`，启动记录见 `validation/gan_qian_suite/launch_20260918.md`。

`seeds/b3/POSCAR`、`seeds/b1/POSCAR` 均为 8 原子 Ga4N4 常规胞，空间群分别 #216/#225。它们是 45.0 GPa 端点优化的**种子**，不是已收敛端点。种子体积从已有 45.7 GPa B4/B1 结果构造，不能直接拿它们计算最终能垒。

已收敛的45.0 GPa端点另存于 `endpoints_45p0/{b3,b1}/CONTCAR`，各目录保留原始优化摘要。两端均在BFGS第2步收敛，通过原子数、广义力、内部压力、几何和多容差晶相审计；详见 `validation/gan_qian_suite/endpoints_45gpa_audit.md`。不要将种子与这些收敛端点混淆。

已完成的 B4→B1 四方机制使用 4 原子 Ga2N2 共格胞，焓垒除以 2；此 B3 常规胞案例除以 4。不同原子数案例不得直接比较每胞总能量。

需要合法的 VASP 可执行程序及 Ga_d、N POTCAR；赝势不随仓库分发。`scripts/prepare_gan_qian_suite.py` 从已成功案例拷贝本地合法输入契约，构造端点种子、六方初始链及优化后的 B3 对角初始链；不执行 DFT。

计算入口：

- `cluster/hf_gan_qian_endpoint.slurm`：BFGS + FrechetCellFilter，45.0 GPa，端点广义力 0.02 eV/Å。
- `cluster/hf_gan_qian_validate.slurm`：预检、29 像初始化及两个端点静态 SCF。
- `cluster/hf_gan_qian_vcneb.slurm`：普通 VCNEB，27 中间像，0.10 eV/Å，9 个 32-MPI 并行 worker；端点固定缓存，无 CI。
- `scripts/submit_gan_qian_suite.py --repo <prepared-remote-repo> --run`：六作业依赖链；已有 launch.json 时拒绝重复提交。

显式链应使用 `examples/run_vcneb_vasp.py --initial-trajectory <initial.traj> --n-images 29`，禁止同时用 MIC、自动匹配或平移对齐重新缩短指定路径。六方锚点只是初始猜测，最终需验证实际机制；B3 原文三峰长链仅用于原图映射诊断，单次转变应该另外审计最短原子对应与分段机制。

文献：[Qian et al., CPC 184 (2013) 2111–2118](https://doi.org/10.1016/j.cpc.2013.04.004)。文献 GaN 使用 QE/PW91/USPP，而此例使用 VASP/PBE/PAW；相近能垒不等同逐参数复现。
