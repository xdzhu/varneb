# 合肥 VASP 执行策略（2026-09-18）

用户指定：后续 VASP 测试、验证与生产路径在 `ssh hf` 可达的合肥超算开展，使用
Slurm `hfacnormal01`。不再默认使用 235/cu17，提交前检查本人队列和可用资源，
不在登录节点执行 DFT，不抢占已有作业。

## 环境与数据集

```bash
source /public/home/iai806/Software/VASP/env.sh 6.3.2
```

VASP 可执行文件：`$VASP_ROOT/bin/vasp_std`。PBE PAW 数据集：
`$VASP_PSEUDO_ROOT/PBE/数据集/POTCAR`；必须区分 `Ga` 和 `Ga_d` 等价电子选择，
按原子种类顺序组合并核对原案例哈希。POTCAR 不进入公开 Git、代码包或论文附件。

当前 GaN 的 `Ga_d + N` 合并 SHA256 为
`f94781ce6cf9b9454c093353320c5e80a2a0aceea6fb8353b33b01ac37b95168`，与 235 案例一致。

## Slurm MPI 启动

本环境为 Intel MPI，使用 Slurm PMI-2，不照搬其他 MPI 环境的 PMIx 参数：

```bash
export OMP_NUM_THREADS=1
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm/lib/libpmi2.so
srun --exclusive --exact --nodes=1 --ntasks=32 --ntasks-per-node=32 \
  --cpus-per-task=1 --mpi=pmi2 "$VASP_ROOT/bin/vasp_std"
```

`--exact` 限定每个 step 的 CPU 数量，`--exclusive` 隔离并发 step，
`--nodes=1` 避免单个 image 跨节点摊开。Intel MPI 的 PMI-2 配置依据
[Intel 官方文档](https://www.intel.com/content/www/us/en/docs/mpi-library/developer-guide-linux/2021-10/job-schedulers-support.html)
与 [Slurm MPI 指南](https://slurm.schedmd.com/mpi_guide.html)。实际是否可运行仍以本平台验证作业为准。

## GaN 迁移管线

- `cluster/hf_gan_vasp_migration_validate.slurm`：1 节点、32 核，验证两个历史失败
  输入的完整 SCF、实平台之间静态能量/最大原子力的一致性、完整断点链的
  29 帧初始化，再生成两个固定端点的新平台静态结果。不进行端点弛豫。
- `cluster/hf_gan_vcneb_vasp_distributed.slurm`：3 节点，每节点 97 核，共 291 核。
  9 个并行 worker，每 worker 32 MPI，DFT 合计 288 ranks；每节点保留 1 核
  余量，包含 batch manager 开销。27 个内部像分三批评估；两个端点使用静态缓存。
- 生产作业必须依赖验证作业 `afterok`，不能仅因为本机几何检查通过就执行 DFT。
- 保持 29 总像、45.7 GPa、Ga2N2、600 eV、Gamma 8×8×6、EDIFF=1e-7、
  ISYM=-1、SYMPREC=1e-4、FIRE、fmax=0.10 eV/Å、普通无 CI。
- 从完整链恢复，不重新插帧。旧平台部分迭代的单点缓存保留在原处，但不直接
  导入新平台生产缓存。新缓存 namespace 绑定合肥 VASP 可执行文件哈希。
- `scripts/validate_vasp_migration_gate.py` 只比较静态能量和最大原子力，明确不是
  全力/应力数组一致性验证，更不是整条路径的物理收敛认证。

本轮迁移的运行目录为
`/public/home/iai806/abacus/agent-runs/20260918-varneb-vasp-hf`；后续案例需显式
设置各自输入与目录，不覆盖这条迁移管线的证据。
