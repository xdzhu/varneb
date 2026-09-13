# VCNEB CPC 论文工作区

论文将参考本机模板包：

`D:\Work\Zstar\zstar-article\submission_packages\ZStar_CPC_pdflatex`

该模板使用 `elsarticle`、pdfLaTeX 和 Computer Physics Communications 的
`Program Summary` 结构。VCNEB 论文暂定沿用以下组织：

1. `Introduction`：晶体相变能垒、固定 cell NEB 的局限、现有 VCNEB 工具和纯 Python calculator-agnostic 需求。
2. `Theory`：扩展构型空间、cell deformation、应力/virial 到 cell force、切线、弹簧、CI、模式与方向约束。
3. `Software`：核心数据模型、ASE-compatible calculator contract、VASP/ABACUS adapter、恢复、并行 image 执行和结果记录。
4. `Examples`：HfO2 T->PO、钙钛矿相变、模式引导与严格约束、能垒和路径结构；有限差分、解析多井势、固定 cell ASE 对照、参数收敛和文献对比作为本节的小节、表格或图展示。
5. `Conclusions/Availability`：可复现性、适用范围、限制、后续扩展、代码和数据获取方式。

## 写作规则

- 理论稿件只使用已经由测试或数据支撑的公式与结论；未完成的内容使用 limitation/future work 标记。
- 每一张图和表都必须能回溯到一个输入 manifest、一个代码 commit 和一个集群运行目录。
- 主文稿、图表生成脚本、BibTeX 和最终复现包在本目录或其明确的子目录中管理。
- 长时间 DFT 只在 `cu17`、`cu22`--`cu26` 执行；论文渲染可在本地完成。

## 当前稿件拆分

- `vcneb_CPC.tex`：按 `Introduction -> Theory -> Software -> Examples -> Conclusions/Availability` 组织的当前主稿件草稿；不单列 Benchmarks。
- `vcneb.bib`：当前稿件引用，正式投稿前需继续补齐软件和材料案例的准确书目信息。
- `zstar-elsarticle-num.bst`：从本地 CPC 模板复制的参考文献样式。
- `figures/`：路径、cell 演化、收敛和 calculator 对比图。
- `data/`：论文使用的汇总 CSV/JSON，不放入大体积 DFT restart 文件。
- `reproduce/`：从 manifest 复现表格和图的脚本。
- `claim_evidence.md`：将主稿件 claim 映射到代码、回归测试、作业和归档结果，并明确当前证据边界。
- 根目录 `scripts/plot_vcneb_metrics.py`：只读取已归档的逐 image CSV，生成焓垒、晶格长度、体积和广义力四联图；不重新调用计算器。

## 当前编译

在仓库根目录执行：

```bash
python C:/Users/zhu/.codex/plugins/cache/openai-bundled/latex/0.2.6/scripts/compile_latex.py paper/vcneb_CPC.tex --output-directory paper/build --json
```

当前草稿已用本机 TeX Live 2023 编译为 `paper/build/vcneb_CPC.pdf`。PDF 和中间文件属于构建产物，不作为论文源文件提交；所有数值结果仍需以 `outputs/` 下的 manifest 和原始集群目录为准。
