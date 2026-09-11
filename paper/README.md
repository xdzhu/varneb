# VCNEB CPC 论文工作区

论文将参考本机模板包：

`D:\Work\Zstar\zstar-article\submission_packages\ZStar_CPC_pdflatex`

该模板使用 `elsarticle`、pdfLaTeX 和 Computer Physics Communications 的
`Program Summary` 结构。VCNEB 论文暂定沿用以下组织：

1. `Introduction`：晶体相变能垒、固定 cell NEB 的局限、现有 VCNEB 工具和纯 Python calculator-agnostic 需求。
2. `Theory`：扩展构型空间、cell deformation、应力/virial 到 cell force、切线、弹簧、CI、模式与方向约束。
3. `Software`：核心数据模型、ASE-compatible calculator contract、VASP/ABACUS adapter、恢复、并行 image 执行和结果记录。
4. `Verification`：有限差分、解析多井势、固定 cell ASE 对照、参数收敛、旧实现和公开方法的差异。
5. `Examples and Benchmarks`：HfO2 T->PO、钙钛矿相变、模式引导与严格约束、能垒和路径结构。
6. `Conclusions`：可复现性、适用范围、限制和后续扩展。

## 写作规则

- 理论稿件只使用已经由测试或数据支撑的公式与结论；未完成的内容使用 limitation/future work 标记。
- 每一张图和表都必须能回溯到一个输入 manifest、一个代码 commit 和一个集群运行目录。
- 主文稿、图表生成脚本、BibTeX 和最终复现包在本目录或其明确的子目录中管理。
- 长时间 DFT 只在 `cu17`、`cu22`--`cu26` 执行；论文渲染可在本地完成。

## 当前稿件拆分

- `vcneb_CPC.tex`：在 P4 验证结果稳定后建立的主稿件。
- `figures/`：路径、cell 演化、收敛和 calculator 对比图。
- `data/`：论文使用的汇总 CSV/JSON，不放入大体积 DFT restart 文件。
- `reproduce/`：从 manifest 复现表格和图的脚本。

