# Four-part narrative and the completed atomic-treatment evidence

2026-10-15. The author explicitly approved four main Results sections with
detailed diagnostics retained in the paper's appendices and repository.
This is a narrative/material-evidence milestone, not a completed JCTC claim.

## English scientific statement

For the five registered preserving-channel local coordinates, most energy
lowering was already obtained by atomic rearrangement at each frozen-Q seed
cell. The original quadratic's shape mismatch was also largely present
before the additional joint-release step. This localizes the measured
treatment difference without identifying a unique microscopic energy
partition, a whole-path mechanism or a stable conditional branch.

## 中文阅读顺序

正文现在只沿四个问题展开：路径竞争如何随应变变化；机械约束如何改变路径
关键结构的能量角色；模式描述覆盖什么、原子与联合放松又改变什么；模型在
已有数据上表现如何、独立预测还缺什么。工作积分、交叉应变、曲率与表示
一致性等细致检查移入附录E，原数值与限制条件保留，不把每次核验写成新发现。

已完成的10点/36SCF原子处理证据写入方法2.7、结果3.3、讨论和附录D.8。
新TableA10的4行同时保留两个中心、两步长及冻结/原子/联合三种处理。
能量统一标明meV/12原子胞，不能混成meV/HfO2。细致数表不是删除的补充材料，
仍属于同一正文文件中的附录。预测部分明确说明：独立验证点尚未计算，
新方法优势与横向稳定性没有证明。没有增加材料或计算任务。

## Evidence and checks

The original complete E194 native proof and comparison are byte-pinned by
check_tables.py. It checks the four TableA10rows, actual complete same-Q
coverage and normal original completion; it does not replace native EFS
validation. Original142table lines remain byte-identical, with their
multiplicities checked independently of the approved new order. Six new
table lines give148in total. All nine figure paths/order/bytes, thirteen
external citation links and the TeX wrapper remain unchanged. Main tables
are1--3; relocated diagnostic tables areA11--A16, with references updated.
Two existing treatment terms were added to the terminology ledger.

First related tests:127pass/11fail. The failures were a new raw-line versus
whitespace-normalized prose assertion, legacy table order/count assumptions,
and eight negative fixtures selecting the first *scalar-matching* table row.
After reordering, the new A10table legitimately precedes the old F/J table
and shares some scalar values. The old negative fixture was duplicating
the wrong table. It now selects the entire original source-derived F/J row;
changed/missing/duplicate source values remain rejected for both tables.
Order checks now preserve all old row contents and multiplicities; the
old9-section expectation becomes the approved4. No physical test tolerance,
input, gradient gate or DFT result was relaxed. The actual first failed XML
is retained as first_related_tests.xml.

Initial compilation rendered36pages. Visual inspection caught two long
appendix captions split across pages. Detailed notes were moved into
preceding prose, leaving short captions attached to the same unchanged
tables. The first PDF/report remain on C:/in render_metadata.json;
render_revision_metadata.json describes the revised rendering. Main paper
and old figure values are not inferred from a compiler success alone.

Reproduce the source-to-table check from the repository root:

```text
python benchmarks/hfo2_channels/20261008/manuscript_atomic_treatment_E195_20261015/check_tables.py --output NEW_TABLE_CHECK.json
```

No new DFT, scheduler operation, physical/private holdout access, parameter,
domain or allocation-budget change, climbing-image run, installation,
release, PyPI publication, automation or journal submission is introduced.
The five protected user edits remain separate. Full goal remains active.

Final scope:138local and138independent related checks across11modules,
including12new checks; zero final failures/errors/skips, not a full suite.
The independent source-to-table CLI reproduces the final report bytewise.
Tested final source treee48fddab65207271559dfcb3187e0da856409e8c/archive
dd93dbc9feac76b1b50fc5feeae89fb5343b1bb6518ecb038c5b399bbd28e23a.
Final manuscript SHA b9e29fca28ed659b87a05e6f26306f993786fb023114462e4916c7e6b540ac91;
36-page PDF SHA6f48c8727e427ee92bff6d4b68754a06aa2e59c4f024bce55b82dde5db3f4b9d.
All36initial pages were visually reviewed as contacts, with detailed checks
of the main treatment/control pages and appendix tables. In the revision,
pages1--31are pixel-identical; changed32--36were all re-inspected. Captions
A13/A15are now on their tables' pages33/35. Final compile has no overfull,
underfull or missing-character warnings. The working-draft/incomplete-tests
notice remains; a clean PDF is not a scientific-completion certificate.
