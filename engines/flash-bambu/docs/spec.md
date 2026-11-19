# `bambu_next` 工程规格与验收契约 (`SPEC`)

> **基线仓库**：[`long_bench/bambu/`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/)（原版 3.13.1 只读源码） & `envs/bambu` 安装包（`bambu 3.12.1`）。  
> **重构目标目录**：[`/gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu_next/`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu_next/)。  
> **参照规范**：[`AGENTS.md`](file:///gpfs/hpcusers/yanhong.hong/.codex/AGENTS.md)、[`long_bench/docs/bambu_full_pipeline_spec.md`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/docs/bambu_full_pipeline_spec.md)、[`isoquant_refactored/docs/spec.md`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/docs/spec.md)。

---

## 1. 核心目标与变更边界

### 1.1 核心目标
1. **端到端全阶段并行 + 低内存深度重写**：打通并重写 Bambu 全流程五大宏观阶段——`S0` 注释准备与索引 (`prepareAnnotations` / `.bbidx`)、`S1` BAM 解码 / 剪接校正 / ReadClass 构建与评分 (`prepareDataFromBam` / `isore.constructJunctionTables` / `isore.constructReadClasses` / `scoreReadClasses`)、`S2` 跨样本候选合并与注释扩展 (`isore.combineTranscriptCandidates` / `isore.extendAnnotations`)、`S3` 距离匹配 / 等价类构造 / EM 定量 (`bambu.quantify` / `emWithL1` / `combineCountSes`)、`S4` 基因汇总与流式导出 (`transcriptToGeneExpression` / `writeBambuOutput` / `writeToGTF`)。
2. **系统性攻克原有 `G8–G17` 性能与内存瓶颈**：通过 **C++17 + `Rhtslib` + 扁平只读 CSR/NCList 索引 + 线程安全 `std::thread`/OpenMP** 脱离 R 分代 GC 与 `MulticoreParam` `fork` CoW 页分裂泥潭，在 16 BAM 规模下大幅压低峰值树 PSS（设计目标较原生 `250+ GiB` 降至数十 GiB 级）与端到端 wall-clock 时间。
3. **公共 API 与产物向后 100% 等价**：公开 R 包函数签名与返回值（`bambu()`, `prepareAnnotations()`, `readFromGTF()`, `writeToGTF()`, `writeBambuOutput()`, `transcriptToGeneExpression()`, `trainBambu()`, `compareTranscripts()`, `plotBambu()`）及输出文件（`extended_annotations.gtf`, `counts_transcript.txt`, `CPM_transcript.txt`, `fullLength_transcript.txt`, `uniqueCounts_transcript.txt`, `counts_gene.txt`, `CPM_gene.txt` 等）保持与原版完全兼容。

### 1.2 只读保护与禁止行为（Out of Scope & Guardrails）
- **严格只读目录**：禁止修改 `long_bench/bambu/`、`bambu_stage1_parallel/`、`bambu_refactor/`、`bambu_production/`、`bambu_lowcopy/`、`bambu_s3core/`、`long_bench/docs/`、`long_bench/results/`、`long_bench/config/`、`envs/*` 及 `/gpfs/.../backup/longpass_benchmark250209/` 下任何文件。
- **不修改既有 Conda 环境**：`bambu_next` R 包通过 `R CMD INSTALL --no-multiarch -l bambu_next/.r_lib` 安装到独立局部库目录；测试及基准运行通过 `R_LIBS_USER` 优先加载 `bambu_next`，绝对对照基准在干净独立进程仅加载系统 `envs/bambu` 安装版运行。
- **禁止写系统 `/tmp`**：全部 Slurm scratch、`.bai` 软链缓存、`TMPDIR`/`TMP`/`TEMP`、中间 checkpoint 与输出全部放在 `/gpfs/group/data_sci_group/weiwu_group/hongyanhong/bambu_next_work/`（用后清理大文件，只留结构化 TSV/JSON/日志/图表）。
- **严禁观测污染基准时间**：精确计时与 `monitor_process_tree.sh` PSS 采样时必须**关闭调试 checkpoint RDS 转储**（与正确性回归完全解耦）；任何未实测的设计预估值不得在最终验收中冒充实测结论。

---

## 2. 集群环境与编译器工具链审计（已实测确认）

| 组件 | 实测版本 / 路径 | 说明与约束 |
| :--- | :--- | :--- |
| **主 R 运行与编译环境** | `conda activate /gpfs/group/data_sci_group/weiwu_group/hongyanhong/miniconda3/miniconda3/envs/bambu` | `R 4.5.3`、`x86_64-conda-linux-gnu-cc` / `c++ -std=gnu++17`、`SHLIB_OPENMP_CXXFLAGS = -fopenmp`、`BLAS_LIBS = -lblas`、`LAPACK_LIBS = -llapack` |
| **Bioconductor / C 接口包** | `bambu 3.12.1`, `BiocParallel 1.44.0`, `Rhtslib 3.6.0`, `GenomicAlignments 1.46.0`, `GenomicRanges 1.62.1`, `GenomicFeatures 1.62.0`, `IRanges 2.44.0`, `S4Vectors 0.48.0`, `Rsamtools 2.26.0`, `SummarizedExperiment 1.40.0`, `BSgenome 1.78.0`, `Biostrings 2.78.0`, `data.table 1.17.8`, `Rcpp 1.1.1`, `RcppArmadillo 15.2.6-1`, `xgboost 3.2.0.1` | 已有全部 C headers：`Rhtslib/usrlib/libhts.a`、`IRanges_interface.h`、`S4Vectors_interface.h`、`Biostrings_interface.h`、`datatableAPI.h`、`RcppArmadillo.h`；**零需新增包** |
| **Python 汇总与绘图环境** | `conda activate /gpfs/group/data_sci_group/weiwu_group/hongyanhong/miniconda3/miniconda3/envs/basicpy` | Python `3.11.15`，用于 Slurm Manifest 生成、严格数值对比、进程树 PSS 聚合与出版级 PDF/图表编译 |
| **Slurm 分区资源** | `chj.cn.1,chj.cn.2,chj.cn.3,chj.cn.4,chj.cp.1`（96–192 CPU / 节点） | 单任务 Array 无上限并发提交（不设 `%K`），按样本规模分配 CPU/内存，严格 `A-B-B-A` 交替配对 |

> **Conda 激活模板**（非交互 Shell / Slurm 必须完整显式执行）：  
> `source /gpfs/group/data_sci_group/weiwu_group/hongyanhong/miniconda3/miniconda3/etc/profile.d/conda.sh && conda activate /gpfs/group/data_sci_group/weiwu_group/hongyanhong/miniconda3/miniconda3/envs/bambu`

---

## 3. 评测数据集与分层测试矩阵

为同时保证**与既有 `G0–G17` 历史基准的严格可比性**以及**对全长读长基准库的全面覆盖**，评测数据分为三层：

### Tier 0：包内置微数据集 & 合成边界套件（秒级单元测试 / 严格 `identical` 回归）
- `bambu/inst/extdata/` 全部夹具：SGNex A549 chr9 (`annotationGranges_txdbGrch38_91_chr9_1_1000000.rds`, `readGrgList_*`, `defaultModels.rds`, `.ubj` 模型与全部期望輸出 `.rds`)；
- `bambu/tests/testthat/` 全部 11 项官方测试集；
- 合成边界用例：unspliced 单外显子 read、0-intron 样本、高重复 readCount 加权中位数（含 `10^7` 计数防展开溢出）、负链微外显子、多基因连通重命名、空矩阵/全兼容/全不兼容分组、`NA` / `factor` levels / tie-breaking 边界。

### Tier 1（`T_core`）：冻结 18 BAM 核心对照集（同源直接对比 `native` / `G8` / `G17`）
来源：[`results/bambu_s3core/g15_20261113/config/samples.tsv`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/results/bambu_s3core/g15_20261113/config/samples.tsv)（16 Human + 2 Mouse BAM）：
1. **单样本代表与扩展**：
   - `ENCFF400BQQ` (`cDNApacbio`, 0.09 GiB, 240,947 alns) — 线程扩展 (`ncore = 1, 4, 8, 16, 32`) 与预建 `.bbidx` 热跑测试；
   - `ENCFF436GKZ` (`cDNApacbio`, 0.48 GiB, 1.29M alns) — `human_single_d1q1` 冻结代表（原生基线 `370.35s`, `4.95 GiB PSS`）；
   - `ENCFF846RID` / `ENCFF207KVO` (`directRNA`, 1.62 GiB / 1.58 GiB, 1.25M / 1.39M alns)；
   - `ENCFF089AMF` (`cappedONT`, 4.78 GiB, 16.42M alns) & `ENCFF887KHY` (`R2C2ONT`, 4.66 GiB, 18.28M alns) — 超大单样本压力测试；
   - `ENCFF598YQO` (`R2C2ONT`, mouse, 3.26 GiB, 10.11M alns) — `mouse_single_d1q1` 冻结代表（原生基线 `943.45s`, `17.43 GiB PSS`）。
2. **多样本全组发现 + 定量**（`discovery=TRUE, quant=TRUE` + 可选参数开关矩阵）：
   - `human_pair_d1q1` (`ENCFF684YOO` + `ENCFF436GKZ`：原生基线 `431.54s`, `12.14 GiB PSS`)；
   - `mouse_all_d1q1` (`ENCFF598YQO` + `ENCFF055REA`：原生基线 `1320.17s`, `37.80 GiB PSS`)；
   - **`human_all_d1q1`**（Human 全部 16 个 BAM 联跑，`>10` 样本临时缓存分支，8 workers：原生基线 `4763.35s`, `250.48 GiB PSS`；`g8_default` `4609.85s`, `254.29 GiB PSS`）。
3. **可选开关覆盖**：`stranded = TRUE`、`trackReads = FALSE`、`returnDistTable = TRUE`、`fusionMode = TRUE`、`lowMemory = TRUE`、`NDR = NULL`（自动 NDR 推荐）、`yieldSize` 分批流式、GTF 原生输入 `prepareAnnotations(gtf)` vs `.prepared_annotation.qs`、`trainBambu()`。

### Tier 2（`G_ext`）：`longpass_benchmark250209/` 全协议全物种扩展验证集
来源：[`longpass_benchmark250209/allalignments_abs.txt`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/backup/longpass_benchmark250209/allalignments_abs.txt)（89 条核心比对 + `SIRV_simulation` + `RBMX_KD` + 海牛 `trichechus_manatus`）：
- 全 5 类人源长读长协议组（`cDNA_pacbio`, `capped_pacbio`, `cDNA_ONT`, `capped_ONT`, `R2C2_ONT`, `directRNA`）；
- `RBMX_KD` (`con_rep1`, `shRBMX_rep2` 7.8 GiB) + **`gencode.v46.annotation.gtf`（1.6 GiB 最新全量 GTF 从头解析 + Discovery + Quant）**；
- `SIRV_simulation/alignments/simulated1_aligned_reads.sorted.bam`（8.77M reads 真值精确率/召回率评估）；
- 非模式生物海牛 `Trichechus_manatus`（含 `de novo` `annotations = NULL` 及缺失 `.bai` 时作业目录只读软链建索引测试）。

---

## 4. 一致性、时间与内存验收口径

1. **科学等价标准（分两级严控）**：
   - **L1（结构与整型对象逐位一致）**：BAM 过滤计数、read class 外显子/内含子坐标、strand、junction motif 与校正映射、read-to-gene 分配、novel transcript 拓扑、`txClassDescription`、排序后 `rowRanges`、`GENEID` 连通集合、`TXNAME`/`GENEID` 命名、`NDR`（浮点误差 `< 1e-12`），所有离散字段要求 **R `identical()`**。
   - **L2（EM 定量数值等价）**：在相同输入等价类矩阵与随机状态下，transcript/gene `counts`、`CPM`、`fullLengthCounts`、`uniqueCounts` 最大绝对误差 `< 1e-6` 且 Pearson `r >= 0.999999`；导出 GTF/TSV 记录集合与数值完全一致。
   - **确定性固定**：默认参数遵循 `seed = 20260821L`、`NDR = 1`（及 `NDR = NULL` 独立验证）、`opt.discovery = list(min.readFractionByGene = 0.001, min.exonDistance = 25)`。
2. **时间与内存采集标准**：
   - 复用与优化独立 [`monitor_process_tree.sh`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/scripts/benchmark/monitor_process_tree.sh) 架构（每秒汇总整个根 PID 进程树的 `/proc/<pid>/smaps_rollup` `Pss`, `Rss`, `Pss_Anon`, `Pss_File` + 实测采样间隔）。
   - 分阶段窗口：`S0` (Annotation Load/Prepare)、`S1` (Per-sample Read Processing & Scoring)、`S2` (Cross-sample Combine & Annotation Extension)、`S3` (Distance Matching, EquiRC & EM Quant)、`S4` (Gene Rollup & GTF/TSV Export) 与总墙钟 `whole_R`。
   - 配对规则：同类实验同节点/Slurm 序列采用 `A-B-B-A` 交替 3 轮（取中位数及 `[min, max]`），监测全关闭调试转储以隔离 I/O 偏差。
