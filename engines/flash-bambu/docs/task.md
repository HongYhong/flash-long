# `bambu_next` 可追踪任务台账 (`TASKS`)

> 状态定义：`TODO`（待开始）｜`IN_PROGRESS`（进行中）｜`DONE`（已完成并验证）｜`BLOCKED`（受阻塞）  
> 优先级：`P0`（关键路径阻塞项）｜`P1`（核心功能与主要优化）｜`P2`（增强/兼容与交付收口）

---

## Phase Pre：调研、架构设计与执行文档定稿

| ID | 阶段/节点 | 任务名称 | 状态 | 优先级 | 依赖 | 交付物与验证标准 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`T00`** | Pre | 复盘 Codex Bambu 优化历史 (`G0–G17`)、审计 15 个 Bioconductor 依赖包官方手册与 C 头文件，新建 [`bambu_next/`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu_next/) 并完成全套执行文档 | `DONE` | `P0` | — | `README.md` 与 `docs/{bioconductor_manuals_audit,spec,optimization_schemes,roadmap,plan,task}.md` 全部落盘并交用户确认执行 |

---

## Phase 0：工程骨架、基准观测底座与官方基线体检 (`M0_PRE`)

| ID | 阶段/节点 | 任务名称 | 状态 | 优先级 | 依赖 | 交付物与验证标准 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`T01`** | Phase 0 (`Node 0`) | 搭建 `bambu_next` R 包骨架（`DESCRIPTION`, `NAMESPACE`, `src/Makevars`, `init.cpp`），验证 `envs/bambu` 下 `Rhtslib` + `IRanges`/`S4Vectors` C stub + `RcppArmadillo` + OpenMP 局部库隔离安装 | `DONE` | `P0` | `T00` | `scripts/install_local.sh` 在 `bambu_next/.r_lib` 编译安装成功，零修改 `envs/bambu` 全局目录 |
| **`T02`** | Phase 0 | 部署 `scripts/benchmark/monitor_process_tree.sh`、`compare_bambu_outputs.R`、`run_paired_ab_ba.py`（全路径隔离 `/tmp`，写入 `/gpfs/.../bambu_next_work/`） | `DONE` | `P0` | `T01` | 合成进程树压力测试通过，准确统计 `/proc/<pid>/smaps_rollup` Tree PSS / RSS 与采样间隔 |
| **`T03`** | Phase 0 (`G0`) | 在 `envs/bambu` 运行官方 `bambu/tests/testthat/` 与 `ENCFF400BQQ` (`92 MiB`) 基线测绘，冻结对照组校验基准 | `DONE` | `P0` | `T02` | 原生 56 项非 ggplot 官方测试 56/56 PASS（跳过 `test_plotBambu.R` 因环境未装 `ggplot2`）；`ENCFF400BQQ` A-B-B-A ×3 (`nc=16`, `5.90x`) 及 7 参数开关全量对齐落盘 |

---

## Phase 1：运行护栏与零拷贝 `Flat-CSR` 注释系统 (`Node 0 ~ Node 1`)

| ID | 阶段/节点 | 任务名称 | 状态 | 优先级 | 依赖 | 交付物与验证标准 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`T04`** | Node 0 | 实现全局线程/缓存护栏 (`setDTthreads(1L)`, XGBoost `nthread=1L`, 主进程安全初始化隔离目录 `BiocFileCache`, C++ 工作线程池 `thread_pool.hpp`) | `DONE` | `P0` | `T01` | 并发线程护栏与 BiocFileCache 多样本隔离单测全通过（56/56 PASS） |
| **`T05`** | Node 1 (`S0`) | 实现 `Flat-CSR` (`flat_csr_index.hpp`) 与单遍 C++ 流式 `GTF` 解析器 (`gtf_stream_parser.cpp`)，兼容 `TxDb` 与预构建 `.prepared_annotation.qs` C 指针导入 | `DONE` | `P0` | `T04` | `human.v39` (`9.62s` vs 原生 `121.42s`)、`mouse.vM28` (`5.42s`)、`gencode.v46` (`16.72s`, `1.28 GiB RSS`) 全通过，与参考 `.qs` `identical()` |
| **`T06`** | Node 1 (`S0`) | 实现一次性参考 `min_eq_class_builder.cpp` 拓扑预计算，固化至 `.bbidx` 与 `FlatCSRAnnotation` 元数据 | `DONE` | `P0` | `T05` | `human.v39` (`245,115` tx) 与 `mouse.vM28` (`142,552` tx) `eqClassById` 逐转录本集合 0 差异 (`diffs=0`) |

---

## Phase 2：Stage 1 `Rhtslib` 并发流式阅读与 ReadClass 构建评分全链路 (`Node 2 ~ Node 4`)

| ID | 阶段/节点 | 任务名称 | 状态 | 优先级 | 依赖 | 交付物与验证标准 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`T07`** | Node 2 (`S1.04–06`) | 实现基于 `Rhtslib` (`sam.h`, `.bai` 分窗并发 + 无索引退路) 的 `bam_exon_extractor.cpp`，单遍解码 `bam1_t` 生成 exon CSR 与线程本地唯一 junction hash | `DONE` | `P0` | `T05` | `extdata`, `ENCFF400BQQ` (`240,940` valid / `144,862` junc), `ENCFF846RID` (`1,247,064` valid / `220,694` junc, `5.82s` vs 原生 `71.70s`) 与原生 100% `identical()` |
| **`T08`** | Node 2 (`S1.07`) | 实现 `htslib/faidx.h` 批量 junction `spliceMotif` 提取 (`faidx_motif_lookup.cpp`)，兼顾内存 `BSgenome` 序列直读回退 | `DONE` | `P0` | `T07` | 全部测试 BAM `uniqueJunctions`、`spliceMotif` 与 `spliceStrand` 与原版 `identical()` |
| **`T09`** | Node 3 (`S1.09–11`) | 实现 C++ 剪接链投票校正 (`junction_corrector.cpp`)、无 `data.frame` 拼接批 XGBoost (`nthread=1L`) 剪接位点校正与高频区间归并 | `DONE` | `P0` | `T08` | `extdata`、`ENCFF400BQQ` (`144,862` junc, 0 diffs) `corrected_strand`/`annotated`/`mergedHighConfJunctionId` 与原版 `identical()` |
| **`T10`** | Node 3 (`S1.12`) | 实现 `read_class_builder.cpp`（64-bit `XXH3` 多键哈希构建 Spliced ReadClass + 双指针区间聚类 Unspliced ReadClass），消除 `dplyr::group_by`/`unstrsplit`/`separate` | `DONE` | `P0` | `T09` | `extdata` (0.056s vs 1.30s) 与 `ENCFF400BQQ` (`59,344` RC / `240,940` reads, `7.81s` vs `29.51s`) 全列及 `readIds` `identical()`；`ENCFF846RID` `1.25M` reads 仅 `11.2s` |
| **`T11`** | Node 4 (`S1.14`) | 实现 C++ Union-Find `assignGeneIds` (`rc_scorer_uf.cpp`)、`isReadClassCompatible`、polyA/T 计数与并行 `trainBambu()`，封装兼容 `ReadClassSe` | `DONE` | `P0` | `T06, T10` | **S1 Gate**：`T_core` 18 BAM `discovery=FALSE, quant=FALSE` 输出 `RangedSummarizedExperiment` 与原版 `identical()` |

---

## Phase 3：Stage 2 跨样本合并与增量注释扩展全链路 (`Node 5 ~ Node 6`)

| ID | 阶段/节点 | 任务名称 | 状态 | 优先级 | 依赖 | 交付物与验证标准 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`T12`** | Node 5 (`S2.01–02`) | 实现 `kway_rc_merger.cpp`（$O(S \cdot N)$ K-way 合并 + `O(U log U)` 无 `rep()` 展开精确 weightedMedian）与 `unspliced_cluster_merge.cpp` | `DONE` | `P0` | `T11` | `human_pair`、`mouse_all` 与 **`human_all_d1q1` 16 样本**合并全部通过（S2 自 `989.26s` 降至 `12.95s`，`76.4x` 加速） |
| **`T13`** | Node 6 (`S2.03–08`) | 实现增量 `addNewSpliced/UnsplicedReadClasses`、`filterTranscripts`、`recommendNDR` 与**增量 `min_eq_class_builder.cpp`** | `DONE` | `P0` | `T06, T12` | 扩展后 `annotations` `GRangesList`、`txClassDescription`、`GENEID`/`TXNAME`、`eqClassById` 与原生 `identical()` (`human_pair`/`mouse_all`/`switches_a/b` 全 PASS) |

---

## Phase 4：Stage 3 零展开匹配与融合无 `dcast` EM 定量全链路 (`Node 7 ~ Node 8`)

| ID | 阶段/节点 | 任务名称 | 状态 | 优先级 | 依赖 | 交付物与验证标准 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`T14`** | Node 7 (`S3.01–02`) | 实现 `csr_splice_dist.cpp`（零展开 `GRangesList`、零 `myGaps()`、零 `tapply` 的多线程双指针匹配核）与 `modifyIncompatibleAssignment` | `DONE` | `P0` | `T13` | `ENCFF089AMF` (`5.60M` hits) 与 `human_all_d1q1` 16 BAM（`nc=8/16`, `S3 356s->13.4s`）彻底消除 `250G` Fork CoW 并全量等价 |
| **`T15`** | Node 8 (`S3.03–08`) | 实现复用一次 `ext_min_eq_csr` 的 `equirc_builder.cpp` + 融合无 `3× dcast` 分组并发 `em_l1_fused.cpp` + 预分配矩阵 `combineCountSes` | `DONE` | `P0` | `T14` | `counts` / `CPM` / `fullLengthCounts` / `uniqueCounts` 最大误差 `< 1e-6`，`r >= 0.999999`；与原版 `_bambu_emWithL1` 数学对齐 (`089/hu/mouse/switches` 全 0 误差/r=1.0 PASS) |

---

## Phase 5：Stage 4 流式导出与完整 R 包封装验证 (`Node 9`)

| ID | 阶段/节点 | 任务名称 | 状态 | 优先级 | 依赖 | 交付物与验证标准 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`T16`** | Node 9 (`S4`) | 实现 C++ 缓冲流式 `fast_export.cpp`（`writeToGTF` / `writeCountsOutput`）、参考引用 `transcriptToGeneExpression`、优化 `compareTranscripts_fast.R` | `DONE` | `P1` | `T15` | 导出 GTF & 5 类 TSV 与原生逐记录对齐（`GTF` 与数值表格等效检验 `089/hu/mouse/extdata` 全通过） |
| **`T17`** | Node 9 | 串接完整 Drop-in [`bambu()`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu_next/docs/plan.md) 入口，跑通全部 `testthat` 11 套官方用例 + 全可选参数开关矩阵 (`stranded`, `trackReads`, `fusionMode`, `lowMemory`, `yieldSize`, `NDR=NULL`, `returnDistTable`) | `DONE` | `P0` | `T16` | 官方非绘图 `testthat` `56/56 PASS` + `extdata 6 模式 PASS` + `switches_a (5/5) & switches_b (8/8) 100% PASS` |

---

## Phase 6：大规模基准对比、全协议 Slurm Array 验收与交付 (`Gates G1 ~ G5`)

| ID | 阶段/节点 | 任务名称 | 状态 | 优先级 | 依赖 | 交付物与验证标准 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`T18`** | Gate G1 | Slurm Array 验证 S1 单样本扩展（`ENCFF400BQQ`, `ENCFF436GKZ`, `ENCFF846RID`, `ENCFF089AMF`, `ENCFF887KHY`, `ENCFF598YQO` 于 `1/4/8/16/32` 线程与 `.bbidx` 热跑） | `DONE` | `P0` | `T17` | 全部 19 项 Slurm 任务 COMPLETED (0:0)，L1 逐项 `identical`；`ENCFF089AMF 19.02x`, `887KHY 18.21x`, `v46 GTF 14.50x` |
| **`T19`** | Gate G2 | Slurm `A-B-B-A` ×3 端到端完整对比：`human_single`, `mouse_single`, `human_pair`, `mouse_all` + 可选参数开关矩阵 | `DONE` | `P0` | `T18` | 全部 7/7 队列与 7/7 开关 COMPLETED (`8.9x–21.06x` 加速，全项 `gtf=TRUE, r=1.0`，`mouse_all 14.42x–15.04x / 2.8G`) |
| **`T20`** | Gate G3 | **`human_all_d1q1` 16 Human BAM 全量联跑核心验收**（`A-B-B-A` ×3 vs `native` / `G8`，评测 `ncore=8` & `ncore=16`） | `DONE` | `P0` | `T19` | `nc=16: 2055.5s->256.9s (8.00x), 79.2->7.0 GiB (11.33x)`；`nc=8: 2490.5s->273.3s (9.11x), 49.4->6.6 GiB`；`G17 NDR=1: 4779s->733.7s (6.51x), 250.5->44.6 GiB`；GTF 结构 100% 零差异 (`r=1.0`) |
| **`T21`** | Gate G4 | `longpass_benchmark250209/` 无上限 Slurm Array 扩展验证：`RBMX_KD/shRBMX_rep2` + `gencode.v46` (1.6G GTF)、`SIRV_simulation` 真值评估、全协议 human/mouse 样本及海牛 `de novo`（只读 `.bai` 软链构建） | `DONE` | `P0` | `T20` | 全 16 项阵列 0 OOM / 0 失败（`8.90x–28.43x`）；SIRV `F1=0.886689` 6 位小数一致（`28.43x`）；`v46`、12 测序协议及海牛 `PB/ONT de novo (18.96x/25.57x)` 均通过 |
| **`T22`** | Gate G5 | 生成结构化汇总 TSV/JSON、绘制多图统计分析（耗时分解、Peak Tree PSS 对比、加速比、一致性散点/相关性、协议分面图）与编译整合单文件 PDF 报告 | `DONE` | `P1` | `T21` | `docs/data/*.{tsv,json}`、`docs/benchmark_report.md`、`docs/figures/*.png`（fig1–6 + 架构图）、`docs/bambu_next_benchmark_report.pdf` (6页) 与 `docs/bambu_next_benchmark_presentation.{pptx,pdf}` (12页) 均已完成并逐页视觉检视零遮挡 |
| **`T23`** | Gate G5 | 彻底清理 `/gpfs/.../bambu_next_work/` 大型 scratch 与临时软链，执行只读目录未篡改审计、`system /tmp` 零写入与 `/dev/shm` 零泄漏核查，更新本台账至全绿完成 | `DONE` | `P0` | `T22` | `bambu_next_work` 缩减至 `25 MiB`（零软链、零 export/rds 残留）；`RECENT_TMP=0`、`SHM=0`；全只读目录零新文件；`bambu 3.13.99` 终态 `testthat 56/56 PASS` (`test_plotBambu.R` 因缺可选依赖 `ggplot2` 跳过) |
