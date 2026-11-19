# `bambu_next` 分节点重构与底层优化技术实施计划 (`PLAN`)

> 本文档在 [bioconductor_manuals_audit.md](bioconductor_manuals_audit.md)、[spec.md](spec.md)、[optimization_schemes.md](optimization_schemes.md) 与 [roadmap.md](roadmap.md) 确定主架构（**方案 B** + `envs/bambu` + `envs/basicpy`）的基础上，给出项目骨架与 `Node 0 ~ Node 9` 的模块设计、C++/R 内存边界、数据结构定义与验收测试路径。

---

## 1. 目标仓库结构骨架 (`bambu_next/`)

`bambu_next/` 既作为可独立安装的 R/C++ 扩展包（通过 `R CMD INSTALL -l bambu_next/.r_lib` 安装隔离，绝不污染 `envs/bambu`），也包含完整的 Slurm/PSS 基准验证与可视化脚本套件：

```text
bambu_next/
├── README.md                              # 总览与快速使用指引
├── DESCRIPTION                            # R 包元数据：Depends/Imports 与 LinkingTo (Rcpp, RcppArmadillo, Rhtslib, S4Vectors, IRanges, XVector, Biostrings, data.table)
├── NAMESPACE                              # 完整导出原有 Bambu API (bambu, prepareAnnotations, readFromGTF, writeToGTF, writeBambuOutput, transcriptToGeneExpression, trainBambu, compareTranscripts, plotBambu)
├── docs/                                  # Phase Pre & 全流程设计/执行文档体系
│   ├── bioconductor_manuals_audit.md      # 依赖包 Bioconductor 官方手册/头文件审计
│   ├── spec.md                            # 契约规格、冻结测试矩阵与度量口径
│   ├── optimization_schemes.md            # G8–G17 瓶颈解剖与三方案对比选型
│   ├── roadmap.md                         # M0_PRE ~ M10 与 Gates G0–G5 路线图
│   ├── plan.md                            # 本文档：Node 0 ~ Node 9 代码级详案
│   └── task.md                            # T00 ~ T23 可追踪任务台账
├── src/                                   # C++17 高性能内核 (Rcpp + Rhtslib + IRanges/S4Vectors C stubs + std::thread/OpenMP)
│   ├── Makevars                           # 安全拼接 Rhtslib::pkgconfig("PKG_CPPFLAGS"/"PKG_LIBS") 与 -fopenmp -O3
│   ├── init.cpp                           # R 原生符号注册 + S4Vectors/IRanges stubs 初始化
│   ├── common/
│   │   ├── flat_csr_index.hpp             # Node 1/6/7 共用：只读 Flat-CSR 参考/扩展注释结构 & 序列化
│   │   ├── fast_hash.hpp                  # 64-bit SplitMix/XXH3-like 整型特征哈希与动态扩容线性探测表
│   │   └── thread_pool.hpp                # Node 0：可控边界 C++17 工作线程池（零接触 R API/SEXP）
│   ├── annotation/
│   │   ├── gtf_stream_parser.cpp          # Node 1：单遍流式 GTF 解析与 exon ranking
│   │   └── min_eq_class_builder.cpp       # Node 1/6：全局一次参考 minEqClass 与增量更新内核
│   ├── stage1/
│   │   ├── bam_exon_extractor.cpp         # Node 2：Rhtslib 直接读 BAM (hts_itr_queryi/.bai + CIGAR + XS/strand + readIds)
│   │   ├── faidx_motif_lookup.cpp         # Node 2：htslib/faidx.h 批量提取唯一 junction 2-bp 碱基判定 spliceMotif
│   │   ├── junction_corrector.cpp         # Node 3：剪接链投票校正、高频 junction 映射归并
│   │   ├── read_class_builder.cpp         # Node 3：XXH3 多键哈希归并 Spliced/Unspliced ReadClass（零 dplyr/splitAsList）
│   │   └── rc_scorer_uf.cpp               # Node 4：Union-Find 连通组件 assignGeneIds、兼容判定、polyA/T
│   ├── stage2/
│   │   ├── kway_rc_merger.cpp             # Node 5：S2 K-way spliced 合并 + O(U) 无展开精确 weightedMedian
│   │   └── unspliced_cluster_merge.cpp    # Node 5：S2 无 separate_rows 整数区间聚类合并
│   ├── stage3/
│   │   ├── csr_splice_dist.cpp            # Node 7：零展开双指针 findSpliceOverlapsByDist + compareTranscripts 内核
│   │   ├── equirc_builder.cpp             # Node 8：直接生成分组 CSR equiRC + 复用参考 minEquiRC（无 3×dcast/unnest）
│   │   └── em_l1_fused.cpp                # Node 8：保持原浮点逻辑的多群组并发 emWithL1 / em_theta
│   └── stage4/
│       └── fast_export.cpp                # Node 9：C++ 缓冲流式 GTF/counts TSV 输出与快速 gene reduction
├── R/                                     # Bioconductor 接口适配与业务流程编排层
│   ├── bambu_next.R                       # 公开 bambu() Drop-in 入口与参数校验
│   ├── config_and_guardrails.R            # Node 0：thread/data.table/xgboost/BiocFileCache 竞争护栏
│   ├── prepareAnnotations_fast.R          # Node 1：S0 GTF/TxDb -> PreparedAnnotation + .bbidx CSR
│   ├── stage1_processReads.R              # Node 2–4：S1 C++ 内核胶水、XGBoost 模型推理/trainBambu、SE 封装
│   ├── stage2_extendAnnotations.R         # Node 5–6：S2 合并、筛选、推荐 NDR、增量扩展注释
│   ├── stage3_quantify.R                  # Node 7–8：S3 C++ 距离匹配、融合 EM 定量、combineCountSes
│   ├── stage4_export_and_gene.R           # Node 9：transcriptToGeneExpression、writeBambuOutput、writeToGTF
│   ├── compareTranscripts_fast.R          # 兼容原版 compareTranscripts()（可选直调 CSR 核）
│   └── plotBambu.R                        # 原位保留兼容可视化函数
├── tests/                                 # testthat 单元测试 + 与原版逐模块 Oracle 一致性对照套件
└── scripts/
    ├── install_local.sh                   # 编译并安装至 bambu_next/.r_lib（零改动 conda envs/bambu）
    └── benchmark/
        ├── monitor_process_tree.sh        # /proc/<pid>/smaps_rollup 进程树 PSS/RSS 精确采样器
        ├── compare_bambu_outputs.R        # L1 结构 identical() + L2 定量精确率/相关性比对器
        ├── run_paired_ab_ba.py            # A-B-B-A 交替基准执行器（strict no-/tmp & 关闭调试 dump）
        └── build_slurm_cohort_and_report.py # 全矩阵 Slurm Array 生成、结果汇总、多图与单文件 PDF 报告
```

---

## 2. 分节点技术实施方案（`Node 0 ~ Node 9`）

### `Node 0`：工程骨架、并发护栏与观测基座（`Phase 0`）
- **Makevars 编译链接修复**：分别获取 `Rhtslib::pkgconfig("PKG_CPPFLAGS")` 与 `Rhtslib::pkgconfig("PKG_LIBS")`（规避连调无空格 `-lz-D_FILE_OFFSET_BITS=64` 缺陷），同时加入 `IRanges`/`S4Vectors`/`XVector`/`Biostrings` C stub 头文件路径与 `-std=gnu++17 -fopenmp`。
- **并发与环境护栏**：
  - 函数入口统一执行 `data.table::setDTthreads(1L)`、锁定 `xgboost` 预测参数 `nthread = 1L`、抑制隐式 BLAS 线程风暴；
  - 针对 `>10` 样本临时 read-class 磁盘分支，主进程预先在当前运行工作区独占子目录完成空 `BiocFileCache` 初始化（解决 `g8` 文档记录的 SQLite 首次建表竞态 `table metadata already exists`，且绝不使用 `/tmp`）。
- **基准工具部署**：编写无 `/tmp` 依赖的 `monitor_process_tree.sh`、`compare_bambu_outputs.R` 与 `install_local.sh`，并先在 `envs/bambu` 跑通官方 `bambu/tests/testthat` 确立 `G0` 基线。

### `Node 1`：流式 `GTF` 编译器、`.bbidx` 零拷贝 CSR 与预计算 `minEquiRC`（`S0`）
- **对应原函数**：[`prepareAnnotations.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/prepareAnnotations.R), [`prepareAnnotations_utilityFunctions.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/prepareAnnotations_utilityFunctions.R), [`getMinimumEqClassByTx`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/bambu-extendAnnotations-utilityExtend.R)。
- **核心设计**：
  1. **单遍 C++ 流式 GTF 解析 (`gtf_stream_parser.cpp`)**：按行流式解析 GTF `exon` 记录、抽提 `transcript_id`/`gene_id`/`strand`，对负链与正链分别按原版 `GenomicFeatures` 一致规则赋予 `exon_rank` / `exon_endRank`，生成标准 `GRangesList` 并**同步写入扁平只读 `FlatCSRAnnotation`**（以 `Rcpp::XPtr` 附着在 `attr(ann, "bbidx_ptr")` 或保存 `.bbidx` 二进制文件）。相比先建 SQLite `TxDb` 再多次 SQLite 查询，解析 `gencode.v46` (1.6G) 预计可从数十秒压至数秒且大幅降低 R 堆膨胀。
  2. **原生 `prepared_annotation.qs` / `TxDb` 零冗余转换**：若用户直接传入已有的 `.qs`（如冻结对比必需的 `human.v39.prepared_annotation.qs`）或 `TxDb` 对象，通过 `IRanges`/`S4Vectors` C 接口一次性线性扫过底层 `int* start/width` 指针直接填充 `FlatCSRAnnotation` 外部 C++ 内存。
  3. **全局一次性构建参考 `minEqClass` 拓扑**：在 `S0` 直接利用 C++ 区间重叠与前缀/后缀匹配构建参考 `tx -> min_eq_tx_set` 拓扑及双向 CSR，并固化为参考句柄的一部分——彻底消除原生流程在 `S2` 重算 `245k × 245k` 与在 `S3` 16 个样本各自重复计算 `16 × processMinEquiRC(annotations)`！

### `Node 2`：`Rhtslib` 并发流式 BAM 解析与参考 Junction Motif 查询（`S1.04–S1.07`）
- **对应原函数**：[`prepareDataFromBam.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/prepareDataFromBam.R), [`bambu-processReads_utilityCreateJunctionTables.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/bambu-processReads_utilityCreateJunctionTables.R)。
- **核心设计**：
  1. **`Rhtslib` 多线程 / 索引分片 BAM 解码 (`bam_exon_extractor.cpp`)**：
     - 对有 `.bai`/`.csi` 索引的 BAM（单样本和多样本场景统一调度！），按染色体或固定参考区间窗口将读取工作切分给 C++ worker 线程（各线程独占自己的 `samFile*` + `hts_idx_t*` 句柄，完全线程安全！）。如果遇到无 `.bai` BAM 且不允许写原目录（如海牛/部分 directRNA），按规格在工作区目录建只读软链并生成临时 `.bai` 或退回单 `samFile` + `hts_set_threads` 多核 BGZF 解压！
     - 精确复现原版语义：过滤 unmapped (`0x4`) 与 secondary (`0x100`)（支持 `supplementary=0x800` 与原版 `readGAlignments` 一致保留规则），由 `bam_get_cigar()` 遍历 `BAM_CMATCH / BAM_CEQUAL / BAM_CDIFF / BAM_CDEL` 合并为相邻 exon blocks，由 `BAM_CREF_SKIP` 识别内含子间隙，精确保留原版“首/末或单 exon ≥1bp 检查、排序及 `seqlevels` 越界裁滤”，支持 `yieldSize` 分块，并将结果按原始 BAM `POS` 与文件序列稳定归并，**绝不构造 `GAlignments` 与千万级 `GRangesList`**！
  2. **零 `GRanges` 唯一 Junction 提取与 `faidx.h` 基序查表 (`faidx_motif_lookup.cpp`)**：
     - 流式解出 exon blocks 时同步用 `(seqid, intron_start, intron_end)` 写入本线程 `unordered/open-addressing` junction hash table（累加 `+`/`-`/unknown 计数，并记录每个 read 所属 splice 边引用 `uint32_t junction_uid`），**完全消除 `unlistIntrons()` 构建 `57.68M` `GRanges` 的巨额分配**！
     - 对去重后的数十万~数百万条 `uniqueJunctions` 按 `(seqid, pos)` 排序后调用 `htslib/faidx`（`faidx_fetch_seq`）顺序扫描读取 5'/3' 各 2-bp 碱基判定 `spliceMotif`（若传入内存 `BSgenome` 对象，则在 R 侧逐条染色体转为原生 `raw` ASCII/2-bit 视图后统一由 C++ 批量索引查找），零构造 `DNAStringSet`！

### `Node 3`：剪接链投票校正、`XGBoost` 批推理、Spliced & Unspliced ReadClass 构建（`S1.09–S1.12`）
- **对应原函数**：[`bambu-processReads_utilityJunctionErrorCorrection.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/bambu-processReads_utilityJunctionErrorCorrection.R), [`bambu-processReads_utilityConstructReadClasses.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/bambu-processReads_utilityConstructReadClasses.R)。
- **核心设计**：
  1. **链特异性 & 剪接链投票 (`updateStrandScoreByRead`)**：利用紧凑 `(read_uid, junction_uid)`CSR 索引直接按 read 与按 junction 双遍整型投票累加（3 步向量扫描），彻底替代 `countMatches` 在 57.68M row 表上来回排序匹配（原 `ENCFF089AMF` `46.26s` 预计缩至 `< 1s`）。
  2. **XGBoost 剪接位点校正 (`findHighConfidenceJunctions` + 特征矩阵)**：用 C++ 一次计算唯一 junction 与邻域支持窗特征矩阵直接喂给 `xgboost`（强制 `nthread = 1L`，避免 R `data.frame` 拼接），并用 C++ 并查集映射高频 junction 坐标校正 (`correctIntronRanges`)。
  3. **`createReadTable` & `constructUnsplicedReadClasses` C++ 化 (`read_class_builder.cpp`)**：
     - 原版 S1 最大单热点（`createReadTable` **`107.85s`**，`constructUnsplicedReadClasses` **`41.78s`**）。
     - `bambu_next` 将读段按内含子链序列 `(seqid, strand_code, intron_start[0..k-1], intron_end[0..k-1])` 直接做 **64-bit 向量哈希 + 坐标逐项比对归并**；单遍流式跟踪最小 `readStart`、最大 `readEnd`、精确或按原版分位数规则记录起止界值、`readCount`（以及按 `trackReads` 选项记录紧凑 `uint32` read 索引）。
     - Unspliced reads 按 `(seqid, strand_code)` 排序后通过 C++ 双指针归并与 `FlatCSR` 参考转录本边界匹配，**彻底告别 `as.data.frame`、`dplyr::group_by(..., nth)`、`unstrsplit(splitAsList(...))`**。

### `Node 4`：Union-Find `assignGeneIds`、ReadClass 评分与并行 `trainBambu`（`S1.14`）
- **对应原函数**：[`bambu-processReads_scoreReadClasses.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/bambu-processReads_scoreReadClasses.R), [`bambu-extendAnnotations-utilityExtend.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/bambu-extendAnnotations-utilityExtend.R) (`assignGeneIds`)。
- **核心设计**：
  1. **C++ Sweep-Line + Union-Find 基因连通重编号 (`rc_scorer_uf.cpp`)**：原 `assignGeneIds` 在 `ENCFF089AMF` 耗时 `24.52s`。新内核基于 `FlatCSRAnnotation` 查找与 ReadClass 外显子重叠的转录本/基因后，直接在 C++ 使用 **Disjoint Set (Union-Find)** 与按原版顺序规则分配 `gene_id` / novel gene IDs，不调用 `igraph`/多轮 `dplyr` join。
  2. **评分向量 (`isReadClassCompatible`, polyA/T, `txScore`)**：直接在 C++ 提取 polyA/T 末位频率特征与兼容标签矩阵，配合单次批量 `predict(xgb_model, mat, nthread = 1L)` 输出 `txScore` 与 `txScore.noFit`；并支持将多样本轻量化 `ReadClassCompact` 直接保留在内存或按紧凑二进制写出（替代耗时的 `648 MB` `saveRDS(readClassSe)` 携带巨大 string 列表）。

### `Node 5`：`S2` 跨样本合并与精准无展开加权分位数（`S2.01–S2.03`）
- **对应原函数**：[`bambu-extendAnnotations-utilityCombine.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/bambu-extendAnnotations-utilityCombine.R)。
- **核心设计**（直接击溃 `human_all_d1q1` 耗时 `1545.2s` 的 S2 阶段，特别是 `S2.01 824.29s` & `S2.02 289.59s`）：
  1. **K-way 哈希合并替代 16 次 `dplyr::full_join` (`kway_rc_merger.cpp`)**：所有样本的 spliced ReadClass（已以紧凑结构或整型 CSR 存储）一次汇入 C++ 多样本归并哈希表，按键合并边界坐标、累加 `readCount`、汇总 `NSampleReadCount`/`maxTxScore` 等指标，复杂度从 $O(S^2 \cdot N)$ 列拷贝降至线性 $O(S \cdot N)$！
  2. **杜绝 `readCountWeightedMedian` 内存爆炸**：原版直接 `rep(values, times = counts)`，当某个高度表达模型跨样本计数高达百万时极度耗内存。新算法将 `(value_i, weight_i)` 按 `value_i` 排序并在累积权重扫至总权一半时直接提取原版精确一致的加权中位数，`O(U log U)` 时间、`O(U)` 空间！
  3. **整数区间聚类替代 `combineUnsplicedTranscriptModels` 字符串切分**：直接扫排好序的 `(seqid, strand, start, end, gene_id)` 整数表合并单外显子候选，**完全移除 `paste(..., collapse="+")` 与 `tidyr::separate_rows()`**！

### `Node 6`：注释筛选、推荐 NDR 与增量式 `combineWithAnnotations`（`S2.04–S2.08`）
- **对应原函数**：[`bambu-extendAnnotations.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/bambu-extendAnnotations.R), [`bambu-extendAnnotations-utilityExtend.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/bambu-extendAnnotations-utilityExtend.R)。
- **核心设计**：
  1. **`addNewSpliced/UnsplicedReadClasses` + `filterTranscriptsByAnnotation` / `recommendNDR`**：借助预先建好的 `FlatCSRAnnotation`，仅需极小开销比对过滤候选并维持与原版一致的 `NDR` 排序和推荐逻辑。
  2. **增量更新 `getMinimumEqClassByTx` (`min_eq_class_builder.cpp`)**：这是原 `isore.extendAnnotations` (`333.75s`) 的元凶——它把 ~245k 参考转录本和新增转录本拼接后，**全量重做一次 $245\text{k} \times 245\text{k}$ 的 `findSpliceOverlapsQuick`**！新引擎直接克隆参考 `minEqClass` CSR 图，**仅对与新发现转录本在同一重叠区间簇的候选转录本执行双向等价类校验并补边**，计算量缩小两个数量级以上！

### `Node 7`：核心战役——C++ Flat-CSR 零展开并发 `calculateDistToAnnotation`（`S3.01–S3.02`）
- **对应原函数**：[`utility_spliceHelper_functions.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/utility_spliceHelper_functions.R), [`bambu-processReads_scoreReadClasses.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/bambu-processReads_scoreReadClasses.R) (`isore.estimateDistanceToAnnotations`), [`bambu-quantify_utilityFunctions.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/bambu-quantify_utilityFunctions.R)。
- **重点目标**：`human_all_d1q1` 峰值 **`250.48–266.57 GiB` PSS** 与 **`2596s` 累计耗时**的正主！
- **核心设计 (`csr_splice_dist.cpp`)**：
  1. **只读共享 `extendedAnnotation CSR` + 多线程 `std::thread`**：主进程唯一一份只读 `FlatCSR` 扩展注释常驻 C++ 内存；针对单样本或多样本队列，直接通过 `std::thread` / OpenMP 并行计算所有 ReadClass 与扩展注释重叠对（**绝不执行 `bplapply` `fork()` 携带 8 份 R 堆，线程内 0 次分配任何 `SEXP`！**）。
  2. **消除 `51.72M` 区间拷贝、`myGaps()` 与 `base::tapply`**：
     - 对第 `i` 条 ReadClass 外显子数组 $Q_i = [(s_0, e_0), \dots, (s_{m-1}, e_{m-1})]$ 与候选第 `j` 条 Subject 转录本 $S_j = [(a_0, b_0), \dots, (a_{n-1}, b_{n-1})]$：
       - 利用有序性 $O(m + n)$ 双指针扫描求得所有 `pintersect` 重叠宽、`uniqueLengthQuery`、`uniqueLengthSubject`、超出 `min.exonDistance` 的起止偏差 `uniqueStartLengthQuery`/`uniqueEndLengthQuery` 与总不匹配外显子标志 (`queryElementsOutsideMaxDist`, `subjectElementsOutsideMaxDist`)，**在 C++ 栈标量上直接累加！**
       - 内含子一致与兼容 (`spliceCompatible`, `equal`) 通过比较内部剪接边界（excl. outer 5'/3'）一次扫完。
     - 计算结果（紧凑整数向量：`query_idx`, `subject_idx`, `dist`, `compatible`, `equal`, `n_match`）直接供后续过滤 (`genFilteredAnTable` 无需分 3 轮过滤重复扫描!) 与 `modifyIncompatibleAssignment` 使用。
  3. 同步暴露该内核给公开工具函数 [`compareTranscripts()`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/compareTranscripts.R) 作为高速后端。

### `Node 8`：复用 `minEquiRC` 的 `genEquiRCs` 与融合无 `dcast` 多群组 `emWithL1`（`S3.03–S3.08`）
- **对应原函数**：[`bambu-quantify.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/bambu-quantify.R), [`bambu-quantify_utilityFunctions.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/bambu-quantify_utilityFunctions.R), [`em.cpp`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/src/em.cpp)。
- **核心设计**：
  1. **`addEmptyRC` 全局仅一次构建**：直接接入 `Node 1/6` 产出的一份扩展注释最小等价类列表 `ext_min_eq_csr`（消除原生在 16 个样本中各自展开 16 次 `tidyr::unnest` + `dplyr::full_join`）。
  2. **告别 `dplyr::group_by(list_col)` 与 `3 × dcast()` (`equirc_builder.cpp` + `em_l1_fused.cpp`)**：
     - 直接以整数序列排序归并统计 `(eq_class_tx_set, unique_count, full_length_count, n_obs)`（无需先转 string key 再调 R `dplyr::group_by`）。
     - `assignGroups` 基于 C++ 并查集将相互共享等价类的转录本合并为独立连通定量簇 (`gene_grp_id`)；直接将各 `gene_grp_id` 的局部 `A`, `A_full`, `A_unique`, `Y`, `K` 在 C++ 中组装为紧凑双精度小矩阵并在多线程间并行运行原版同公式 `emWithL1`（若单个超大基因簇可继续借用 OpenMP/BLAS，绝不经 R `for` 循环调 3 次 `data.table::dcast()`！）。
  3. **`combineCountSes` 低拷贝汇总**：所有样本直接往预分配好的最终 `num_tx × num_samples` 双精度数值矩阵填列，告别先构造 16 份大 `RangedSummarizedExperiment` 再调 `do.call(cbind, ...)`！

### `Node 9`：C++ 流式 GTF/TSV 导出、基因表达汇总与完整 CLI/Drop-in 验证（`S4`）
- **对应原函数**：[`readWrite.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/readWrite.R), [`transcriptToGeneExpression.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/transcriptToGeneExpression.R), [`bambu.R`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/R/bambu.R)。
- **核心设计**：
  1. **C++ 缓冲流式 `writeToGTF` & `writeCountsOutput` (`fast_export.cpp`)**：直接遍历 `rowRanges`/`FlatCSR` 与最终计数矩阵按 4 MiB chunk 写出规范 GTF 与 TSV（完全绕开 R 几十万行 `data.frame` + `paste0`）。
  2. **`transcriptToGeneExpression` 加速**：用原位矩阵按 `GENEID` 行组累加（`rowsum()` 或 C++ 索引累加）替代 `melt` + 双 `dcast`，保留精确一致的 `incompatibleCounts` 并入与 `reducedRangesByGenes`。
  3. **严格向后兼容边界回退策略**：当收到极端非常规输入（如罕见第三方自定义非标准对象无法匹配 C 指针契约）时记录明确状态并支持优雅回退，确保公开 API 零破坏；全面跑通 `testthat` 与 `G0–G5` 全矩阵 Slurm Array 对比与交付。
