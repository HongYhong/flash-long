# Bambu 核心依赖包 Bioconductor 官方手册与 C 接口深度审计 (`Bioconductor Manuals Audit`)

> **环境版本基准**（实测于 `/gpfs/group/data_sci_group/weiwu_group/hongyanhong/miniconda3/miniconda3/envs/bambu`）：  
> `R 4.5.3` | `bambu 3.12.1` | `BiocParallel 1.44.0` | `GenomicAlignments 1.46.0` | `GenomicRanges 1.62.1` | `GenomicFeatures 1.62.0` | `IRanges 2.44.0` | `S4Vectors 0.48.0` | `Rsamtools 2.26.0` | `Rhtslib 3.6.0` | `SummarizedExperiment 1.40.0` | `GenomeInfoDb 1.46.2` | `BSgenome 1.78.0` | `Biostrings 2.78.0` | `data.table 1.17.8` | `Rcpp 1.1.1` | `RcppArmadillo 15.2.6-1` | `xgboost 3.2.0.1`。

---

## 1. 为什么既往 `G8–G17` 屡试屡滞：从 Bioconductor 官方手册找答案

通读 `envs/bambu/lib/R/library/*/doc/` 与 `*/include/` 中的一线文档后，可以明确看出原版 Bambu（以及沿用原版对象图的前几代候选）在架构上**违背了 Bioconductor 官方手册反复强调的多条可扩展性原则**：

---

## 2. 逐包手册精要与面向 `bambu_next` 的工程映射

### 2.1 `BiocParallel 1.44.0`（`Introduction_To_BiocParallel.Rmd`）
- **官方手册核心警告**（L162–179, L364–463, L771–785）：
  1. **Fork CoW 陷阱**：`MulticoreParam` 基于 POSIX `fork()`。虽然文档称非修改状态可共享只读页，但只要 R 子进程内分配对象并触发 **R 分代垃圾回收器（Generational GC）标记共享对象头部的 `sxpinfo` GC mark bit**，内核会瞬间对整个大对象的内存页触发 **Copy-on-Write (CoW) 页分裂**，外加 I/O / 数据库句柄共享竞态。
  2. **官方推荐模式**：引用 Lawrence & Morgan《*Scalable Genomics with R and Bioconductor*》——四大支柱：**限制查询范围 (restricting queries)**、**压缩数据 (compressing data)**、**块迭代 (`bpiterate` / streaming chunks)**、**无巨型 R 对象依赖的并行计算**。
  3. **多进程 vs 线程池隔离**：对可重用固定对象集群优先 `bpstart()` / `bpstop()` 或将大只读索引留在进程外/原生堆（Non-R-GC C++ heap / `mmap` 只读匿名共享映射）。
- **对 `bambu_next` 的直接影响**：
  - **`human_all_d1q1` 峰值 `250.48–266.57 GiB` PSS 的根本推手就是 `bplapply(..., bambu.quantify)` 在 8 个 `fork` R 进程中各自持有并遍历巨型 `annotations`/`extendedAnnotations` S4 堆且同时爆发数百 GB S4 临时分配引发的全盘 CoW！**
  - **对策**：重计算与巨型基因组索引全部脱离 R GC 堆，存放于 C++ `std::vector` / 只读扁平 CSR 结构中，内部用 **C++17 `std::thread` 或 OpenMP** 并行扫描，保证并行线程零触碰 R `SEXP`、零 `MulticoreParam` 8 份 R 堆复制。

### 2.2 `Rhtslib 3.6.0` & `Rsamtools 2.26.0`（`Rhtslib.Rmd`, `Rsamtools-Overview.Rmd`, `include/*.h`）
- **官方能力实测确认**：
  1. `envs/bambu` 内 `Rhtslib::pkgconfig("PKG_LIBS")` 与 `Rhtslib::pkgconfig("PKG_CPPFLAGS")` 可直接链接静态库 `libhts.a`（`-lcurl -lbz2 -llzma -lz`），且包含完整 `sam.h`、`hts.h`、`faidx.h`、`bgzf.h`、`thread_pool.h`。（注意实测 quirk：连在一起调用会缺少空格导致 `-lz-D_FILE_OFFSET_BITS=64`，`Makevars` 中需分开换行调用。）
  2. `Rsamtools::scanBam` / `scanFa` / `getSeq,FaFile` 本质上也是单线程封装 `.Call()` 进入 htslib，然后构造沉重的 R `list` / `DNAStringSet`。而原版 `bambu::prepareDataFromBam` 还额外多包一层 `GenomicAlignments::readGAlignments()` + `grglist()`！
- **对 `bambu_next` 的直接影响**：
  - **Node 2 (`S1.04`) BAM 解析革新**：直接基于 `Rhtslib` (`sam_open`, `hts_set_threads` / `bgzf_mt` 或按 `.bai` 分区索引并发解码) 扫描 `bam1_t`，单遍在 C++ 中解析 CIGAR (`bam_get_cigar`) 生成 exon blocks（精确复现原版过滤 `secondary=0x100`、1bp block 合并/裁剪、`flag` 链向与 `XS` tag 提取逻辑），直接填入连续 `int32` CSR 缓冲。彻底跳过构建千万级 `GAlignments` 与 `GRangesList`！
  - **Node 2 (`S1.07`) 剪接基序查询 (`faidx.h`)**：当参考为 FASTA 路径（或 `.fa.fai`）时，直接在 C++ 使用 `fai_load` / `faidx_fetch_seq` 对唯一 junction 批按 `(chr, start)` 顺序缓存 5' 与 3' 2-bp 碱基判定 `GT-AG`/`CT-AC`/`GC-AG`/`CT-GC`/`GT-AT`/`AT-AC` 并统计 motif 支持度，省去 `BSgenome::getSeq` + `Biostrings::reverseComplement` 构建数百万 `DNAString` 对象的耗时与堆膨胀（当用户输入是 R 包形式 `BSgenome` 对象时，通过 `Biostrings` / 提取染色质字节序列缓存保留 100% 兼容）。

### 2.3 `GenomicAlignments 1.46.0`（`GenomicAlignmentsIntroduction.Rnw`, `OverlapEncodings.Rnw`, `.findSpliceOverlaps`）
- **官方源码审计关键发现**：
  - 检查 `GenomicAlignments:::.findSpliceOverlaps` 定义可见：它执行 `olap <- findOverlaps(query, subject)` 后，**物理执行了 `query <- query[queryHits(olap)]`, `subject <- subject[subjectHits(olap)]`, `splice <- .gaps(query)`**。
  - 原版 Bambu 的 `utility_spliceHelper_functions.R`（`findSpliceOverlapsQuick`, `findSpliceOverlapsByDist`, `compareTranscripts`）**几乎照搬了这个扩展模式**，而且在 `ENCFF089AMF` 上有 **`5,604,315` 条重叠对**，每对又展开平均 ~9 个外显子 = **`51.72M` 个区间元素**，紧接着调用 `myGaps()`、`setdiff()`、`intersect()`、`pintersect()` 和 **`base::tapply(..., myId, sum)`**，单样本累计瞬时分配就高达数百 GB！
  - `OverlapEncodings.Rnw` 手册的核心思想本是利用外显子/内含子位置编码避免完整基因组对象展开，但 Bambu 实际落成了全量 `GRangesList` 复制。
- **对 `bambu_next` 的直接影响**：
  - **Node 7 (`S3.01` `calculateDistToAnnotation` / `findSpliceOverlapsByDist`) 与 Node 5 (`S2.04–S2.05` `findSpliceOverlapsQuick`)**：在 C++ 中以紧凑 CSR 数组直接表示 ReadClass（query exons）和 Transcripts（subject exons），对每条 `(query_i, subject_j)` 候选对直接**原位双指针扫描**计算：
    - exon-by-exon 重叠长度 `overlap_bp`；
    - `queryElementsOutsideMaxDist` / `subjectElementsOutsideMaxDist`（等价 `rangesDist`）；
    - 首/末外显子匹配与起止偏移 `startMatch` / `endMatch` / `dist`（等价 `calculateFirstLastExonsDist`，**零 `tapply` 字符串分组，零 `myGaps()` 对象生成**）；
    - `spliceCompatible` / `equal` 标志（复现原版所有 `maxDist`/单外显子/multi-exon 边界细节）。

### 2.4 `GenomicRanges 1.62.1` & `IRanges 2.44.0`（`GNCList-class.Rd`, `NCList-class.Rd`, `Ten_things_slides.Rnw`, `IRangesOverview.Rnw`, `IRanges_interface.h`）
- **官方手册最佳实践**：
  1. `GNCList` / `NCList` 明确指出：Nested Containment List 构建开销不低，但若**预先构建一次 `GNCList(subject)`** 并在多次查询间复用，可彻底避免每次 `findOverlaps` 重复建树。
  2. 官方明确建议：不需要具体命中的场景用 `countOverlaps()`（或直接在 C 级遍历 `NCList`）以防止在 R 堆中实例化海量 `Hits` 对象。
  3. `IRanges` 与 `S4Vectors` 公开了标准 C ABI (`IRanges_interface.h`, `S4Vectors_interface.h`, `_IRanges_stubs.c`, `_S4Vectors_stubs.c`)，可直接零转换从 R `IRanges`/`CompressedIRangesList`/`Rle` 读取 `int*` 指针与 `PartitioningByEnd` 边界，也可直接构造高速 `CompressedList` 返回给 R。
  4. `Ten_things_slides.Rnw` (#1, #3, #4) 强调 `mcols(grl)` (outer) vs `mcols(unlist(grl))` (inner) 与 `extractList()` / `with.revmap=TRUE`，原版 `createExonsByReadClass` / `reducedRangesByGenes` 多处用了低效的 `splitAsList` + factor/string `match()` 重排，均可改为整数索引 + `relist(unlistData, PartitioningByEnd)` 或 C 级构造。

### 2.5 `S4Vectors 0.48.0` & `GenomeInfoDb 1.46.2`（`S4QuickOverview.Rmd`, `RleTricks.Rmd`, `HTS_core_package_stack.txt`）
- **核心发现**：
  - `S4QuickOverview.Rmd` 与 Bambu 运行剖面共同证明：**高频 S4 方法分派 (`start()`, `end()`, `strand()`, `[`, `mcols<-`) 放在按 read-class 或按 exon 循环体里，耗时远超实际整型计算本身**。
  - 原版 Bambu 对 `GRangesList` 的 `as.character(strand(grl))`、`unlist(..., use.names=TRUE)` 会强制生成数百万条全局 R `CHARSXP` 字符串哈希缓存条目（即 `data.table-benchmarking.Rmd` 也警告的 **R global string cache 污染与持久驻留**！）。一旦插入 R global string hash table，这些字符串即使对象被释放也永不缩小哈希表并拖慢所有后续比较。
- **对 `bambu_next` 的直接影响**：
  - 在中间阶段（S1/S2/S3 核心匹配与合并）全程使用**整型编码**（`int32_t seqid`, `int8_t strand_code` `1=+ / -1=- / 0=*`, `int32_t exon_start/end`），绝不产生数百万级 `paste0("rc_", ...)`、`paste(exons, collapse="+")` 或 `factor(seqnames)` 中间字符串。只在最终 S4/输出接口层才按需物化 ID 字符串。

### 2.6 `GenomicFeatures 1.62.0`（`GenomicFeatures.Rmd`）
- **核心发现**：
  - 官方文档明确指出：`TxDb` 的底层是 **SQLite 数据库**，`prepareAnnotations()` 调用 `exonsBy(txdb, by="tx")` 等函数需要多层 SQLite JOIN、Dplyr/AnnotationDbi 查询与 S4 转换；处理 `gencode.v46.annotation.gtf`（`1.6 GiB`）等大注释往往需要数十秒到上百秒。
  - 更致命的是，`bambu-extendAnnotations` 里的 `getMinimumEqClassByTx(extendedAnnotationRanges)` 在 S2 末端把 `245k` 参考转录本加上几百条新发现转录本合在一起，再做全量 `245k × 245k` `cutStartEndFromGrangesList()` + `findSpliceOverlapsQuick()`！同时 S3 的 `addEmptyRC()` 在 16 个样本内又各自重复执行 `processMinEquiRC(annotations)`！
- **对 `bambu_next` 的直接影响**：
  - **Node 1**：实现 C++ 流式 `prepareAnnotationsFromGTF` 解析器（保留与 `GenomicFeatures::makeTxDbFromGFF` + `bambu::prepareAnnotations` 完全一致的 exon ranking、gene_id 过滤、`txClassDescription`/`minEqClass` 计算契约），且直接在预处理对象上附带并**缓存不可变的参考 `minEqClass` 拓扑与 Flat-CSR 索引**；在 S2 扩展阶段**增量计算**仅受新转录本/新重叠基因影响的等价类，在 S3 所有样本直接复用一次结果，消除 `16 ×` 重复计算！

### 2.7 `SummarizedExperiment 1.40.0`（`SummarizedExperiment.Rmd`, `Extensions.Rmd`）
- **核心发现**：
  - 官方 `cbind()`、`rowRanges()`、`metadata()` 契约要求严格维护行名、列名与 `DataFrame` 类型对齐。
  - 原版 `bambu()` 在多样本 `combineCountSes()` 中重复传递/合并含庞大 `metadata(se)$readToTranscriptMap` 的样本级 `SE`，并在全流程末尾耗时数百秒做 `saveRDS`。
- **对 `bambu_next` 的直接影响**：
  - 节点流转与内部传递采用轻量化整数矩阵 + 结构化 `metadata`；仅在 `bambu()` 最终返回与 `writeBambuOutput()` 前按官方契约组装标准的 `RangedSummarizedExperiment` 对象，保证与下游 `transcriptToGeneExpression()`、`plotBambu()`、`DESeq2`/`DEXSeq` 100% 接口与内容一致。

### 2.8 `data.table 1.17.8` & `RcppArmadillo 15.2.6-1` & `xgboost 3.2.0.1`
- **`data.table-benchmarking.Rmd` 与 `datatable-reference-semantics.Rmd`**：
  - 强烈建议使用原位引用更新 `set()`（尤其循环内）与 `:=`、`setDT()`、`setkey()` / binary search join (`X[Y, on=...]`)，显式避开 `dplyr::full_join` 全表拷贝、`dplyr::group_by` 对 `list-column`（如 `eqClassById`）的高开销迭代、以及 `tidyr::separate_rows` / `tidyr::unnest` 的爆表展开。
  - `setDTthreads(1L)` 防范多核嵌套线程抢占。
- **`RcppArmadillo 15.2.6-1` (`src/em.cpp`)**：
  - 原版 `bambu.quantDT` 在 R 循环中对每个 `gene_grp_id` 调用 **3 次 `data.table::dcast()`** 来把稀疏 `equiRC` 转成稠密矩阵 `A, A_full, A_unique` 再喂给 `arma::mat` 版 `emWithL1`。
  - **Node 8 重构**：直接在 C++ 中按 `gene_grp_id` 一次排序分桶构建局部紧凑数组并直接运行 `emWithL1`（无需先建全库长 table 再 3 次 `dcast`，并严格保持与原 `em.cpp` 一致的浮点运算顺序与 `conv`/`minvalue`/`maxiter`）。
- **`xgboost 3.2.0.1`**：
  - 注意 `envs/bambu` 为 `xgboost 3.2.0.1`，默认加载 `read_class_ME.ubj`/`read_class_SE.ubj`/`model_spliceSitePrediction*.ubj` 模型。必须显式传 `nthread = 1L`（避免多线程超配与潜在 fork 竞争）并一次性对整块矩阵 `predict()`；同时调研其内部 `UBJ` 树节点评估以确认是否可兼顾纯 C++ 无锁评估与原生 `predict()` 双模式精确对齐（误差 `< 1e-12` 或严格 `identical`）。
