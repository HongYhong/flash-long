# `bambu_next` 低内存与高并发架构优化方案对比 (`Optimization Schemes`)

> **目标**：彻底解决原版 Bambu 及既往 `G8–G17` 重构版本在多样本全流程中“单节点优化被其他阶段吞没、8 核 `fork` CoW 撑爆 `250+ GiB` PSS、以及百万级 S4 对象与 `dplyr/tidyr/tapply` 热路径”导致的停滞，设计兼顾**极端内存压缩**、**真并发加速**与**Bioconductor R 生态 100% 兼容**的下一代架构。

---

## 1. 原版 Bambu & 既往 `G8–G17` 真实瓶颈解剖（基于实测剖面）

以下所有数字均直接来自 `long_bench` 在 `envs/bambu` 上的真实测量记录（[`g8_20261104`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/results/bambu_refactor/g8_20261104/REPORT.md)、[`g12_20261110`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/results/bambu_lowcopy/g12_20261110/analysis/full49_prod45_g14_v1/paired/paired_summary.tsv)、[`g15_20261113`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/results/bambu_s3core/g15_20261113/analysis/full49_prod45_g17_v1/paired/run_medians_spread.tsv)）：

### 1.1 `human_all_d1q1`（16 Human BAM, 8 Workers）与 `ENCFF089AMF` (`16.42M` alns) 宏观瓶颈

| 阶段 | 实测耗时（`human_all` 中位数） | 实测 Peak Tree PSS | 既往尝试 (`G8/G14/G17`) 为何失效 | 根因源码与数据结构陷阱 |
| :--- | :--- | :--- | :--- | :--- |
| **`S1` 单样本构建与评分** | `1495.6s` (native) / `1503.6s` (`G8`) | `58.0–58.9 GiB` | 多 BAM 模式为防 8×8 超配降级回 `sample8 × reader1`，等同原地跑原生 `readGAlignments`；仅优化 `createReadTable` 的 `CRT` 全局几何均值仅 `1.0044×` | 1. 全量加载 16M+ 比对入 `GAlignments -> grglist()` (`ENCFF089AMF` 原生 `338.82s`)。<br>2. `unlistIntrons` 构建 **`57.68M` `GRanges`**；`BSgenome::getSeq` + `dplyr` + `updateStrandScoreByRead` (`82s`)。<br>3. `createReadTable` (`107.85s`) + `constructUnspliced` (`41.78s`) + `assignGeneIds` (`24.52s`) + `io.save_read_class_cache` (`557.4s` 累计；单样本 `648 MB` RDS 含全量 `readIds` 列表)。 |
| **`S2` 跨样本合并与扩展** | `1545.2s` (native) / `1415.5s` (`G8`) | `33.6–38.8 GiB` | `G8`–`G17` 完全未优化 `S2` 算法结构 | 1. `combineSplicedTranscriptModels` (`824.3s`) 执行 16 次 `dplyr::full_join` 持续复制增宽后缀列。<br>2. **`readCountWeightedMedian` 物理展开 `rep(val, times = readCount)`** 暴占内存。<br>3. `combineUnsplicedTranscriptModels` (`289.6s`) 用 `+` 拼字符串后调 `tidyr::separate_rows()`。<br>4. `isore.extendAnnotations` (`333.8s`) 在 `combineWithAnnotations` 内对 `245k × 245k` 全部扩展注释重新 `cutStartEnd` + `findSpliceOverlapsQuick` 算 `minEqClass`！ |
| **`S3` 距离匹配与 EM 定量** | `1422.2s` (native) / `1388–1414s` (`G8/G17`) | **`250.48–268.39 GiB`（全程最高峰！）** | `G14 lowcopy` / `G17 s3core_v3` 仅降 PSS `-5.1% ~ +4.9%`，耗时 `+0.8% ~ +0.01%` | 1. **8 核 `fork` CoW 灾难**：每进程持有巨型 S4 对象且 S4 临时分配频繁触发 R GC 标记写入，`254 GiB` PSS / `294–311 GiB` RSS！<br>2. **`calculateDistToAnnotation` (`2596s` 累计)**：在 `ENCFF089AMF` 展开 `5.60M` 重叠对的 **`51.72M` 个外显子 `GRangesList` 双份拷贝** + `myGaps` + `pintersect` + **`tapply(..., myId, sum)` 分配 5.6M 字符串分组**（单步累计 `404.91 GiB` 瞬时分配）。<br>3. **`genEquiRCs` (`886.5s` 累计)**：`as.data.frame` + `dplyr::group_by(list_col)` + **16 个样本独立重算 `processMinEquiRC(annotations)` (`tidyr::unnest` + `full_join`)**。<br>4. **`bambu.quantDT` (`997.3s` 累计)**：每个 `gene_grp_id` 调 3 次 `data.table::dcast()` 构建稠密矩阵喂给 `emWithL1`。 |
| **`S4` 基因与文件导出** | `55–69s` (clean) / `~1144s` (full RDS) | `33.7–42.3 GiB` | 未触动 | 1. `writeToGTF` 将全量 `GRangesList` 展成表再用 R `paste0` 逐行拼字符串。<br>2. `transcriptToGeneExpression` 用 `melt` + 双 `dcast` + `table`/`relist`。 |

---

## 2. 三套候选重构方案对比

| 维度 | 方案 A：纯 R 向量化微调 (`data.table` 替代 `dplyr/tidyr` + `GNCList`) | **方案 B（主推）：C++17 + `Rhtslib` + 扁平只读 CSR/NCList 引擎 + Drop-in S4 Adapter** | 方案 C：独立 Rust / 外部 CLI 二进制重构并由 R 解析文本回包 |
| :--- | :--- | :--- | :--- |
| **核心技术栈** | `data.table` `:=`/`set()`, `GenomicRanges::GNCList`, 保持 R `MulticoreParam` 与 `GRangesList` | **`Rhtslib` (`sam.h`/`faidx.h`) + C++17 `std::thread`/OpenMP + Flat CSR (`int32`/`uint32`) + `IRanges/S4Vectors` C 接口 + `RcppArmadillo` + `data.table` + 兼容 `SummarizedExperiment` 组装** | 外部独立可执行程序 + HDF5/GTF 落盘后再读回 R 构造 `SE` |
| **消除 `GAlignments`/57.7M 展开** | ❌ BAM 仍需 `readGAlignments -> grglist()` | ✅ C++ 直接解码 `bam1_t` CIGAR 到 `int32` CSR 缓冲，零 S4 中间对象 | ✅ |
| **消除 S3 `51.7M` 展开与 `tapply`** | ❌ `findSpliceOverlapsByDist` 若保留 `GRangesList` 依旧要多次索引/向量拆分 | ✅ C++ 原位双指针扫描 `(query_i, subject_j)`，零中间区间对象、零 `tapply` | ✅ |
| **解决 8 核 `fork` `250 GiB` CoW** | ❌ R GC 每次标记仍会分裂 `annotations`/`readClass` 页 | ✅ 只读扁平 CSR 索引与样本块留存在 C++堆/只读 `mmap` 区，**`std::thread` 零触碰 R `SEXP`**，多线程共享内存不分裂 | ✅ |
| **R 接口 & `Bioconductor` 语义保真** | ✅ 极高，但峰值改善天花板低（预计耗时仅 `1.5x~2x`、PSS 仍超 `100 GiB`） | ✅ **完全兼容**：`envs/bambu` 零新包即 compile & load；原生返回 `RangedSummarizedExperiment`、完整支持 `BSgenome`/`TxDb`/`trainBambu` 全部参数 | ❌ 外部二进制无法直接接收/返回内存级 `TxDb`/`BSgenome`/`trainBambu` R 对象，IPC 序列化损失大 |
| **预计收益（设计目标，待 Gate 实测）** | 耗时 `1.5x–2x`；PSS `-25%~-40%` | **全链路耗时预计 `4x–12x+` 提速；`human_all_d1q1` 峰值 PSS 预计压缩 `75%~90%`（从 `250 GiB` 级压至 `~20–50 GiB` 级）** | 同等性能但打破 R 包契约 |

---

## 3. 方案 B 五大架构支柱详解

1. **支柱一：单遍 C++ `GTF -> Flat-CSR`（`.bbidx`）+ 预计算可复用 `minEqClass` 拓扑**
   - 直接流式解析 GTF（或从 `TxDb`/`prepared_annotation GRangesList` 零转换提取 C 指针），生成连续整型数组：`chrom_id`, `strand_code`, `tx_exon_ptr[]`, `exon_start[]`, `exon_end[]`, `intron_start[]`, `intron_end[]`, `gene_id[]` 及按 `(seqid, strand)` 排序的边界包围盒索引。
   - **消除 `16 × processMinEquiRC`**：参考注释的最小等价类矩阵 `ref_min_eq_csr` 在索引准备期一次性算出并固化在 C++ 结构中；S2 `combineWithAnnotations` 仅对有新发现转录本重叠的基因连通域执行局部增量更新！
2. **支柱二：`Rhtslib` 并发流式 BAM 解析 + 按基因座/区块流式聚合（`S1` 零 `GRangesList`）**
   - 直接调 `Rhtslib`（`sam_hdr_read` + `hts_itr_queryi` 按 `.bai`/染色体切片，`hts_tpool` / worker 级 `std::thread`）解析 `bam1_t`：一次遍历即同时完成 secondary 过滤、1bp 边界清理、`CIGAR` exon-block 切分、splice intron 抽提与 `QS`/`strand`/`XS` 解析。
   - Spliced read 在读取结束后立刻通过 `(chr, strand, intron_starts..., intron_ends...)` **64-bit XXH3/FNV1a + 精确校验多键哈希归并**为唯一 read class 拓扑并流式维护 5'/3' 端点分位数/边界范围（完全取代 `unlistIntrons` 的 `57.68M` 对象 + `createReadTable` 的 `dplyr::group_by` + `unstrsplit(splitAsList)`！）；若 `trackReads=TRUE`，仅保留紧凑 `uint32` `read_id` 偏移而非巨型字符串嵌套列表！
   - 使用 `htslib/faidx.h` 批量提取唯一 junction 上下各 2bp `dinucleotide` 完成 `spliceMotif` 与链校正投票；Unspliced reads 直接按 `(seqid, strand)` 双指针区间聚类构建。
3. **支柱三：单写多读零 CoW 并发模型（告别 `MulticoreParam` `fork` R 堆拷贝）**
   - **`bambu_next` 核心原则**：多样本与多 BAM 并发不再把整个 `annotations` + 全部 16 个 S4 `RangedSummarizedExperiment` 拉入 8 个互不共享的 R `fork` 子进程执行重计算。
   - 改为：**主 R 进程维护唯一轻量句柄**，C++ 引擎通过 `std::thread` + `libhts` 对各样本 / 染色体分片并发运算（只读写各自的连续扁平 C++ 缓冲区与共用 `const` 只读 `.bbidx` CSR 数组），**线程内零分配 `SEXP`、零触发 R GC**。16 份样本同时跑定量时，公共参考注释永远只有 **1 份物理内存页**！
4. **支柱四：`S3` 零展开 C++ 双指针匹配（`findSpliceOverlapsByDist`）+ 融合无 `dcast` `emWithL1`**
   - 为参考与扩展后转录本构建按 `(seqid, strand)` 排列的 **Flat-CSR 区间树 / Sweep-Line 索引**；每个 `readClass` 查找到候选转录本后，在 C++ 中用原始 `int32` 指针单遍算出全部距离与兼容性指标，直接写入紧凑列向量返回/缓存，彻底拔除 `5.6M` `GRangesList` 复制、`myGaps()`、`pintersect` 与 `tapply`。
   - 直接由兼容关系数组构建 per-`gene_grp_id` CSR 计数簇并调用多线程 `C++ emWithL1`（移除 R 循环内 `3 × dcast` 和 `dplyr` `list-column` 分组），同时保持与原生 [`em.cpp`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/src/em.cpp) 完全一致的数值迭代逻辑与收敛阈值。
5. **支柱五：`S2` $O(N \log N)$ K-way 合并、精准分位数（无 `rep` 展开）与 `S4` C++ 流式导出**
   - `combineSplicedTranscriptModels`：由 C++ K-way 归并哈希表一次性汇总所有样本的特征支持、边界起止及权重（`O(U)` Exact Weighted Median 直接在 `(value, weight)` 排序数组扫前缀和，彻底删去 `rep(x, times = count)` 展开）；
   - `combineUnsplicedTranscriptModels`：以整数归并聚类替代 `paste(row_id, collapse="+")` + `tidyr::separate_rows`；
   - `writeToGTF` & `writeBambuOutput`：直接在 C++ 中以 buffered I/O（如 4 MiB 写缓冲）流式生成标准 GTF 与 counts/CPM TSV；仅在 R 层返回要求完整 `RangedSummarizedExperiment` 对象时，通过 `IRanges`/`S4Vectors` C stub（`_IRanges_stubs.c` / `_S4Vectors_stubs.c`）低拷贝组装 S4 壳返回。
