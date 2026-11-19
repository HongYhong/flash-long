# IsoQuant 分节点重构技术实施计划 (`PLAN`)

> 本文档在 [spec.md](spec.md)、[optimization_schemes.md](optimization_schemes.md) 与 [roadmap.md](roadmap.md) 确定主架构（**方案 B** + `sqanti3` Conda 环境）的基础上，给出每个重构节点的代码级模块划分、底层数据结构、Cython 内核边界、Bug 修复落点与验证命令模板。

---

## 1. 目标代码库骨架（`isoquant_refactored/`）

```text
isoquant_refactored/
├── README.md                              # 导航入口
├── REFACTORING_PLAN.md                    # 历史汇总与新文档体系指针
├── pyproject.toml                         # 打包规范 (setuptools + Cython, Python >=3.10)
├── setup.py                               # Cython Extensions 编译脚本 (gcc -O3, inplace)
├── docs/
│   ├── spec.md                            # 契约规格、Conda 环境、Benchmark 数据与口径
│   ├── optimization_schemes.md            # 低内存/高并发架构方案对比与包手册规范
│   ├── roadmap.md                         # M0~M10 与 G0~G4 里程碑
│   ├── plan.md                            # 本文档：Node 0 ~ Node 10 实施详案
│   └── task.md                            # T00 ~ T31 可追踪任务表
├── scripts/
│   └── benchmark/
│       ├── monitor_process_tree.sh        # 周期采样 /proc/<pid>/smaps_rollup 的 PSS/RSS 监控器
│       ├── run_tiered_benchmark.py        # A-B-B-A 交替基准测试执行器
│       └── compare_outputs.py             # 原版 vs 重构版 GTF/TSV/MTX/read_info 自动差分校验
├── src/
│   └── isoquant_next/
│       ├── __init__.py
│       ├── cli/                           # Node 0 & 9: isoquant-next, detect-barcodes, visualize 入口
│       ├── config/                        # Node 0: 不可变 dataclass 配置、run_manifest.json、策略表
│       ├── core/                          # Node 0 & 1: 日志/ExitCode/命名/ID池/二进制IO/WorkerPool/Arena
│       │   └── _cy_hash_pool.pyx          # 可选高速字节池/哈希辅助
│       ├── annotation/                    # Node 1: GTF/GFF3 流式校验、.iqidx 编译与 mmap 只读引擎、FASTA
│       ├── mapping/                       # Node 1: minimap2 / STARlong 子进程管道与 blake2b 缓存键
│       ├── barcodes/                      # Node 2: FastxFile 批读、2-bit SharedMemory、Cython SW & Myers UMI
│       │   └── _cy_barcode_umi.pyx        # Cython nogil Smith-Waterman + Bit-Parallel Edit Distance
│       ├── alignment/                     # Node 3: BAM 堆归并、Locus 边界分块、SoA 缓冲、Illumina 校正、ReadGroup
│       │   └── _cy_cigar_blocks.pyx       # Cython CIGAR block & polyA tail 修剪内核
│       ├── assignment/                    # Node 4: Cython 区间轮廓匹配、轻量化分配结果、精简多映射归约、IO
│       │   └── _cy_profile_assign.pyx     # Cython nogil 排序双指针/位图轮廓匹配与分类内核
│       ├── terminal/                      # Node 5: 稀疏直方图向量矩、零 XGBoost JSON 树评估、轻量 CAGE 索引
│       │   ├── _cy_hist_stats.pyx         # Cython 稀疏权重矩计算与峰值特征提取
│       │   └── _cy_tree_eval.pyx          # Cython L1-cache float32 gbtree 评估器
│       ├── discovery/                     # Node 6: 显式 ChromosomeDiscoveryState、IntronGraph、末端精修与过滤
│       ├── quantification/                # Node 7: 策略权重、区域计数、Polars 聚合、SciPy COO MTX 转换、辅助脚本
│       ├── fusion/                        # Node 8: 断点扫描、有界实例缓存、minimap2/Cython 校验器、静态区间索引
│       ├── pipeline/                      # Node 9: Stage Graph、Locus Streaming Worker、分片合并、原子 Resume
│       ├── visualization/                 # Node 10: OutputConfig (JSON 优先) 与 Agg 安全 Matplotlib 绘图
│       └── data/                          # 引用/复制 model_polya.json 与 model_tss.json
└── tests/                                 # Node 0~10 单元测试、Cython vs Python 核等价测试、G0~G4 回归套件
```

---

## 2. Phase 0：基线摸底与工程底座建立 (`M0_PRE`)

1. **构建包骨架与 Cython 原地编译测试**：
   - 编写最小 `pyproject.toml` 与 `setup.py`，在 `sqanti3` 环境下通过 `CC=gcc python setup.py build_ext --inplace` 编译并导入占位 `.pyx` 扩展，验证无头文件或权限阻碍。
2. **编写 `scripts/benchmark/monitor_process_tree.sh`**：
   - 借鉴 `long_bench/scripts/bambu_exon_block_parallel/monitor_process_tree.sh` 逻辑，按固定采样周期追踪测试根 PID 及其所有后代 PID 的 `/proc/<pid>/smaps_rollup` (`Pss`, `Rss`, `Pss_Anon`, `Pss_File`) 与精确运行时间。
3. **原版基线可行性与基线摸底体检（`G0` 前置）**：
   - 在 `sqanti3` 下分别实测原版 `isoquant.py` 哪些子模式可直接运行（预期：`simple_data` / `toy_data` / `ENCFF400BQQ` 不开启需 `xgboost` 的端模型预测或触发 `--complete_genedb` 时的行为），精确记录触发 `ImportError` (`xgboost` / `mappy` / `ssw` / `editdistance`) 的具体 CLI 选项。
   - 对原版可跑通模式采集 Tier 0 与 `ENCFF400BQQ` (`92 MiB`)、`ENCFF436GKZ` (`489 MiB`) 的真实基线耗时与 Peak PSS，作为后续对照基准表。

---

## 3. 分节点技术实施方案（Node 0 ~ Node 10）

### Node 0：配置契约、日志基础设施与安全并发底座 (`config/`, `core/`)
- **覆盖原模块**：`modes.py`, `error_codes.py`, `file_naming.py`, `file_utils.py`, `id_policy.py`, `common.py`（基础集）, `isoquant.py` 参数解析部分。
- **修复明确 Bug 与技术债**：
  1. 用不可变 `@dataclass(frozen=True, slots=True)` 替代向 `argparse.Namespace` 动态塞属性，提供 JSON Schema 版本化持久化 `run_manifest.json`（淘汰 `pickle` `.params`）。
  2. 修复 `parse_args(None)` 异常、修复 `get_path_to_program()` 仅比对字符串主版本号的缺陷（改用 `shutil.which` + 语义版本元组对比）。
  3. 统一唯一 `reverse_complement()`（`str.maketrans`，完整支持大小写及 IUPAC 简并碱基），清查并关闭所有裸 `open()`。
- **性能与并发防护重点**：
  - 统一 `WorkerPoolExecutor` 抽象，显式设置 `mp.get_context("forkserver")`（Linux）或 `"spawn"`，同时强制设置单线程环境变量 `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `NUMEXPR_NUM_THREADS=1`，杜绝 BLAS/OpenMP 与多进程抢核导致的吞吐骤降与死锁。
- **验收标准**：全部 CLI 默认策略表（ON/PacBio/Single-Cell/CAPTrap/DirectRNA）、`IsoQuantExitCode`、输出命名规则 100% 对齐；`run_manifest.json` 往返校验通过。

### Node 1：流式 GTF 编译器、`.iqidx` 零拷贝 `mmap` 索引、序列化与比对封装 (`core/`, `annotation/`, `mapping/`)
- **覆盖原模块**：`serialization.py`, `string_pools.py`, `gtf2db.py`, `gene_info.py`, `input_data_storage.py`, `read_mapper.py`。
- **修复明确 Bug**：
  1. `serialization.py`：`read_dict()` 错调 `read_int()` 读带符号整数的问题，统一使用 `struct.Struct` 定长编解码并补全边界单测。
  2. `gtf2db.py`：修复 `inner_ext.lower() == '.gff3'`（补全前缀点号）；流式逐行校验写入取代 `readlines()` + `+=` 字符串拼接，支持属性字段内嵌空格（正则状态机/合规词法切分）。
  3. `read_mapper.py`：用 `hashlib.blake2b` 替代内置随机 `hash()` 做比对缓存指纹；用 `contextlib.ExitStack` 与明确等待顺序重写 `minimap2 | samtools sort` 管道，消除 `SIGPIPE` 与句柄泄漏；大染色体 自动 `.csi` 索引回退。
- **核心低内存 & 高速架构设计**：
  1. **`.iqidx` 编译型 `mmap` 零拷贝注释索引（消除 `gffutils` SQLite N+1 与 `GeneInfo` 内存膨胀）**：
     - 单遍解析 GTF/GFF3（或读取既有 `.db` 后一次性导出），对每个基因内部采用**单遍聚合子树 + 排序扫描线算法**消除原版 `GeneInfo` 的 5 次 `db.children()` 查询、`set_gene_attributes()` $O(T \cdot E)$ 与 `set_feature_properties()` $O(N^2)$ 双重循环。
     - 输出紧凑定长二进制索引目录/文件 `annotation.iqidx`：
       - `chrom_table`（固定记录：chrom_id, name_offset, gene_start, gene_end）
       - `gene_table`（升序 `(start, end, strand, chrom_id, tx_start_idx, tx_end_idx, exon_start_idx, exon_end_idx, intron_start_idx, intron_end_idx, name_id)`）
       - `tx_table` / `known_exons` / `known_introns` / `split_exons` 全部以连续 `int32` CSR 数组紧凑存放；
       - 字符串表（Gene ID, Transcript ID, Gene Name, Biotype）与 `StringPool` 紧致编码对齐，按 `int32` 下标 O(1) 反查。
     - 运行期各 worker 通过 Python `mmap` + `np.frombuffer(..., dtype=...)` 直接获得零拷贝只读视图：1.6G `gencode.v46.annotation.gtf` 在多达 32 个 worker 并发时物理 `PSS` 仅占一份 OS 页缓存！
  2. **FASTA 访问统一**：全部收口至 worker 内惰性只读 `pysam.FastaFile`（原生直读 `.fa` 及 `.fa.gz` bgzip `.fai+.gzi`），消除 `pyfaidx` 与 `pysam.FastaFile` 双库混用与跨进程非安全持有。
- **验收标准**：对比 `test_gene_info.py` / `test_serialization.py` / `test_string_pools.py` / `gencode.v46.annotation.gtf` 解析输出，基因/转录本/外显子/内含子坐标及特征属性 100% 一致，统计建库耗时与多核 PSS。

### Node 2：FastxFile 流式 Barcode Calling 与 Cython Bit-Parallel UMI 去重 (`barcodes/`)
- **覆盖原模块**：`isoquant_detect_barcodes.py`, `barcode_calling/*`, `table_splitter.py`。
- **修复明确 Bug**：
  - 修复 `isoquant.py` 中 `call_barcodes()` 的 `with ProcessPoolExecutor` 误缩进循环体内导致 N 倍重复提交的致命缺陷；
  - 清除不可达的 `load_h5_barcodes_bit` (`h5py`) 死代码；修复 `BarcodeCaller.process` 与 `load_barcodes`、`table_splitter.py` 文件句柄泄漏；
  - `SharedMemoryArray2BitKmerIndexer` 强制全部位掩码使用 `np.uint64`（满足 NumPy 2.0 NEP 50），k-mer 表大小使用 `1 << (2 * k)`，并在 Python 3.13+ 子进程挂载传 `track=False` 加显式清理句柄。
- **核心性能架构**：
  1. 使用 **`pysam.FastxFile`** 替代重量级 `Bio.SeqIO.parse`，按批（如 10,000 reads/chunk）连续投喂工作进程；采用 `str.translate` 替代手工字典反向互补。
  2. **`_cy_barcode_umi.pyx` (`nogil`)**：
     - 在核对 `barcode_calling/common.py` `AlignmentMgr` 确切评分矩阵后，用 Cython 实现可复用缓冲区的 Smith-Waterman 短模式匹配（替代未安装的 `ssw-py`，单核批量打分无 Python 对象开销）；
     - 实现 Myers 64-bit 位并行编辑距离 + 前置长度差 $O(1)$ 剪枝与短串等长 Hamming 判定（替代未安装的 `editdistance`），消除 Python 回调盒装开销。
- **验收标准**：对照 `test_barcode_*.py`、`test_kmer_indexer.py`、`test_shared_mem_index.py`、`test_umi_filtering.py`、`test_universal_extraction.py` 全量通过。

### Node 3：BAM 在线堆归并、Locus 边界流式分片与外显子校正 (`alignment/`)
- **覆盖原模块**：`alignment_info.py`, `alignment_processor.py`, `exon_corrector.py`, `illumina_exon_corrector.py`, `read_groups.py`。
- **修复明确 Bug**：
  - 修复 `illumina_exon_corrector.py` 中 `type(i) != "int"` 与 `logger.debug("old:", introns)` 崩溃漏洞；统一 `with` 管理短读 BAM；整合 `read_groups.py` 重复解析分支。
- **核心低内存架构**：
  1. **`BAMStreamMerger` + `LocusChunker`**：按参考染色体有序流式扫描 BAM，仅在**基因组上读段与基因注释均无跨越重叠的空隙处（Locus Gap）**触发 Chunk 提交；极端超长覆盖密集区支持安全最大读段水位上限切片。
  2. **`_cy_cigar_blocks.pyx`**：从 pysam `cigartuples` 直接在 Cython 填充连续 `int32` 块缓冲与 polyA 裁剪坐标，不生成临时冗余小元组；使用 `SoA`（read_id, start, end, strand, flags, blocks_offset, blocks_len, polya_pos...）存储单一活跃 Locus 窗口内读段。
- **验收标准**：对齐 `test_alignment_info.py`, `test_illumina_exon_corrector.py`, `test_read_groups.py`。

### Node 4：Cython `nogil` 轮廓匹配、Read-to-Isoform 分配与零对象 IPC 归约 (`assignment/`)
- **覆盖原模块**：`isoform_assignment.py`, `long_read_profiles.py`, `junction_comparator.py`, `long_read_assigner.py`, `multimap_resolver.py`, `processed_read_manager.py`, `assignment_loader.py`, `assignment_io.py`, `assignment_aggregator.py`。
- **核心低内存 & 高速架构**：
  1. **`_cy_profile_assign.pyx` (`nogil`)**：直接操作只读 `mmap` `.iqidx` 的已知外显子/内含子 `int32` 视图，通过二分起始定位 + 双指针扫掠生成紧凑 `int8` 匹配剖面，彻底淘汰 `OverlappingFeaturesProfileConstructor` 每条 read 构造 3 个 `MappedReadProfile` Python `list[int]` 的分配风暴。
  2. **职责拆解与紧致分配记录**：将 830 行 `long_read_assigner.py` 拆为 `CandidateFilter`、`SingleExonAssigner`、`MultiExonAssigner` 与 `InconsistencyClassifier`，严格保持全部 SQANTI 分类与 mismatch 事件代码等价。`ReadAssignment` 在工作阶段尽量保持紧凑数值元组/`@dataclass(slots=True)` + `StringPool` ID。
  3. **取消大型对象 IPC 传递与双重 `.save` I/O**：若单核/按基因座闭环完成全处理，分配结果直接喂给同 Locus 的下游 `Node 5 / 6 / 7` 消费后就地释放；仅在显式需要 `--read_assignments` `.save` 或跨样本联合建模时写出紧凑分片二进制流（废弃旧版 `ProcessedReadsManager*` 把整染色体百万对象列表 pickle 塞给主进程管道的做法）。Multimap Resolver 通过轻量定长数组 `(read_hash_or_id, gene_idx, match_score, assignment_type)` 在最终归约期一次解决。
- **验收标准**：对照 `test_long_read_profile.py`, `test_junction_comparator.py`, `test_long_read_assigner.py`, `test_output_formats.py`, `test_assess_assignment_quality.py` 及 `ENCFF400BQQ` 逐 read `read_assignments.tsv` 100% 匹配。

### Node 5：稀疏直方图向量矩与零 XGBoost L1-Cache 树推理引擎 (`terminal/`)
- **覆盖原模块**：`terminal_prediction/*` (`polya_finder.py`, `polya_verification.py`, `terminal_peaks.py`, `terminal_counter.py`, `cage_finder.py`)。
- **核心低内存 & 高速架构**：
  1. **合并重复模块并废弃列表展开**：统一 `terminal_peaks.py` 与 `terminal_counter.py` 两套重复直方图逻辑；`_cy_hist_stats.pyx` 直接利用排序后的 `(pos_int32, count_int32)` 数组一次扫描算完 8 个特征 (`var`, `skew`, `peak_count`, `peak_width`, `entropy`, `mean_height`, `peak_heights`, `relative_height`)，**完全消灭 `data.extend([pos] * count)` OOM 风险**。
  2. **`_cy_tree_eval.pyx`（替代 `XGBClassifier` + `pd.DataFrame`）**：
     - 实测验证 `model_polya.json` 与 `model_tss.json` 均为标准 `gbtree`、`binary:logistic`、`split_type={0}`（共 417 树、<=2,637 节点）。
     - 一次解析写入紧致 `float32` / `int32` 数组（`< 60 KB`），批处理预测直接计算 IEEE-754 `float32` `<` 判断与 `logit(base_score)` 累加 sigmoid，`nogil` 运行且彻底摆脱 `xgboost` 包与 OpenMP fork 隐患。
  3. **矢量化辅助接口**：`TerminalCounter._write_grouped` 用 `np.searchsorted` 替换 $O(N_{\text{peaks}} \cdot N_{\text{groups}})$ `iterrows()`；`CagePeakFinder` 改为轻量纯 Python 流式读取 BED 写入静态有序区间数组，**移除未声明包 `pybedtools` 依赖**。
- **验收标准**：对照 `test_polya_prediction.py`, `test_terminal_peaks.py`, `test_polya_cage_finder.py`，并与原始模型离线推理结果对比误差 `< 1e-6`。

### Node 6：无状态污染的 IntronGraph 与局部转录本模型重建 (`discovery/`)
- **覆盖原模块**：`model_construction/*`（`intron_graph.py`, `intron_path.py`, `fl_graph_constructor.py`, `assignment_based_constructor.py`, `end_processor.py`, `model_filter.py`, `read_assigner.py`, `model_read_counter.py`, `transcript_printer.py`, `transcript_splice_site_corrector.py`, `transcript_to_gene_joiner.py`, `context.py`, `model_construction.py`）。
- **重点修复与优化**：
  1. **消除 `ModelStore` 类级静态可变状态** (`detected_known_isoforms`, `extended_transcript_ids`)，改为每个基因座/染色体显式传递 `ChromosomeDiscoveryState` 对象，确保跨样本/多轮调用及并行安全。
  2. 将 Node 5 统一峰值特征与树推理无缝接入 `end_processor.py` 的 5'/3' alt-end 判定，减少循环中直方图重建开销。
  3. 明确五阶段流式契约，并在基因座结束即写出 GTF / `novel_vs_known.SQANTI-like.tsv` 碎片，释放该基因座图结构内存。
- **验收标准**：对照 `test_intron_graph*.py`, `test_graph_alt_ends.py`, `test_transcript_splice_site_corrector.py` 与真实样本重建 `transcript_models.gtf`。

### Node 7：Polars / SciPy Out-of-Core 定量计数与 MTX 导出 (`quantification/`)
- **覆盖原模块**：`long_read_counter.py`, `convert_grouped_counts.py`, `stats.py`, `scripts/*`。
- **重点优化**：
  1. 将 `long_read_counter.py` (1,129 行) 拆解为策略权重、`AssignedFeatureCounter`、`JointExonCounter`、`ExonSpliceSiteCounter`、`IntronRetentionCounter` 与 `CompositeCounter`；底层累加全面使用连续 `float64/int64` NumPy 向量，在最后汇总输出基因/转录本矩阵。
  2. `convert_grouped_counts.py`：彻底重写 MTX 导出，采用 Polars / NumPy 列映射将 `(feature, group, count)` 转成 `scipy.sparse.coo_matrix` 并由 `scipy.io.mmwrite` 原生 Cython/SciPy 写出；修正所有未关闭句柄。
  3. 多分组 / 亿级单细胞 `grouped_counts` 合并使用 `Polars 0.20` lazy streaming scan，将聚合内存保持有界。
- **验收标准**：全部定量计数单测通过；Tier 0/1/2/3 与原版各种 `--counting_strategy`（`fsm_only`, `unique_only`, `default` 等）计数 100% 等价，MTX 稀疏内容完全一致。

### Node 8：无全局污染融合基因检测与多后端重比对 (`fusion/`)
- **覆盖原模块**：`fusion_detector.py`, `fusion_validator.py`, `fusion_metadata.py`, `genomic_interval_index.py`。
- **重点修复与优化**：
  1. 消除 `fusion_detector.py` 4 个模块级全局可变字典（`_CONTEXT_CACHE`, `_GENE_ASSIGNMENT_CACHE`, `_ALIGNER_MAP_CACHE`, `_CIGAR_CACHE`），替换为绑定到单次 `FusionDetector` 实例的有界 `lru_cache`，并废除外部强行注入 `fd.bam_path = ...`。
  2. `genomic_interval_index.py` 由支持静态排序数组的 `bisect`/`np.searchsorted` 索引替代高内存开销的纯 Python `IntervalTree`（或保持一致封装 API）。
  3. **`mappy` 无依赖适配**：针对局部软剪切断点校验，支持“批收集候选融合序列 -> 单次管道调用 `sqanti3/bin/minimap2`”与 Cython 局部 Smith-Waterman 双模式（存在 `mappy` 时亦支持原 API 校验对照），清理 `q_st/qstart` 猜属性分支。
- **验收标准**：`test_fusion_*.py` & `test_genomic_interval_index.py` 全通过。

### Node 9 & Node 10：全流程流式编排、原子 Resume 与现代可视化 (`pipeline/`, `visualization/`)
- **Node 9 实施要点**：
  - 将 `isoquant.py` + `dataset_processor.py` + `parallel_workers.py` 重写为清晰 Stage Graph：
    `Preflight -> Barcode Calling (可选) -> .iqidx 准备 -> Mapping (可选) -> 染色体/Locus 分片并发流水线 (Node 3->4->5->6->7 单遍闭环，最小化中间落盘) -> 多样本合并 & MTX -> Fusion (可选) -> 可视化 (可选)`。
  - 实现基于临时写入 + `os.replace()` 的原子级 Stage Manifest (`--resume` 状态机)，清理陈旧锁文件竞赛隐患，并修复 `filter_umis` 缺失 worker 日志初始化问题。
- **Node 10 实施要点**：
  - `OutputConfig` 优先读取 Node 0 的 `run_manifest.json`，仅对历史运行目录只读降级兼容 `.params`；
  - 统一 `matplotlib.use("Agg")` 无头绘制并在单基因图绘制后强制 `plt.close(fig)`，修复 Pandas 3.0 CoW 下 `gene_model.py` 等文件可能存在的链式赋值告警。
  - 运行全量基准测试矩阵（G0~G4），并在 `basicpy` (`NumPy 2.4.6` + `Pandas 3.0.2`) 校验核心矩阵/表处理层的前瞻兼容性。

---

## 4. 关键对照策略说明（当原版因缺失包受限时）

在 `sqanti3` 环境实测已确认缺少 `xgboost`、`mappy`、`ssw-py`、`editdistance`。为保证验收客观严谨：
1. **可直接与原版全量闭环对照的基准**：常规 `--bam` + `--gtf`/`--genedb` bulk/全长定量及模型重建流程（需先实测原版 `isoquant.py` 在 `--no_polya` 或未开启融合/barcode 模式下是否可在 `sqanti3` 零缺包崩溃跑完）。凡原版能跑通的路径，一律进行 `isoquant_py vs isoquant-next` 同环境 `A-B-B-A` 实测壁钟时间与进程树 Peak PSS 严格对比。
2. **受缺失包阻断路径的模块级 Oracle 验证**：针对依赖 `xgboost`（Node 5 PolyA/TSS 树模型）、`ssw`/`editdistance`（Node 2 单细胞 barcode/UMI）、`mappy`（Node 8 融合校验），通过：
   - 数学形式完全等价检验（如 JSON 树 `float32` 逐样本对比参考数学公式、`editdistance` vs Myers 暴力 DP 对比）；
   - `SIRV_simulation` 已知真值转录本集上的灵敏度（Sensitivity / Recall）、精确度（Precision）和定量相关性校验；
   - 保持清晰接口预留可选验证开关（若在含这 4 个包的节点或容器测试，自动开启与原生库逐条对照）。
