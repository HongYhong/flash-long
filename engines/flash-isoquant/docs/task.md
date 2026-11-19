# IsoQuant 重构可执行任务清单 (`TASKS`)

> 状态定义：`TODO`（待开始）｜`IN_PROGRESS`（进行中）｜`DONE`（已完成并验证）｜`BLOCKED`（受阻塞）
> 优先级：`P0`（关键路径阻塞项）｜`P1`（核心功能与主要优化）｜`P2`（增强/辅助/兼容）

---

## Phase 0：规划、工程骨架与原版基线摸底 (`M0_PRE`)

| ID | 阶段/节点 | 任务名称 | 状态 | 优先级 | 依赖 | 交付物与验证标准 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`T00`** | Pre | 多文档规划体系（`README` / `spec` / `optimization_schemes` / `roadmap` / `plan` / `task`）与环境/数据调研 | `DONE` | `P0` | — | `isoquant_refactored/README.md` 与 `docs/*.md` 落盘并经人工确认 |
| **`T01`** | Phase 0 | 新建 `isoquant_refactored/` 包结构、`pyproject.toml`、`setup.py`（Cython `sqanti3` inplace 编译验证） | `DONE` | `P0` | `T00` | `setup.py build_ext --inplace` 与 `tests/test_smoke_build.py` 通过 |
| **`T02`** | Phase 0 | 部署 `scripts/benchmark/monitor_process_tree.sh`、`compare_outputs.py` 与 `run_tiered_benchmark.py` | `DONE` | `P0` | `T01` | 已实测采集 `/proc/<pid>/smaps_rollup` Tree PSS/RSS 与 `tests/test_benchmark_tools.py` 通过 |
| **`T03`** | Phase 0 (G0) | 在 `sqanti3` 实测原版 `isoquant.py` 可跑通路径 vs 缺失 4 包阻塞点，采集 Tier 0 & `ENCFF400BQQ` 基线指标 | `DONE` | `P0` | `T02` | 已归档至 `docs/g0_baseline.md` |

---

## Phase 1：基础底座与零拷贝共享注释系统 (`Node 0 ~ Node 1`)

| ID | 阶段/节点 | 任务名称 | 状态 | 优先级 | 依赖 | 交付物与验证标准 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`T04`** | Node 0 | 重写 `config/` 强类型 `@dataclass`、`run_manifest.json`、策略默认参数表与 `IsoQuantExitCode` | `DONE` | `P0` | `T01` | `tests/test_node0_config_core.py` 全通过 |
| **`T05`** | Node 0 | 重写 `core/` 基础工具（唯一 `reverse_complement`、`shutil.which` 元组版本检查、`WorkerPoolExecutor` 强制单线程 OpenMP 防护） | `DONE` | `P0` | `T04` | `tests/test_node0_config_core.py` 全通过 |
| **`T06`** | Node 1 | 修复二进制序列化（`read_dict` 有符号 int 与 `struct`）与紧凑 `StringPool` | `DONE` | `P0` | `T05` | `tests/test_node1_annotation_io.py` 全通过 |
| **`T07`** | Node 1 | 流式 GTF/GFF3 校验器（修复 `.gff3` 与属性内空格，淘汰 `readlines()`）+ `.db` 兼容导入 | `DONE` | `P0` | `T06` | `gencode.v46` (`1.6G`) 流式去重校验 `12.497s` 通过；`tests/test_node1_annotation_io.py` 全通过 |
| **`T08`** | Node 1 | 设计与实现 `GTF/GFF3 -> .iqidx` 零拷贝 `mmap` CSR 索引，合并 5x `children()` 与 $O(N^2)$ 扫描线替换 | `DONE` | `P0` | `T06, T07` | 与 `GeneInfo` Oracle & `gencode.v46` (`63,086` 基因/`254,070` 转录本/`1,668,627` 外显子) 100% 一致；`v46` `.iqidx` `141.24 MiB`，8 核扫描 `4.085s`，全页 `mmap` 8 核共享 `PSS` 仅 `304,647 KiB`（文件页 `1.01x` 共享） |
| **`T09`** | Node 1 | 统一参考基因组 FASTA（worker 内 `pysam.FastaFile` 支持 `.fa`/bgzip）与 `read_mapper.py`（`blake2b` + `ExitStack` + `.csi`） | `DONE` | `P1` | `T05` | `tests/test_node1_annotation_io.py` 全通过 |

---

## Phase 2：Barcode/UMI 与流式 Locus 比对分块 (`Node 2 ~ Node 3`)

| ID | 阶段/节点 | 任务名称 | 状态 | 优先级 | 依赖 | 交付物与验证标准 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`T10`** | Node 2 | `pysam.FastxFile` 批流式读 FASTA/Q 替换 `Bio.SeqIO` + `SharedMemory` (`np.uint64` + 3.13 `track=False`) + `table_splitter` 修复 | `DONE` | `P1` | `T06` | `tests/test_node2_barcodes_umi.py` 与 `basicpy` NEP 50 `uint64` 校验全通过 |
| **`T11`** | Node 2 | 编写 `_cy_barcode_umi.pyx`（`nogil` Smith-Waterman 替代 `ssw` + Myers bit-parallel 替代 `editdistance`），重写 Caller 与 `umi_filtering` | `DONE` | `P1` | `T10` | 25/25 单元测试与 `isoquant_lib` Oracle 逐序列输出一致；`/dev/shm` 零泄漏 |
| **`T12`** | Node 3 | 实现 `BAMStreamMerger`、`LocusChunker`（基因组间隙与超密区安全水位）与 `_cy_cigar_blocks.pyx`（SoA 块缓冲） | `DONE` | `P0` | `T08, T09` | `_cy_cigar_blocks` 编译通过；`chr9.4M` 与 `ENCFF400BQQ` `chr22` `.iqidx` 流式分块无漏读；`AlignmentInfo` Oracle 100% 一致 |
| **`T13`** | Node 3 | 重写 `exon_corrector`、`illumina_exon_corrector`（修复 `type(i) != "int"` 与 `logger.debug` bug）与 `read_groups` | `DONE` | `P1` | `T12` | 29/29 单元测试与 `short_toy.bam` / `chr9.4M.ont.sim.read_groups.tsv` Oracle 全通过 |

---

## Phase 3：核心匹配、零 XGBoost 树引擎与模型重建 (`Node 4 ~ Node 6`)

| ID | 阶段/节点 | 任务名称 | 状态 | 优先级 | 依赖 | 交付物与验证标准 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`T14`** | Node 4 | 编写 `_cy_profile_assign.pyx`（读取 `.iqidx` `mmap` 视图生成紧凑轮廓与剪接匹配），替代 Python `list[int]` 轮廓 | `DONE` | `P0` | `T08, T12` | `test_long_read_profile.py` & `test_junction_comparator.py` 100% 等价 |
| **`T15`** | Node 4 | 拆分 `long_read_assigner`、轻量 `ReadAssignment`、`assignment_io`（上下文管理）与轻量整数数组多映射归约（废弃跨进程对象 pickling） | `DONE` | `P0` | `T14` | `test_long_read_assigner.py`, `test_output_formats.py`, `test_assess_assignment_quality.py` 全通过 |
| **`T16`** | Node 5 | 实现统一 `TerminalHistogram`（`_cy_hist_stats.pyx` 稀疏 `(pos, count)` 矩，废弃 `[pos]*count` 展开）、`np.searchsorted` 分组输出与纯 Python `CagePeakFinder` | `DONE` | `P0` | `T06, T12` | `test_terminal_peaks.py`, `test_polya_cage_finder.py` 全通过，移除 `pybedtools` 依赖 |
| **`T17`** | Node 5 | 实现 `_cy_tree_eval.pyx`（`model_polya.json` & `model_tss.json` L1-cache `float32 gbtree` 无锁推理）与 `polya_finder`/`verification` | `DONE` | `P0` | `T16` | `test_polya_prediction.py` 全通过，模型推理与原生树定义数学误差 `< 1e-6` |
| **`T18`** | Node 6 | 实现显式 `ChromosomeDiscoveryState` 替换 `ModelStore` 类级全局状态，重写流式 `IntronGraph` 与 5 阶段局部转录本重建 | `DONE` | `P0` | `T15, T17` | `test_intron_graph*.py`, `test_graph_alt_ends.py`, `test_transcript_splice_site_corrector.py` 全通过 |

---

## Phase 4：Out-of-Core 定量、融合检测、全管道与可视化 (`Node 7 ~ Node 10`)

| ID | 阶段/节点 | 任务名称 | 状态 | 优先级 | 依赖 | 交付物与验证标准 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`T19`** | Node 7 | 拆分 `long_read_counter`、Polars Out-of-Core 聚合、`scipy.sparse.coo_matrix + mmwrite` 重写 MTX 转换与清理辅助脚本 | `DONE` | `P0` | `T15, T18` | `tests/test_node7_quantification.py` 全通过，各计数策略 TSV (`atol <= 1e-6`) 与 MTX 稀疏矩阵 100% 等价 |
| **`T20`** | Node 8 | 重写无全局缓存字典的 `FusionDetector`、静态排序 `GenomicIntervalIndex` 与 `minimap2` CLI / Cython 局部校验后端 | `DONE` | `P2` | `T08, T09` | `tests/test_node8_fusion.py` 全通过，清除 4 个全局缓存字典并用有界 `lru_cache` 与 `.iqidx` 重现 Oracle 产出 |
| **`T21`** | Node 9 | 搭建端到端 Stage Graph 与染色体/Locus 单遍流式 Worker 编排器（Node 3->4->5->6->7 在线消费即释放），支持多样本分片归并 | `DONE` | `P0` | `T11, T13, T18, T19, T20` | `tests/test_node9_pipeline.py` 全通过（`simple_data/` 两阶段与单遍全文件、`MAPT` FASTQ + CAGE、`--yaml` + `--illumina_bam` 与 Oracle 一致） |
| **`T22`** | Node 9 | 实现基于 `os.replace()` 的原子 Stage Manifest `--resume` 状态机、YAML 多样本集与 `filter_umis` 日志初始化修复 | `DONE` | `P1` | `T21` | `tests/test_node9_pipeline.py` 验证 `tenX_v3` SC/UMI `detect_barcodes` 逐行等价及中断 `--resume` 状态恢复 |
| **`T23`** | Node 10 | 重写 `visualizer/`（优先读取 `run_manifest.json`、非交互 `Agg` 后端、显式 `plt.close(fig)` 防 OOM、Pandas 3.0 CoW 兼容） | `DONE` | `P2` | `T04, T21` | `tests/test_node10_visualization_and_compat.py` 全通过（AST 零禁导、50 次渲染 RSS 稳定且 `get_fignums() == []`） |
| **`T24`** | Node 10 | 在 `basicpy` (`NumPy 2.4.6 + Pandas 3.0.2`) 跑通核心数学/矩阵/表格前瞻兼容性验证测试集 | `DONE` | `P1` | `T19, T23` | `basicpy` `-W error::DeprecationWarning -W error::FutureWarning` 下全绿通过 |

---

## Phase 5：大规模基准验证与回归验收 (`Gate G1 ~ G4`)

| ID | 阶段/节点 | 任务名称 | 状态 | 优先级 | 依赖 | 交付物与验证标准 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`T25`** | Phase 5 (G1) | Tier 1 轻量多样本基准对比（`ENCFF400BQQ`, `ENCFF684YOO`, `ENCFF339FMQ`，`-t 1/8/16`，交替 `A-B-B-A` >=3 轮，YAML 多样本） | `DONE` | `P0` | `T03, T21` | `ENCFF400BQQ` `-t 8` 3 轮 `A-B-B-A`：`1840.98s -> 130.43s`（`14.11x`）、Peak PSS `-73.9%`、计数/`read_info` `100%` 精确；预编译 `.iqidx` `31.02x` |
| **`T26`** | Phase 5 (G2) | Tier 2 多协议对比测试（`ENCFF436GKZ`, `ENCFF846RID`, `ENCFF207KVO`, `RBMX_KD/con_rep1` + `gencode.v46` 1.6G GTF + `CAGE`） | `DONE` | `P0` | `T25` | `8.43x–14.18x` 加速，Peak PSS `-62.0% ~ -69.4%`（`Peak RSS -71.7% ~ -75.7%`），`read_info >= 0.999990` |
| **`T27`** | Phase 5 (G3) | Tier 3 千万 reads 全规模压力测绘（`simulated1` 8.77M, `ENCFF598YQO` 10.1M, `ENCFF089AMF` 16.4M, `shRBMX_rep2` 7.8G BAM） | `DONE` | `P0` | `T26` | 零失败/零 `/dev/shm` 残留；`simulated1` 精确率 `0.966097` vs 原版 `0.966096`、`ext_gtf Jaccard=1.0`；`shRBMX_rep2` `1.72h -> 21.9m`（`4.71x`）、RSS `55.68 -> 19.41 GiB`（`-65.1%`） |
| **`T28`** | Phase 5 (G4) | `ENCFF008HMY.fastq` 从头 `minimap2` 端到端全链路测绘 & 固化 Slurm 复现脚本与阶段基准报告 | `DONE` | `P1` | `T27` | `ENCFF008HMY.fastq` `61.4m -> 6.34m`（`9.68x`），`read_info` 一致率 `0.999998` |

---

## Phase 6：全样本・全协议・全物种无上限 Slurm Array 大闭环 (`G0 ~ G4` 全量 125 项)

| ID | 阶段/节点 | 任务名称 | 状态 | 优先级 | 依赖 | 交付物与验证标准 |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **`T29`** | Phase 6 | 扩展 `summarize_pair.py`（增加 `discovered_transcript_counts` 与外显子链对齐对比）与自动化校验全目录队列生成器 `build_all_samples_manifest.py` | `DONE` | `P0` | `T28` | `test_benchmark_tools.py` 全通过；产出零缺失 **125 项**全队列 [all_samples_manifest.json](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/benchmarks/all_samples_manifest.json) |
| **`T30`** | Phase 6 | 通过无并发上限 Slurm Array（`all_samples_slurm_array.slurm`：`2785890`, `2785902`, `2785903` + `2786427`）完成全部人类/小鼠/海牛、`directRNA`/`cDNA`/`capped`/`R2C2`/`rmdup`/`RBMX_KD`+`CAGE`/`SIRV_simulation`/`Sequin`/`YAML`/`FASTQ`/`--illumina_bam`/`--resume`+`isoquant-next-visualize` | `DONE` | `P0` | `T29` | **125 / 125 完成（0 失败）**：121 组配对累计 `144.89h -> 23.51h`（中位 `6.3x`，最高 `37.55x`，PSS `-56.8%`，RSS `-66.6%`），全队列 `tx_r >= 0.9999818`、`gene_r >= 0.99998199`、`read_info >= 0.999057`；海牛无 GTF 从头发现 `-80.9%` PSS、`r>=0.99999987` |
| **`T31`** | Phase 6 | 运行 `aggregate_cohort.py --strict`、`plot_all_cohort.py` 渲染 6 张出版级 PNG 与 `build_all_cohort_pdf.py` 编译 12 页整合单文件 PDF 并逐页视觉验收 | `DONE` | `P0` | `T30` | `cohort_all_results.{json,tsv}`、`cohort_by_{protocol,gate}.tsv`、`docs/figures/fig1..6.png` 与 `docs/isoquant_next_benchmark_report.pdf`（12 页，`5.35 MiB`）全部零遮挡、零裁切 |
| **`T32`** | Phase 6 | 更新 `docs/benchmark_report.md`、`docs/task.md`、彻底清理 `/gpfs/.../iq*` 临时目录并完成 AST 禁导、`sqanti3`/`basicpy` 全测、`/dev/shm` 与只读 `git status` 最终验收 | `DONE` | `P0` | `T31` | 所有门禁校验全绿，无 `/tmp` 写入或 `/dev/shm` 残留，原版代码与只读数据零修改 |
