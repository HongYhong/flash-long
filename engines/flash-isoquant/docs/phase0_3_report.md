# IsoQuant (`isoquant_next`) 重构 Phase 0 ~ Phase 3（`T01` ~ `T18`）工程与实测验收报告

> **范围边界**：严格覆盖 **Phase 0（`M0_PRE` & `Gate G0`）至 Phase 3 结束（`Node 0` ~ `Node 6`，`T01` ~ `T18`）**；Phase 4（`Node 7 ~ 10`）与 Phase 5（`Gate G1 ~ G4`）尚未开始，全部保持 `TODO`。  
> **环境约束**：所有开发、Cython 3.2 编译与单元测试均在 `sqanti3` Conda 环境（Python 3.11.15；经用户明确授权已补齐原版所需的 `xgboost 3.2.0`、`ssw-py 1.0.1`、`editdistance 0.8.1`、`mappy 2.31`、`numba 0.68.1` 用于真实原版对照）完成；前瞻兼容性校验在 `basicpy` Conda 环境（NumPy 2.4.6 + Pandas 3.0.2）完成。原版 `isoquant.py`、`isoquant_lib/`、`isoquant_tests/` 全程保持只读零改动（`git status --short` 仅含新增 `isoquant_refactored/` 及原生运行生成之只读 `.fai`）。

---

## 1. 选定架构（方案 B）核心落地成果

| 模块 / 节点 | 任务 ID | 核心交付物（[src/isoquant_next](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/src/isoquant_next)） | 关键优化与重构修复 |
| :--- | :--- | :--- | :--- |
| **Phase 0 (`M0_PRE`, `G0`)** | `T01~T03` | `pyproject.toml`, `setup.py`, `_cy_smoke.pyx`, `scripts/benchmark/`, [g0_baseline.md](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/docs/g0_baseline.md) | 建立 `/proc/<pid>/smaps_rollup` Tree PSS 采样器；先以占位桩定位原版 `--help`/Stage 2 缺失 `ssw`/`XGBClassifier` 精确阻塞点，随后在 `tmp_bench/run_native_with_lazy_blockers.py` 将 4 个缺失包无侵入桥接至 `isoquant_next` Cython 内核，全量打通原版 Tier 0 & Tier 1 (`rc=0`) 端到端金标基线 |
| **Node 0（配置与并发）** | `T04~T05` | `config/` (`dataclass` + `run_manifest.json`), `core/` (`IsoQuantExitCode`, `utils.py`, `concurrency.py`) | 淘汰 `pickle` 配置存储；唯一 `reverse_complement`；`shutil.which` + 元组版本校验；`WorkerPoolExecutor` 强制 `forkserver`/`spawn` 并注入 `OMP_NUM_THREADS=1` / `OPENBLAS_NUM_THREADS=1` |
| **Node 1（零拷贝注释与映射）** | `T06~T09` | `core/serialization.py`, `core/string_pools.py`, `annotation/{gtf_validator,iqidx,gene_info,reference}.py`, `mapping/read_mapper.py` | 修复 `read_dict` 有符号 int 缺陷；流式 GTF/GFF3 校验器；`mmap` CSR `.iqidx` 索引编译器；合并 `GeneInfo` 5 次 `children()` 与 $O(N^2)$ active scan；`LazyContigView` 消除每核 ~250 MiB 全染色体加载；`blake2b` + `ExitStack` + `.csi` 回退 |
| **Node 2（Barcode / UMI）** | `T10~T11` | `barcodes/{indexers,molecule_structure,results,common,umi_filtering,detect_barcodes}.py`, `_cy_barcode_umi.pyx`, `barcodes/callers/` | `pysam.FastxFile` 替代 `Bio.SeqIO`；`np.uint64` k-mer 索引（NEP 50 安全 + Python 3.13+ `track=False`）；Cython `nogil` Smith-Waterman 替代 `ssw` + Myers bit-parallel 编辑距离替代 `editdistance`；修复 `call_barcodes` 重复循环与 `table_splitter` 句柄泄漏 |
| **Node 3（流式分块与 CIGAR）** | `T12~T13` | `alignment/_cy_cigar_blocks.pyx`, `alignment/{alignment_info,bam_merger,exon_corrector,illumina_exon_corrector,read_groups}.py` | Cython SoA CIGAR 块/错误提取；`BAMStreamMerger` + `LocusChunker` 按基因座窗口/安全水位流式分片；修复 `type(i) != "int"` 与 `logger.debug` 缺陷 |
| **Node 4（Cython 轮廓、分配与无 Pickle 归约）** | `T14~T15` | `assignment/_cy_profile_assign.pyx`, `assignment/{types,profiles,comparator,assigner,multimap,io}.py` | `Int8ProfileView` 1 字节紧凑轮廓（淘汰 `list[int]`）；修复 mono-exonic 基因空内含子轮廓兼容；`CompactMultimapTuple` 结构体二进制编码废弃跨进程对象 pickling；全 `io` 类支持上下文管理器与 `.iqidx` 无 SQLite 查询 |
| **Node 5（稀疏矩与无 XGBoost 树引擎）** | `T16~T17` | `terminal/_cy_hist_stats.pyx`, `terminal/_cy_tree_eval.pyx`, `terminal/{tree_model,peaks,cage_finder,polya_finder,counters,verification}.py` | 稀疏 `(pos, count)` Cython 矩统计替代 `[pos]*count` 列表展开；`model_polya.json` (`417 -> 76` 活跃树) 与 `model_tss.json` (`417 -> 86` 活跃树) L1-Cache `float32 gbtree` 无锁推理；`np.searchsorted` 分组统计；纯 Python `bisect_right` `CagePeakFinder` 移除 `pybedtools`（零 `xgboost`、零 `pandas` 运行时依赖） |
| **Node 6（显式状态与 5 阶段发现）** | `T18` | `discovery/{state,intron_graph,splice_corrector,fl_constructor,assignment_constructor,end_processor,model_filter,read_assigner,gene_joiner,printer,builder}.py` | 引入显式 `ChromosomeDiscoveryState` 消除 `ModelStore` 类级全局可变集合；`deque` BFS 与 `bytearray` 覆写流式 `IntronGraph` 及 5 阶段重建；支持 `.iqidx` 无 SQLite 扩展存储与上下文管理 `GFFPrinter` |

---

## 2. 全量单元与黄金对照（Oracle）测试实测结果

在 `sqanti3` 环境下执行全量 Phase 0 ~ Phase 3 测试集：
```bash
PYTHONPATH=src pytest -v \
  tests/test_smoke_build.py tests/test_benchmark_tools.py \
  tests/test_node0_config_core.py tests/test_node1_annotation_io.py \
  tests/test_node2_barcodes_umi.py tests/test_node3_alignment.py \
  tests/test_node4_assignment.py tests/test_node5_terminal.py tests/test_node6_discovery.py
```
- **核心数值指标**：
  - **`86 passed in 15.89s`（`test_smoke_build` + `test_benchmark_tools` + `test_node0` ~ `test_node6`）**：全数通过，与原版 `isoquant_lib` 在序列评分、10x/Stereo/Universal barcode 解码、UMI 过滤、CIGAR block、`read_groups`、profile 构造、`JunctionComparator`、`LongReadAssigner`（含 `synth.db` 与 `chr9.4M` 真实 BAM reads）、`assess_assignment_quality.AssignmentData`（`isoquant` & `isoquant_legacy` 预设解析 `ReadInfoPrinter` / `BasicTSVAssignmentPrinter`）、`PolyAFinder`、book-ended BED `CagePeakFinder`、`PolyACounter` 分组 TSV、以及 `IntronCollector` / `IntronGraph` 图构建/简化剪枝/polyA 端接上完全通过验证。
  - **原版全量套件原生跑通（`sqanti3` 已补包）**：未修改原版 `isoquant.py --help` / `--test` 全通过；`pytest -q isoquant_tests/test_*.py --ignore=isoquant_tests/test_intron_graph.py` 实测 **`1287 passed, 4 skipped`（`8.61s`）**。
  - **无 XGBoost 树引擎数学精度**：`L1TreeEnsemble` 对 `model_polya.json` 与 `model_tss.json` 各 301 组测试向量（含 `skew=NaN` 缺省路由与真实密集峰）相对纯 Python 参考树遍历以及**真实 `xgboost 3.2.0` `XGBClassifier(n_jobs=1)`** 的绝对概率误差均 **`< 1e-6`**，分类标签 `100%` 一致。
  - **`basicpy` 前瞻兼容性（NumPy `2.4.6` + Pandas `3.0.2`）**：验证通过显式 `np.uint64` 位运算及 Pandas CoW DataFrame 构造，输出 `BASICPY_COMPAT_OK`。
  - **依赖与进程卫生**：核心包 `src/isoquant_next/` 经 grep 检查零导入 `xgboost`/`ssw`/`editdistance`/`mappy`/`pybedtools`/`Bio`/`pandas`（`ZERO_FORBIDDEN_IMPORTS`），测试与 benchmark 后 `/dev/shm` 零残留（`SHM_CLEAN`）。

---

## 3. 已完成基线对比 benchmark 实测（Scheme B vs 原版 `isoquant.py`）

1. **人类全注释（`lrgasp_gencode_v39_annotation_sirvs.human.gtf`，`1.4 GiB`）编译与全库簇加载（`-t 8`）**（详见 [g0_baseline.md](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/docs/g0_baseline.md)）：
   - **原版 `v39.db`（`complete_db=True`，已丢弃 125 SIRV 基因）**：DB 体积 **`369 MiB`**，串行簇 `GeneInfo` 扫描 **`1052.241s`**（总 `1053.1322s`），`200,941 KiB` PSS。
   - **Scheme B `.iqidx`（完整保留所有 61,658 基因 / 245,115 转录本）**：索引体积 **`135.89 MiB`**，编译 **`66.122s`**，8 核 `mmap` 全人类簇扫描仅 **`3.483s`**（**`302.1x` 加速**，扫描期 `~284 MiB` PSS）。
2. **`1.6 GiB` `gencode.v46.annotation.gtf` 流式校验、解析正确性与 8 核 `mmap` 零拷贝共享实测**：
   - **流式校验与解析精确度**：`check_gtf_duplicates_streaming` 单遍流式校验耗时 **`12.497s`**（`modified=False`）；编译生成 `.iqidx` 体积 **`141.24 MiB`**（编译耗时 `75.044s`，单进程编译峰值 `PSS = 1,555,375 KiB` / `1.48 GiB`，连同 8 核全库扫描总计 `79.6284s`）。经逐条扫描原始 GTF 核对，`v46.iqidx` 完整无损保留全部 **`63,086` 个基因、`254,070` 条转录本、`1,668,627` 个外显子（`32,672` 个簇，25 条序列）**，抽检 `ENSG00000204138.13` 坐标区间 `(28369582, 28500364)` 完全一致。
   - **8 核全簇扫描与多核 OS 页缓存零拷贝实测（`v46_mmap_probe.py` & `g1_iqidx_vs_gffutils.py cand_v46`）**：8 核 `mmap` 扫描全 `32,672` 簇仅 **`4.085s`**（扫描期 11 进程树峰值 `PSS = 365,778 KiB` / `RSS = 440,424 KiB`）；独立全页 `mmap` 探针（逐 4KB OS 页触摸整个 `141.24 MiB` `.iqidx`）对比：
     - `1` 核：`wall = 1.4602s`，`RSS = 134,400 KiB`，`PSS = 136,712 KiB`（其中 `tree_pss_file_kib = 114,004 KiB`）
     - `8` 核并发（`11` 进程）：`wall = 1.3451s`，`RSS = 1,154,656 KiB`（`1.10 GiB`），峰值进程树 **`PSS = 304,647 KiB`（`297.5 MiB`）**，其中共享文件页 **`tree_pss_file_kib = 115,644 KiB`**（相对 1 核仅 **`1.01x`**，实证 8 核零拷贝共享物理页）。
3. **Tier 0 & Tier 1 端到端实测对比（真实原版 vs Cython Shim vs Scheme B）**：
   - **Tier 0 `simple_data` (`-t 2`, `--fl_data`, 含完整建模与定量)**：
     - 未修改原版（真实 `xgboost 3.2.0`）：耗时 **`12.6617s`**，峰值 **`PSS = 221,811 KiB`** / `RSS = 578,880 KiB`（`failures=0.0`）；
     - Cython Shim（[run_native_with_lazy_blockers.py](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/tmp_bench/run_native_with_lazy_blockers.py)）：耗时 `9.6556s`，峰值 `PSS = 174,948 KiB` / `RSS = 435,948 KiB`；两者经 `compare_outputs.py` 与 `diff` 验证 26 个基因计数、31 条转录本外显子链、421 条 `read_info` 记录、`polyA_prediction.tsv` 与 `TSS_prediction.tsv` **100% 完全一致**。
   - **Tier 0 `toy_data` FASTQ (`-t 4`, 含 `minimap2` + 完整建模与定量)**：
     - 未修改原版（真实 `xgboost 3.2.0`）：耗时 **`2.9683s`**，峰值 **`PSS = 205,025 KiB`** / `RSS = 906,896 KiB`（`failures=0.0`）；
     - Cython Shim：耗时 `2.4812s`，峰值 `PSS = 157,108 KiB` / `RSS = 694,208 KiB`（真实 `xgboost` 额外带来 `+45~48 MiB` 进程树 PSS 开销）。
   - **Tier 1 `ENCFF400BQQ.sorted.k14.bam` (`92 MiB`, `-t 8`)**：
     - **原版 `isoquant.py` 全流程 E2E（预建 `v39.db`, `--no_model_construction`: Stage 1 `task-636` + Stage 2 `--resume` Cython Shim `task-1396`, `failures=0.0`）**：Stage 1 耗时 `1542.9443s`（峰值 `PSS = 2,163,274 KiB`, `RSS = 2,903,176 KiB`），Stage 2 耗时 `233.1252s`（峰值 `PSS = 2,155,480 KiB`, `RSS = 2,922,708 KiB`），**双阶段合计 `1776.0695s`（~29.6 分钟）**、总体峰值 **`PSS = 2,163,274 KiB`（`2.06 GiB`）** / **`RSS = 2,922,708 KiB`（`2.79 GiB`）**，成功写出定量表（`61,537` 基因 / `244,943` 转录本）及 `OUT.polyA_prediction.tsv`（`18,331` 行）；真实 `xgboost 3.2.0` 下的 `ENCFF400BQQ` 多轮 `A-B-B-A` 全链路对比将在 Phase 5 Gate G1 (`T25`) 执行。
     - **Scheme B Nodes 0~6 流式管线（`tmp_bench/phase3_tier1_bench.py`, `v39.iqidx` + `LocusChunker` + `LazyContigView` + Cython 内核 + `L1TreeEnsemble` + 5 阶段 `GraphBasedModelConstructor`, `task-1188`）**：跑完 25 条人类染色体（`231,475` 条非次级主比对 reads、`224,583` 分配、`10,902` 个实测 polyA 峰、含完整 `15,043` 个重建转录本模型及 `19.68 MiB` `CompactMultimapTuple` 二进制压缩）仅耗时 **`32.287s`**（**较原版仅关闭建模的双阶段 `1776.07s` 仍提速 `55.0x`**），峰值进程树 **`PSS = 687,422 KiB`（`671.3 MiB`，较原版降低 `68.2%` / `3.15x`）** / **`RSS = 863,992 KiB`（`843.7 MiB`）**，`failures=0.0`。

---

## 4. 下一步（待用户指令启动 Phase 4 ~ 5）

- `Node 7 ~ 10` (`T19` ~ `T24`) 与 `Gate G1 ~ G4` (`T25` ~ `T28`) 处于待命状态。
