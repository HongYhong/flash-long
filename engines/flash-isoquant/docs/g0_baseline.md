# Gate G0 实测基线与原版可运行边界报告（`sqanti3` 环境）

> 测量环境：`/gpfs/group/data_sci_group/weiwu_group/hongyanhong/miniconda3/miniconda3/envs/sqanti3`（Python 3.11.15）。
> 环境变更说明：Phase 0 初始体检时 `sqanti3` 缺少 `ssw` / `editdistance` / `xgboost` / `mappy` / `numba`，曾以 [run_native_with_lazy_blockers.py](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/tmp_bench/run_native_with_lazy_blockers.py) 定位阻塞点并在 Phase 3 后以内存 Cython Shim 跑通原版对照；随后经用户明确授权已在 `sqanti3` 安装补齐 `xgboost 3.2.0`、`ssw-py 1.0.1`、`editdistance 0.8.1`、`mappy 2.31`、`numba 0.68.1`，现已可直接运行未修改的原版 `python isoquant.py` 全套对照（`--help` / `--test` / `isoquant_tests/` / 端到端流程）。
> 数据来源：所有指标均由 `scripts/benchmark/run_tiered_benchmark.py`（`monitor_process_tree.sh` 采样 `/proc/<pid>/smaps_rollup`）与 `scripts/benchmark/compare_outputs.py` 实测生成。

---

## 1. 原版 `isoquant.py` 可运行性边界（装包前后对比）

| 检查项 | 命令/入口 | 结果 | 阻塞位置与实测结论 |
| :--- | :--- | :--- | :--- |
| Phase 0（未装 5 包）：CLI `--help` / `--test` | `python isoquant.py --help` / `--test` | FAIL (`rc=1`) | [isoquant.py:55](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant.py#L55) 导入 `ssw` 报错 |
| **已补包（真实原版）：CLI `--help` / `--test`** | `python isoquant.py --help` & `python isoquant.py --test` | **PASS (`rc=0`)** | `--help` 正常退出；`--test` 完整跑完 `TEST_DATA` 并输出 `=== TEST PASSED CORRECTLY ===` |
| **已补包（真实原版）：`isoquant_tests/` 全集** | `pytest -q isoquant_tests/test_*.py --ignore=isoquant_tests/test_intron_graph.py` | **1287 passed, 4 skipped** (`8.61s`) | Phase 0 原本因缺包导致 13 组 collect error，补包后零报错全量跑通 |
| **已补包（真实原版）：Tier 0 `simple_data/`** | `python isoquant.py` (`-t 2`, `--fl_data`, 含建模+真实 `XGBClassifier 3.2.0`) | **PASS (`rc=0`, 完整产出 14 类结果)** | **`12.6617s`**，`PSS = 221,811 KiB`（`216.6 MiB`），`RSS = 578,880 KiB`；经 `compare_outputs.py` 与 `diff` 校验与 Cython Shim 输出（26 特征计数、31 转录本外显子链、421 `read_info` 记录、`polyA_prediction.tsv`、`TSS_prediction.tsv`）**100% 完全一致** |
| **已补包（真实原版）：Tier 0 `toy_data/`** | `python isoquant.py` (`MAPT.Mouse.ONT.simulated.fastq`, `-t 4`, `minimap2` + 真实 `XGBClassifier 3.2.0`) | **PASS (`rc=0`, 完整产出 13 类结果)** | **`2.9683s`**，`PSS = 205,025 KiB`（`200.2 MiB`），`RSS = 906,896 KiB` |
| Cython Shim `--resume`：Tier 1 `ENCFF400BQQ` | `Stage 1` (`1542.9443s`) + `Stage 2 --resume` (`233.1252s`) (`-t 8`, `--no_model_construction`) | **PASS (`rc=0`, Stage 1 + Stage 2 完成)** | **总耗时 `1776.0695s`**（~29.6 分钟），总体峰值 **`PSS = 2,163,274 KiB` (`2.06 GiB`)** / **`RSS = 2,922,708 KiB` (`2.79 GiB`)**；产出 `OUT.gene_counts.tsv` (`61,537` 行)、`OUT.transcript_counts.tsv` (`244,943` 行)、`OUT.polyA_prediction.tsv` (`18,331` 行)；真实 `xgboost 3.2.0` 的 `ENCFF400BQQ` 多轮 `A-B-B-A` 全链路对比将在 Phase 5 Gate G1 (`T25`) 执行 |

---

## 2. 基线耗时与进程树峰值内存实测表

| 评测项 | 样本与配置 | 重复 | Exit | Wall Median (s) | Peak Tree PSS | Peak Tree RSS | 关键输出计数 |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- |
| Native `gtf2db + GeneInfo` | `chr9.4M.gtf.gz` | 3 | 0 | `1.3822` (`gtf2db=0.286s`, `gi=0.213s`) | `64,011 KiB` (`62.5 MiB`) | `68,992 KiB` | 23 genes / 54 tx |
| Native Tier 0 (`simple_data` Cython Shim) | `chr9.4M.ont.sim.polya.bam` (`-t 2`, `--fl_data`) | 1 | 0 | `9.6556` | `174,948 KiB` (`170.8 MiB`) | `435,948 KiB` | `unique=137`, `known=25`, `novel=6`, 全量 TSV/GTF |
| **Native Tier 0 E2E (`simple_data` 真实原版)** | `chr9.4M.ont.sim.polya.bam` (`-t 2`, `--fl_data`, `xgboost 3.2.0`) | 1 | **0** | **`12.6617`** | **`221,811 KiB`** (`216.6 MiB`) | `578,880 KiB` | 26 基因计数 / 31 转录本 / 421 `read_info` 与 Shim **100% 一致** |
| Native Tier 0 (`toy_data` FASTQ Cython Shim) | `MAPT.Mouse.ONT.simulated.fastq` (`-t 4`) | 1 | 0 | `2.4812` | `157,108 KiB` (`153.4 MiB`) | `694,208 KiB` | `unique=117`, `noninformative=15`, `known=10` |
| **Native Tier 0 E2E (`toy_data` FASTQ 真实原版)** | `MAPT.Mouse.ONT.simulated.fastq` (`-t 4`, `xgboost 3.2.0`) | 1 | **0** | **`2.9683`** | **`205,025 KiB`** (`200.2 MiB`) | `906,896 KiB` | 全量 13 类 TSV/GTF 完备 |
| Native Tier 1 (`ENCFF400BQQ` Stage 1) | `ENCFF400BQQ.sorted.k14.bam` + `v39.db` (`-t 8`, `--no_model_construction`) | 1 | 1 (`xgboost` 占位阻断) | `1542.9443` | `2,163,274 KiB` (`2.06 GiB`) | `2,903,176 KiB` (`2.77 GiB`) | Stage 1 全部染色体 (`OUT.save_chr*`) |
| Native Tier 1 (`ENCFF400BQQ` Stage 2 `--resume` Shim) | `encff400bqq_quant --resume -t 8` | 1 | **0** | `233.1252` | `2,155,480 KiB` (`2.06 GiB`) | `2,922,708 KiB` (`2.79 GiB`) | `unique=87,369`, `unique_minor=22,249`, `polyA=18,331` 行 |
| **Native Tier 1 E2E 合计 (`ENCFF400BQQ` Stage 1+2)** | `Stage 1 + Stage 2 (--resume)` (`-t 8`, `--no_model_construction`) | 1 | **0** | **`1776.0695`** (~29.6 min) | **`2,163,274 KiB`** (`2.06 GiB`) | **`2,922,708 KiB`** (`2.79 GiB`) | 完整产出定量矩阵、`polyA_prediction.tsv` 与 `read_info.tsv.gz` |
| Native `gtf2db(complete=False)` | `lrgasp_gencode_v39_sirvs` (1.4 GiB) | 1 | 取消（超 443s 仅 8%） | `> 443.74` | `1,041,016+ KiB` (`~1.01 GiB`) | `~1.01 GiB` | `~2.7 GiB` DB，卡在 `980,108` 关系推理 |
| Native `gtf2db(complete=True)` | 同上 | 1 | 0 | N/A (峰值 `> 3.1 GiB` PSS) | `3,191,115 KiB` (`3.04 GiB`) | `3,192,912 KiB` | `v39.db` (`369 MiB`)；丢弃 125 SIRV 基因 |
| Native Cluster `GeneInfo` Scan | `v39.db` 全库串行簇构造 (`task-511`) | 1 | 0 | **`1053.1322`** (`1052.241s` scan) | `200,941 KiB` (`196.2 MiB`) | `204,696 KiB` | 33,028 clusters / 61,533 genes / 244,939 tx |
| **Candidate `.iqidx` (Scheme B)** | **直接编译 `v39.gtf` -> `.iqidx` (`-t 8` 全库簇构造)** | 1 | 0 | **`70.5608`** (`compile=66.122s`, `scan=3.483s`) | `1,277,737 KiB`（含编译峰值；多核扫描期 `~284 MiB`） | `1,261,308 KiB` | `v39.iqidx` (**`135.89 MiB`**)，完整保留 SIRV：**33,142 clusters / 61,658 genes / 245,115 tx** |

> 关键发现：`v39.iqidx` 将基因库体积由默认 DB 的 `~2.7 GiB`（`complete_db=True` 为 `369 MiB` 且丢失 SIRV 注释）压缩至 **`135.89 MiB`**，并通过 8 核零拷贝 `mmap` 将全人类 `GeneInfo` 集群加载时间从原版 SQLite **`1052.241s`** 压缩至 **`3.483s`（302.1x 加速）**。
