# IsoQuant Scheme B (`isoquant-next`) 全样本全协议全物种基准验收报告

> **对比基线**：未修改的原版 `python isoquant.py`（IsoQuant `4.0.0`，使用真实 `xgboost 3.2.0`、`ssw-py 1.0.1`、`editdistance 0.8.1`、`mappy 2.31`、`numba 0.68.1`、`gffutils`）。  
> **优化候选**：重构版 `isoquant-next`（Scheme B：`.iqidx` 零拷贝 `mmap`、`LocusChunker`、`LazyContigView`、`_cy_*` `nogil` Cython 内核、`L1TreeEnsemble`、`CompactMultimapTuple`、`ChromosomeDiscoveryState`、Polars Out-of-Core + `scipy.sparse.coo_matrix` 定量）。  
> **测量基础设施**：`scripts/benchmark/run_tiered_benchmark.py`（底层调用 `scripts/benchmark/monitor_process_tree.sh` 1 Hz 采样 `/proc/<pid>/smaps_rollup` 进程树 `Peak PSS` / `Peak RSS` 与 `wall/user/sys time`）、`scripts/benchmark/summarize_pair.py`（含 `gene_counts`、`transcript_counts`、`discovered_transcript_counts` 及外显子链对齐对比）及 `scripts/benchmark/eval_simulated_sirv.py`。  
> **全队列 Slurm 复现作业**：
> - 先前 Phase 5 锚点：`scripts/benchmark/phase5_slurm_array.slurm`（Slurm Job Array `2784905` / `2784995`，17/17 完成）；
> - 全量无上限 Slurm Array：`scripts/benchmark/build_all_samples_manifest.py` -> [all_samples_manifest.json](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/benchmarks/all_samples_manifest.json)（**125 项队列**），经 `scripts/benchmark/all_samples_slurm_array.slurm`（无 `%K` 并发限制 `--array=0-N`：Job IDs `2785890`、`2785902`、`2785903`，以及 `ENCFF400BQQ` 全基因组散点采样作业 `2786427`）全部执行完毕：**125 / 125 任务完成（`missing_count = 0`, `failed_count = 0`）**。

---

## 1. Gate G0：核心注释引擎微基准与原版基线

- 原版全套单元测试 `pytest -q isoquant_tests/test_*.py --ignore=isoquant_tests/test_intron_graph.py` 在 `sqanti3` 实测 `1287 passed, 4 skipped`。
- Tier 0 `simple_data/`（`-t 2`, `--fl_data`）：原版 `12.6617s`、`Peak PSS = 221,811 KiB`（`216.6 MiB`），详见 [g0_baseline.md](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/docs/g0_baseline.md)。
- 人类全注释库（`lrgasp_gencode_v39_annotation_sirvs.human.gtf`，`1.4 GiB`）：
  - 原版 `gtf2db(complete=False)` 超过 `443.74s` 仍停留在 ~8%（`~2.7 GiB` DB）；`gtf2db(complete=True)` 峰值 `PSS = 3,191,115 KiB`（`3.04 GiB`）且丢失 `125` 个 SIRV 基因，后续 33,028 cluster 加载耗时 `1052.241s`。
  - 重构版 `.iqidx` 直接编译 GTF 生成 `135.89 MiB` 索引（`66.122s`），保留全量 `33,142` clusters / `61,658` genes / `245,115` transcripts，8 核 `mmap` 全库 `GeneInfo` 构造仅耗时 **`3.483s`（302.1x 加速）**；`gencode.v46`（`1.6 GiB`）`.iqidx` 体积 `141.24 MiB`，8 核扫描 `4.085s`、多进程共享 `PSS` 仅 `304,647 KiB`。
- `ENCFF400BQQ`（`-t 1 / -t 8 交替 A-B-B-A 3 轮 / -t 8 热索引 / -t 16`）：
  - 原版从 `-t 8 -> -t 16` 出现负扩展（`1840.98s -> 2090.46s`），`Peak PSS` 从 `9.13 GiB` 升至 `14.54 GiB`（`Peak RSS = 21.17 GiB`）。
  - `isoquant-next` 在 `-t 1 -> -t 8 -> -t 16` 保持正扩展（`414.40s -> 130.43s -> 111.06s`），热 `.iqidx` 预编译对比可降至 **`64.69s`（`31.02x` 加速）**，且 4 组配置均实现定量计数与 `read_info` **`100%` 逐行完全一致**。

---

## 2. 全队列（125 任务 / 121 配对双引擎）宏观性能、内存与精度总览

结构化汇总文件：[cohort_all_results.json](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/benchmarks/cohort_all_results.json) ｜ [cohort_all_results.tsv](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/benchmarks/cohort_all_results.tsv) ｜ [cohort_by_protocol.tsv](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/benchmarks/cohort_by_protocol.tsv) ｜ [cohort_by_gate.tsv](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/benchmarks/cohort_by_gate.tsv)。

- **累计墙钟时间（121 组配对引擎基准）**：原版累计 **`144.89 小时`** vs 重构版 **`23.51 小时`**；
- **全队列中位加速比**：**`6.3x`**（区间 **`0.94x – 37.55x`**）；
- **进程树峰值内存压缩（中位数）**：`Peak Tree PSS` **`-56.8%`**，`Peak Tree RSS` **`-66.6%`**；
- **全队列精度底线（Min Floor across all 121 paired tasks）**：`gene_pearson_r >= 0.99998199`，`transcript_pearson_r >= 0.9999818`，`read_info_full_agreement >= 0.999057`。

### 2.1 按测序协议与运行模式分层汇总（125 项全覆盖）

| Protocol / Mode | 任务数 / 配对数 | 原版中位 Wall | 新版中位 Wall | 中位加速比 `[Min–Max]` | 原版 -> 新版 Peak PSS (中位) | PSS 降幅 | RSS 降幅 | 最低 Tx `r` | 最低 `read_info` | 最低模型 Jaccard |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `cDNA_PacBio`（含人/鼠/海牛/`-t1/8/16/warm`） | `17 / 17` | `3888.8s` (`1.08h`) | `261.7s` (`4.4m`) | **`16.08x`** `[1.04–31.02]` | `15.42 -> 4.41 GiB` | **`-69.5%`** | **`-74.7%`** | `0.99999998` | `0.999907` | `0.997672` |
| `cDNA_PacBio_v46`（`ENCFF370NFS` + `v46`） | `1 / 1` | `4199.5s` (`1.17h`) | `276.2s` (`4.6m`) | **`15.20x`** | `16.35 -> 4.70 GiB` | **`-71.2%`** | **`-76.6%`** | **`1.0`** | `0.999988` | `0.999122` |
| `cDNA_PacBio_vM36`（`ENCFF313VYZ` + `vM36`） | `1 / 1` | `4126.6s` (`1.15h`) | `280.2s` (`4.7m`) | **`14.73x`** | `12.71 -> 4.32 GiB` | **`-66.0%`** | **`-71.2%`** | **`1.0`** | `0.999969` | `0.998378` |
| `capped_PacBio`（全人/鼠 12 复本） | `12 / 12` | `4007.3s` (`1.11h`) | `184.1s` (`3.1m`) | **`21.33x`** `[12.04–37.55]` | `14.71 -> 3.89 GiB` | **`-72.7%`** | **`-73.0%`** | `0.99999993` | `0.999960` | `0.998043` |
| `MultiSample_YAML`（2rep / 跨系 / RBMX 3cond） | `3 / 3` | `6429.2s` (`1.79h`) | `455.7s` (`7.6m`) | **`14.11x`** `[7.10–19.69]` | `15.62 -> 4.84 GiB` | **`-69.0%`** | **`-73.3%`** | `0.9999818` | `0.999057` | `0.997926` |
| `directRNA`（全人/鼠 12 复本） | `12 / 12` | `2915.9s` (`48.6m`) | `139.0s` (`2.3m`) | **`16.85x`** `[9.62–31.59]` | `14.64 -> 4.31 GiB` | **`-68.6%`** | **`-73.2%`** | `0.99999998` | `0.999976` | `0.996032` |
| `cDNA_ONT`（全人/鼠/海牛 15 复本） | `15 / 15` | `7113.0s` (`1.98h`) | `932.4s` (`15.5m`) | **`6.20x`** `[1.60–12.17]` | `20.12 -> 8.40 GiB` | **`-52.9%`** | **`-61.3%`** | `0.99999987` | `0.999829` | `0.996500` |
| `RBMX_KD`（3 条件 x 2 复本 + `v46`） | `6 / 6` | `5802.3s` (`1.61h`) | `815.5s` (`13.6m`) | **`8.30x`** `[4.71–9.85]` | `17.60 -> 10.12 GiB` | **`-42.5%`** | **`-70.1%`** | `0.9999818` | `0.999057` | `0.997926` |
| `RBMX_KD_CAGE`（`--cage` + `v46`） | `3 / 0` | —（CAGE 特有验证） | `489.9–851.7s` | — | `-> 6.31–11.77 GiB` | — | — | — | — | `5,539–7,689` 模型 |
| `capped_ONT`（全人/鼠 12 复本） | `12 / 12` | `4411.2s` (`1.23h`) | `1175.8s` (`19.6m`) | **`3.47x`** `[2.93–5.69]` | `21.06 -> 13.78 GiB` | **`-39.6%`** | **`-62.3%`** | `0.99999992` | `0.999953` | `0.996324` |
| `R2C2_ONT`（全人/鼠 24 复本） | `24 / 24` | `4015.1s` (`1.12h`) | `990.6s` (`16.5m`) | **`5.08x`** `[3.00–9.88]` | `16.70 -> 8.00 GiB` | **`-51.0%`** | **`-64.1%`** | `0.99999992` | `0.999936` | `0.996928` |
| `R2C2_ONT_rmdup`（`seqkit rmdup` 12 复本） | `12 / 12` | `5034.4s` (`1.40h`) | `1188.1s` (`19.8m`) | **`4.17x`** `[2.82–5.92]` | `18.76 -> 10.08 GiB` | **`-47.6%`** | **`-61.5%`** | `0.99999998` | `0.999959` | `0.996928` |
| `SIRV_simulation`（`simulated1..3` + `v46`） | `3 / 3` | `3998.0s` (`1.11h`) | `703.5s` (`11.7m`) | **`5.64x`** `[3.98–5.71]` | `15.58 -> 6.83 GiB` | **`-56.1%`** | **`-65.9%`** | **`1.0`** | `0.999992` | `0.997932` |
| `Sequin_SIRV4`（`ont_combined_no_supp`） | `1 / 1` | `878.3s` (`14.6m`) | `937.0s` (`15.6m`) | `0.94x`（超高深度少数 contig） | `17.74 -> 17.23 GiB` | `-2.9%` | `-11.6%` | **`1.0`** | `0.999998` | `0.985401` |
| `FASTQ_directRNA`（从头 `minimap2`） | `1 / 1` | `3683.8s` (`1.02h`) | `380.6s` (`6.3m`) | **`9.68x`** | `21.56 -> 21.60 GiB` | `-0.2%` | **`-3.4%`** | **`1.0`** | `0.999998` | `0.996763` |
| `Illumina_Hybrid`（`chr9` `--illumina_bam`） | `1 / 1` | `4.74s` | `2.94s` | **`1.61x`** | `0.22 -> 0.05 GiB` | **`-75.8%`** | **`-95.8%`** | **`1.0`** | **`1.0`** | **`1.0`** |
| `Resume_Visualize`（`ENCFF400BQQ` 断点+出图） | `1 / 0` | `110.28s`（完整运行） | `40.65s`（`resume`） | `2.71x` | `3.83 -> 3.68 GiB` | `-3.8%` | `-2.4%` | **`1.0`** | **`1.0`** | **`1.0`** |

### 2.2 按验证阶段 (`G0_G1 ~ G4`) 与物种分层汇总

| 分组维度 (`Gate / Species`) | 任务数 / 配对数 | 原版总耗时 | 新版总耗时 | 中位加速比 | 原版 -> 新版中位 PSS | 中位 PSS 降幅 | 中位 RSS 降幅 | 最低 Tx / Disc `r` | 最低模型 Jaccard |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`G0_G1`**（Tier 1 并发扩展 & 热索引） | `3 / 3` | `11,988.6s` (`3.33h`) | `609.5s` (`10.2m`) | **`19.64x`** | `6.27 -> 2.24 GiB` | `-64.3%` | `-57.2%` | **`1.0`** | `0.998406` |
| **`G1`**（PacBio 全量 + 海牛 + v46/vM36 + YAML） | `30 / 30` | `110,606.4s` (`30.72h`) | `7,903.0s` (`2.20h`) | **`17.67x`** | `15.01 -> 4.15 GiB` | `-70.8%` | `-73.7%` | `0.99999993` | `0.997672` |
| **`G2`**（directRNA + cDNA_ONT + RBMX + CAGE + YAML） | `36 / 33` | `161,087.1s` (`44.75h`) | `21,332.4s` (`5.93h`) | **`9.25x`** | `16.09 -> 6.64 GiB` | `-61.0%` | `-71.7%` | `0.9999818` | `0.996032` |
| **`G3`**（capped_ONT + R2C2 + rmdup + SIRV + Sequin） | `53 / 53` | `234,231.5s` (`65.06h`) | `54,413.2s` (`15.11h`) | **`4.52x`** | `18.56 -> 9.76 GiB` | `-47.6%` | `-63.1%` | `0.99999985` | `0.985401` |
| **`G4`**（FASTQ + `--illumina_bam` + `--resume` + Vis） | `3 / 2` | `3,688.6s` (`1.02h`) | `383.5s` (`6.4m`) | **`5.65x`** | `10.89 -> 10.83 GiB` | `-37.8%` | `-49.6%` | **`1.0`** | `0.996763` |
| **Human**（`H1_DE`, `H1_ES`, `WTC11`, `RBMX_KD`） | `93 / 89` | `425,637.6s` (`118.23h`) | `67,149.6s` (`18.65h`) | **`7.10x`** | `17.25 -> 7.89 GiB` | `-55.0%` | `-66.6%` | `0.9999818` | `0.996324` |
| **Mouse**（`F121_cell` + `vM28` / `vM36`） | `22 / 22` | `77,489.4s` (`21.52h`) | `10,442.3s` (`2.90h`) | **`7.61x`** | `12.78 -> 4.81 GiB` | `-59.6%` | `-67.7%` | `0.99999988` | `0.996032` |
| **Manatee**（无 GTF 从头发现，15,608 contigs） | `5 / 5` | `6,815.0s` (`1.89h`) | `4,003.7s` (`1.11h`) | `1.83x` | `22.26 -> 4.25 GiB` | **`-80.9%`** | **`-81.8%`** | `0.99999987` | `0.999443` |
| **Human_SIRV**（NanoSim `simulated1..3` + `v46`） | `3 / 3` | `10,777.1s` (`2.99h`) | `2,106.0s` (`35.1m`) | **`5.64x`** | `15.58 -> 6.83 GiB` | `-56.1%` | `-65.9%` | **`1.0`** | `0.997932` |
| **Sequin_SIRV4**（高深度外源混合体） | `1 / 1` | `878.3s` (`14.6m`) | `937.0s` (`15.6m`) | `0.94x` | `17.74 -> 17.23 GiB` | `-2.9%` | `-11.6%` | **`1.0`** | `0.985401` |
| **Human_chr9**（`chr9.4M` + `--illumina_bam`） | `1 / 1` | `4.7s` | `2.9s` | **`1.61x`** | `0.22 -> 0.05 GiB` | **`-75.8%`** | **`-95.8%`** | **`1.0`** | **`1.0`** |

> **关键全队列洞察**：
> 1. **ENCFF400BQQ 全基因组逐行零误差**：[scatter_encff400bqq](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/benchmarks/scatter_encff400bqq) 实测显示全基因组共 `61,533` 条基因记录与 `244,939` 条转录本记录中，原版与重构版相关性 **`Pearson r = 1.000000`**、`exact match = 100.00%`（`max |diff| = 0`）。
> 2. **海牛无参考从头发现模式（15,608 contigs, `-80.9%` PSS / `-92.1%` Peak PSS）**：`Trichechus manatus` 5 个样本均产出 **`1.0`** 或 **`0.9999999`** 发现转录本相关度与 **`0.999443–1.0`** 模型结构 Jaccard。其中 `ENCFF589LDN` CPU User Time 由 `3581.9s -> 48.9s`（**-98.6%**）、Peak PSS 从 `22.06 -> 1.75 GiB`（**-92.1%**），墙钟加速 `1.04x`（`433.9s -> 415.6s`）受 GPFS 对 ~10,988 条空 contig 进行元数据 merge/unlink 限制；较高通量的 `ENCFF835YMN` / `ENCFF090GGQ` / `ENCFF629AWL` 取得 **`1.60x–1.87x`** 端到端加速。
> 3. **Sequin_SIRV4（仅 exon-only GTF）**：`lrgasp_sirvs4.gtf.gz` 仅含 `464` 条 `exon` 行（无 `gene`/`transcript` 父行），故不传 `--complete_genedb` 实测 `ont_combined_no_supp`：`gene_r = 1.0`, `tx_r = 1.0`, `ext_gtf_jaccard = 1.0`, `read_info = 0.999998`。因 `3.82M` reads 极度集中在少数人工 spike-in contig 上受单 contig 分区并发度上限影响，耗时 `878.3s -> 937.0s`（`0.94x`）。

---

## 3. Gate G3 `SIRV_simulation`（`simulated1..3` + `gencode.v46`）NanoSim 真值评估

真值标签直接解析 read ID 中携带的 `ENST...` 来源标注：

| Replicate | 引擎 | Wall (s) | Speedup | Peak PSS (GiB) | Read Precision | Read Recall | Model Known Prec. | Model Recall (`>=2r`) | Quant Precision | Quant Recall (`>=2r`) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`simulated1`** | 原版 `isoquant.py` | `2759.9` | — | `15.58` | `0.966096` | `0.275456` | `0.970258` | `0.647117` | `0.879445` | `0.869042` |
| **`simulated1`** | **`isoquant-next`** | **`693.6`** | **`3.98x`** | **`6.83`** (`-56.1%`) | **`0.966097`** | **`0.275457`** | **`0.970628`** | `0.647087` | **`0.879458`** | `0.869017` |
| **`simulated2`** | 原版 `isoquant.py` | `4019.1` | — | `15.42` | `0.966106` | `0.275584` | `0.970663` | `0.648927` | `0.881198` | `0.869536` |
| **`simulated2`** | **`isoquant-next`** | **`703.5`** | **`5.71x`** | **`6.91`** (`-55.2%`) | **`0.966108`** | **`0.275585`** | **`0.970880`** | **`0.648986`** | `0.881195` | `0.869511` |
| **`simulated3`** | 原版 `isoquant.py` | `3998.0` | — | `15.78` | `0.965965` | `0.275630` | `0.970185` | `0.646309` | `0.880636` | `0.869965` |
| **`simulated3`** | **`isoquant-next`** | **`708.9`** | **`5.64x`** | **`6.79`** (`-57.0%`) | **`0.965966`** | **`0.275629`** | **`0.970431`** | **`0.646315`** | `0.880634` | `0.869937` |

- 全部 3 个 `simulated1..3`（每组 ~18.98M read records + `gencode.v46`）上，`isoquant-next` 的 `Read Precision` 与 `Model Known Precision` **全部略优于或精确追平原版**，且 `extended_annotation.gtf` Jaccard 全部达到 **`1.000000`**。

---

## 4. 高级特性、多样本 YAML、CAGE、Hybrid 短读段校正、断点续跑与可视化验收

1. **多样本 `--yaml` 分片归并（3 组全覆盖）**：
   - `g1_multisample_yaml_t16`（`ENCFF684YOO + ENCFF339FMQ`）：`4419.5s -> 224.5s`（**`19.69x`**），`PSS -73.6%`，`tx_r = 1.0`, `models_jacc = 0.998448`。
   - `g1_cross_cell_yaml_t16`（`H1_DE + H1_ES + WTC11` 跨细胞系）：`6429.2s -> 455.7s`（**`14.11x`**），`PSS -69.0%`，`tx_r = 0.99999998`, `read_info = 0.999988`。
   - `g2_rbmx_yaml_t32`（`con_rep1 + ama_rep1 + shRBMX_rep1` 3 条件 + `v46`, `-t 32`）：`12036.6s (3.34h) -> 1694.9s (28.2m)`（**`7.10x`**），`Peak RSS` `48.95 -> 13.08 GiB`（**-73.3%**），`tx_r = 0.9999818`, `ext_gtf = 1.0`。
2. **`--cage` 辅助 TSS 发现（`RBMX_KD` 3 条件全覆盖 + `gencode.v46` + `CAGE_WTC11.all_reps_v2.bed`）**：
   - `con_rep1`：`489.9s` (`8.2m`)，`19,767` 条 `TSS_prediction.tsv` 记录、`6,052` 转录本模型；
   - `ama_rep1`：`625.2s` (`10.4m`)，`17,725` 条 `TSS_prediction.tsv` 记录、`5,539` 转录本模型；
   - `shRBMX_rep1`：`851.7s` (`14.2m`)，`24,417` 条 `TSS_prediction.tsv` 记录、`7,689` 转录本模型。
3. **`--illumina_bam` 与 `shortreadSJ` 格式审计（`Gate G4`）**：
   - `g4_chr9_illumina_t8`（`chr9.4M.ont.sim.bam` + `--illumina_bam short_toy.bam`）：`4.74s -> 2.94s`（`1.61x`），`Peak RSS` `1.62 -> 0.07 GiB`（**-95.8%**），所有输出表与模型 **`1.0` 完全一致**。
   - **数据格式核查备注**：`backup/longpass_benchmark250209/genome/shortreadSJ/es_rep1SJ.out.tab` 为 9 列 STAR `SJ.out.tab` 剪接位点表（`chrom, start, end, strand_code, motif, ...`），而非 `isoquant` `--illumina_bam` 或 `--junc_bed_file` 所需的 BAM/BED 格式，因此脚本未做违背原语义的强行输入。
4. **`--resume` 状态机与 `isoquant-next-visualize`（`g4_resume_visualize_t16`）**：
   - 使用完整运行完成后的 `.stage_state.json`（保留 `reference_and_annotation` 已完成阶段并重跑下游阶段）触发 `--resume`：`110.28s -> 40.65s`（**`2.71x`**），日志与状态检验 `resume_skipped_reference_and_annotation = true`，最终产出与从头运行相比 `gene_r = 1.0, tx_r = 1.0, tx_models_jacc = 1.0, ext_gtf_jacc = 1.0, read_info = 1.0`。
   - `isoquant-next-visualize` 直接消费 `run_manifest.json` 在非交互 `Agg` 后端渲染 **6 张 PNG 与 6 张 PDF** 基因转录本丰度与外显子结构图（`ENSG00000074800.16` 等），已在 Figure 6D 中嵌图确认可读。

---

## 5. 前瞻兼容性与运行期约束验收

1. **`basicpy` 前瞻环境（`Python 3.11.15`, `NumPy 2.4.6`, `Pandas 3.0.2`, `SciPy 1.17.1`, `Matplotlib 3.10.8`, `pysam 0.23.3`）**：启用 `-W error::DeprecationWarning -W error::FutureWarning` 执行 `test_node10_visualization_and_compat.py` + `test_node0_config_core.py` + `test_node2_barcodes_umi.py`，覆盖 NEP 50 `uint64` 位运算、Pandas 3.0 CoW 与 `coo_matrix`/`mmwrite` 稀疏写出及 50 次循环绘图无内存泄露。
2. **依赖与只读铁律校验**：
   - AST 全扫描 `src/isoquant_next/` 确认零 `xgboost`、`ssw`、`editdistance`、`mappy`、`pybedtools`、`Bio`、`pandas` 运行时依赖；
   - 原版 `isoquant.py`、`isoquant_lib/`、`isoquant_tests/` 及只读数据目录 `backup/longpass_benchmark250209/` 零改动；
   - 全程禁用系统 `/tmp`，Slurm 与测试临时目录均收敛于 `/gpfs/group/data_sci_group/weiwu_group/hongyanhong/` 工作区并在收口时彻底删除，`/dev/shm` 与节点共享内存零残留。

---

## 6. 整合 PDF 报告与 6 张全队列出版级图表

- **单文件整合 PDF 报告（12 页：全量协议/Gate/物种总表 + `simulated1..3` + 特殊功能汇总 + Fig 1–6 + 125 行全量明细附录）**：  
  [isoquant_next_benchmark_report.pdf](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/docs/isoquant_next_benchmark_report.pdf)
- [fig1_speedup_and_walltime.png](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/docs/figures/fig1_speedup_and_walltime.png)：121 组配对 Slurm 基准跨协议墙钟耗时（log 刻度+单样本散点）与加速比分布（中位数 `6.3x`，最高 `37.55x`）。
- [fig2_memory_and_scaling.png](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/docs/figures/fig2_memory_and_scaling.png)：全协议 Peak Tree PSS（中位 `-56.8%`）/ RSS（中位 `-66.6%`）压缩与 `ENCFF400BQQ` `-t 1 / 8 / 16` + Warm `.iqidx` 扩展曲线。
- [fig3_correlation_and_accuracy.png](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/docs/figures/fig3_correlation_and_accuracy.png)：`ENCFF400BQQ` 全基因组 61,533 gene / 244,939 transcript 散点图（`Pearson r = 1.000000`, `100.00%` exact match）、121 组一致率分布与 `simulated1..3` NanoSim 真值对比。
- [fig4_pss_memory_profiles.png](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/docs/figures/fig4_pss_memory_profiles.png)：Tier 1 (`ENCFF400BQQ`)、Tier 3 (`shRBMX_rep2` 20.44M ONT + `v46`, `-t 32`)、海牛 (`ENCFF589LDN` 15,608 contigs, `-92% PSS`) 与 `ENCFF089AMF` (`16.42M`) 实时进程树内存轨迹。
- [fig5_species_cellline_heatmap.png](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/docs/figures/fig5_species_cellline_heatmap.png)：`Protocol x Species x Cell Line` 二维加速比、PSS 降幅与最小相关系数热图（Human / Mouse / Manatee）。
- [fig6_features_cage_resume_vis.png](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/docs/figures/fig6_features_cage_resume_vis.png)：`RBMX_KD` `--cage` TSS 预测、`--yaml`/`--illumina_bam`/`--resume` 一致率、`v46`/`vM36`/`Sequin`/`FASTQ` 对比及原生 `isoquant-next-visualize` 样图。

---

## 7. 复现方式

```bash
source /gpfs/group/data_sci_group/weiwu_group/hongyanhong/miniconda3/miniconda3/etc/profile.d/conda.sh && conda activate sqanti3
# 1. 生成并校验 125 项全样本全协议队列清单
python isoquant_refactored/scripts/benchmark/build_all_samples_manifest.py
# 2. 按需无并发上限提交（标准 -t 16 / 重负载 -t 32）
# sbatch --array=0-<N_STD-1> ... isoquant_refactored/scripts/benchmark/all_samples_slurm_array.slurm
# 3. 严格 125/125 汇聚、重绘 6 张总图与生成 12 页整合 PDF
python isoquant_refactored/scripts/benchmark/aggregate_cohort.py --manifest isoquant_refactored/benchmarks/all_samples_manifest.json --bench-root isoquant_refactored/benchmarks --strict
python isoquant_refactored/scripts/benchmark/plot_all_cohort.py --cohort-json isoquant_refactored/benchmarks/cohort_all_results.json --bench-root isoquant_refactored/benchmarks --scatter-dir isoquant_refactored/benchmarks/scatter_encff400bqq --out-dir isoquant_refactored/docs/figures
python isoquant_refactored/scripts/benchmark/build_all_cohort_pdf.py --cohort-json isoquant_refactored/benchmarks/cohort_all_results.json --fig-dir isoquant_refactored/docs/figures --out-pdf isoquant_refactored/docs/isoquant_next_benchmark_report.pdf
```
