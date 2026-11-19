# IsoQuant 重构工程规格与验收契约 (`SPEC`)

> **基线仓库**：IsoQuant `4.0.0`（`/gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/`，只读保留）
> **重构目录**：`/gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/`
> **参照规范**：`~/.codex/AGENTS.md`、`long_bench/docs/{spec,bambu_g4_performance,bambu_full_pipeline_spec}.md`

---

## 1. 目标、非目标与变更边界

### 1.1 核心目标
1. **分节点彻底重构**：在 `isoquant_refactored/` 中新建包 `isoquant_next`，严格按依赖拓扑（Node 0 ~ Node 10）拆分重整耦合状态、全局缓存、`pickle` 参数传递、上下文泄露与已知语义/并发缺陷。
2. **大幅降低运行内存与运行时间**：通过零拷贝共享注释索引、按基因座（Locus）流式处理、消除跨进程 Python 对象序列化、以及 Cython `nogil` 原生内核，系统性压低进程树峰值 `PSS` 与端到端 `wall_time`。
3. **零 Conda 环境改动可运行**：完全依托集群现有 Conda 环境编译本地 Cython 扩展与运行流水线，不在运行时强制依赖集群缺失的包（如 `xgboost`、`mappy`、`ssw-py`、`editdistance`）。
4. **输出向后严格兼容**：所有最终交付产物（GTF/GFF3、基因与转录本定量 TSV、Matrix Market MTX、`read_info.tsv.gz`、`read_assignments.tsv`、`*.polyA_prediction.tsv`、`*.TSS_prediction.tsv`、`*.fusions.tsv`、`*.allinfo` 等）格式及生物学判定结果与原版对齐。

### 1.2 明确非目标（Out of Scope）
- **禁止改动原目录**：不修改 `isoquant.py`、`isoquant_lib/`、`isoquant_tests/`、集群 Conda 环境或 `~/.codex/` 配置。
- **不承诺中间缓存二进制格式向后兼容**：`.params`（原版 `pickle`）与 `.save`（原版未修复有符号整型的自定义二进制缓存）在 `isoquant_next` 中升级为带 `schema_version` 的 JSON Manifest 与自研定长/分片二进制格式；仅保留 `visualizer` 对旧 `.params` 的可选只读降级。
- **不将未观测的估测值当作实测结论**：在运行真实基线与 Candidate 对比前，文档中所有加速比与内存降幅均标注为**设计目标**，须经 Gate G0~G4 实测验收。

---

## 2. 集群运行环境审计与依赖策略

### 2.1 现有 Conda 环境实测盘点（遵循 `~/.codex/AGENTS.md`）

| 环境名称 | 路径 | Python | 已验证关键包 | CLI 工具 | 角色定位 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **`sqanti3`** (**主选**) | `.../miniconda3/miniconda3/envs/sqanti3` | `3.11.15` | `pysam 0.22.1`, `gffutils 0.14`, `pyfaidx 0.9.0.4`, `Cython 3.2.4`, `polars 0.20.31`, `pyarrow 14.0.2`, `intervaltree`, `numpy 1.26.4`, `pandas 2.2.3`, `scipy 1.11.4`, `scikit-learn 1.5.2`, `Bio 1.81`, `pytest 7.4.4`, `PyYAML 6.0.3` | `minimap2`, `STAR`, `samtools`, `bedtools`, 系统 `/usr/bin/gcc 13.3.0` | **主力开发、Cython 编译与性能评测环境** |
| **`chorus`** | `.../miniforge3/envs/chorus` | `3.10.20` | `pysam 0.24.0`, `gffutils 0.14`, `pyfaidx 0.9.0.4`, `polars 1.43.2`, `pyarrow 25.0.0`, `intervaltree`, `numpy 1.26.4`, `pandas 2.2.3`, `scipy 1.15.2`, `Bio 1.88`, `pytest 9.1.1` | `samtools`, `bedtools` | 纯 Python/Polars 次选校验环境（缺 `Cython`/`minimap2`） |
| **`basicpy`** | `.../miniconda3/miniconda3/envs/basicpy` | `3.11.15` | `pysam 0.23.3`, `numpy 2.4.6`, `pandas 3.0.2`, `scipy 1.17.1`, `sklearn 1.8.0`, `Bio 1.87`, `pytest 9.1.1` | 无生信 CLI | **NumPy 2.0 (NEP 50) & Pandas 3.0 (CoW) 前瞻兼容性回归环境** |
| **`scanpy`** | `.../miniconda3/miniconda3/envs/scanpy` | `3.12.13` | `numba 0.64.0`, `pyarrow 23.0.1`, `numpy 2.4.6`, `pandas 2.3.3`, `scipy 1.18.0`, `sklearn 1.9.0` | 无 `pysam`/`gffutils` | 仅辅助校验 Numba 数学/矩阵微内核（不作主环境） |
| **`PTSDsnDLPFC`** | `.../miniforge3/envs/PTSDsnDLPFC` | `3.9.16` | `pysam 0.23.3`, `numba 0.57.1`, `pyarrow 16.1.0`, `numpy 1.24.3`, `pandas 1.5.3` | `samtools`, `bedtools` | Python 3.9 偏旧，不作为目标主环境 |

> **集群未安装工具与包审计（实测确认）**：
> - 扫描上述全部 Conda `envs/*/lib/python*/site-packages/` 后，确认 **没有任何环境安装 `xgboost`、`mappy`、`ssw-py` 或 `editdistance`**；同时系统不存在 `cargo` / `rustc`。
> - 规范命令激活方式（每次运行命令必须显式执行）：
>   `source /gpfs/group/data_sci_group/weiwu_group/hongyanhong/miniconda3/miniconda3/etc/profile.d/conda.sh && conda activate sqanti3`

### 2.2 四项缺失外部依赖的无缺包重构契约

| 原版缺失依赖 | 原版源码精确调用点 | 重构内建替代机制（`sqanti3` 零新包可运行） | 兼容性与对照 Oracle 说明 |
| :--- | :--- | :--- | :--- |
| **`xgboost`** (`XGBClassifier`) | `terminal_counter.py`, `terminal_peaks.py`, `end_processor.py` | **L1-Cache `float32` JSON 树评估器**：经查 `model_polya.json`（417 树，最多 79 节点/树，总 1,937 节点）与 `model_tss.json`（417 树，最多 81 节点/树，总 2,637 节点）均为 `gbtree`、`binary:logistic`、`split_type={0}`（纯数值 `<` 比较）。使用 Cython/NumPy 按 IEEE-754 `float32` 逐树累加叶权重加 `logit(base_score)` 后 sigmoid，**零 OpenMP 死锁且零 `xgboost`/`pandas` 运行时依赖**。 | 若检测到环境变量或测试机存在 `xgboost`，自动化回归测试额外校验与 `XGBClassifier` 逐样本差值 `< 1e-6`。 |
| **`ssw-py`** (`AlignmentMgr`) | `barcode_calling/common.py`（3 个函数，被 `tenx.py`/`curio.py`/`universal_extraction.py`/`stereo.py` 调用） | **Cython `nogil` Striped/Affine Smith-Waterman 内核**：针对 12~30bp barcode/linker 短序列匹配，直接复用对齐评分矩阵与回溯起止点逻辑；若存在 `ssw` 则保留等价一致性测试。 | Node 2 需核对 `common.py` 中 `AlignmentMgr` 默认 match/mismatch/gap 分值以保证最高分选链 100% 一致。 |
| **`editdistance`** (`eval`) | `umi_filtering.py:203`（唯一调用点：`ed = editdistance.eval(occ, umi)`） | **Cython Myers Bit-Parallel / 带界带状 Levenshtein 内核**：前置长度差剪枝 + 短串等长 Hamming 快路径 + 64-bit 位并行编辑距离（`nogil`）。 | 与标准 Levenshtein 距离数学完全一致，可在任意数据集逐 UMI 验证 0 差异。 |
| **`mappy`** (`mp.Aligner`) | `fusion_detector.py:7, 941` | **双模式局部重比对后端**：优先调用 `minimap2` CLI（`sqanti3/bin/minimap2`）批处理融合候选软剪切序列，小批量短片段回退至 Cython 局部 Smith-Waterman / k-mer seed 校验；若环境存在 `mappy` 则自动启用可选后端对照。 | Node 8 保留对齐坐标与 MAPQ/mlen 过滤契约测试。 |

---

## 3. 评测数据集规范（内置微套件 + `longpass_benchmark250209/`）

评测分为 4 个递增规模层级，覆盖从单模块单元测试到千万级 reads 全转录组压测：

### Tier 0：仓库内置微数据集（单元测试 & 快速逻辑等效验收）
- `/gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_tests/simple_data/`：
  `chr9.4M.gtf.gz`（9.9K）、`chr9.4M.fa.gz`（299K）、`chr9.4M.ont.sim.polya.bam`（135K）、`chr9.4M.ont.sim.fq.gz`（89K）、`chr9.4M.Illumina.bam`（1.1M）、`chr9.4M.ont.sim.read_groups.tsv`、`chr9.4M.yaml`
- `/gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_tests/toy_data/`：
  `MAPT.Mouse.genedb.gtf`（89K）、`MAPT.Mouse.reference.fasta`（100K）、`MAPT.Mouse.ONT.simulated.fastq`（331K）、`MAPT.Mouse.CAGE.bed`（1.9K）、`synth.db`（80K）
- `/gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_tests/short_reads_toy/`：`short_toy.bam`（3.8M）

### Tier 1：轻量真实长读长烟测与快速迭代集（< 500 MB BAM）
来源：`/gpfs/group/data_sci_group/weiwu_group/hongyanhong/backup/longpass_benchmark250209/`（及 `long_bench/config/bambu_refactor_samples.tsv` 清单校验数据）：

| 样本 ID | 协议 | 物种 | Reads 数 / BAM 大小 | 路径与参考文件 | 验证覆盖节点 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **`ENCFF400BQQ`** | `cDNApacbio` | human | 240,947 reads / 92 MiB | `cDNA_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF400BQQ.sorted.k14.bam` + `genome/lrgasp_gencode_v39_annotation_sirvs.human.gtf` + `genome/GRCh38.p14.genome.fa` | Node 1, 3, 4, 5, 6, 7（单核与多核基准快跑） |
| **`ENCFF684YOO`** | `cappedpacbio` | human | 887,281 reads / 149 MiB | `capped_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF684YOO.bam` + Gencode v39 + GRCh38 | 多样本定量 + CapTrap TSS 模式 |
| **`ENCFF339FMQ`** | `cappedpacbio` | human | 1,279,509 reads / 197 MiB | `capped_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF339FMQ.bam` + Gencode v39 + GRCh38 | 多样本联跑 & 合并 GTF |

### Tier 2：跨协议中型代表数据集（1M ~ 3M reads / 0.5 ~ 3.0 GiB BAM）
| 样本 ID | 协议 / 场景 | Reads 数 / 大小 | 路径与附加参考 | 验证目标 |
| :--- | :--- | :--- | :--- | :--- |
| **`ENCFF436GKZ`** | `cDNApacbio` (human H1) | 1,292,572 reads / 489 MiB | `cDNA_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF436GKZ.sorted.k14.bam` | 标准全节点基准 |
| **`ENCFF846RID`** | `directRNA` (human H1) | 1,247,419 reads / 1.6 GiB | `directRNA_fulpass/human/H1_ES_cell/alignments/ENCFF846RID.UtoT.fastq.sorted.k14.bam` | 链特异性 directRNA + polyA |
| **`ENCFF207KVO`** | `directRNA` (human H1) | 1,385,786 reads / 1.6 GiB | `directRNA_fulpass/human/H1_ES_cell/alignments/ENCFF207KVO.UtoT.fastq.sorted.k14.bam` | directRNA 多重复验证 |
| **`con_rep1`** | `RBMX_KD` ONT | 2.4 GiB BAM | `RBMX_KD/alignments/con_rep1.fastq.sorted.nosecondary.bam` + `genome/gencode.v46.annotation.gtf` (1.6G) | 大规模最新 Gencode v46 GTF 建库与分配内存 |
| **`CAGE & shortreadSJ`** | 多组学辅助输入 | 1.9G CAGE / 102M SJ | `genome/CAGEPAS/CAGE_WTC11.all_reps_v2.bed` + `genome/shortreadSJ/wtc11_rep1SJ.out.tab` + `genome/human.refTSS_v3.1.hg38.bed` | Node 3/5 CAGE & 短读 splice junction 支持 |

### Tier 3：高通量压力与百万读段全系统验收集（3 GB ~ 8 GB BAM & FASTQ）
| 样本 ID | 协议 / 场景 | Reads 数 / 大小 | 路径 | 压力测试指标 |
| :--- | :--- | :--- | :--- | :--- |
| **`simulated1`** | `SIRV_simulation` | 8,772,606 mapped reads / 1.8 GiB | `SIRV_simulation/alignments/simulated1_aligned_reads.sorted.bam` + `genome/GRCh38.p14.genome.fa` | 精确 Known/Novel 灵敏度、精度与 8.77M reads 吞吐量 |
| **`ENCFF598YQO`** | `R2C2ONT` (mouse F121) | 10,111,651 reads / 3.3 GiB | `R2C2_ONT_fulpass/mouse/F121_cell/alignments/ENCFF598YQO.bam` + `genome/lrgasp_gencode_vM28_sirvs.mouse.gtf` + `genome/GRCm39.genome.fa` | 千万 reads 小鼠全流程峰 PSS 与 wall time |
| **`ENCFF089AMF`** | `cappedONT` (human H1) | 16,424,133 reads / 4.8 GiB | `capped_ONT_fulpass/human/H1_ES_cell/alignments/ENCFF089AMF.bam` + Gencode v39 + GRCh38 | 1.64 千万 reads ONT 超高并发防 OOM 压测 |
| **`shRBMX_rep2`** | `RBMX_KD` ONT | 7.8 GiB BAM | `RBMX_KD/alignments/shRBMX_rep2.fastq.sorted.nosecondary.bam` + `gencode.v46` (1.6G) | 最大单 BAM 文件流式峰值 PSS 检验 |
| **`ENCFF008HMY`** | FASTQ 从头比对入口 | 1.6 GiB FASTQ | `/gpfs/group/data_sci_group/weiwu_group/hongyanhong/backup/longpass_benchmark250209/ENCFF008HMY.fastq` | Node 1 `minimap2` 子进程流式比对 & 缓存路径指纹 |

---

## 4. 系统结构与 Node 0 ~ Node 10 数据契约

### 4.1 全局数据契约原则
1. **显式不可变配置（Node 0）**：禁止向 `argparse.Namespace` 注入衍生字段；配置拆分为 `RunConfig`、`InputSpec`、`MatchingConfig`、`DiscoveryConfig`、`QuantConfig`、`BarcodeConfig`（均 `@dataclass(slots=True, frozen=True)`），输出 `run_manifest.json`（含 `schema_version`, `isoquant_next_version`, `command_fingerprint`, `stage_states`）。
2. **只读共享注释索引（Node 1）**：将 `GTF/GFF3`（及兼容性 `.db`）单遍编译为 **`IsoQuantMmapIndex` (`.iqidx`)**，内部以紧凑 NumPy/Struct 定长数组存储（CSR 布局：`chroms -> genes -> transcripts -> exons/introns`，所有字符串入池为 `int32` ID），通过 `mmap(..., MAP_SHARED)` 在主进程与所有 worker 间零拷贝共享。
3. **按 Gene Locus 流式分片（Node 3-4-6-7）**：
   - 废除整条染色体同时持有数百万条 `ReadAssignment` + `MappedReadProfile` Python 对象的模式；
   - 以“不与任何下游读段/基因区间重叠的基因组间隙（Locus Boundary）”为释放点：每个基因座计算完成后，立刻执行多映射统计累加、末端直方图累加、转录本模型构建与计数，仅将轻量二进制/Arrow 结果块写入磁盘分片，随即释放该基因座全部对象。
4. **数值与类型安全**：
   - 位掩码、k-mer 编码、哈希统一使用显式 `np.uint64`（严守 NumPy 2.0 NEP 50）；
   - Pandas 数据表操作仅限导出/汇总边界，禁用 `df[col][mask] = ...` 链式赋值与 `iterrows()`；
   - 所有文件/HTSlib 句柄强制 `with` 上下文生命周期管理，移除 `__del__` 隐式关闭。

---

## 5. 等价性判定与性能/内存评测口径

遵循 `long_bench/docs/bambu_g4_performance.md` 工程严谨标准：
1. **对照运行隔离（Oracle vs Candidate）**：
   - 所有计时基准必须处于**正式非调试模式**（关闭逐 read debug/trace checkpoint 转储，避免测量序列化开销污染纯算法耗时）。
   - Baseline 与 Candidate 测试同一输入时**交替运行顺序**（`A-B-B-A` 或轮替样本序列），每项核心基准至少 `>= 3` 次重复，汇报 `median` 与 `min-max` 区间，而非单次 luck-run 最小值。
   - 对缺失 4 个非核心包的原版对照组（如涉及 polyA XGBoost/barcode SSW 的模式）：基准对比分为**可直接在 `sqanti3` 跑通的原版功能子集**（如不开启 `polya` XGBoost 与 barcode calling 的常规 `--bam + --genedb/--gtf` BULK 定量/建模基准，待实测确认原版哪些参数路径受缺失包阻断）与**模块级合成/已知真值参考集**（如 XGBoost 离线金标推理表 / `SIRV_simulation` 真值评估），绝不虚构未跑通的对照基准。
2. **内存采集标准（进程树 `PSS`，而非单一 `VmRSS`）**：
   - 多进程 Python 程序中，简单对 `RSS` 求和会将只读共享页（如 `mmap` 页缓存、Python 共享 `.so` 代码段）重复统计 $N$ 次，造成严重失真；
   - 性能基准工具 `scripts/benchmark/monitor_process_tree.sh` 参照 `long_bench/scripts/bambu_exon_block_parallel/monitor_process_tree.sh` 设计，周期性遍历主进程及全部子/孙进程 `/proc/<pid>/smaps_rollup` 累加 `Pss`，同步记录 `RSS_sum`、`Peak_PSS`、`wall_time`、`user_time`、`sys_time` 与 `/usr/bin/time -v` 最大驻留集。
3. **等价性验收门槛（Correctness Gate）**：
   - **定量矩阵（Node 7）**：`gene_counts.tsv`、`transcript_counts.tsv` 在所有测试数据集上基因/转录本集合及整数/EM 计数达到**100% 精确一致**（浮点 EM 允许 `atol <= 1e-6`）；MTX 非零条目及行列坐标一致。
   - **Read 分配（Node 4）**：`read_assignments.tsv` / `read_info.tsv.gz` 逐条 read 的 `isoform_id`、`assignment_type`、`events` 一致率 `100%`（或仅当修复原版明确 bug 如 `read_dict` 负数/未定义行为时明确记录差异条目）。
   - **末端预测与模型重建（Node 5-6）**：`model_polya.json` / `model_tss.json` 推理得分最大绝对误差 `< 1e-6`；`transcript_models.gtf` 与 `extended_annotation.gtf` 转录本外显子链（intron chain）与原版完全一致，边界坐标偏差为 0。
   - **设计性能目标（待 G1~G4 实测验证）**：相比原版在同线程数与同 `sqanti3` 环境下，目标达成 **峰值进程树 `PSS` 降低 >= 60%（高并发大样本目标 70%~85%）**，**核心计算与端到端 `wall_time` 加速 >= 3x~8x**。
