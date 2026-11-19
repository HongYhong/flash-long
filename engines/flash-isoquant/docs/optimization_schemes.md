# 低内存与高性能重构架构方案对比 (`Optimization Schemes`)

> **目标**：在保证生物学结果与接口向后兼容、且**不改动集群现有 Conda 环境**的前提下，大幅压缩 IsoQuant 在亿级/千万级长读长样本（如 `ENCFF089AMF` 16.4M reads、`RBMX_KD/shRBMX_rep2` 7.8G BAM + `gencode.v46` 1.6G GTF）上的**运行时间（Wall-clock Time）**与**多进程峰值内存（Process-Tree PSS）**。

---

## 1. 原版 IsoQuant 4.0.0 内存与时间瓶颈解剖

基于对 `isoquant.py` 与 `isoquant_lib/` 共 80 个 Python 文件的源码审计，当前高内存、长耗时的根因集中在六个结构性问题：

| 瓶颈维度 | 原版关键源码位置 | 具体低效机理 | 实测风险影响 |
| :--- | :--- | :--- | :--- |
| **1. GTF/SQLite N+1 与 `GeneInfo` 膨胀** | `gtf2db.py`, `gene_info.py`, `dataset_processor.py` | 1. 建库全量 `readlines()` + `+=` 字符串拼接，且判断 `.gff3` 缺前缀点号。<br>2. 启动全表扫描 `all_features()`。<br>3. 每个基因多次调用 `db.children()`（至少 5 轮冗余 SQLite 查询）；`set_gene_attributes()` 对每个转录本再查全部外显子（`O(T*E)`）；`set_feature_properties()` 标注存在 `O(N^2)` 双重循环。<br>4. 每个基因创建完整 Python `GeneInfo` + 嵌套列表/字典，在多进程中被复制/反序列化多份。 | Gencode v46 (1.6G GTF) 加载阶段慢且各 worker 堆内存随基因数线性倍增。 |
| **2. 染色体级 Python 对象图与双重 IPC 序列化** | `dataset_processor.py`, `parallel_workers.py`, `processed_read_manager.py`, `long_read_profiles.py`, `isoform_assignment.py` | 1. 每条 read 生成 `CombinedReadProfiles` + 3 个 `MappedReadProfile`（Python `list[int]`）+ `BasicReadAssignment` + `IsoformMatch` + `MatchEvent` 等多层动态属性对象。<br>2. 按整条染色体或大区域累积对象列表，`ProcessPoolExecutor` 经由 pipe 将巨大 Python 对象图 pickle 回传主进程；部分模式先写成 `.save` 自定义二进制，下一阶段再完整反序列化回 Python 对象。 | 单条染色体（如 human chr1/chr2）上千万 reads 的对象包装与 pickle 缓冲区极易撑爆十数 GB~数十 GB PSS，且 GC（垃圾回收）停顿严重。 |
| **3. Fork Copy-on-Write (CoW) 失效与 OpenMP 并发雷区** | `isoquant.py`, `dataset_processor.py`, `detect_barcodes.py` | 1. Python 引用计数（`Py_INCREF`/`Py_DECREF`）与垃圾回收头字段修改会在读取共享 Python 对象时立即触发内核 **CoW 页分裂**，导致原本希望通过 `fork` 共享的注释数据在 16~32 worker 中变成 $16\sim32\times$ 物理内存。<br>2. `fork` 与 `xgboost`/`scipy`/BLAS OpenMP 线程层混用存在子进程死锁风险。 | 多核并发（`-t 16/32`）下实际物理 `PSS` 急剧攀升而非线性摊薄。 |
| **4. 末端特征展开分配与 XGBoost + Pandas 推理开销** | `terminal_peaks.py`, `terminal_counter.py`, `end_processor.py` | 1. 将稀疏 `{pos: count}` 末端直方图通过 `data.extend([pos] * count)` **物理展开成百万元素整数列表**再算均值/方差/偏度/熵。<br>2. 对仅有 8 个标量特征、417 棵浅树（总计 <=2,637 节点）的 `model_polya.json` / `model_tss.json`，每次候选窗构造 `pd.DataFrame` 调用 `XGBClassifier`，OpenMP 调度开销远大于树遍历本身。<br>3. `TerminalCounter._write_grouped` 使用 `df.iterrows()` 嵌套匹配。 | 大量无谓堆分配、Pandas DataFrame 构造消耗以及缺失 `xgboost` 无法运行。 |
| **5. 单细胞 FASTQ / Barcode / UMI 纯解释器热循环** | `detect_barcodes.py`, `barcode_calling/common.py`, `umi_filtering.py` | 1. 使用 `Bio.SeqIO.parse` 对亿级 FASTQ reads 为每条序列建立 `SeqRecord` 与质量字典。<br>2. Barcode 候选打分逐条调用 `AlignmentMgr`（缺乏向量/无 GIL 批处理）。<br>3. `umi_filtering.py` 两两比较缺少 Hamming/长度剪枝。 | FASTQ 解析及 Python 回调占单细胞流程绝大部分耗时。 |
| **6. 定量统计与 MTX 导出逐行迭代** | `long_read_counter.py`, `convert_grouped_counts.py` | 1. 计数聚合依赖高层字典/StringPool 查找循环。<br>2. `convert_grouped_counts.py` 用 `df.iterrows()` 逐行拼接文本写 Matrix Market，耗时长且文件句柄缺少上下文关闭。 | 多分组/百万细胞-基因稀疏矩阵导出极为缓慢。 |

---

## 2. 四套候选优化方案设计对比

围绕用户提出的“**大幅降低运行内存与运行时间**”与“**复用系统既有 Conda 环境（优先 `sqanti3`）**”双约束，评估以下四种架构方案：

### 方案 A：零编译纯列式 Python 重构（`pysam` + `NumPy` + `Polars` + `PyArrow` + `mmap`）
- **核心思路**：不写 Cython/C 扩展，纯用 `sqanti3` 现有向量化库（`NumPy`, `Polars 0.20`, `PyArrow 14`, `SciPy`, `pysam 0.22`）重写核心：
  - `GTF` 预解析为 Parquet / Arrow IPC 只读表并通过 `pyarrow.memory_map` 按需切片读取；
  - `ReadAssignment` 以列式 Arrow/Struct numpy array 存储；
  - 末端直方图改 `np.bincount` 加权矩；MTX 用 `scipy.sparse.coo_matrix` + `mmwrite`；表格拆分/聚合用 `Polars` lazy scan。
- **优点**：零编译门槛，可同时直接在 `sqanti3` 与 `chorus` 运行。
- **短板**：
  1. 外显子/内含子区间剖面构建、IntronGraph 路径搜索、reads CIGAR block 解析仍是 Python `for` 循环，解释器调度开销依旧显著；
  2. 对集群缺失的 `ssw`、`editdistance`、本地局部比对无法通过 `Polars` 加速（纯 Python Smith-Waterman 速度较 C/SIMD 慢 50~200x）。
- **预计收益（设计估算）**：运行时间 2x~4x 提速；峰值 PSS 降低 ~40%~55%。

### 方案 B（主推方案）：`mmap` 零拷贝 CSR 注释索引 + `O(Locus)` 流式分片 + Cython `nogil` 原生内核 + 无 OpenMP 树引擎
- **核心思路**：充分利用 `sqanti3` 环境原生自带的 **`Cython 3.2.4` + `/usr/bin/gcc 13.3.0` + `pysam 0.22.1` + `Polars 0.20.31` + `PyArrow 14.0.2` + `NumPy 1.26.4` + `SciPy 1.11.4` + `minimap2`/`samtools`**，将全部热路径与系统缺失的 4 个底层算法下沉至编译级 `.so` 模块 + OS 内核共享页缓存，同时在架构上用**基因座流式窗口**替代整染色体全量对象留存。
- **五大核心支柱**：
  1. **支柱一：单遍流式 `GTF -> .iqidx` 零拷贝共享内存索引（消灭 `GeneInfo` 与 CoW 页分裂）**
     - 主进程（或离线构建步）单遍扫描 GTF/GFF3，直接生成扁平定长 C 布局二进制索引 `.iqidx`（Chrom -> Gene CSR -> Transcript CSR -> Exon/Intron/SpliceSite 排序 `int32` 区间数组 + 字符串表）。
     - **零 CoW 复制**：各 worker 通过 `mmap(PROT_READ, MAP_SHARED)` 直接将 `.iqidx` 映射为只读 C `memoryview` / NumPy 视图。因为内部不含任何 `PyObject*`（无 `Py_INCREF` 修改），Linux 内核物理内存页**在 32 个并发 worker 间绝对共享（PSS 按 1 份分摊！）**。
     - 彻底消除 `gffutils` 运行时 SQLite 查询，并内置兼容性 `.db` 读取转换器（仍支持传入预建 `.db` 转 `.iqidx`）。
  2. **支柱二：`O(Locus)` 流式窗口 + SoA 紧凑内存池（消灭亿级对象与跨进程 Pickling）**
     - 对有序 BAM（支持自动检测 `.bai`/`.csi` 并发区间切分），每个 worker 直接负责指定基因座群组（Locus Chunk），在 Cython 结构数组（Structure-of-Arrays, SoA：连续 `int32`/`uint8`/`uint32` 缓冲）中完成：CIGAR block 解析 -> 基因区间二分定位 -> 外显子/内含子位图剖面（Bitset/Array 替代 Python `list[int]`）-> Read-to-Isoform 分配与事件编码。
     - 一旦穿过读取间隔间隙（Locus 完结），立即就地执行末端直方图累计、局部 IntronGraph 转录本重建与基因计数归约，**将高度压缩的最终结果块（计数增量、GTF 记录、`read_info` 列分片）流式刷盘，立刻重用预分配 Arena 缓冲区**！
     - **内存复杂度**：峰值单进程内存从 $O(\text{全染色体读段对象数})$ 降为 $O(\text{单个基因座读段窗口 + 固定 I/O 缓冲})$，跨进程零大型对象 `pickle` 回传！
  3. **支柱三：Cython 3.2 `nogil` 热点微内核（同时消除缺失依赖与百倍加速内循环）**
     - `cy_fasta_cigar.pyx`：HTSlib/pysam 对接或高速 CIGAR 区间切片 + polyA/TSS 末端修剪；
     - `cy_profile_assign.pyx`：排序双指针/二分定位实现外显子/内含子/单外显子特征匹配（替换 `OverlappingFeaturesProfileConstructor`、`JunctionComparator` 和 `GeneInfo` 的 $O(N^2)$ 与重复查询）；
     - `cy_tree_eval.pyx`：`model_polya.json` / `model_tss.json` 在初始化时解析为连续扁平数组（总计 <=2,637 个节点，占用 < 60 KB L1/L2 Cache），`nogil` 状态下执行 IEEE-754 `float32` 精确 `<` 决策树遍历 + sigmoid（替代 `xgboost` + `pd.DataFrame`）；
     - `cy_hist_stats.pyx`：对 `(positions, counts)` 稀疏直方图单遍计算加权均值、方差、偏度、众数、熵与平滑寻峰（彻底取消 `[pos] * count` 列表构造）；
     - `cy_barcode_umi.pyx`：`pysam.FastxFile` 批读字符串转为 ASCII `const uint8_t[::1]` 视图后在 `nogil` 下运行 Myers/Affine-Gap Smith-Waterman Barcode 对齐（替代 `ssw-py`）及带界 bit-parallel UMI 编辑距离（替代 `editdistance`），结合 Python 3.13/3.11 兼容 `SharedMemory`（含 3.13+ `track=False` 保护与 `np.uint64` k-mer 位操作）。
  4. **支柱四：并发安全与无锁归约（Forkserver/Spawn + Thread-Local Shards）**
     - 显式禁用第三方 C 库的隐式 OpenMP 多线程（`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`），消除 `fork`/并发超卖与死锁。
     - 若需全局去重/跨基因多映射解析（Multimap Resolver）：只通过共享 Compact Int IDs (`uint32 read_uid, uint32 gene_uid, int32 score`) 做轻量基数/排序归约，避免把整个 `BasicReadAssignment` 送回主进程。
  5. **支柱五：Polars / SciPy Out-of-Core 定量与矩阵导出引擎**
     - 多样本合并、`TerminalCounter._write_grouped`（用 `np.searchsorted` 替换 `iterrows()`）、`convert_grouped_counts`（`scipy.sparse.coo_matrix` + `mmwrite`）及 `.allinfo`/`read_info` 分表切分交由 `Polars` 惰性流式引擎或 SciPy 稀疏库处理，防止百万级行表爆内存。
- **预计收益（设计目标，待 G1~G4 实测校验）**：
  - **峰值进程树 `PSS`**：降低 **60%~85%**（尤其是 `gencode.v46` 注释内存由每进程各自持有变为 1 份内核只读页；千万 reads 不再全量驻留堆对象）；
  - **端到端 `wall_time`**：**3x~10x** 加速（单细胞 barcode/UMI、PolyA XGBoost 峰预测、大规模 MTX/表格导出环节可达更高的局部加速）。
- **可落地性**：**100% 在现有 `sqanti3` 环境开箱即用，零需新增系统软件或包安装！**

### 方案 C：Numba `@njit` JIT + Arrow IPC（`scanpy` / `PTSDsnDLPFC` 思路）
- **核心思路**：用 Numba `@njit` / `@njit(parallel=True)` 替代 Cython 编写数值热循环。
- **针对本集群的阻塞分析**：
  - 实测 `scanpy` 环境（含 `numba 0.64.0`）**缺少 `pysam`、`gffutils`、`pyfaidx`、`minimap2`、`samtools`**；
  - `PTSDsnDLPFC` 环境虽同时有 `pysam` 与 `numba 0.57.1`，但 Python 版本为 **`3.9.16`** 且缺少 `gffutils`、`minimap2`、`polars`；
  - 此外，Numba 启动存在 LLVM JIT 首次编译冷启动开销，且 `@njit(parallel=True)` 引入 TBB/OpenMP 并发层，易与多进程 worker 发生线程过载竞争；Cython `.so` 编译后启动 <1ms，且能直接无损耗操纵 `mmap` 字节指针与 C 结构体数组。
- **结论**：不作为主干架构；仅保持数学函数契约纯正，必要时不阻断纯 NumPy 回退。

### 方案 D：Rust (`PyO3`/`maturin`) 或全 C++20 重写底层引擎 + Python 编排层
- **核心思路**：完全用 Rust (`rust-htslib`/`noodles` + `rayon`) 或 C++20 (`htslib` + `pybind11`) 重写从 BAM 解码到 IntronGraph/定量的全管线，仅暴露 Python CLI 包装。
- **针对本集群的客观评估**：
  - 实测本集群 **`which cargo rustc` 不可用（无 Rust 工具链）**，且新建 `rust`/自定义 `pybind11` Conda 环境违反最小环境干预原则。
  - 若未来需突破 Python `IntronGraph` 或 HTSlib BGZF 解压极限，Scheme B 下的 Cython 3.2 原生支持直接链接 C++17 (`# distutils: language = c++` + `/usr/bin/g++ 13.3` + `libcpp` 向量/哈希表) 与 HTSlib C 头文件，即可平滑演进至近似全 C++ 性能，而无需引入跨语言 FFI 工程复杂性或新工具链！

---

## 3. 方案综合横向评分与最终选型

| 维度 | 方案 A (纯列式 Python) | **方案 B (`mmap` CSR + `O(Locus)` 流式 + Cython `nogil`)** | 方案 C (Numba + Arrow) | 方案 D (Rust/C++ 全重写) |
| :--- | :--- | :--- | :--- | :--- |
| **集群现有环境兼容性** | 可用 (`sqanti3`) | **极佳（`sqanti3` 全包+`Cython 3.2.4`+`gcc 13.3` 就绪）** | 差（缺统一 env 覆盖全依赖） | 受阻（无 `rustc/cargo`） |
| **缺失 4 包解决度** | 部分受阻 (`ssw`/`editdistance` 慢) | **彻底解决（原生 Cython 内核 + `minimap2` CLI）** | 部分受阻（缺 `pysam` 联合 env） | 需全部手写 |
| **多核内存共享（降 PSS）** | 中等（Arrow 零拷贝，但业务对象仍在堆内） | **极高（`mmap .iqidx` 零 CoW + `O(Locus)` 有界内存池）** | 中等（JIT 缓存占每进程内存） | 极高（单进程多线程模型） |
| **运行时间提速潜力** | 中（2~4x） | **极高（3~10x 端到端，微内核 10~50x）** | 高（但冷启动与线程冲突） | 极限（5~15x） |
| **生物逻辑精确等价性风险** | 低 | **低~中（逐节点对比 + Cython float32 位级复现）** | 中（Numba 类型推导限制） | 高（整库重造难逐节点等价验证） |

> **选型结论**：正式采用 **方案 B（核心）+ 方案 A 生态组合（Polars/PyArrow/SciPy 用于下游大规模表与 MTX 导出）** 作为默认重构架构。
> 所有 Cython 模块必须同步维护清晰的纯 NumPy/Python 参照实现或回退路径（用于单元测试验证“向量化/Cython 内核与原版逻辑逐位一致”）。

---

## 4. 外部依赖包官方规范在低内存架构中的落点

整合重构前对所有依赖包/工具官方手册与版本特性的深度调研，各环节刚性约束如下：

1. **`pysam` (`0.22+` / HTSlib)**：
   - 全部 `AlignmentFile`、`FastaFile`、`FastxFile` 强制 `with` 上下文；
   - 超越 $2^{29}-1$ (~536Mb) 染色体自动回退 `.csi` 索引 (`samtools index -c`)；
   - 使用 **`pysam.FastxFile`** 替代 `Bio.SeqIO.parse` 流式解压 `.fq.gz`/`.fa.gz`，省去每条 read 的 `SeqRecord` 与注解字典创建。
2. **`gffutils` (`0.14`)**：
   - 仅用于兼容用户提供的旧 `.db` 文件导入与转换为 `.iqidx`，及可选构建 `--complete_genedb` SQLite；不再作为运行期高并发 worker 热点查询源；
   - 使用时启用 `db.seqids()`（废弃 `all_features()` 全扫描）、WAL 只读独立连接（绝不跨进程传连接对象）与 `disable_infer_genes=True, disable_infer_transcripts=True` 优化建库参数。
3. **`pyfaidx` (`0.9`) vs `pysam.FastaFile`**：
   - 统一封装基因序列读取器，优先以 worker 内惰性只读打开 `pysam.FastaFile`（原生直读 `.fa.gz` bgzip `.fai`+`.gzi`，线程/进程隔离），清理 `fusion_detector.py` 与其他模块混用两种 FASTA 库的问题，并移除未使用的重复索引状态。
4. **`NumPy` (`1.26` -> `2.4+` NEP 50)**：
   - 严格遵从 NumPy 2.0 标量提升规范，k-mer 位移 (`np.uint64(1) << np.uint64(2 * k)`) 与掩码运算禁止混用 Python 有符号 `int`，提前通过 `basicpy` (`numpy 2.4.6`) 验证零溢出或退化。
5. **`Pandas` (`2.2` -> `3.0+` Copy-on-Write)**：
   - 严格遵循 CoW（废弃链式赋值 `df[col][mask] = ...`，统一使用 `.loc`）；从 Node 5/7 的逐行热路径（`iterrows()`）彻底移除 Pandas，仅在最终导出边界使用向量化 Pandas / Polars。
6. **`SciPy` (`1.11`~`1.17`)**：
   - `scipy.stats.mode` 对齐现代 `keepdims=False` 语义；MTX 导出直接使用 `scipy.sparse.coo_matrix` + `scipy.io.mmwrite`，替换嵌套字符串拼接循环；Peak 统计完全由稀疏加权直方图在 Cython/NumPy 中直接算完。
7. **`XGBoost JSON` (`model_polya.json`, `model_tss.json`)**：
   - 调研证实两模型均为标准 `gbtree`（无 DART），417 棵树共 <=2,637 节点，8 个 `float32` 特征；采用 Cython/NumPy 无锁只读树阵直接预测并设置主进程环境变量禁用 OpenMP/BLAS 多线程超卖。
8. **`IntervalTree` -> 静态有序区间与扫描线**：
   - 运行期已知外显子/内含子/融合基因区间查询全部转为静态升序 C 数组 + 二分搜索/扫描线（$O(\log N)$ 或均摊 $O(1)$，零 Python 对象指针头开销），替换 `GeneInfo` $O(N^2)$ 与重量级 `IntervalTree` 内存树。
9. **Python `multiprocessing` 与 `SharedMemory`**：
   - 并行运行器默认安全使用 `forkserver` / `spawn`（避开裸 `fork` 与原生库状态冲突）；
   - 共享 k-mer 索引块（Node 2）按 Python 版本适配 `SharedMemory(..., track=False)`（在 Python 3.13+ 防止 `resource_tracker` 在子进程退出时误释放或报错）并用显式生命周期管理器替代 `__del__`；
   - 注释索引（Node 1）使用只读文件 `mmap`（由内核按需加载页并多进程物理共享），从根本上免除 POSIX `shm` 配额溢出（`/dev/shm` 空间受限导致的大集群崩溃问题）。
10. **外部 CLI 工具（`minimap2` / `samtools` / 死代码清除）**：
    - `shutil.which` + 语义版本元组校验；`hashlib.blake2b` 替代随机化 `hash()` 做比对缓存指纹；`contextlib.ExitStack` 管理双向管道防 `SIGPIPE`/死锁；
    - 彻底移除未声明且冗余的 `pybedtools`（Node 5 改内置轻量 BED + 二分搜索区间表）与不可达的 `h5py` 死代码（Node 2）。
