# `bambu_next`: 并行与低内存全链路优化版 Bambu 工程文档与实施方案

> **状态**：**Phase 0–6 (`T00–T23` / `Gates G0–G5`) 全链路重构与无上限 Slurm Array 验收已全部完成（56/56 `testthat` PASS，29 组主流队列 `r=1.000000` / 100% GTF 结构一致，详见 [docs/benchmark_report.md](docs/benchmark_report.md) 与 [docs/bambu_next_benchmark_report.pdf](docs/bambu_next_benchmark_report.pdf)）。**  
> **基准基线**：原版 Bioconductor `bambu` (`3.12.1` 安装包 / `3.13.1` 只读源码树 [`bambu/`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu/)) 与历史最佳可用版本 [`bambu_refactor/`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu_refactor/) (`G8`)。  
> **运行环境**：`envs/bambu` (`R 4.5.3`, `bambu 3.12.1`, `BiocParallel 1.44.0`, `Rhtslib 3.6.0`, `GenomicAlignments 1.46.0`, `data.table 1.17.8`, `xgboost 3.2.0.1`) 与 `envs/basicpy`（统计汇总与绘图）。

---

## 1. 为什么在 `G8` 之后新开独立目录 `bambu_next/`

复盘先前的重构与优化轨迹（`019fe711-da4e-7b40-9f87-9e5a8c33c1ba`）：
1. **`G0–G8`（[`bambu_refactor/`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu_refactor/)）** 完成了模块化拆分与可观测性，但在 `human_all_d1q1`（16 个 human BAM，约 `96.2M` mapped alignments，`8 workers`）上仅取得 `4763s -> 4610s`、`Peak PSS 250.48 GiB -> 254.29 GiB`（略增）的微弱变化；单 BAM `S1.04` 并行虽然将 `ENCFF089AMF` 从 `338.82s` 压至 `92.44s`，但在 16 BAM 并发场景因切换回 `8 sample workers × 1 reader` 而退化回原生速度。
2. **`G9–G11`（[`bambu_production/`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu_production/)）**、`G12–G14`（[`bambu_lowcopy/`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu_lowcopy/)）与 `G15–G17`（[`bambu_s3core/`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/bambu_s3core/)）均因保留了 R S4 `GRangesList`/`IRangesList`/`Hits` 巨型对象扩展、保留了 `MulticoreParam` `fork()` + R 分代 GC 触发 Copy-on-Write (CoW) 页分裂、且未能协同优化 S1/S2/S3/S4 全阶段，未能达标（例如 `G17 candidate` 在 `human_all_d1q1` 上耗时 `4585.6s`，Peak PSS 反而升至 `266.57 GiB`）。
3. 参照 [`isoquant_refactored/`](file:///gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/) 的工程范式，结合对 Bambu 全部核心 Bioconductor 依赖包官方手册与 C 头文件的系统审计，本目录从**数据结构与并行模型根部**重新设计：
   - **C++17 + `Rhtslib` + `IRanges/S4Vectors` C 接口 + `RcppArmadillo` 零拷贝扁平 CSR 引擎**；
   - 消灭 `MulticoreParam` 多进程持有巨型 R S4 堆导致的 `250+ GiB` CoW 爆炸；
   - 全链路重构 **S0 注释编译 (`GTF/TxDb`) → S1 (`Rhtslib` BAM + 剪接矩阵 + ReadClass 构建/评分) → S2 (跨样本合并 + 增量注释扩展) → S3 (`CSR` 距离匹配 + Fused `equiRC/EM`) → S4 (流式 GTF/counts/SummarizedExperiment 导出)**。

---

## 2. 文档导航

| 文档 | 内容说明 |
| :--- | :--- |
| [docs/bioconductor_manuals_audit.md](docs/bioconductor_manuals_audit.md) | `envs/bambu` 下 15 个 Bioconductor/R 依赖包手册（`.Rmd`/`.Rnw`/`.Rd`）与 C 头文件深度审计及针对 Bambu 的工程落点 |
| [docs/spec.md](docs/spec.md) | 工程规格、只读保护边界、Conda/本地库隔离策略、冻结 `T_core` (18 BAM) + 扩展 `G_ext` 测试矩阵、等价性与性能/PSS 验收契约 |
| [docs/optimization_schemes.md](docs/optimization_schemes.md) | `G8–G17` 停滞根因量化解剖、三套新架构候选方案对比，以及选定的 **方案 B（C++17/`Rhtslib`/Flat-CSR 零 CoW 引擎 + Bambu S4 Drop-in 适配器）** 详细剖析 |
| [docs/roadmap.md](docs/roadmap.md) | `Phase 0 ~ Phase 6`、`Milestones M0_PRE ~ M10` 与端到端性能/精度验收门禁 `Gate G0 ~ G5` |
| [docs/plan.md](docs/plan.md) | 目标包结构骨架（`bambu_next/`）与逐节点 `Node 0 ~ Node 9` 代码级实施详案 |
| [docs/task.md](docs/task.md) | `T00 ~ T23` 可追踪任务台账、优先级与依赖拓扑 |

---

## 3. 只读与安全铁律

- **严禁改动**：`long_bench/bambu/`、`bambu_stage1_parallel/`、`bambu_refactor/`、`bambu_production/`、`bambu_lowcopy/`、`bambu_s3core/`、`long_bench/docs/`、`long_bench/results/`、`miniconda3/envs/*` 以及 `/gpfs/.../backup/longpass_benchmark250209/` 下任何文件。
- **禁止写 `/tmp`**：所有编译临时文件、R/Slurm `TMPDIR`、scratch 与大型中间产出仅放在 `/gpfs/group/data_sci_group/weiwu_group/hongyanhong/` 下专用工作目录内，验收完即清大临时目录。
- **本地安装隔离**：`bambu_next` 通过 `R CMD INSTALL --no-multiarch -l <project_local_lib>` 安装至项目自有隔离库，绝不覆盖共享 Conda 环境。
