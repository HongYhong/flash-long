# IsoQuant Refactored (`isoquant_next`)

> **基线来源**：IsoQuant `4.0.0`（上游目录 `/gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/`，保持完全只读作为行为与结果对照基线）
> **重构目标目录**：`isoquant_refactored/`（`src/isoquant_next/` + `tests/` + `scripts/benchmark/` + `docs/`）
> **核心目标**：按 **Node 0 ~ Node 10** 分节点彻底重写架构、消除已知缺陷与多进程并发隐患，并**大幅降低峰值运行内存（PSS）与端到端运行时间**，同时严格复用集群既有 Conda 环境（不新建或修改 Conda 环境）。

---

## 文档导航（Documentation Suite）

| 文档 | 用途 | 核心内容 |
| :--- | :--- | :--- |
| **[docs/spec.md](docs/spec.md)** | 工程规格与验收契约 | 重构范围、现有 Conda 环境与缺失包审计、PSS/Wall-time 观测规范、`longpass_benchmark250209/` 分层测试集、Node 0~10 输入输出契约与等价性门禁 |
| **[docs/optimization_schemes.md](docs/optimization_schemes.md)** | 低内存与高性能方案对比 | 原版六大性能/内存瓶颈剖析、4 套架构方案（Scheme A~D）对比、推荐方案设计细节（`mmap` 零拷贝 CSR 索引 + `O(Locus)` 流式窗口 + Cython `nogil` 内核 + 零 XGBoost 树引擎）与依赖包官方规范落点 |
| **[docs/roadmap.md](docs/roadmap.md)** | 分阶段路线图与放行门禁 | Phase 0~5 里程碑（M0~M10）及性能基准门禁（G0~G4） |
| **[docs/plan.md](docs/plan.md)** | 分节点技术实施计划 | 包骨架设计、Node 0 ~ Node 10 逐节点重写方案、Bug 修复清单、测试集与基准评测执行设计 |
| **[docs/task.md](docs/task.md)** | 可追踪任务清单 | `T00` ~ `T28` 细粒度任务表、状态跟踪、前置依赖与验收命令 |
| **[REFACTORING_PLAN.md](REFACTORING_PLAN.md)** | 执行摘要与索引入口 | 指向上述多文档规划体系的索引入口 |

---

## 运行环境与数据基准速览

- **主开发与基准测试环境**：`sqanti3`（`Python 3.11.15`, `pysam 0.22.1`, `gffutils 0.14`, `pyfaidx 0.9.0.4`, `Cython 3.2.4`, `polars 0.20.31`, `pyarrow 14.0.2`, `numpy 1.26.4`, `pandas 2.2.3`, `scipy 1.11.4`, `scikit-learn 1.5.2`, `Bio 1.81`, 系统 `gcc 13.3.0`，含 CLI `minimap2`/`STAR`/`samtools`/`bedtools`）。
- **前瞻兼容验证环境**：`basicpy`（`Python 3.11.15`, `numpy 2.4.6` NEP 50, `pandas 3.0.2` Copy-on-Write, `pysam 0.23.3`, `scipy 1.17.1`）。
- **集群环境纪律（遵循 `~/.codex/AGENTS.md`）**：
  - 所有命令通过显式激活进入现有环境：`source /gpfs/group/data_sci_group/weiwu_group/hongyanhong/miniconda3/miniconda3/etc/profile.d/conda.sh && conda activate sqanti3`；
  - 严禁使用裸系统 `python3`，严禁擅自新建、修改、安装或升级 Conda 环境。
- **数据来源**：
  - 单元测试与快速功能回归：`/gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_tests/{simple_data,toy_data,short_reads_toy}/`
  - 大规模真实与模拟基准评测：`/gpfs/group/data_sci_group/weiwu_group/hongyanhong/backup/longpass_benchmark250209/`（覆盖 6 类长读长协议、GRCh38/GRCm39/SIRV/RBMX_KD/CAGEPAS/shortreadSJ）。

---

## 下一步行动边界

当前处于**规划文档确认阶段**。待您确认**主推优化方案（推荐 `Scheme B`）**、**基准运行环境（`sqanti3`）**及**Node 0 前架构约束**后，再依次从 `T01`（工程骨架与基线摸底）-> `Node 0` 开始编码实现。
