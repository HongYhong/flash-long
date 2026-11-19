# IsoQuant 重构计划总览（Legacy Index -> Multi-Document Suite）

> **说明**：根据多文档规划要求与集群环境、基准数据集、`~/.codex/AGENTS.md` 规范深化调研结果，原单文件重构计划已升级并拆分为结构化多文档体系（位于 [isoquant_refactored/docs/](docs/)），请直接查阅对应专项文档：

1. **[README.md](README.md)** — 文档导航、环境与数据入口总览
2. **[docs/spec.md](docs/spec.md)** — 工程规格与验收契约（含现有 Conda 环境 `sqanti3`/`chorus`/`basicpy` 实测审计、缺失 4 包无缺包替代方案、`/gpfs/group/data_sci_group/weiwu_group/hongyanhong/backup/longpass_benchmark250209/` Tier 0~3 测试集清单与进程树 `PSS` 评测口径）
3. **[docs/optimization_schemes.md](docs/optimization_schemes.md)** — 低内存与高性能架构方案对比（剖析原版 6 大性能/内存瓶颈，深入对比 4 套方案并论证主推 **方案 B：`mmap` 零拷贝 CSR 注释索引 + `O(Locus)` 流式分片 + Cython `nogil` 原生内核 + 零 XGBoost 树引擎**）
4. **[docs/roadmap.md](docs/roadmap.md)** — 分阶段路线图（Phase 0~5、里程碑 M0~M10 与性能验收门禁 G0~G4）
5. **[docs/plan.md](docs/plan.md)** — 分节点技术实施详案（Node 0 ~ Node 10 模块结构、Bug 修复落点、核心内核设计与对照策略）
6. **[docs/task.md](docs/task.md)** — 可追踪任务清单（`T00` ~ `T28`）
