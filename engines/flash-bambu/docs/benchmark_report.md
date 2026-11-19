# `bambu_next` (`bambu 3.13.99`) 全链路重构与 Slurm 验收测试报告

> 生成时间：2026-11-19 07:50:38 | 方法：Strict `A-B-B-A` ×3 配对集群直测（`monitor_process_tree.sh` 1 Hz 进程树 PSS 采样）+ G17 官方生产基线对照

## 1. 核心摘要与关键成果

- **队列总体覆盖**：完成 Gate G1–G4 共 **29** 个 ABBA ×3 主对比任务（每项 6 轮独立运行）、6 个样本 `ncore=1..32` 扩展曲线、6 个 Flat-CSR 冷/热缓存对比、`gencode.v46` (1.6 GB) S0 GTF 基准、以及 7 个 Bambu 可选 API 参数开关专项校验，全部零 OOM、零进程死锁。
- **端到端加速区间**：全队列核心运行加速比 **5.90x – 28.43x**（中位数 **16.23x**）。
- **16-BAM `human_all_d1q1` (`nc=16`, ~96M alignments)**：核心耗时 `2055.46s -> 256.93s` (**`8.00x`**)，全进程树峰值 PSS **`79.16 GiB -> 6.99 GiB`** (**`11.33x`** 压缩，仅为原生 `8.8%`)。
- **16-BAM `human_all_d1q1` (`nc=8`)**：核心耗时 `2490.54s -> 273.33s` (**`9.11x`**)，峰值 PSS **`49.39 GiB -> 6.61 GiB`** (**`7.47x`** 压缩)。
- **G17 严格同口径 (`NDR=1, trackReads=TRUE`, 468,933 tx)**：`G17 Native` `4779s / 250.48 GiB` (`G8` `4626s / 254.29 GiB`) $\rightarrow$ `bambu_next` **`733.68s` core (`751.08s` tot) / `44.60 GiB`** (**`6.51x` 提速, `5.62x` 内存压缩**)。
- **`SIRV_simulation` 真值验证（8,795,698 真值模拟读段）**：`902.39s -> 31.74s` (**`28.43x`**)，Recall/Precision/F1（`0.879261 / 0.894244 / 0.886689`）与原生 `bambu 3.12.1` **6 位小数完全一致**，5 份导出表格 `max_abs_diff = 0`。
- **L1/L2 等价性**：**23/29** 个主队列实现 **100% GTF 坐标与基因归属一致**，全部 **29/29** 队列转录本定量 Pearson $r \ge 0.999996$（其中单/双样本及协议测试 **`25/29` 达到 `max_diff <= 1e-5` 严格全通**；`human_pair` / `RBMX_KD` 仅 2–3 格因 `conv=0.01` 第 99/100 步 EM 分界出现 `0.00019` 波动，`r=1.000000`；`>=3` 样本 `16-BAM` 队列因原生内部调用随机 `sample()` 重排样本合并顺序导致 `BambuTx<序号>` 置换，剔除序号后 3.52M / 1.85M 行 GTF 坐标结构 100% 零差异）。

## 2. 全队列汇总表 (G1–G4 ABBA ×3 中位数)

| Gate | Case ID | ncore | Ref Core | Next Core | 加速比 | Ref PSS | Next PSS | 内存压缩 | 结构GTF一致 | Pearson r |
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | :---: | ---: |
| G1 | `g1_s1_ENCFF089AMF_nc16` | 16 | 546.00s | 28.71s | **19.02x** | 12.30 GiB | 2.85 GiB | 4.31x | S1 L1 | S1 L1 exact |
| G1 | `g1_s1_ENCFF400BQQ_nc16` | 16 | 36.46s | 6.18s | **5.90x** | 1.63 GiB | 1.24 GiB | 1.31x | S1 L1 | S1 L1 exact |
| G1 | `g1_s1_ENCFF436GKZ_nc16` | 16 | 98.93s | 12.10s | **8.18x** | 2.65 GiB | 1.67 GiB | 1.58x | S1 L1 | S1 L1 exact |
| G1 | `g1_s1_ENCFF598YQO_nc16` | 16 | 391.91s | 24.39s | **16.07x** | 10.10 GiB | 2.26 GiB | 4.46x | S1 L1 | S1 L1 exact |
| G1 | `g1_s1_ENCFF846RID_nc16` | 16 | 103.04s | 11.34s | **9.09x** | 2.64 GiB | 1.51 GiB | 1.75x | S1 L1 | S1 L1 exact |
| G1 | `g1_s1_ENCFF887KHY_nc16` | 16 | 575.05s | 31.57s | **18.21x** | 12.84 GiB | 2.96 GiB | 4.35x | S1 L1 | S1 L1 exact |
| G2 | `human_pair_d1q1` | 8 | 332.74s | 21.60s | **15.40x** | 9.57 GiB | 2.18 GiB | 4.40x | YES | 1.0000000 |
| G2 | `human_single_d1q1` | 8 | 319.70s | 19.34s | **16.53x** | 4.24 GiB | 2.05 GiB | 2.07x | YES | 1.0000000 |
| G2 | `mouse_all_d1q1` | 8 | 806.37s | 55.91s | **14.42x** | 18.28 GiB | 2.76 GiB | 6.63x | YES | 1.0000000 |
| G2 | `mouse_all_d1q1_nc16` | 16 | 803.91s | 53.45s | **15.04x** | 18.22 GiB | 3.10 GiB | 5.88x | YES | 1.0000000 |
| G2 | `mouse_single_d1q1` | 8 | 681.49s | 32.36s | **21.06x** | 10.10 GiB | 2.30 GiB | 4.39x | YES | 1.0000000 |
| G3 | `human_all_d1q1_nc16` | 16 | 2055.46s | 256.93s | **8.00x** | 79.16 GiB | 6.99 GiB | 11.33x | YES | 1.0000000 |
| G3 | `human_all_d1q1_nc8` | 8 | 2490.54s | 273.33s | **9.11x** | 49.39 GiB | 6.61 GiB | 7.47x | YES | 1.0000000 |
| G4 | `g4_RBMX_KD_v46_pair` | 16 | 479.34s | 42.06s | **11.40x** | 12.60 GiB | 2.84 GiB | 4.44x | YES | 1.0000000 |
| G4 | `g4_SIRV_simulated1` | 16 | 902.39s | 31.74s | **28.43x** | 14.01 GiB | 2.97 GiB | 4.71x | YES | 1.0000000 |
| G4 | `g4_hum_R2C2_ONT_ENCFF887KHY` | 8 | 1009.03s | 42.58s | **23.70x** | 13.77 GiB | 2.81 GiB | 4.90x | YES | 1.0000000 |
| G4 | `g4_hum_cDNA_ONT_ENCFF714YOZ` | 8 | 783.59s | 34.10s | **22.98x** | 11.35 GiB | 2.79 GiB | 4.07x | YES | 1.0000000 |
| G4 | `g4_hum_cDNA_pacbio_ENCFF400BQQ` | 8 | 140.15s | 13.74s | **10.20x** | 3.32 GiB | 1.72 GiB | 1.93x | YES | 1.0000000 |
| G4 | `g4_hum_capped_ONT_ENCFF089AMF` | 8 | 909.13s | 39.16s | **23.22x** | 14.52 GiB | 2.68 GiB | 5.42x | YES | 1.0000000 |
| G4 | `g4_hum_capped_pacbio_ENCFF684YOO` | 8 | 164.69s | 12.92s | **12.74x** | 3.57 GiB | 1.76 GiB | 2.03x | YES | 1.0000000 |
| G4 | `g4_hum_directRNA_ENCFF373LMV` | 8 | 78.71s | 8.85s | **8.90x** | 2.99 GiB | 1.68 GiB | 1.78x | YES | 1.0000000 |
| G4 | `g4_manatee_cDNA_ONT_denovo` | 16 | 1502.17s | 58.75s | **25.57x** | 14.92 GiB | 4.08 GiB | 3.66x | YES | 1.0000000 |
| G4 | `g4_manatee_pacbio_denovo` | 16 | 187.83s | 9.91s | **18.96x** | 3.23 GiB | 1.48 GiB | 2.18x | YES | 1.0000000 |
| G4 | `g4_mus_R2C2_ONT_ENCFF598YQO` | 8 | 690.07s | 32.84s | **21.01x** | 10.08 GiB | 2.26 GiB | 4.46x | YES | 1.0000000 |
| G4 | `g4_mus_cDNA_ONT_ENCFF232YSU` | 8 | 667.90s | 31.03s | **21.53x** | 7.95 GiB | 2.44 GiB | 3.26x | YES | 1.0000000 |
| G4 | `g4_mus_cDNA_pacbio_ENCFF313VYZ` | 8 | 282.06s | 17.38s | **16.23x** | 4.62 GiB | 1.86 GiB | 2.49x | YES | 1.0000000 |
| G4 | `g4_mus_capped_ONT_ENCFF056EOI` | 8 | 802.34s | 36.44s | **22.02x** | 14.32 GiB | 2.60 GiB | 5.52x | YES | 1.0000000 |
| G4 | `g4_mus_capped_pacbio_ENCFF310IPO` | 8 | 180.84s | 11.87s | **15.23x** | 3.42 GiB | 1.48 GiB | 2.31x | YES | 1.0000000 |
| G4 | `g4_mus_directRNA_ENCFF349BIN` | 8 | 293.43s | 17.00s | **17.26x** | 3.67 GiB | 1.66 GiB | 2.22x | YES | 1.0000000 |
