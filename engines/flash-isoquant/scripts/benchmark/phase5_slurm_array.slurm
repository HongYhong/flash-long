#!/usr/bin/env bash
#SBATCH -J iq_p5_bench
#SBATCH -p chj.cn.1,chj.cn.2,chj.cn.3,chj.cn.4,chj.cp.1
#SBATCH -N 1
#SBATCH --cpus-per-task=48
#SBATCH --mem=320G
#SBATCH --array=0-16
#SBATCH -o /gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/benchmarks/slurm_logs/p5_%A_%a.out
#SBATCH -e /gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/benchmarks/slurm_logs/p5_%A_%a.err

set -euo pipefail

CONDA_SH="/gpfs/group/data_sci_group/weiwu_group/hongyanhong/miniconda3/miniconda3/etc/profile.d/conda.sh"
source "${CONDA_SH}"
conda activate sqanti3

REPO_ROOT="/gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant"
BENCH_ROOT="${REPO_ROOT}/isoquant_refactored"
DATA_ROOT="/gpfs/group/data_sci_group/weiwu_group/hongyanhong/backup/longpass_benchmark250209"
GENOME_DIR="${DATA_ROOT}/genome"
V39_GTF="${GENOME_DIR}/lrgasp_gencode_v39_annotation_sirvs.human.gtf"
V46_GTF="${GENOME_DIR}/gencode.v46.annotation.gtf"
VM28_GTF="${GENOME_DIR}/lrgasp_gencode_vM28_sirvs.mouse.gtf"
GRCH38_FA="${GENOME_DIR}/GRCh38.p14.genome.fa"
GRCM39_FA="${GENOME_DIR}/GRCm39.genome.fa"
CAGE_BED="${GENOME_DIR}/CAGEPAS/CAGE_WTC11.all_reps_v2.bed"

TASK_ID="${SLURM_ARRAY_TASK_ID:-${1:-0}}"
BASE_SCRATCH="/gpfs/group/data_sci_group/weiwu_group/hongyanhong/iq_p5"
SCRATCH="${BASE_SCRATCH}/j${SLURM_JOB_ID:-0}_${TASK_ID}"
SHORT_TMP="${BASE_SCRATCH}/t${TASK_ID}"
mkdir -p "${SCRATCH}/mpl" "${SHORT_TMP}"
export TMPDIR="${SHORT_TMP}"
export TMP="${SHORT_TMP}"
export TEMP="${SHORT_TMP}"
export MPLCONFIGDIR="${SCRATCH}/mpl"

cleanup() {
    rm -rf "${SCRATCH}" "${SHORT_TMP}"
    rmdir "${BASE_SCRATCH}" 2>/dev/null || true
}
trap cleanup EXIT

run_pair() {
    local slug="$1"
    local repeats="$2"
    local prefix="$3"
    local common_orig="$4"
    local common_next="$5"
    local log_dir="${BENCH_ROOT}/benchmarks/${slug}"
    mkdir -p "${log_dir}"

    local orig_cmd="source '${CONDA_SH}' && conda activate sqanti3 && cd '${REPO_ROOT}' && rm -rf '${SCRATCH}/orig' && python isoquant.py ${common_orig} -o '${SCRATCH}/orig'"
    local next_cmd="source '${CONDA_SH}' && conda activate sqanti3 && cd '${REPO_ROOT}' && rm -rf '${SCRATCH}/next' && isoquant-next ${common_next} -o '${SCRATCH}/next'"

    python "${BENCH_ROOT}/scripts/benchmark/run_tiered_benchmark.py" \
        --repeats "${repeats}" \
        --interval 1.0 \
        --log-dir "${log_dir}" \
        --baseline-cmd "${orig_cmd}" \
        --candidate-cmd "${next_cmd}"

    cp -f "${SCRATCH}/orig/isoquant.log" "${log_dir}/orig.isoquant.log" 2>/dev/null || true
    cp -f "${SCRATCH}/next/isoquant.log" "${log_dir}/next.isoquant.log" 2>/dev/null || true
    cp -f "${SCRATCH}/next/run_manifest.json" "${log_dir}/next.run_manifest.json" 2>/dev/null || true

    if [[ "${prefix}" != "NONE" ]]; then
        python "${BENCH_ROOT}/scripts/benchmark/summarize_pair.py" \
            --ref-dir "${SCRATCH}/orig" \
            --cand-dir "${SCRATCH}/next" \
            --prefix "${prefix}" \
            --json-out "${log_dir}/accuracy_summary.json" || true
    fi
}

echo "[$(date -Is)] Starting Phase 5 task ${TASK_ID} on $(hostname)"

case "${TASK_ID}" in
    0)
        # Gate G1: Tier 1 ENCFF400BQQ alternating A-B-B-A (3 repeats) @ -t 8
        COMMON="--clean_start -d pacbio_ccs --fl_data --bam ${DATA_ROOT}/cDNA_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF400BQQ.sorted.k14.bam --genedb ${V39_GTF} --complete_genedb -r ${GRCH38_FA} -t 8 -p ENCFF400BQQ"
        run_pair "g1_encff400bqq_t8_abba3" 3 "ENCFF400BQQ" "${COMMON}" "${COMMON}"
        ;;
    1)
        # Gate G1: ENCFF400BQQ single-thread scaling @ -t 1
        COMMON="--clean_start -d pacbio_ccs --fl_data --bam ${DATA_ROOT}/cDNA_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF400BQQ.sorted.k14.bam --genedb ${V39_GTF} --complete_genedb -r ${GRCH38_FA} -t 1 -p ENCFF400BQQ"
        run_pair "g1_encff400bqq_t1" 1 "ENCFF400BQQ" "${COMMON}" "${COMMON}"
        ;;
    2)
        # Gate G1: ENCFF400BQQ 16-thread scaling @ -t 16
        COMMON="--clean_start -d pacbio_ccs --fl_data --bam ${DATA_ROOT}/cDNA_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF400BQQ.sorted.k14.bam --genedb ${V39_GTF} --complete_genedb -r ${GRCH38_FA} -t 16 -p ENCFF400BQQ"
        run_pair "g1_encff400bqq_t16" 1 "ENCFF400BQQ" "${COMMON}" "${COMMON}"
        ;;
    3)
        # Gate G1: ENCFF684YOO (capped PacBio, 887k reads) @ -t 16
        COMMON="--clean_start -d pacbio_ccs --fl_data --bam ${DATA_ROOT}/capped_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF684YOO.bam --genedb ${V39_GTF} --complete_genedb -r ${GRCH38_FA} -t 16 -p ENCFF684YOO"
        run_pair "g1_encff684yoo_t16" 1 "ENCFF684YOO" "${COMMON}" "${COMMON}"
        ;;
    4)
        # Gate G1: ENCFF339FMQ (capped PacBio, 1.28M reads) @ -t 16
        COMMON="--clean_start -d pacbio_ccs --fl_data --bam ${DATA_ROOT}/capped_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF339FMQ.bam --genedb ${V39_GTF} --complete_genedb -r ${GRCH38_FA} -t 16 -p ENCFF339FMQ"
        run_pair "g1_encff339fmq_t16" 1 "ENCFF339FMQ" "${COMMON}" "${COMMON}"
        ;;
    5)
        # Gate G1: Multi-sample YAML combine_counts (ENCFF684YOO + ENCFF339FMQ) @ -t 16
        YAML_PATH="${SCRATCH}/g1_multi.yaml"
        cat > "${YAML_PATH}" <<EOF
- data format: "bam"
- name: ENCFF684YOO
  long read files:
  - ${DATA_ROOT}/capped_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF684YOO.bam
- name: ENCFF339FMQ
  long read files:
  - ${DATA_ROOT}/capped_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF339FMQ.bam
EOF
        COMMON="--clean_start -d pacbio_ccs --fl_data --yaml ${YAML_PATH} --genedb ${V39_GTF} --complete_genedb -r ${GRCH38_FA} -t 16"
        run_pair "g1_multisample_yaml_t16" 1 "NONE" "${COMMON}" "${COMMON}"
        for P in ENCFF684YOO ENCFF339FMQ; do
            python "${BENCH_ROOT}/scripts/benchmark/summarize_pair.py" \
                --ref-dir "${SCRATCH}/orig" \
                --cand-dir "${SCRATCH}/next" \
                --prefix "${P}" \
                --json-out "${BENCH_ROOT}/benchmarks/g1_multisample_yaml_t16/accuracy_${P}.json" || true
        done
        ;;
    6)
        # Gate G2: Tier 2 ENCFF436GKZ (cDNA PacBio, 1.29M reads) @ -t 16
        COMMON="--clean_start -d pacbio_ccs --fl_data --bam ${DATA_ROOT}/cDNA_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF436GKZ.sorted.k14.bam --genedb ${V39_GTF} --complete_genedb -r ${GRCH38_FA} -t 16 -p ENCFF436GKZ"
        run_pair "g2_encff436gkz_t16" 1 "ENCFF436GKZ" "${COMMON}" "${COMMON}"
        ;;
    7)
        # Gate G2: Tier 2 ENCFF846RID (directRNA stranded, 1.25M reads; temp symlink+.bai outside backup/) @ -t 16
        mkdir -p "${SCRATCH}/links"
        ln -s "${DATA_ROOT}/directRNA_fulpass/human/H1_ES_cell/alignments/ENCFF846RID.UtoT.fastq.sorted.k14.bam" "${SCRATCH}/links/ENCFF846RID.bam"
        samtools index -@ 8 "${SCRATCH}/links/ENCFF846RID.bam"
        COMMON="--clean_start -d nanopore --stranded forward --fl_data --bam ${SCRATCH}/links/ENCFF846RID.bam --genedb ${V39_GTF} --complete_genedb -r ${GRCH38_FA} -t 16 -p ENCFF846RID"
        run_pair "g2_encff846rid_drna_t16" 1 "ENCFF846RID" "${COMMON}" "${COMMON}"
        ;;
    8)
        # Gate G2: Tier 2 ENCFF207KVO (directRNA replicate, 1.39M reads; temp symlink+.bai outside backup/) @ -t 16
        mkdir -p "${SCRATCH}/links"
        ln -s "${DATA_ROOT}/directRNA_fulpass/human/H1_ES_cell/alignments/ENCFF207KVO.UtoT.fastq.sorted.k14.bam" "${SCRATCH}/links/ENCFF207KVO.bam"
        samtools index -@ 8 "${SCRATCH}/links/ENCFF207KVO.bam"
        COMMON="--clean_start -d nanopore --stranded forward --fl_data --bam ${SCRATCH}/links/ENCFF207KVO.bam --genedb ${V39_GTF} --complete_genedb -r ${GRCH38_FA} -t 16 -p ENCFF207KVO"
        run_pair "g2_encff207kvo_drna_t16" 1 "ENCFF207KVO" "${COMMON}" "${COMMON}"
        ;;
    9)
        # Gate G2: Tier 2 RBMX_KD con_rep1 (6.74M ONT reads) + gencode.v46 (1.6G GTF) @ -t 16
        COMMON="--clean_start -d nanopore --bam ${DATA_ROOT}/RBMX_KD/alignments/con_rep1.fastq.sorted.nosecondary.bam --genedb ${V46_GTF} --complete_genedb -r ${GRCH38_FA} -t 16 -p con_rep1"
        run_pair "g2_con_rep1_v46_t16" 1 "con_rep1" "${COMMON}" "${COMMON}"
        ;;
    10)
        # Gate G2: isoquant-next CAGE-assisted discovery on con_rep1 + gencode.v46 + CAGE_WTC11 @ -t 16
        LOG_DIR="${BENCH_ROOT}/benchmarks/g2_con_rep1_cage_t16"
        mkdir -p "${LOG_DIR}"
        CAND_CMD="source '${CONDA_SH}' && conda activate sqanti3 && cd '${REPO_ROOT}' && rm -rf '${SCRATCH}/next' && isoquant-next --clean_start -d nanopore --fl_data --bam ${DATA_ROOT}/RBMX_KD/alignments/con_rep1.fastq.sorted.nosecondary.bam --genedb ${V46_GTF} --complete_genedb -r ${GRCH38_FA} --cage ${CAGE_BED} -t 16 -p con_rep1 -o '${SCRATCH}/next'"
        python "${BENCH_ROOT}/scripts/benchmark/run_tiered_benchmark.py" \
            --repeats 1 --interval 1.0 --log-dir "${LOG_DIR}" --candidate-cmd "${CAND_CMD}"
        cp -f "${SCRATCH}/next/isoquant.log" "${LOG_DIR}/next.isoquant.log" 2>/dev/null || true
        python -c "
import json
from pathlib import Path
s = Path('${SCRATCH}/next/con_rep1')
stats = {
    'tss_prediction_lines': sum(1 for _ in open(s / 'con_rep1.TSS_prediction.tsv')) if (s / 'con_rep1.TSS_prediction.tsv').exists() else 0,
    'transcript_models_count': sum(1 for l in open(s / 'con_rep1.transcript_models.gtf') if '\ttranscript\t' in l) if (s / 'con_rep1.transcript_models.gtf').exists() else 0,
    'gene_counts_lines': sum(1 for _ in open(s / 'con_rep1.gene_counts.tsv')) if (s / 'con_rep1.gene_counts.tsv').exists() else 0,
}
Path('${LOG_DIR}/cage_summary.json').write_text(json.dumps(stats, indent=2) + '\n')
print(json.dumps(stats, indent=2))
"
        ;;
    11)
        # Gate G3: Tier 3 simulated1 (8.77M NanoSim reads) + gencode.v46 + Recall/Precision eval @ -t 16
        COMMON="--clean_start -d nanopore --bam ${DATA_ROOT}/SIRV_simulation/alignments/simulated1_aligned_reads.sorted.bam --genedb ${V46_GTF} --complete_genedb -r ${GRCH38_FA} -t 16 -p simulated1"
        run_pair "g3_simulated1_t16" 1 "simulated1" "${COMMON}" "${COMMON}"
        python "${BENCH_ROOT}/scripts/benchmark/eval_simulated_sirv.py" \
            --ref-dir "${SCRATCH}/orig" \
            --cand-dir "${SCRATCH}/next" \
            --prefix "simulated1" \
            --json-out "${BENCH_ROOT}/benchmarks/g3_simulated1_t16/eval_simulated1.json" || true
        ;;
    12)
        # Gate G3: Tier 3 ENCFF598YQO (mouse F121 R2C2 ONT, 10.11M reads) + vM28 + GRCm39 @ -t 16
        COMMON="--clean_start -d nanopore --fl_data --bam ${DATA_ROOT}/R2C2_ONT_fulpass/mouse/F121_cell/alignments/ENCFF598YQO.bam --genedb ${VM28_GTF} --complete_genedb -r ${GRCM39_FA} -t 16 -p ENCFF598YQO"
        run_pair "g3_encff598yqo_mouse_t16" 1 "ENCFF598YQO" "${COMMON}" "${COMMON}"
        ;;
    13)
        # Gate G3: Tier 3 ENCFF089AMF (capped ONT human, 16.42M reads) + v39 + GRCh38 @ -t 32
        COMMON="--clean_start -d nanopore --fl_data --bam ${DATA_ROOT}/capped_ONT_fulpass/human/H1_ES_cell/alignments/ENCFF089AMF.bam --genedb ${V39_GTF} --complete_genedb -r ${GRCH38_FA} -t 32 -p ENCFF089AMF"
        run_pair "g3_encff089amf_t32" 1 "ENCFF089AMF" "${COMMON}" "${COMMON}"
        ;;
    14)
        # Gate G3: Tier 3 shRBMX_rep2 (7.8G BAM, 20.44M mapped reads) + gencode.v46 (1.6G GTF) @ -t 32
        COMMON="--clean_start -d nanopore --bam ${DATA_ROOT}/RBMX_KD/alignments/shRBMX_rep2.fastq.sorted.nosecondary.bam --genedb ${V46_GTF} --complete_genedb -r ${GRCH38_FA} -t 32 -p shRBMX_rep2"
        run_pair "g3_shrbmx_rep2_v46_t32" 1 "shRBMX_rep2" "${COMMON}" "${COMMON}"
        ;;
    15)
        # Gate G4: ENCFF008HMY.fastq (1.6G, 667k directRNA reads) from-scratch minimap2 -> full pipeline @ -t 16
        COMMON="--clean_start -d nanopore --stranded forward --fl_data --fastq ${DATA_ROOT}/ENCFF008HMY.fastq --genedb ${V39_GTF} --complete_genedb -r ${GRCH38_FA} -t 16 -p ENCFF008HMY"
        run_pair "g4_encff008hmy_fastq_t16" 1 "ENCFF008HMY" "${COMMON}" "${COMMON}"
        ;;
    16)
        # Gate G1/G4: Warm precompiled annotation comparison (.db vs .iqidx) on ENCFF400BQQ @ -t 8 A-B-B-A 3 reps
        LOG_DIR="${BENCH_ROOT}/benchmarks/g1_encff400bqq_t8_warm_abba3"
        mkdir -p "${LOG_DIR}"
        cd "${REPO_ROOT}"
        python -c "
from isoquant_lib.gtf2db import gtf2db
gtf2db('${V39_GTF}', '${SCRATCH}/pre_v39.db', True, False)
"
        python -c "
from isoquant_next.annotation.iqidx import compile_gtf_to_iqidx
compile_gtf_to_iqidx('${V39_GTF}', '${SCRATCH}/pre_v39.iqidx')
"
        BASE_FLAGS="-d pacbio_ccs --fl_data --bam ${DATA_ROOT}/cDNA_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF400BQQ.sorted.k14.bam --complete_genedb -r ${GRCH38_FA} -t 8 -p ENCFF400BQQ"
        ORIG_CMD="source '${CONDA_SH}' && conda activate sqanti3 && cd '${REPO_ROOT}' && rm -rf '${SCRATCH}/orig' && python isoquant.py ${BASE_FLAGS} --genedb '${SCRATCH}/pre_v39.db' -o '${SCRATCH}/orig'"
        NEXT_CMD="source '${CONDA_SH}' && conda activate sqanti3 && cd '${REPO_ROOT}' && rm -rf '${SCRATCH}/next' && isoquant-next ${BASE_FLAGS} --genedb '${SCRATCH}/pre_v39.iqidx' -o '${SCRATCH}/next'"
        python "${BENCH_ROOT}/scripts/benchmark/run_tiered_benchmark.py" \
            --repeats 3 --interval 1.0 --log-dir "${LOG_DIR}" \
            --baseline-cmd "${ORIG_CMD}" --candidate-cmd "${NEXT_CMD}"
        python "${BENCH_ROOT}/scripts/benchmark/summarize_pair.py" \
            --ref-dir "${SCRATCH}/orig" --cand-dir "${SCRATCH}/next" \
            --prefix "ENCFF400BQQ" --json-out "${LOG_DIR}/accuracy_summary.json" || true
        ;;
    *)
        echo "Unknown SLURM_ARRAY_TASK_ID=${TASK_ID}" >&2
        exit 2
        ;;
esac

echo "[$(date -Is)] Finished Phase 5 task ${TASK_ID}; /dev/shm check:"
ls -la /dev/shm || true
