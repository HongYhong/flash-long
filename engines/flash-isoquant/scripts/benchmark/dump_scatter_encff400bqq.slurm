#!/usr/bin/env bash
#SBATCH -J iq_scatter
#SBATCH -p chj.cn.1,chj.cn.2,chj.cn.3,chj.cn.4,chj.cp.1
#SBATCH -N 1
#SBATCH --cpus-per-task=20
#SBATCH --mem=160G
#SBATCH -o /gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/benchmarks/slurm_logs/scatter_%j.out
#SBATCH -e /gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/benchmarks/slurm_logs/scatter_%j.err

set -euo pipefail
source "/gpfs/group/data_sci_group/weiwu_group/hongyanhong/miniconda3/miniconda3/etc/profile.d/conda.sh"
conda activate sqanti3

REPO_ROOT="/gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant"
DATA_ROOT="/gpfs/group/data_sci_group/weiwu_group/hongyanhong/backup/longpass_benchmark250209"
OUT_SCATTER="${REPO_ROOT}/isoquant_refactored/benchmarks/scatter_encff400bqq"
SHORT_TMP="/gpfs/group/data_sci_group/weiwu_group/hongyanhong/iq/sc${SLURM_JOB_ID}"
SCRATCH="/gpfs/group/data_sci_group/weiwu_group/hongyanhong/iq/j${SLURM_JOB_ID}_sc"
mkdir -p "${SHORT_TMP}" "${SCRATCH}/mpl" "${OUT_SCATTER}"
export TMPDIR="${SHORT_TMP}" TMP="${SHORT_TMP}" TEMP="${SHORT_TMP}" MPLCONFIGDIR="${SCRATCH}/mpl"

cleanup() { rm -rf "${SCRATCH}" "${SHORT_TMP}"; }
trap cleanup EXIT

COMMON="--clean_start -d pacbio_ccs --fl_data --bam ${DATA_ROOT}/cDNA_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF400BQQ.sorted.k14.bam --genedb ${DATA_ROOT}/genome/lrgasp_gencode_v39_annotation_sirvs.human.gtf --complete_genedb -r ${DATA_ROOT}/genome/GRCh38.p14.genome.fa -t 16 -p ENCFF400BQQ"
cd "${REPO_ROOT}"
python isoquant.py ${COMMON} -o "${SCRATCH}/orig"
isoquant-next ${COMMON} -o "${SCRATCH}/next"

cp -f "${SCRATCH}/orig/ENCFF400BQQ/ENCFF400BQQ.gene_counts.tsv" "${OUT_SCATTER}/orig_gene.tsv"
cp -f "${SCRATCH}/next/ENCFF400BQQ/ENCFF400BQQ.gene_counts.tsv" "${OUT_SCATTER}/next_gene.tsv"
cp -f "${SCRATCH}/orig/ENCFF400BQQ/ENCFF400BQQ.transcript_counts.tsv" "${OUT_SCATTER}/orig_tx.tsv"
cp -f "${SCRATCH}/next/ENCFF400BQQ/ENCFF400BQQ.transcript_counts.tsv" "${OUT_SCATTER}/next_tx.tsv"
echo "Scatter dump completed cleanly."
