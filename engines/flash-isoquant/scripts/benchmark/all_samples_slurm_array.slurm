#!/usr/bin/env bash
#SBATCH -J iq_all_cohort
#SBATCH -p chj.cn.1,chj.cn.2,chj.cn.3,chj.cn.4,chj.cp.1
#SBATCH -N 1
#SBATCH --cpus-per-task=20
#SBATCH --mem=160G
#SBATCH -o /gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/benchmarks/slurm_logs/all_%A_%a.out
#SBATCH -e /gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored/benchmarks/slurm_logs/all_%A_%a.err

set -euo pipefail

CONDA_SH="/gpfs/group/data_sci_group/weiwu_group/hongyanhong/miniconda3/miniconda3/etc/profile.d/conda.sh"
source "${CONDA_SH}"
conda activate sqanti3

BENCH_ROOT="/gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant/isoquant_refactored"
MANIFEST="${MANIFEST_PATH:-${BENCH_ROOT}/benchmarks/pending_standard_manifest.json}"
TAG="${TIER_TAG:-s}"
TASK_ID="${SLURM_ARRAY_TASK_ID:-${1:-0}}"

BASE_SCRATCH="/gpfs/group/data_sci_group/weiwu_group/hongyanhong/iq"
# Keep TMPDIR <= 62 bytes so forkserver AF_UNIX listener binds succeed cleanly without writing /tmp
SHORT_TMP="${BASE_SCRATCH}/${TAG}${TASK_ID}"
SCRATCH="${BASE_SCRATCH}/j${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-0}}_${TAG}${TASK_ID}"
mkdir -p "${SHORT_TMP}" "${SCRATCH}/mpl"
export TMPDIR="${SHORT_TMP}"
export TMP="${SHORT_TMP}"
export TEMP="${SHORT_TMP}"
export MPLCONFIGDIR="${SCRATCH}/mpl"

cleanup() {
    rm -rf "${SCRATCH}" "${SHORT_TMP}"
}
trap cleanup EXIT

echo "[$(date -Is)] Task ${TASK_ID} (tag=${TAG}, manifest=${MANIFEST}) on $(hostname) TMPDIR=${TMPDIR} (${#TMPDIR} chars)"
python "${BENCH_ROOT}/scripts/benchmark/run_manifest_task.py" \
    --manifest "${MANIFEST}" \
    --index "${TASK_ID}" \
    --scratch "${SCRATCH}" \
    --skip-existing

echo "[$(date -Is)] Finished task ${TASK_ID}; /dev/shm check:"
ls -la /dev/shm || true
