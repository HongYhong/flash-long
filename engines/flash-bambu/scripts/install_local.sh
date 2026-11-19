#!/usr/bin/env bash
set -euo pipefail

PKG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
LIB_DIR="${PKG_DIR}/.r_lib"
WORK_TMP="/gpfs/group/data_sci_group/weiwu_group/hongyanhong/bambu_next_work/tmp"

mkdir -p "${LIB_DIR}" "${WORK_TMP}"
export TMPDIR="${WORK_TMP}"
export TMP="${WORK_TMP}"
export TEMP="${WORK_TMP}"
export R_SESSION_TMPDIR="${WORK_TMP}"

source /gpfs/group/data_sci_group/weiwu_group/hongyanhong/miniconda3/miniconda3/etc/profile.d/conda.sh
conda activate /gpfs/group/data_sci_group/weiwu_group/hongyanhong/miniconda3/miniconda3/envs/bambu

rm -f "${PKG_DIR}"/src/*.o "${PKG_DIR}"/src/*/*.o "${PKG_DIR}/src/bambu.so"
R CMD INSTALL --no-multiarch --preclean -l "${LIB_DIR}" "${PKG_DIR}"
echo "Installed bambu_next cleanly into ${LIB_DIR}"
