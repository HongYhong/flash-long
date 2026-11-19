#!/usr/bin/env bash
# Unified installer helper for Flash-Long (Flash-IsoQuant + Flash-Bambu)
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TARGET="${1:-all}"

install_isoquant() {
  echo "==> Installing Flash-IsoQuant (engines/flash-isoquant)..."
  (
    cd "${ROOT_DIR}/engines/flash-isoquant"
    python3 -m pip install --no-build-isolation --no-deps -e .
  )
  echo "==> Flash-IsoQuant installed. CLI entrypoints: isoquant-next, isoquant-next-detect-barcodes, isoquant-next-visualize"
}

install_bambu() {
  local RLIB="${BAMBU_RLIB:-${ROOT_DIR}/engines/flash-bambu/.r_lib}"
  echo "==> Installing Flash-Bambu (engines/flash-bambu) into isolated library: ${RLIB}"
  mkdir -p "${RLIB}"
  R CMD INSTALL --no-multiarch -l "${RLIB}" "${ROOT_DIR}/engines/flash-bambu"
  echo "==> Flash-Bambu installed. Load via: R_LIBS_USER=\"${RLIB}\" Rscript -e 'library(bambu); packageVersion(\"bambu\")'"
}

case "${TARGET}" in
  isoquant) install_isoquant ;;
  bambu)    install_bambu ;;
  all)      install_isoquant; install_bambu ;;
  *)
    echo "Usage: $0 [all|isoquant|bambu]" >&2
    exit 1
    ;;
esac
