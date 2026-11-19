"""Unified logging configuration and OpenMP/BLAS thread-pinning guards."""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path
from typing import Optional

MATH_THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
)


def pin_single_thread_math_env(force: bool = True) -> None:
    """Pin C/BLAS/OpenMP libraries to 1 thread per process to prevent oversubscription & fork deadlocks."""
    for key in MATH_THREAD_ENV_VARS:
        if force or key not in os.environ:
            os.environ[key] = "1"


def configure_logging(
    *,
    debug: bool = False,
    log_file: Optional[str | Path] = None,
    logger_name: str = "IsoQuant",
) -> logging.Logger:
    """Configure logger without leaking unclosed file handlers across repeated calls."""
    logger = logging.getLogger(logger_name)
    for h in list(logger.handlers):
        logger.removeHandler(h)
        try:
            h.close()
        except OSError:
            pass

    logger.setLevel(logging.DEBUG if debug else logging.INFO)
    logger.propagate = False
    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setLevel(logging.DEBUG if debug else logging.INFO)
    stream_handler.setFormatter(fmt)
    logger.addHandler(stream_handler)

    if log_file is not None:
        p = Path(log_file)
        p.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(p, mode="a", encoding="utf-8")
        file_handler.setLevel(logging.DEBUG if debug else logging.INFO)
        file_handler.setFormatter(fmt)
        logger.addHandler(file_handler)

    return logger
