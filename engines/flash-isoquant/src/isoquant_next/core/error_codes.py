"""Exit codes and structured runtime exception types for isoquant_next."""

from __future__ import annotations

import logging
import sys
from enum import IntEnum
from typing import Optional

logger = logging.getLogger("IsoQuant")


class IsoQuantExitCode(IntEnum):
    """Standardized exit codes matching IsoQuant CLI contract."""

    SUCCESS = 0
    INPUT_FILE_NOT_FOUND = 1
    OUTPUT_ALREADY_EXISTS = 2
    BAM_NOT_INDEXED = 3
    GENE_DB_NOT_FOUND = 4
    REFERENCE_NOT_FOUND = 5
    NO_INPUT_DATA = 6
    FILE_FORMAT_ERROR = 7
    YAML_PARSING_ERROR = 8
    INVALID_FILE_FORMAT = 9

    INVALID_PARAMETER = 20
    MISSING_REQUIRED_OPTION = 21
    INCOMPATIBLE_OPTIONS = 22
    RESUME_CONFIG_NOT_FOUND = 23
    CORRUPTED_GTF = 24
    INVALID_CONVERSION_MODE = 25

    LABEL_COUNT_MISMATCH = 40
    DUPLICATE_FILES = 41
    CHROMOSOME_MISMATCH = 42
    REPLICA_COUNT_MISMATCH = 43

    ALIGNER_NOT_FOUND = 60
    INDEXING_FAILED = 61
    ALIGNMENT_FAILED = 62
    SAMTOOLS_FAILED = 63
    SUBPROCESS_FAILED = 64

    BARCODE_FILE_COUNT_MISMATCH = 80
    BARCODE_WHITELIST_MISSING = 81

    TEST_FAILED = 90
    UNCAUGHT_EXCEPTION = 99


class IsoQuantError(RuntimeError):
    """Structured exception carrying an IsoQuantExitCode."""

    def __init__(self, code: IsoQuantExitCode, message: str = "") -> None:
        self.code = IsoQuantExitCode(code)
        super().__init__(message or self.code.name)


def exit_with_code(code: IsoQuantExitCode, message: Optional[str] = None) -> None:
    """Log optional error message and exit with canonical exit code."""
    if message:
        logger.error(message)
    sys.exit(int(code))
