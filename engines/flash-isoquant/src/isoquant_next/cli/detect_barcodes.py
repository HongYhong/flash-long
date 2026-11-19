"""Standalone barcode detection CLI entrypoint (`isoquant-next-detect-barcodes`)."""

from __future__ import annotations

import argparse
import logging
import os
import sys
from typing import List, Optional

from isoquant_next.barcodes.detect_barcodes import (
    BARCODE_CALLING_MODES,
    process_in_parallel,
    process_single_thread,
)
from isoquant_next.config.enums import IsoQuantMode
from isoquant_next.core.error_codes import IsoQuantExitCode
from isoquant_next.core.logging_setup import configure_logging, pin_single_thread_math_env

logger = logging.getLogger("IsoQuant")


def parse_args(sys_argv: List[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="isoquant-next-detect-barcodes")
    parser.add_argument("--output", "-o", type=str, required=True, help="output prefix name")
    parser.add_argument("--barcodes", "-b", nargs="+", type=str, required=False, help="barcode whitelist(s)")
    parser.add_argument("--mode", type=str, required=True, choices=[x.name for x in BARCODE_CALLING_MODES.keys()])
    parser.add_argument("--molecule", type=str, default=None, help="MDF file for custom_sc mode")
    parser.add_argument("--input", "-i", nargs="+", type=str, required=True, help="input reads (FASTA/FASTQ/BAM)")
    parser.add_argument("--threads", "-t", type=int, default=16, help="threads to use [16]")
    parser.add_argument("--tmp_dir", type=str, default=None)
    parser.add_argument("--debug", action="store_true", default=False)

    args = parser.parse_args(sys_argv)
    args.mode = IsoQuantMode[args.mode]
    args.out_fasta = None
    args.output_tsv = None
    return args


def check_args(args: argparse.Namespace) -> None:
    n = len(args.input)
    if args.output_tsv is None:
        args.output_tsv = (
            [f"{args.output}.barcoded_reads.tsv"]
            if n == 1
            else [f"{args.output}_{i}.barcoded_reads.tsv" for i in range(n)]
        )
    if args.out_fasta is None and args.mode.produces_new_fasta():
        args.out_fasta = (
            [f"{args.output}.split_reads.fasta"]
            if n == 1
            else [f"{args.output}_{i}.split_reads.fasta" for i in range(n)]
        )


def main(cmd_args: Optional[List[str]] = None) -> int:
    pin_single_thread_math_env(force=True)
    argv = list(sys.argv[1:] if cmd_args is None else cmd_args)
    args = parse_args(argv)
    configure_logging(debug=args.debug)
    check_args(args)
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    if not args.tmp_dir:
        args.tmp_dir = out_dir or "."

    if args.threads == 1 or args.mode.enforces_single_thread():
        process_single_thread(args)
    else:
        process_in_parallel(args)
    return int(IsoQuantExitCode.SUCCESS)


if __name__ == "__main__":
    raise SystemExit(main())
