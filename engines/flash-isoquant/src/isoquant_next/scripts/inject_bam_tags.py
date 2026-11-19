"""Inject barcode (CB) and UMI (UB) tags into a BAM file from barcoded_reads.tsv."""

from __future__ import annotations

import argparse
import logging
import os
import sys

import pysam

logger = logging.getLogger(__name__)


def load_barcoded_reads(barcoded_reads_path: str) -> dict[str, tuple[str, str]]:
    barcode_dict: dict[str, tuple[str, str]] = {}
    with open(barcoded_reads_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or line.startswith("read_id\t"):
                continue
            tokens = line.strip().split("\t")
            if len(tokens) < 3 or tokens[1] == "*":
                continue
            barcode_dict[tokens[0]] = (tokens[1], tokens[2])
    return barcode_dict


def inject_tags(
    input_bam: str,
    output_bam: str,
    barcode_dict: dict[str, tuple[str, str]],
    barcode_tag: str = "CB",
    umi_tag: str = "UB",
    suffix: str = "",
) -> None:
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        with pysam.AlignmentFile(output_bam, "wb", header=infile.header) as outfile:
            for read in infile:
                read_id = read.query_name
                if read_id in barcode_dict:
                    barcode, umi = barcode_dict[read_id]
                    read.set_tag(barcode_tag, barcode + suffix)
                    if umi and umi != "*":
                        read.set_tag(umi_tag, umi)
                outfile.write(read)
    pysam.index(output_bam)


def main() -> None:
    parser = argparse.ArgumentParser(description="Inject barcode/UMI tags into BAM")
    parser.add_argument("--bam", required=True)
    parser.add_argument("--barcoded_reads", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--barcode_tag", default="CB")
    parser.add_argument("--umi_tag", default="UB")
    parser.add_argument("--suffix", default="")
    args = parser.parse_args()
    if not os.path.exists(args.bam) or not os.path.exists(args.barcoded_reads):
        sys.exit(1)
    inject_tags(args.bam, args.output, load_barcoded_reads(args.barcoded_reads), args.barcode_tag, args.umi_tag, args.suffix)
