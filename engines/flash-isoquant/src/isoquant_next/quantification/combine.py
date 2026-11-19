"""Multi-sample count/TPM table combination powered by Polars (zero Pandas runtime dependency)."""

from __future__ import annotations

import logging
import os
from typing import Callable

import polars as pl

from isoquant_next.core.utils import proper_plural_form

logger = logging.getLogger("IsoQuant")


def transform_counts(
    path_to_tsv: str, label: str, column_name: str = "count", full: bool = False
) -> pl.DataFrame:
    df = pl.read_csv(
        path_to_tsv,
        separator="\t",
        schema_overrides={"feature_id": pl.String, column_name: pl.Float64},
        null_values=[],
    )
    if not full:
        df = df.filter(~pl.col("feature_id").str.starts_with("__"))
    return df.select(["feature_id", pl.col(column_name).alias(label)])


def combine_table(
    input_data,
    output: str,
    get_file_name: Callable,
    output_file_name: str,
    column_name: str = "count",
    full: bool = False,
) -> None:
    sample_0 = input_data.samples[0]
    combined = transform_counts(get_file_name(sample_0), sample_0.prefix, column_name, full)
    for sample in input_data.samples[1:]:
        nxt = transform_counts(get_file_name(sample), sample.prefix, column_name, full)
        combined = combined.join(nxt, on="feature_id", how="full", coalesce=True)
    combined.sort("feature_id").write_csv(os.path.join(output, output_file_name), separator="\t", null_value="")


def combine_counts(input_data, output: str) -> None:
    logger.info("Aggregating counts for %s", proper_plural_form("sample", len(input_data.samples)))
    combine_table(
        input_data,
        output,
        lambda x: f"{x.out_gene_counts_tsv}_counts.tsv",
        "combined_gene_counts.tsv",
    )
    combine_table(
        input_data,
        output,
        lambda x: f"{x.out_gene_counts_tsv}_tpm.tsv",
        "combined_gene_tpm.tsv",
        column_name="TPM",
        full=True,
    )
    combine_table(
        input_data,
        output,
        lambda x: f"{x.out_transcript_counts_tsv}_counts.tsv",
        "combined_transcript_counts.tsv",
    )
    combine_table(
        input_data,
        output,
        lambda x: f"{x.out_transcript_counts_tsv}_tpm.tsv",
        "combined_transcript_tpm.tsv",
        column_name="TPM",
        full=True,
    )
    logger.info("Aggregation finished")


__all__ = ["combine_counts", "combine_table", "transform_counts"]
