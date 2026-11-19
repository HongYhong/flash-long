"""Polars Out-of-Core matrix aggregation, SciPy COO MTX export, and zero-copy .iqidx/DB annotation lookup."""

from __future__ import annotations

import gzip
import logging
import re
from contextlib import ExitStack
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import polars as pl
from scipy.io import mmwrite
from scipy.sparse import coo_matrix

from isoquant_next.annotation.gtf_validator import detect_annotation_format, parse_gtf_attributes
from isoquant_next.annotation.iqidx import MmapAnnotationIndex
from isoquant_next.core.file_naming import (
    mtx_barcodes_file,
    mtx_exclude_matrix_file,
    mtx_features_file,
    mtx_include_matrix_file,
    mtx_matrix_file,
)

GROUP_COUNT_CUTOFF = 100
PROFILE_LINEAR_COLS = [
    "chr",
    "start",
    "end",
    "strand",
    "flags",
    "gene_ids",
    "group_id",
    "include_counts",
    "exclude_counts",
]

logger = logging.getLogger("IsoQuant")


def _open_write(path: str, gzipped: bool):
    if gzipped:
        return gzip.open(path + ".gz", "wt", encoding="utf-8")
    return open(path, "w", encoding="utf-8")


def _open_write_bin(path: str, gzipped: bool):
    if gzipped:
        return gzip.open(path + ".gz", "wb")
    return open(path, "wb")


def _count_header_cols(input_linear_counts: str) -> int:
    with open(input_linear_counts, "r", encoding="utf-8") as fh:
        first_line = fh.readline().rstrip("\r\n")
    return len(first_line.split("\t")) if first_line else 0


def _compute_group_norm_factors(
    df: pl.DataFrame, usable_reads_per_group: Optional[Dict[str, float]] = None
) -> Dict[str, float]:
    if usable_reads_per_group:
        return {
            gid: (1000000.0 / float(cnt)) if float(cnt) > 0.0 else 0.0
            for gid, cnt in usable_reads_per_group.items()
        }
    totals = df.group_by("group_id").agg(pl.col("count").sum().alias("total"))
    return {
        row[0]: (1000000.0 / float(row[1])) if float(row[1]) > 0.0 else 0.0
        for row in totals.iter_rows()
    }


def _scan_linear_counts(input_linear_counts: str) -> pl.DataFrame:
    df = (
        pl.scan_csv(
            input_linear_counts,
            separator="\t",
            has_header=False,
            skip_rows=1,
            schema_overrides={"column_1": pl.String, "column_2": pl.String, "column_3": pl.Float64},
            null_values=[],
            raise_if_empty=False,
        )
        .rename({"column_1": "gene_id", "column_2": "group_id", "column_3": "count"})
        .collect(streaming=True)
    )
    if df.width == 0:
        return pl.DataFrame(schema={"gene_id": pl.String, "group_id": pl.String, "count": pl.Float64})
    return df


def _scan_profile_linear(input_linear_counts: str) -> Optional[pl.DataFrame]:
    cols = _count_header_cols(input_linear_counts)
    if cols != 9:
        logger.error("Unexpected number of columns in %s: %d (expected 9)", input_linear_counts, cols)
        return None
    schema_overrides = {
        "column_1": pl.String,
        "column_2": pl.Int64,
        "column_3": pl.Int64,
        "column_4": pl.String,
        "column_5": pl.String,
        "column_6": pl.String,
        "column_7": pl.String,
        "column_8": pl.Float64,
        "column_9": pl.Float64,
    }
    rename_map = {f"column_{i + 1}": col for i, col in enumerate(PROFILE_LINEAR_COLS)}
    df = (
        pl.scan_csv(
            input_linear_counts,
            separator="\t",
            has_header=False,
            skip_rows=1,
            schema_overrides=schema_overrides,
            null_values=[],
            raise_if_empty=False,
        )
        .rename(rename_map)
        .with_columns(
            (
                pl.col("chr")
                + pl.lit(":")
                + pl.col("start").cast(pl.String)
                + pl.lit("-")
                + pl.col("end").cast(pl.String)
                + pl.lit(":")
                + pl.col("strand")
            ).alias("feature_id")
        )
        .collect(streaming=True)
    )
    if df.width == 0:
        schema = {c: schema_overrides[f"column_{i + 1}"] for i, c in enumerate(PROFILE_LINEAR_COLS)}
        schema["feature_id"] = pl.String
        return pl.DataFrame(schema=schema)
    return df


def _fmt_num(val: float) -> str:
    if val == 0.0:
        return "0"
    if float(int(val)) == val and abs(val) < 1e15:
        return f"{val:.1f}" if isinstance(val, float) and "." in str(val) else str(val)
    return str(val)


def convert_to_matrix(
    input_linear_counts: str,
    output_file_path: str,
    feature_id_to_name_dict: Optional[Dict[str, str]] = None,
    gzipped: bool = False,
    feature_type: str = "gene",
    convert_to_tpm: bool = False,
    usable_reads_per_group: Optional[Dict[str, float]] = None,
    max_groups: int = 0,
) -> int:
    """Convert 3-column linear grouped counts into full dense TSV matrix via Polars."""
    logger.info("Converting %s to full matrix", input_linear_counts)
    cols = _count_header_cols(input_linear_counts)
    if cols != 3:
        logger.error("Unexpected number of columns in %s: %d", input_linear_counts, cols)
        return 0

    if max_groups > 0:
        num_groups_df = (
            pl.scan_csv(
                input_linear_counts,
                separator="\t",
                has_header=False,
                skip_rows=1,
                schema_overrides={"column_1": pl.String, "column_2": pl.String, "column_3": pl.Float64},
                raise_if_empty=False,
            )
            .select(pl.col("column_2").n_unique())
            .collect(streaming=True)
        )
        num_groups = int(num_groups_df.item()) if num_groups_df.height > 0 else 0
        if num_groups > max_groups:
            logger.info(
                "Skipping full matrix conversion: %d groups exceeds limit of %d. "
                "Use '--counts_format matrix' to force conversion.",
                num_groups,
                max_groups,
            )
            return int(num_groups)

    df = _scan_linear_counts(input_linear_counts)
    if df.height == 0:
        out_path = output_file_path + ".tsv"
        with _open_write(out_path, gzipped) as outfile:
            outfile.write(f"{feature_type}_id\t\n")
        return 0

    unique_genes = sorted(df.get_column("gene_id").unique().to_list())
    unique_groups = sorted(df.get_column("group_id").unique().to_list())
    num_groups = len(unique_groups)
    gene_to_idx = {g: i for i, g in enumerate(unique_genes)}
    group_to_idx = {g: i for i, g in enumerate(unique_groups)}

    norm_factors = (
        {g: 1.0 for g in unique_groups}
        if not convert_to_tpm
        else _compute_group_norm_factors(df, usable_reads_per_group)
    )

    matrix = np.zeros((len(unique_genes), num_groups), dtype=np.float64)
    for g_id, grp_id, cnt in df.iter_rows():
        matrix[gene_to_idx[g_id], group_to_idx[grp_id]] += float(cnt) * norm_factors.get(grp_id, 1.0)

    out_path = output_file_path + ".tsv"
    if num_groups > GROUP_COUNT_CUTOFF:
        logger.warning(
            "You have %d groups in your matrix, conversion might take a lot of time and the output file can be very large",
            num_groups,
        )

    with _open_write(out_path, gzipped) as outfile:
        outfile.write(f"{feature_type}_id\t" + "\t".join(unique_groups) + "\n")
        for r_idx, gene_id in enumerate(unique_genes):
            gene_name = (
                gene_id
                if feature_id_to_name_dict is None
                else feature_id_to_name_dict.get(gene_id, gene_id)
            )
            row_vals = matrix[r_idx]
            cells = ["0" if v == 0.0 else str(v) for v in row_vals]
            outfile.write(gene_name + "\t" + "\t".join(cells) + "\n")

    logger.info("Matrix was saved to %s", out_path + (".gz" if gzipped else ""))
    return num_groups


def convert_to_mtx(
    input_linear_counts: str,
    output_file_prefix: str,
    feature_id_to_name: Optional[Dict[str, str]] = None,
    gzipped: bool = False,
    convert_to_tpm: bool = False,
    usable_reads_per_group: Optional[Dict[str, float]] = None,
) -> None:
    """Convert 3-column linear counts to Matrix Market MTX using Polars + scipy.sparse.coo_matrix."""
    logger.info("Converting %s to MTX format", input_linear_counts)
    cols = _count_header_cols(input_linear_counts)
    if cols != 3:
        logger.error("Unexpected number of columns in %s: %d", input_linear_counts, cols)
        return

    df = _scan_linear_counts(input_linear_counts)
    unique_genes = sorted(df.get_column("gene_id").unique().to_list()) if df.height > 0 else []
    unique_groups = sorted(df.get_column("group_id").unique().to_list()) if df.height > 0 else []

    mtx_file = mtx_matrix_file(output_file_prefix)
    features_file = mtx_features_file(output_file_prefix)
    barcodes_file = mtx_barcodes_file(output_file_prefix)

    with open(barcodes_file, "w", encoding="utf-8") as bc_out:
        for group in unique_groups:
            bc_out.write(f"{group}\n")

    with open(features_file, "w", encoding="utf-8") as ft_out:
        for gene_id in unique_genes:
            gene_name = (
                feature_id_to_name.get(gene_id, gene_id) if feature_id_to_name is not None else gene_id
            )
            ft_out.write(f"{gene_id}\t{gene_name}\n")

    norm_factors = (
        {g: 1.0 for g in unique_groups}
        if not convert_to_tpm
        else _compute_group_norm_factors(df, usable_reads_per_group)
    )

    if df.height == 0:
        with _open_write_bin(mtx_file, gzipped) as out_bin:
            mat = coo_matrix((len(unique_genes), len(unique_groups)), dtype=np.float64)
            mmwrite(out_bin, mat)
        logger.info("Matrix was saved to 3 files with prefix %s", output_file_prefix)
        return

    gene_df = pl.DataFrame({"gene_id": unique_genes, "_r": np.arange(len(unique_genes), dtype=np.int32)})
    group_df = pl.DataFrame(
        {
            "group_id": unique_groups,
            "_c": np.arange(len(unique_groups), dtype=np.int32),
            "_scale": np.asarray([norm_factors.get(g, 1.0) for g in unique_groups], dtype=np.float64),
        }
    )

    encoded = (
        df.join(gene_df, on="gene_id", how="inner")
        .join(group_df, on="group_id", how="inner")
        .with_columns((pl.col("count") * pl.col("_scale")).alias("_val"))
        .group_by(["_r", "_c"], maintain_order=False)
        .agg(pl.col("_val").sum())
        .filter(pl.col("_val") != 0.0)
        .sort(["_c", "_r"])
    )

    rows = encoded.get_column("_r").to_numpy()
    cols_idx = encoded.get_column("_c").to_numpy()
    vals = encoded.get_column("_val").to_numpy().astype(np.float64, copy=False)

    mat = coo_matrix((vals, (rows, cols_idx)), shape=(len(unique_genes), len(unique_groups)), dtype=np.float64)
    with _open_write_bin(mtx_file, gzipped) as out_bin:
        mmwrite(out_bin, mat)

    logger.info("Matrix was saved to 3 files with prefix %s", output_file_prefix)


def _ordered_unique_features(df: pl.DataFrame) -> Tuple[List[str], Dict[str, str]]:
    uniq = df.select(["feature_id", "gene_ids"]).unique(subset=["feature_id"], keep="first", maintain_order=True)
    fids = uniq.get_column("feature_id").to_list()
    gids = uniq.get_column("gene_ids").to_list()
    return fids, dict(zip(fids, gids))


def convert_profile_to_matrix(
    input_linear_counts: str,
    output_file_path: str,
    gzipped: bool = False,
    feature_type: str = "exon",
    max_groups: int = 0,
) -> int:
    """Convert 9-column exon/intron linear counts to dense 'include,exclude' matrix TSV."""
    logger.info("Converting %s to full matrix", input_linear_counts)
    df = _scan_profile_linear(input_linear_counts)
    if df is None:
        return 0

    if max_groups > 0:
        num_groups = df.get_column("group_id").n_unique()
        if num_groups > max_groups:
            logger.info(
                "Skipping full matrix conversion: %d groups exceeds limit of %d. "
                "Use '--counts_format matrix' to force conversion.",
                num_groups,
                max_groups,
            )
            return int(num_groups)

    feature_order, _ = _ordered_unique_features(df)
    columns = sorted(df.get_column("group_id").unique().to_list()) if df.height > 0 else []
    num_groups = len(columns)

    f_map = {fid: i for i, fid in enumerate(feature_order)}
    c_map = {grp: i for i, grp in enumerate(columns)}
    incl_mat = np.zeros((len(feature_order), num_groups), dtype=np.int64)
    excl_mat = np.zeros((len(feature_order), num_groups), dtype=np.int64)

    for fid, grp, inc_v, exc_v in df.select(
        ["feature_id", "group_id", "include_counts", "exclude_counts"]
    ).iter_rows():
        r = f_map[fid]
        c = c_map[grp]
        incl_mat[r, c] += int(inc_v)
        excl_mat[r, c] += int(exc_v)

    out_path = output_file_path + ".tsv"
    if num_groups > GROUP_COUNT_CUTOFF:
        logger.warning(
            "You have %d groups in your matrix, conversion might take a lot of time and the output file can be very large",
            num_groups,
        )

    with _open_write(out_path, gzipped) as outfile:
        outfile.write(f"{feature_type}_id\t" + "\t".join(columns) + "\n")
        for r, feature_id in enumerate(feature_order):
            inc_row = incl_mat[r]
            exc_row = excl_mat[r]
            cells = [f"{inc_row[c]},{exc_row[c]}" for c in range(num_groups)]
            outfile.write(feature_id + "\t" + "\t".join(cells) + "\n")

    logger.info("Matrix was saved to %s", out_path + (".gz" if gzipped else ""))
    return num_groups


def convert_profile_to_mtx(
    input_linear_counts: str,
    output_file_prefix: str,
    gzipped: bool = False,
) -> None:
    """Convert 9-column exon/intron linear counts to paired include & exclude MTX files via scipy coo_matrix."""
    logger.info("Converting %s to MTX format", input_linear_counts)
    df = _scan_profile_linear(input_linear_counts)
    if df is None:
        return

    feature_order, feature_to_gene = _ordered_unique_features(df)
    unique_groups = sorted(df.get_column("group_id").unique().to_list()) if df.height > 0 else []

    features_file = mtx_features_file(output_file_prefix)
    barcodes_file = mtx_barcodes_file(output_file_prefix)
    include_file = mtx_include_matrix_file(output_file_prefix)
    exclude_file = mtx_exclude_matrix_file(output_file_prefix)

    with open(barcodes_file, "w", encoding="utf-8") as bc_out:
        for group in unique_groups:
            bc_out.write(f"{group}\n")

    with open(features_file, "w", encoding="utf-8") as ft_out:
        for feature_id in feature_order:
            ft_out.write(f"{feature_id}\t{feature_to_gene.get(feature_id, '')}\n")

    feat_df = pl.DataFrame({"feature_id": feature_order, "_r": np.arange(len(feature_order), dtype=np.int32)})
    grp_df = pl.DataFrame({"group_id": unique_groups, "_c": np.arange(len(unique_groups), dtype=np.int32)})

    agg = (
        df.join(feat_df, on="feature_id", how="inner")
        .join(grp_df, on="group_id", how="inner")
        .group_by(["_r", "_c"], maintain_order=False)
        .agg(
            [
                pl.col("include_counts").sum().alias("_inc"),
                pl.col("exclude_counts").sum().alias("_exc"),
            ]
        )
        .sort(["_c", "_r"])
    )

    shape = (len(feature_order), len(unique_groups))

    def _emit_sparse(path: str, col_name: str) -> None:
        nz = agg.filter(pl.col(col_name) != 0)
        if nz.height == 0:
            mat = coo_matrix(shape, dtype=np.int64)
        else:
            r = nz.get_column("_r").to_numpy()
            c = nz.get_column("_c").to_numpy()
            v = nz.get_column(col_name).to_numpy().astype(np.int64, copy=False)
            mat = coo_matrix((v, (r, c)), shape=shape, dtype=np.int64)
        with _open_write_bin(path, gzipped) as out_bin:
            mmwrite(out_bin, mat)

    _emit_sparse(include_file, "_inc")
    _emit_sparse(exclude_file, "_exc")
    logger.info("MTX output saved with prefix %s (include + exclude matrices)", output_file_prefix)


def load_gene_name_map(genedb_path: str | Path) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, Tuple[str, str]]]:
    """Load gene_name, transcript_name, and transcript2gene from `.iqidx` or legacy `.db`."""
    path = Path(genedb_path)
    fmt, _ = detect_annotation_format(path)
    gene_names: Dict[str, str] = {}
    transcript_names: Dict[str, str] = {}
    transcript2gene: Dict[str, Tuple[str, str]] = {}

    if fmt == "iqidx":
        with MmapAnnotationIndex(path) as idx:
            n_genes = len(idx.gene_id_str)
            gene_id_by_idx: List[str] = [""] * n_genes
            for gi in range(n_genes):
                gid = idx.get_str(int(idx.gene_id_str[gi]))
                attrs = parse_gtf_attributes(idx.get_str(int(idx.gene_attr_str[gi])))
                gene_names[gid] = attrs.get("gene_name", gid)
                gene_id_by_idx[gi] = gid

            n_tx = len(idx.tx_id_str)
            for ti in range(n_tx):
                tid = idx.get_str(int(idx.tx_id_str[ti]))
                gidx = int(idx.tx_gene[ti])
                gid = gene_id_by_idx[gidx] if 0 <= gidx < n_genes else "unknown_gene"
                attrs = parse_gtf_attributes(idx.get_str(int(idx.tx_attr_str[ti])))
                transcript2gene[tid] = (gid, gene_names.get(gid, gid))
                transcript_names[tid] = attrs.get("transcript_name", tid)
        return gene_names, transcript_names, transcript2gene

    import gffutils

    db = gffutils.FeatureDB(str(path), keep_order=True)
    for g in db.all_features(featuretype="gene"):
        gene_names[g.id] = g["gene_name"][0] if "gene_name" in g.attributes else g.id
    for t in db.all_features(featuretype=("transcript", "mRNA")):
        if "gene_id" in t.attributes:
            gid = t["gene_id"][0]
        elif "Parent" in t.attributes:
            gid = t["Parent"][0]
        else:
            gid = "unknown_gene"
        transcript2gene[t.id] = (gid, gene_names.get(gid, gid))
        transcript_names[t.id] = t["transcript_name"][0] if "transcript_name" in t.attributes else t.id
    return gene_names, transcript_names, transcript2gene


__all__ = [
    "GROUP_COUNT_CUTOFF",
    "PROFILE_LINEAR_COLS",
    "convert_profile_to_matrix",
    "convert_profile_to_mtx",
    "convert_to_matrix",
    "convert_to_mtx",
    "load_gene_name_map",
]
