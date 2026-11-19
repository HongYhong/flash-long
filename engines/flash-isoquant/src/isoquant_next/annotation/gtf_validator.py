"""Streaming GTF/GFF3 format detection, quote-aware attribute lexer, and duplicate deduplication."""

from __future__ import annotations

import gzip
import logging
import os
import re
import tempfile
from pathlib import Path
from typing import Dict, Iterable, Optional, TextIO, Tuple

logger = logging.getLogger("IsoQuant")

# Matches GTF key "val" or key val pairs safely even when values contain whitespace
_GTF_ATTR_RE = re.compile(r'([A-Za-z0-9_.\-]+)\s+(?:"((?:[^"\\]|\\.)*)"|([^;\s]+))')


def detect_annotation_format(file_path: str | Path) -> Tuple[str, bool]:
    """Return (format in {'gtf','gff3','db'}, is_gzipped) fixing missing-dot '.gff3' bug."""
    p = Path(file_path)
    name_lower = p.name.lower()
    if name_lower.endswith(".db") or name_lower.endswith(".sqlite"):
        return "db", False
    is_gz = name_lower.endswith(".gz")
    stem_lower = name_lower[:-3] if is_gz else name_lower
    ext = os.path.splitext(stem_lower)[1]
    if ext in (".gff3", ".gff"):
        return "gff3", is_gz
    if ext == ".iqidx":
        return "iqidx", False
    return "gtf", is_gz


def open_text_maybe_gz(file_path: str | Path) -> TextIO:
    _, is_gz = detect_annotation_format(file_path)
    if is_gz:
        return gzip.open(file_path, "rt", encoding="utf-8", errors="replace")
    return open(file_path, "r", encoding="utf-8", errors="replace")


def parse_gtf_attributes(attr_field: str) -> Dict[str, str]:
    """Parse column 9 GTF attributes supporting quoted spaces and unquoted tokens."""
    attrs: Dict[str, str] = {}
    for match in _GTF_ATTR_RE.finditer(attr_field):
        key = match.group(1)
        val = match.group(2) if match.group(2) is not None else match.group(3)
        if key not in attrs:
            attrs[key] = val
    return attrs


def parse_gff3_attributes(attr_field: str) -> Dict[str, str]:
    attrs: Dict[str, str] = {}
    for item in attr_field.strip().split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def _format_gtf_attributes(attrs: Dict[str, str]) -> str:
    return " ".join(f'{k} "{v}";' for k, v in attrs.items())


def check_gtf_duplicates_streaming(
    gtf_path: str | Path, out_clean_path: Optional[str | Path] = None
) -> Tuple[bool, Optional[str]]:
    """
    Single-pass streaming GTF/GFF3 cleaner replacing full-file readlines() + O(N^2) '+=' concat.
    Returns (was_modified, cleaned_file_path_if_modified).
    """
    fmt, _ = detect_annotation_format(gtf_path)
    if fmt == "gff3":
        return _check_gff3_duplicates_streaming(gtf_path, out_clean_path)

    seen_gene_ids: Dict[str, int] = {}
    seen_tx_ids: Dict[str, int] = {}
    modified = False

    if out_clean_path is None:
        fd, tmp_name = tempfile.mkstemp(prefix="isoquant_gtf_", suffix=".gtf")
        os.close(fd)
        target_path = Path(tmp_name)
    else:
        target_path = Path(out_clean_path)

    with open_text_maybe_gz(gtf_path) as inf, open(target_path, "w", encoding="utf-8") as outf:
        for line in inf:
            if not line or line[0] == "#":
                outf.write(line)
                continue
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 9:
                outf.write(line)
                continue
            ftype = cols[2]
            if ftype in ("gene", "transcript", "mRNA"):
                attrs = parse_gtf_attributes(cols[8])
                gid = attrs.get("gene_id")
                tid = attrs.get("transcript_id")
                if ftype == "gene" and gid is not None:
                    cnt = seen_gene_ids.get(gid, 0) + 1
                    seen_gene_ids[gid] = cnt
                    if cnt > 1:
                        attrs["gene_id"] = f"{gid}.dup{cnt}"
                        cols[8] = _format_gtf_attributes(attrs)
                        modified = True
                        outf.write("\t".join(cols) + "\n")
                        continue
                elif ftype in ("transcript", "mRNA") and tid is not None:
                    cnt = seen_tx_ids.get(tid, 0) + 1
                    seen_tx_ids[tid] = cnt
                    if cnt > 1:
                        attrs["transcript_id"] = f"{tid}.dup{cnt}"
                        cols[8] = _format_gtf_attributes(attrs)
                        modified = True
                        outf.write("\t".join(cols) + "\n")
                        continue
            outf.write(line)

    if not modified and out_clean_path is None:
        try:
            target_path.unlink()
        except OSError:
            pass
        return False, None
    return modified, str(target_path)


def _check_gff3_duplicates_streaming(
    gff_path: str | Path, out_clean_path: Optional[str | Path] = None
) -> Tuple[bool, Optional[str]]:
    seen_ids: Dict[str, int] = {}
    modified = False
    if out_clean_path is None:
        fd, tmp_name = tempfile.mkstemp(prefix="isoquant_gff3_", suffix=".gff3")
        os.close(fd)
        target_path = Path(tmp_name)
    else:
        target_path = Path(out_clean_path)

    with open_text_maybe_gz(gff_path) as inf, open(target_path, "w", encoding="utf-8") as outf:
        for line in inf:
            if not line or line[0] == "#":
                outf.write(line)
                continue
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 9:
                outf.write(line)
                continue
            attrs = parse_gff3_attributes(cols[8])
            fid = attrs.get("ID")
            if fid is not None:
                cnt = seen_ids.get(fid, 0) + 1
                seen_ids[fid] = cnt
                if cnt > 1:
                    attrs["ID"] = f"{fid}.dup{cnt}"
                    cols[8] = ";".join(f"{k}={v}" for k, v in attrs.items())
                    modified = True
                    outf.write("\t".join(cols) + "\n")
                    continue
            outf.write(line)

    if not modified and out_clean_path is None:
        try:
            target_path.unlink()
        except OSError:
            pass
        return False, None
    return modified, str(target_path)
