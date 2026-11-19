"""Streaming GTF/GFF3 and mmap `.iqidx` annotation dictionary builder (zero gffutils)."""

from __future__ import annotations

from dataclasses import dataclass, field
import gzip
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from isoquant_next.annotation.gtf_validator import parse_gff3_attributes, parse_gtf_attributes
from isoquant_next.annotation.iqidx import MmapAnnotationIndex


@dataclass
class TranscriptAnnotation:
    """Lightweight transcript exon structure for gene model plotting."""

    transcript_id: str
    gene_id: str
    chrom: str
    strand: str
    exons: List[Tuple[int, int]] = field(default_factory=list)
    transcript_type: str = ""
    is_discovered: bool = False


@dataclass
class GeneModelRecord:
    """Per-gene reference and discovered transcript structures."""

    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    transcripts: Dict[str, TranscriptAnnotation] = field(default_factory=dict)

    @property
    def span(self) -> Tuple[int, int]:
        starts = [ex[0] for tx in self.transcripts.values() for ex in tx.exons]
        ends = [ex[1] for tx in self.transcripts.values() for ex in tx.exons]
        if not starts:
            return (0, 0)
        return (min(starts), max(ends))


def _parse_attrs(raw: str, is_gff3: bool) -> Dict[str, str]:
    return parse_gff3_attributes(raw) if is_gff3 else parse_gtf_attributes(raw)


def load_gtf_models_streaming(
    gtf_path: str | Path,
    *,
    gene_filter: Optional[Iterable[str]] = None,
    is_discovered: bool = False,
) -> Dict[str, GeneModelRecord]:
    """Stream `.gtf` / `.gtf.gz` / `.gff3` without creating any `.db` SQLite database."""
    target_genes = set(gene_filter) if gene_filter is not None else None
    path = Path(gtf_path)
    is_gff3 = ".gff" in path.name.lower()
    opener = gzip.open if path.suffix.lower() == ".gz" else open

    genes: Dict[str, GeneModelRecord] = {}
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            ftype = cols[2]
            if ftype != "exon":
                continue
            attrs = _parse_attrs(cols[8], is_gff3)
            gid = attrs.get("gene_id") or attrs.get("Parent", "")
            tid = attrs.get("transcript_id") or attrs.get("Parent", "")
            if not gid or not tid:
                continue
            if target_genes is not None and gid not in target_genes:
                continue

            gname = attrs.get("gene_name") or gid
            rec = genes.get(gid)
            if rec is None:
                rec = GeneModelRecord(gene_id=gid, gene_name=gname, chrom=cols[0], strand=cols[6])
                genes[gid] = rec
            tx = rec.transcripts.get(tid)
            if tx is None:
                tx = TranscriptAnnotation(
                    transcript_id=tid,
                    gene_id=gid,
                    chrom=cols[0],
                    strand=cols[6],
                    transcript_type=attrs.get("transcript_type", ""),
                    is_discovered=is_discovered,
                )
                rec.transcripts[tid] = tx
            tx.exons.append((int(cols[3]), int(cols[4])))

    for rec in genes.values():
        for tx in rec.transcripts.values():
            tx.exons.sort()
    return genes


def load_iqidx_models(
    iqidx_path: str | Path, *, gene_filter: Optional[Iterable[str]] = None
) -> Dict[str, GeneModelRecord]:
    """Extract GeneModelRecord map directly from zero-copy `.iqidx` index."""
    target_genes = set(gene_filter) if gene_filter is not None else None
    _strand = {1: "+", -1: "-", 0: "."}
    genes: Dict[str, GeneModelRecord] = {}

    with MmapAnnotationIndex(iqidx_path) as idx:
        g_tx = idx.gene_tx_ptr
        g_id_s = idx.gene_id_str
        g_attr_s = idx.gene_attr_str
        t_id_s = idx.tx_id_str
        t_str = idx.tx_strand
        t_ex = idx.tx_exon_ptr
        e_start = idx.exon_start
        e_end = idx.exon_end

        for chrom in idx.seqids():
            g0, g1 = idx.chrom_gene_slice(chrom)
            for gi in range(g0, g1):
                gid = idx.get_str(int(g_id_s[gi]))
                if target_genes is not None and gid not in target_genes:
                    continue
                gattr = idx.get_str(int(g_attr_s[gi]))
                parsed = parse_gtf_attributes(gattr) if '"' in gattr else parse_gff3_attributes(gattr)
                gname = parsed.get("gene_name") or gid
                tx_dict: Dict[str, TranscriptAnnotation] = {}
                strand_val = "."
                for ti in range(int(g_tx[gi]), int(g_tx[gi + 1])):
                    tid = idx.get_str(int(t_id_s[ti]))
                    strand_val = _strand.get(int(t_str[ti]), ".")
                    exons = sorted((int(e_start[ei]), int(e_end[ei])) for ei in range(int(t_ex[ti]), int(t_ex[ti + 1])))
                    if exons:
                        tx_dict[tid] = TranscriptAnnotation(
                            transcript_id=tid,
                            gene_id=gid,
                            chrom=chrom,
                            strand=strand_val,
                            exons=exons,
                            is_discovered=False,
                        )
                if tx_dict:
                    genes[gid] = GeneModelRecord(
                        gene_id=gid, gene_name=gname, chrom=chrom, strand=strand_val, transcripts=tx_dict
                    )
    return genes


class DictionaryBuilder:
    """Unified builder merging reference transcripts (`.iqidx`/GTF) and discovered models."""

    def __init__(self, config) -> None:
        self.config = config

    def build(
        self,
        *,
        gene_ids: Optional[Iterable[str]] = None,
    ) -> Dict[str, GeneModelRecord]:
        merged: Dict[str, GeneModelRecord] = {}
        target = set(gene_ids) if gene_ids is not None else None

        if getattr(self.config, "iqidx_path", None):
            merged = load_iqidx_models(self.config.iqidx_path, gene_filter=target)
        elif getattr(self.config, "reference_gtf", None):
            merged = load_gtf_models_streaming(self.config.reference_gtf, gene_filter=target, is_discovered=False)

        if not getattr(self.config, "ref_only", False):
            for sample in getattr(self.config, "samples", []):
                disc_gtf = getattr(sample, "transcript_models_gtf", None)
                if disc_gtf is None:
                    continue
                disc_map = load_gtf_models_streaming(disc_gtf, gene_filter=target, is_discovered=True)
                for gid, grec in disc_map.items():
                    if gid not in merged:
                        merged[gid] = grec
                    else:
                        for tid, tx in grec.transcripts.items():
                            key = tid if tid not in merged[gid].transcripts else f"{sample.prefix}:{tid}"
                            merged[gid].transcripts[key] = tx
        return merged

    @staticmethod
    def load_transcript_counts(tsv_path: str | Path) -> Dict[str, Dict[str, float]]:
        """Parse transcript counts / TPM table into {tx_id: {group: count}}."""
        counts: Dict[str, Dict[str, float]] = {}
        with open(tsv_path, "r", encoding="utf-8") as fh:
            header_cols: Optional[List[str]] = None
            for raw in fh:
                line = raw.rstrip("\n")
                if not line:
                    continue
                if line.startswith("#"):
                    if header_cols is None:
                        header_cols = line.lstrip("#").split("\t")
                    continue
                cols = line.split("\t")
                if header_cols is None:
                    header_cols = cols
                    continue
                if cols[0].startswith("__"):
                    continue
                row: Dict[str, float] = {}
                for col_name, val in zip(header_cols[1:], cols[1:]):
                    try:
                        row[col_name] = float(val)
                    except ValueError:
                        row[col_name] = 0.0
                counts[cols[0]] = row
        return counts
