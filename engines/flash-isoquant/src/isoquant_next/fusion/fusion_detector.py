"""Stateless fusion detector backed by `.iqidx`, sorted intervals, minimap2 CLI, and Cython SW."""

from __future__ import annotations

from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from functools import lru_cache
import logging
import os
import re
import shutil
import subprocess
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pysam

from isoquant_next.annotation.iqidx import MmapAnnotationIndex, load_or_build_iqidx
from isoquant_next.barcodes._cy_barcode_umi import sw_align_bytes
from isoquant_next.fusion.fusion_metadata import FusionMetadata
from isoquant_next.fusion.fusion_validator import FusionValidator
from isoquant_next.fusion.genomic_interval_index import GenomicIntervalIndex

logger = logging.getLogger("IsoQuant")

ANTISENSE_SUFFIX_RE = re.compile(r"-(AS\d+|DT|DIVERGENT|NAT)$", re.IGNORECASE)
_CIGAR_TOKEN_RE = re.compile(r"(\d+)([MIDNSHP=XB])")
_MITO_CHRS = frozenset({"chrM", "MT", "M", "chrMT", "mitochondrion"})
_RC_TRANS = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass(frozen=True, slots=True)
class AlignmentHit:
    """PAF / Smith-Waterman alignment record matching mappy's hit fields."""

    ctg: str
    r_st: int
    r_en: int
    q_st: int
    q_en: int
    mapq: int
    mlen: int
    alen: int


def _revcomp(seq: str) -> str:
    return seq.translate(_RC_TRANS)[::-1]


def _encode_ascii(seq: str) -> np.ndarray:
    return np.frombuffer(seq.upper().encode("ascii", errors="replace"), dtype=np.uint8)


class _BoundedLRUDict:
    """Per-instance bounded LRU cache (no shared module-level global mutable state)."""

    __slots__ = ("_data", "_maxsize")

    def __init__(self, maxsize: int) -> None:
        self._data: OrderedDict = OrderedDict()
        self._maxsize = max(1, int(maxsize))

    def get(self, key):
        if key in self._data:
            self._data.move_to_end(key)
            return self._data[key]
        return None

    def __contains__(self, key) -> bool:
        return key in self._data

    def __getitem__(self, key):
        val = self._data[key]
        self._data.move_to_end(key)
        return val

    def __setitem__(self, key, value) -> None:
        if key in self._data:
            self._data.move_to_end(key)
        self._data[key] = value
        while len(self._data) > self._maxsize:
            self._data.popitem(last=False)

    def clear(self) -> None:
        self._data.clear()

    def __len__(self) -> int:
        return len(self._data)


@lru_cache(maxsize=65536)
def _aligned_len_from_cigarstring_cached(cigar: str) -> int:
    total = 0
    for m in _CIGAR_TOKEN_RE.finditer(cigar):
        if m.group(2) in ("M", "=", "X"):
            total += int(m.group(1))
    return total


class _Minimap2SWAligner:
    """Zero-mappy realignment backend using minimap2 CLI (PAF) and Cython SW for local windows."""

    def __init__(
        self,
        reference_fasta: str,
        minimap2_bin: Optional[str] = None,
        threads: int = 1,
    ) -> None:
        self.reference_fasta = str(reference_fasta)
        self.minimap2_bin = minimap2_bin or shutil.which("minimap2")
        self.threads = max(1, int(threads))
        self._ref_fasta: Optional[pysam.FastaFile] = None

    def close(self) -> None:
        if self._ref_fasta is not None:
            self._ref_fasta.close()
            self._ref_fasta = None

    def _get_fasta(self) -> pysam.FastaFile:
        if self._ref_fasta is None:
            self._ref_fasta = pysam.FastaFile(self.reference_fasta)
        return self._ref_fasta

    def map_batch(self, seqs: Dict[str, str]) -> Dict[str, Tuple[AlignmentHit, ...]]:
        """Map named sequences against `reference_fasta` in a single minimap2 CLI invocation."""
        if not seqs or not self.minimap2_bin or not os.path.isfile(self.reference_fasta):
            return {k: () for k in seqs}

        fasta_payload = "".join(f">{qname}\n{seq}\n" for qname, seq in seqs.items() if seq)
        if not fasta_payload:
            return {k: () for k in seqs}

        cmd = [
            self.minimap2_bin,
            "-c",
            "-t", str(self.threads),
            "-N", "5",
            "--secondary=yes",
            self.reference_fasta,
            "-",
        ]
        try:
            proc = subprocess.run(
                cmd,
                input=fasta_payload,
                text=True,
                capture_output=True,
                check=False,
            )
        except OSError:
            return {k: () for k in seqs}

        if proc.returncode != 0:
            logger.debug("minimap2 realignment failed rc=%s: %s", proc.returncode, proc.stderr.strip())
            return {k: () for k in seqs}

        grouped: Dict[str, List[AlignmentHit]] = {k: [] for k in seqs}
        for line in proc.stdout.splitlines():
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 12:
                continue
            qname = cols[0]
            if qname not in grouped:
                continue
            try:
                grouped[qname].append(
                    AlignmentHit(
                        ctg=cols[5],
                        r_st=int(cols[7]),
                        r_en=int(cols[8]),
                        q_st=int(cols[2]),
                        q_en=int(cols[3]),
                        mapq=int(cols[11]),
                        mlen=int(cols[9]),
                        alen=int(cols[10]),
                    )
                )
            except ValueError:
                continue
        return {k: tuple(v) for k, v in grouped.items()}

    def map(self, seq: str) -> Iterable[AlignmentHit]:
        if not seq:
            return ()
        return self.map_batch({"q0": seq}).get("q0", ())

    def align_to_breakpoint_windows(
        self,
        seq: str,
        windows: List[Tuple[str, int, int]],
        min_match_len: int = 30,
    ) -> List[AlignmentHit]:
        """Targeted Cython `sw_align_bytes` fallback around candidate loci."""
        if not seq or len(seq) < min_match_len or not os.path.isfile(self.reference_fasta):
            return []
        try:
            ref = self._get_fasta()
        except Exception:
            return []
        q_fwd = _encode_ascii(seq)
        q_rev = _encode_ascii(_revcomp(seq))
        hits: List[AlignmentHit] = []
        for chrom, win_start, win_end in windows:
            try:
                ref_len = ref.get_reference_length(chrom)
            except (KeyError, ValueError):
                continue
            s0 = max(0, int(win_start) - 1)
            e0 = min(ref_len, int(win_end))
            if e0 <= s0:
                continue
            ref_seq = ref.fetch(chrom, s0, e0)
            if not ref_seq:
                continue
            r_bytes = _encode_ascii(ref_seq)
            best_score = 0
            best_res = None
            for qb in (q_fwd, q_rev):
                score, r_st, r_en, q_st, q_en = sw_align_bytes(qb, r_bytes, 2, 2, 3, 1)
                if score > best_score:
                    best_score = score
                    best_res = (r_st, r_en, q_st, q_en, score)
            if best_res is None:
                continue
            r_st, r_en, q_st, q_en, score = best_res
            mlen = max(0, q_en - q_st + 1)
            alen = max(mlen, r_en - r_st + 1)
            if mlen >= min_match_len:
                mapq = min(60, max(10, int((score * 60) // max(1, 2 * mlen))))
                hits.append(
                    AlignmentHit(
                        ctg=chrom,
                        r_st=s0 + r_st,
                        r_en=s0 + r_en + 1,
                        q_st=q_st,
                        q_en=q_en + 1,
                        mapq=mapq,
                        mlen=mlen,
                        alen=alen,
                    )
                )
        return hits


class FusionDetector:
    """Stateless fusion detector backed by `.iqidx`, sorted intervals, minimap2 CLI, and Cython SW."""

    def __init__(
        self,
        bam_path: str,
        gene_db_path: Any,
        reference_fasta: Optional[str] = None,
        *,
        aligner: Optional[Any] = None,
        cache_dir: Optional[str] = None,
    ) -> None:
        self.bam_path = bam_path
        self.genedb_path = gene_db_path
        self._owns_db = False

        if isinstance(gene_db_path, MmapAnnotationIndex):
            self.db: Any = gene_db_path
        elif isinstance(gene_db_path, str):
            iqidx_path = load_or_build_iqidx(gene_db_path, cache_dir=cache_dir)
            self.db = MmapAnnotationIndex(iqidx_path)
            self._owns_db = True
        else:
            self.db = gene_db_path

        self.fusion_candidates: Dict[str, Set[str]] = defaultdict(set)
        self.fusion_breakpoints: Dict[str, Dict[Tuple, int]] = defaultdict(lambda: defaultdict(int))
        self.reference_fasta = reference_fasta
        if aligner is not None:
            self.aligner = aligner
            self._owns_aligner = False
        elif reference_fasta:
            self.aligner = _Minimap2SWAligner(reference_fasta)
            self._owns_aligner = True
        else:
            self.aligner = None
            self._owns_aligner = False

        self.fusion_metadata: Dict[str, dict] = {}
        self.fusion_assigned_pairs: Dict[str, Dict[str, Tuple[str, str]]] = defaultdict(dict)
        self.fusion_read_scores: Dict[str, Dict[str, Tuple[float, float]]] = defaultdict(dict)
        self._resolved_name_cache: Dict[str, str] = {}
        self._symbol_biotype_cache: Dict[str, str] = {}
        self.exon_cache: Dict[str, List[Tuple[int, int]]] = {}

        # Per-instance bounded caches (zero module-level global mutable state)
        self._context_cache = _BoundedLRUDict(200_000)
        self._gene_assignment_cache = _BoundedLRUDict(500_000)
        self._aligner_map_cache = _BoundedLRUDict(10_000)

        self.interval_index = GenomicIntervalIndex(self.db)
        self.debug = False
        self._ref_fasta: Optional[pysam.FastaFile] = None

    def close(self) -> None:
        if getattr(self, "_ref_fasta", None) is not None:
            self._ref_fasta.close()
            self._ref_fasta = None
        if getattr(self, "_owns_aligner", False) and self.aligner is not None and hasattr(self.aligner, "close"):
            self.aligner.close()
        if getattr(self, "_owns_db", False) and hasattr(self.db, "close"):
            self.db.close()

    def __enter__(self) -> "FusionDetector":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def safe_reference_start(self, read) -> Optional[int]:
        try:
            return read.reference_start + 1 if read.reference_start is not None else None
        except Exception:
            return None

    def compute_aligned_length(self, read) -> int:
        try:
            if getattr(read, "cigartuples", None):
                return self.aligned_len_from_cigartuples(read.cigartuples)
            return self.aligned_len_from_cigarstring(getattr(read, "cigarstring", None))
        except Exception:
            return 0

    def clear_state(self) -> None:
        self.fusion_candidates.clear()
        self.fusion_breakpoints.clear()
        self.fusion_metadata.clear()
        self.fusion_assigned_pairs.clear()
        self.fusion_read_scores.clear()

    def clear_caches(self) -> None:
        self._context_cache.clear()
        self._gene_assignment_cache.clear()
        self._aligner_map_cache.clear()
        self._resolved_name_cache.clear()
        self._symbol_biotype_cache.clear()
        self.exon_cache.clear()

    def _lookup_db_feature(self, feature_id: str) -> Any:
        if self.interval_index is not None and isinstance(self.db, MmapAnnotationIndex):
            return self.interval_index.lookup_feature(feature_id)
        return self.db[feature_id]

    def _build_exon_cache(self) -> None:
        if self.exon_cache:
            return
        try:
            genes = (
                self.interval_index.get_all_genes()
                if self.interval_index is not None
                else list(self.db.features_of_type("gene"))
            )
            for gene in genes:
                gene_id = getattr(gene, "id", None)
                if not gene_id:
                    continue
                self.exon_cache[gene_id] = [
                    (int(ex.start), int(ex.end))
                    for ex in self.interval_index.get_exons_of_gene(gene)
                ]
        except Exception as exc:
            logger.warning("Failed to build exon cache: %s", exc)

    def _get_cached_exons(self, gene) -> List[Tuple[int, int]]:
        gene_id = getattr(gene, "id", None)
        if not gene_id:
            return []
        if gene_id in self.exon_cache:
            return self.exon_cache[gene_id]
        try:
            if self.interval_index is not None:
                exons = self.interval_index.get_exons_of_gene(gene)
            else:
                exons = list(self.db.children(gene, featuretype="exon", order_by="start"))
            exon_list = [(int(ex.start), int(ex.end)) for ex in exons]
            self.exon_cache[gene_id] = exon_list
            return exon_list
        except Exception:
            return []

    def _context_query(self, chrom: str, pos: int) -> Tuple[Optional[str], str]:
        cache_key = (chrom, pos)
        cached = self._context_cache.get(cache_key)
        if cached is not None:
            return cached
        try:
            if self.interval_index is not None:
                genes = self.interval_index.get_genes_at(chrom, pos)
                exons = self.interval_index.get_exons_at(chrom, pos)
            else:
                genes = list(self.db.region(region=(chrom, pos, pos), featuretype="gene"))
                exons = list(self.db.region(region=(chrom, pos, pos), featuretype="exon"))
            if genes:
                attrs = getattr(genes[0], "attributes", {}) or {}
                gene_name = attrs.get("gene_name", [genes[0].id])[0]
                result = (gene_name, "exonic" if exons else "intronic")
            else:
                result = (None, "exonic" if exons else "intergenic")
        except Exception:
            result = (None, "unknown")
        self._context_cache[cache_key] = result
        return result

    def _get_parent_gene_symbol(self, feature) -> str:
        attrs = getattr(feature, "attributes", {}) or {}
        if getattr(feature, "featuretype", None) == "gene":
            return attrs.get("gene_name", [feature.id])[0]
        if "gene_id" in attrs:
            gid = attrs["gene_id"][0]
            try:
                parent = self._lookup_db_feature(gid)
                pattrs = getattr(parent, "attributes", {}) or {}
                return pattrs.get("gene_name", [parent.id])[0]
            except Exception:
                return gid
        if "Parent" in attrs:
            pid = attrs["Parent"][0]
            try:
                parent = self._lookup_db_feature(pid)
                pattrs = getattr(parent, "attributes", {}) or {}
                return pattrs.get("gene_name", [parent.id])[0]
            except Exception:
                return pid
        return attrs.get("gene_name", [feature.id])[0] if attrs else feature.id

    def _lookup_feature_by_id(self, name_or_id: str):
        try:
            feat = self._lookup_db_feature(name_or_id)
            return feat, self._get_parent_gene_symbol(feat)
        except Exception:
            return None, name_or_id

    def _fallback_gene_lookup(self, name_or_id: str) -> str:
        if isinstance(self.db, MmapAnnotationIndex):
            matches = self.interval_index.find_genes_by_name(name_or_id)
            if matches:
                attrs = getattr(matches[0], "attributes", {}) or {}
                return attrs.get("gene_name", [matches[0].id])[0]
            return name_or_id
        try:
            genes = list(self.db.features_of_type("gene", id=name_or_id))
            if not genes:
                return name_or_id
            gene = genes[0]
            attrs = getattr(gene, "attributes", {}) or {}
            return attrs.get("gene_name", [gene.id])[0]
        except (AttributeError, KeyError, TypeError):
            return name_or_id

    def resolve_gene_name(self, name_or_id: Optional[str]) -> Optional[str]:
        if not name_or_id:
            return name_or_id
        cached = self._resolved_name_cache.get(name_or_id)
        if cached is not None:
            return cached
        feat, resolved = self._lookup_feature_by_id(name_or_id)
        if feat is None:
            resolved = self._fallback_gene_lookup(name_or_id)
        self._resolved_name_cache[name_or_id] = resolved
        return resolved

    def has_antisense_suffix(self, gene_name: Optional[str]) -> bool:
        if not gene_name or not isinstance(gene_name, str):
            return False
        return ANTISENSE_SUFFIX_RE.search(gene_name) is not None

    def strip_antisense_suffix(self, gene_name: Optional[str]) -> Optional[str]:
        if not gene_name or not isinstance(gene_name, str):
            return gene_name
        return ANTISENSE_SUFFIX_RE.sub("", gene_name)

    def _lookup_gene_by_name(self, gene_name: str, chrom: str, pos: int):
        try:
            gene_features = (
                self.interval_index.get_genes_at(chrom, pos, window=500)
                if self.interval_index is not None
                else list(self.db.region(region=(chrom, max(1, pos - 500), pos + 500), featuretype="gene"))
            )
            if not gene_features:
                return None, None, None
            selected = None
            for gf in gene_features:
                attrs = getattr(gf, "attributes", {}) or {}
                gname = attrs.get("gene_name", [getattr(gf, "id", None)])[0]
                if gname == gene_name:
                    selected = gf
                    break
            if selected is None:
                selected = gene_features[0]
            return selected, int(selected.start), int(selected.end)
        except Exception:
            return None, None, None

    def _compute_gene_score(self, target, pos: int, gstart=None, gend=None) -> float:
        if gstart is None or gend is None:
            raise ValueError("gstart and gend required")
        return 1.0 if gstart <= pos <= gend else 0.0

    def _get_genes_at(self, chrom: str, pos: int, window: int):
        if self.interval_index is None:
            return None
        try:
            return self.interval_index.get_genes_at(chrom, pos, window=window)
        except Exception:
            return None

    def _collect_exonic_genes(self, genes, pos: int):
        exonic_genes = []
        for gene in genes:
            attrs = getattr(gene, "attributes", {}) or {}
            gene_name = attrs.get("gene_name", [getattr(gene, "id", None)])[0]
            gene_type = (
                attrs.get("gene_type", [None])[0]
                or attrs.get("gene_biotype", [None])[0]
                or attrs.get("transcript_biotype", [None])[0]
            )
            try:
                gene_start = int(gene.start)
                gene_end = int(gene.end)
                exons = self._get_cached_exons(gene)
            except (AttributeError, TypeError, ValueError):
                continue
            if any(start <= pos <= end for start, end in exons):
                exonic_genes.append((gene, gene_name, gene_type, gene_start, gene_end, exons))
        return exonic_genes

    def _compute_gene_distances(self, g, pos: int):
        exon_min_dist = None
        boundary_min_dist = None
        exonic_hit = False
        try:
            exons = self._get_cached_exons(g)
            for es, ee in exons:
                if es <= pos <= ee:
                    d_span = 0
                    exonic_hit = True
                else:
                    d_span = min(abs(pos - es), abs(pos - ee))
                exon_min_dist = d_span if exon_min_dist is None else min(exon_min_dist, d_span)
                d_bound = min(abs(pos - es), abs(pos - ee))
                boundary_min_dist = d_bound if boundary_min_dist is None else min(boundary_min_dist, d_bound)
        except Exception:
            exons = []
        gstart = int(getattr(g, "start", 0))
        gend = int(getattr(g, "end", 0))
        body_dist = 0 if gstart <= pos <= gend else min(abs(pos - gstart), abs(pos - gend))
        return exon_min_dist, boundary_min_dist, exonic_hit, body_dist, exons

    def _build_gene_key(self, gtype, score, boundary_min_dist, exon_min_dist, body_dist, glen, gname):
        inf = 10**9
        bnd = boundary_min_dist if boundary_min_dist is not None else inf
        exd = exon_min_dist if exon_min_dist is not None else inf
        bod = body_dist if body_dist else 0
        return (
            1 if gtype == "protein_coding" else 0,
            score,
            -bnd,
            -exd,
            -bod,
            -min(glen, 500000),
            gname or "",
        )

    def _score_exonic_genes(self, exonic_genes, pos: int):
        if len(exonic_genes) == 1:
            g, gname, _gtype, gstart, gend, _exons = exonic_genes[0]
            return gname, self._compute_gene_score(g, pos, gstart=gstart, gend=gend)
        best_gene = None
        best_score = None
        best_key = None
        for g, gname, gtype, gstart, gend, exons in exonic_genes:
            boundary_distances = [min(abs(pos - start), abs(pos - end)) for start, end in exons]
            boundary_min_dist = min(boundary_distances) if boundary_distances else 0
            score = self._compute_gene_score(g, pos, gstart=gstart, gend=gend)
            glen = max(1, gend - gstart + 1)
            key = self._build_gene_key(gtype, score, boundary_min_dist, 0, 0, glen, gname)
            if best_key is None or key > best_key:
                best_key = key
                best_gene = gname or getattr(g, "id", None)
                best_score = score
        return best_gene, best_score if best_score is not None else 0.0

    def _score_intronic_genes(self, genes, pos: int):
        best_gene = None
        best_score = None
        best_key = None
        for g in genes:
            attrs = getattr(g, "attributes", {}) or {}
            gname = attrs.get("gene_name", [getattr(g, "id", None)])[0]
            gtype = (
                attrs.get("gene_type", [None])[0]
                or attrs.get("gene_biotype", [None])[0]
                or attrs.get("transcript_biotype", [None])[0]
            )
            gstart = int(getattr(g, "start", 0))
            gend = int(getattr(g, "end", 0))
            exon_min_dist, boundary_min_dist, _exhit, body_dist, _exons = self._compute_gene_distances(g, pos)
            score = self._compute_gene_score(g, pos, gstart=gstart, gend=gend)
            glen = max(1, gend - gstart + 1)
            key = self._build_gene_key(gtype, score, boundary_min_dist, exon_min_dist, body_dist, glen, gname)
            if best_key is None or key > best_key:
                best_key = key
                best_gene = gname or getattr(g, "id", None)
                best_score = score
        return best_gene, best_score if best_score is not None else 0.0

    def assign_fusion_gene(self, chrom: str, pos: int, window: int = 1000) -> Tuple[Optional[str], float]:
        genes = self._get_genes_at(chrom, pos, window)
        if not genes:
            return None, 0.0
        exonic_genes = self._collect_exonic_genes(genes, pos)
        if exonic_genes:
            return self._score_exonic_genes(exonic_genes, pos)
        return self._score_intronic_genes(genes, pos)

    def assign_fusion_gene_cached(self, chrom: str, pos: int) -> Tuple[Optional[str], float]:
        cache_key = (chrom, pos)
        cached = self._gene_assignment_cache.get(cache_key)
        if cached is not None:
            return cached
        result = self.assign_fusion_gene(chrom, pos)
        self._gene_assignment_cache[cache_key] = result
        return result

    def _passes_read_filters(self, read, min_al_len_primary: int, min_sa_mapq: int):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            return False, None, 0
        if self.compute_aligned_length(read) < min_al_len_primary:
            return False, None, 0
        if getattr(read, "mapping_quality", 0) < min_sa_mapq:
            return False, None, 0
        clip_side, clip_len = self.detect_softclip(read)
        return True, clip_side, clip_len

    def _filter_sa_entries(self, sa_entries, min_al_len_sa: int, min_sa_mapq: int):
        filtered = []
        for sa_chr, sa_pos, sa_strand, sa_cigar, sa_mapq, sa_nm in sa_entries:
            if sa_pos is None or sa_cigar is None:
                continue
            sa_al_len = self.aligned_len_from_cigarstring(sa_cigar)
            if sa_al_len < min_al_len_sa or sa_mapq < min_sa_mapq:
                continue
            filtered.append((sa_chr, sa_pos, sa_strand, sa_cigar, sa_mapq, sa_nm, sa_al_len))
        if len(filtered) > 3:
            filtered.sort(key=lambda x: (x[4], x[6]), reverse=True)
            filtered = filtered[:3]
        return filtered

    def _process_sa_entries(self, read, sa_entries, clip_side, jitter_window: int, min_sa_mapq: int) -> bool:
        seen_pairs = set()
        sa_used = False
        for sa_chr, sa_pos, _sa_strand, sa_cigar, _sa_mapq, _sa_nm, _sa_al_len in sa_entries:
            left_chr, left_pos, _r_chr, right_pos = self.estimate_breakpoint(
                read, sa_pos, clip_side=clip_side, sa_cigar=sa_cigar
            )
            if left_pos is None or right_pos is None:
                continue
            key = (left_chr, left_pos // jitter_window, sa_chr, right_pos // jitter_window)
            if key in seen_pairs:
                continue
            seen_pairs.add(key)
            left_gene, left_score = self.assign_fusion_gene_cached(left_chr, left_pos)
            right_gene, right_score = self.assign_fusion_gene_cached(sa_chr, right_pos)
            if left_gene is None or right_gene is None:
                continue
            if left_gene != right_gene:
                self.record_fusion(
                    left_gene,
                    right_gene,
                    read.query_name,
                    left_chr,
                    left_pos,
                    sa_chr,
                    right_pos,
                    left_score=left_score,
                    right_score=right_score,
                )
                sa_used = True
        return sa_used

    def detect_softclip(self, read, min_len: int = 50) -> Tuple[Optional[str], int]:
        cigartuples = getattr(read, "cigartuples", None)
        if cigartuples:
            left_clip = cigartuples[0][0] == 4 and cigartuples[0][1] >= min_len
            right_clip = cigartuples[-1][0] == 4 and cigartuples[-1][1] >= min_len
            if left_clip and right_clip:
                return None, 0
            if left_clip:
                return "left", cigartuples[0][1]
            if right_clip:
                return "right", cigartuples[-1][1]
            return None, 0
        cigar = getattr(read, "cigarstring", "") or ""
        m_left = re.match(r"^(\d+)S", cigar)
        m_right = re.search(r"(\d+)S$", cigar)
        left_len = int(m_left.group(1)) if m_left else 0
        right_len = int(m_right.group(1)) if m_right else 0
        if left_len >= min_len and right_len >= min_len:
            return None, 0
        if left_len >= min_len:
            return "left", left_len
        if right_len >= min_len:
            return "right", right_len
        return None, 0

    def aligned_len_from_cigartuples(self, cigartuples) -> int:
        if not cigartuples:
            return 0
        return sum(length for op, length in cigartuples if op in (0, 7, 8))

    def aligned_len_from_cigarstring(self, cigar: Optional[str]) -> int:
        if not cigar:
            return 0
        return _aligned_len_from_cigarstring_cached(cigar)

    def _cached_aligner_map(self, seq: str) -> Tuple[Any, ...]:
        cached = self._aligner_map_cache.get(seq)
        if cached is not None:
            return cached
        try:
            result = tuple(self.aligner.map(seq))
        except Exception:
            result = ()
        self._aligner_map_cache[seq] = result
        return result

    def _is_softclip_orientation_valid(self, read, clip_side, sa_chr, sa_pos, tol=10) -> bool:
        if clip_side is None or read.reference_name != sa_chr:
            return True
        prim_start = self.safe_reference_start(read)
        prim_end = read.reference_end
        if prim_start is None or prim_end is None:
            return True
        if not read.is_reverse:
            if clip_side == "right":
                return sa_pos > prim_end - tol
            if clip_side == "left":
                return sa_pos < prim_start + tol
        else:
            if clip_side == "left":
                return sa_pos > prim_end - tol
            if clip_side == "right":
                return sa_pos < prim_start + tol
        return True

    def realign_clipped_seq(self, seq: str):
        if not seq or not self.aligner:
            return None
        try:
            for hit in self._cached_aligner_map(seq):
                r_en = getattr(hit, "r_en", None)
                r_st = getattr(hit, "r_st", None)
                if r_en is not None and r_st is not None:
                    hit_len = r_en - r_st
                else:
                    hit_len = getattr(hit, "mlen", None) or getattr(hit, "alen", None) or 0
                if hit_len and hit_len >= 50:
                    return hit
            return None
        except Exception:
            return None

    def _safe_gene_token(self, g) -> str:
        if not g or (isinstance(g, str) and not g.strip()):
            return "intergenic"
        return str(g)

    def record_fusion(
        self,
        context1: Optional[str],
        context2: Optional[str],
        read_name: str,
        chrom1: str,
        pos1: int,
        chrom2: str,
        pos2: int,
        left_score: Optional[float] = None,
        right_score: Optional[float] = None,
    ) -> None:
        if self._is_mitochondrial_candidate(chrom1, chrom2, context1, context2):
            return
        left_gene = self.normalize_gene_label(context1)
        right_gene = self.normalize_gene_label(context2)
        sorted_genes = sorted([left_gene, right_gene])
        fusion_key = f"{sorted_genes[0]}--{sorted_genes[1]}"
        self.fusion_assigned_pairs[fusion_key][read_name] = (left_gene, right_gene)
        self.fusion_candidates[fusion_key].add(read_name)
        bp = (chrom1, int(pos1), chrom2, int(pos2))
        self.fusion_breakpoints[fusion_key][bp] += 1
        if left_score is not None and right_score is not None:
            self.fusion_read_scores[fusion_key][read_name] = (left_score, right_score)
        meta = self.fusion_metadata.setdefault(
            fusion_key,
            {"supporting_reads": set(), "consensus_bp": None, "left_gene": None, "right_gene": None, "support": 0},
        )
        meta["supporting_reads"].add(read_name)
        meta["support"] = len(meta["supporting_reads"])

    def normalize_gene_label(self, gene: Optional[str]) -> str:
        if not gene:
            return "intergenic"
        gene = self.resolve_gene_name(gene)
        gene = self.canonical_locus_name(gene)
        if gene is None or gene.startswith(("ENSG", "RP11")):
            return "intergenic"
        return gene

    def build_metadata(self, min_support: int = 1) -> None:
        FusionValidator(self).filter_early_non_coding_genes()
        try:
            FusionMetadata(self).process_all(min_support=min_support)
        except Exception as exc:
            logger.error("Fusion metadata processing failed: %s", exc)
            logger.debug("Traceback:", exc_info=True)

    def cluster_breakpoints(self, bp_counts: dict, window: int):
        if not bp_counts:
            return None
        by_pair = defaultdict(list)
        for (c1, p1, c2, p2), w in bp_counts.items():
            by_pair[(c1, c2)].append((p1, p2, w))

        best_pair = None
        best = (None, None, 0)
        for pair, items in by_pair.items():
            best_total = 0
            best_p1 = 0
            best_p2 = 0
            for p1_i, p2_i, _ in items:
                sum_w = sum_p1w = sum_p2w = 0
                for p1_j, p2_j, wj in items:
                    if abs(p1_i - p1_j) <= window and abs(p2_i - p2_j) <= window:
                        sum_w += wj
                        sum_p1w += p1_j * wj
                        sum_p2w += p2_j * wj
                if sum_w > best_total:
                    best_total = sum_w
                    best_p1 = sum_p1w // sum_w
                    best_p2 = sum_p2w // sum_w
            if best_total > best[2]:
                best_pair = pair
                best = (best_p1, best_p2, best_total)

        if best_pair is None:
            return None
        c1, c2 = best_pair
        p1, p2, total = best
        return (c1, p1, c2, p2), total

    def parse_sa_entries(self, sa_tag: Optional[str]) -> List[Tuple]:
        entries: List[Tuple] = []
        if not sa_tag:
            return entries
        for sa in sa_tag.split(";"):
            if not sa:
                continue
            fields = sa.split(",")
            if len(fields) < 6:
                chrom = fields[0]
                pos = int(fields[1]) if len(fields) > 1 and fields[1].isdigit() else None
                strand = fields[2] if len(fields) > 2 else "+"
                cigar = fields[3] if len(fields) > 3 else None
                mapq = int(fields[4]) if len(fields) > 4 and fields[4].isdigit() else 0
                nm = int(fields[5]) if len(fields) > 5 and fields[5].isdigit() else 0
            else:
                chrom = fields[0]
                pos = int(fields[1])
                strand = fields[2]
                cigar = fields[3]
                mapq = int(fields[4])
                nm = int(fields[5])
            entries.append((chrom, pos, strand, cigar, mapq, nm))
        return entries

    def estimate_breakpoint(
        self,
        read,
        sa_pos: Optional[int],
        clip_side: Optional[str] = None,
        sa_cigar: Optional[str] = None,
    ) -> tuple:
        left_chr = read.reference_name
        prim_start = self.safe_reference_start(read)
        prim_end = int(read.reference_end) if read.reference_end is not None else None
        if clip_side == "right" and prim_end is not None:
            left_pos = prim_end
        elif clip_side == "left" and prim_start is not None:
            left_pos = prim_start
        else:
            left_pos = prim_end if prim_end is not None else prim_start
        right_chr = read.reference_name
        right_pos = sa_pos
        if sa_cigar:
            sa_al_len = self.aligned_len_from_cigarstring(sa_cigar)
            if clip_side == "left" and sa_pos is not None and sa_al_len:
                right_pos = sa_pos + sa_al_len
            else:
                right_pos = sa_pos
        return (
            left_chr,
            int(left_pos) if left_pos is not None else None,
            right_chr,
            int(right_pos) if right_pos is not None else None,
        )

    def realign_softclip(
        self,
        read,
        clip_side: str,
        clip_len: int,
        seen_pairs: set,
        jitter_window: int = 25,
        min_sa_mapq: int = 20,
    ) -> None:
        seq = read.query_sequence
        if not seq:
            return
        clipped_seq = seq[:clip_len] if clip_side == "left" else seq[-clip_len:]
        best_hit = self.realign_clipped_seq(clipped_seq)
        if not best_hit:
            return
        if getattr(best_hit, "mapq", 0) < min_sa_mapq:
            return
        sa_chr = getattr(best_hit, "ctg", None)
        sa_pos = getattr(best_hit, "r_st", None)
        if sa_chr is None or sa_pos is None:
            return
        sa_pos = int(sa_pos) + 1
        left_chr, left_pos, _rchr, right_pos = self.estimate_breakpoint(
            read, sa_pos, clip_side=clip_side, sa_cigar=None
        )
        if left_pos is None or right_pos is None:
            return
        key = (left_chr, left_pos // jitter_window, sa_chr, right_pos // jitter_window)
        if key in seen_pairs:
            return
        left_gene, left_score = self.assign_fusion_gene_cached(left_chr, left_pos)
        right_gene, right_score = self.assign_fusion_gene_cached(sa_chr, right_pos)
        if left_gene is not None and right_gene is not None and left_gene != right_gene:
            self.record_fusion(
                left_gene,
                right_gene,
                read.query_name,
                left_chr,
                left_pos,
                sa_chr,
                right_pos,
                left_score=left_score,
                right_score=right_score,
            )

    def detect_fusions(
        self,
        min_al_len_primary: int = 50,
        min_al_len_sa: int = 40,
        min_sa_mapq: int = 10,
        min_softclip_len: int = 30,
        jitter_window: int = 50,
    ) -> None:
        logger.info("Starting fusion detection on %s", self.bam_path)
        pending_softclips: List[Tuple[Any, str, int, str]] = []
        unique_clipped: Dict[str, str] = {}

        with pysam.AlignmentFile(self.bam_path, "rb") as bam:
            for read in bam:
                passes, clip_side, clip_len = self._passes_read_filters(read, min_al_len_primary, min_sa_mapq)
                if not passes:
                    continue
                if read.has_tag("SA"):
                    sa_entries = self.parse_sa_entries(read.get_tag("SA"))
                    if sa_entries:
                        filtered = self._filter_sa_entries(sa_entries, min_al_len_sa, min_sa_mapq)
                        self._process_sa_entries(read, filtered, clip_side, jitter_window, min_sa_mapq)
                if clip_side and clip_len >= min_softclip_len and self.aligner is not None:
                    seq = getattr(read, "query_sequence", None)
                    if seq:
                        cseq = seq[:clip_len] if clip_side == "left" else seq[-clip_len:]
                        if hasattr(self.aligner, "map_batch") and cseq not in self._aligner_map_cache:
                            unique_clipped.setdefault(f"sc_{len(unique_clipped)}", cseq)
                            pending_softclips.append((read, clip_side, clip_len, cseq))
                            continue
                    self.realign_softclip(
                        read, clip_side, clip_len, set(), jitter_window=jitter_window, min_sa_mapq=min_sa_mapq
                    )

        if unique_clipped and self.aligner is not None and hasattr(self.aligner, "map_batch"):
            batch_hits = self.aligner.map_batch(unique_clipped)
            for qname, cseq in unique_clipped.items():
                self._aligner_map_cache[cseq] = batch_hits.get(qname, ())
            for read, clip_side, clip_len, _cseq in pending_softclips:
                self.realign_softclip(
                    read, clip_side, clip_len, set(), jitter_window=jitter_window, min_sa_mapq=min_sa_mapq
                )

    def reconstruct_fusion_transcript(
        self, c1: str, p1: int, c2: str, p2: int, exon_padding: int = 100
    ):
        try:
            if self.interval_index is not None:
                left_exons = self.interval_index.get_exons_at(c1, p1, window=exon_padding)
                right_exons = self.interval_index.get_exons_at(c2, p2, window=exon_padding)
            else:
                left_exons = list(self.db.region(region=(c1, max(1, p1 - exon_padding), p1 + exon_padding), featuretype="exon"))
                right_exons = list(self.db.region(region=(c2, max(1, p2 - exon_padding), p2 + exon_padding), featuretype="exon"))
            if not left_exons or not right_exons:
                return None, None, None
            left_coords = sorted((int(e.start), int(e.end)) for e in left_exons)
            right_coords = sorted((int(e.start), int(e.end)) for e in right_exons)
            if self._ref_fasta is None:
                self._ref_fasta = pysam.FastaFile(self.reference_fasta)
            left_seq = "".join(self._ref_fasta.fetch(c1, s - 1, e) for s, e in left_coords)
            right_seq = "".join(self._ref_fasta.fetch(c2, s - 1, e) for s, e in right_coords)
            if not left_seq or not right_seq:
                return None, None, None
            return left_seq + right_seq, left_coords, right_coords
        except (ValueError, KeyError, AttributeError, OSError, TypeError):
            return None, None, None

    def realign_fusion_transcript(
        self,
        fusion_seq: str,
        min_match_len: int = 30,
        min_mapq: int = 10,
        breakpoint_windows: Optional[List[Tuple[str, int, int]]] = None,
    ) -> list:
        if not self.aligner or not fusion_seq or len(fusion_seq) < min_match_len:
            return []
        try:
            hits = list(self.aligner.map(fusion_seq))
            if not hits and breakpoint_windows and hasattr(self.aligner, "align_to_breakpoint_windows"):
                hits = self.aligner.align_to_breakpoint_windows(
                    fusion_seq, breakpoint_windows, min_match_len=min_match_len
                )
            good_hits = []
            for h in hits:
                q_st = getattr(h, "q_st", None)
                q_en = getattr(h, "q_en", None)
                if q_st is None or q_en is None:
                    q_st = getattr(h, "qstart", None) or getattr(h, "query_start", None)
                    q_en = getattr(h, "qend", None) or getattr(h, "query_end", None)
                match_len = None
                if q_st is not None and q_en is not None:
                    try:
                        match_len = int(q_en) - int(q_st)
                    except Exception:
                        match_len = None
                if not match_len:
                    match_len = getattr(h, "mlen", None) or getattr(h, "alen", None) or 0
                mapq = getattr(h, "mapq", 0)
                if match_len >= min_match_len and mapq >= min_mapq:
                    good_hits.append(h)
            return good_hits
        except Exception:
            return []

    def _is_mitochondrial_candidate(
        self, c1: str, c2: str, left_gene: Optional[str], right_gene: Optional[str]
    ) -> bool:
        if c1 in _MITO_CHRS or c2 in _MITO_CHRS:
            return True
        for g in (left_gene, right_gene):
            if not g:
                continue
            gu = str(g).upper()
            if gu.startswith("MT-") or gu in ("MT", "MTDNA"):
                return True
        return False

    def _apply_classification_and_filters(
        self,
        meta,
        flags,
        c1,
        p1,
        c2,
        p2,
        g1_name,
        r1,
        g2_name,
        r2,
        require_gene_names,
        max_intra_chr_distance,
    ) -> None:
        if g1_name and g2_name:
            g1 = meta["left_gene"]
            g2 = meta["right_gene"]
            if g1 == g2:
                flags["class"] = "intragenic"
                flags["is_valid"] = False
                flags["reasons"].append(f"Intralocus (antisense/divergent) event at {g1}")
            elif c1 == c2 and max_intra_chr_distance is not None and abs(p2 - p1) <= max_intra_chr_distance:
                flags["class"] = "cis-SAGe"
        elif r1 == "intergenic" or r2 == "intergenic":
            flags["class"] = "intergenic"
        if require_gene_names and not (g1_name and g2_name):
            flags["is_valid"] = False
            flags["reasons"].append("Missing gene name on one side")

    def _attempt_reconstruction_and_realignment(self, meta, flags, c1, p1, c2, p2) -> None:
        if not (self.reference_fasta and self.aligner):
            meta.setdefault("notes", []).append("Realignment skipped: reference or aligner not available")
            return
        fusion_seq, left_coords, right_coords = self.reconstruct_fusion_transcript(c1, p1, c2, p2, exon_padding=100)
        if not fusion_seq:
            flags["is_valid"] = False
            flags["reasons"].append("Failed to reconstruct fusion transcript")
            meta["reconstruction_ok"] = False
            return
        meta["reconstruction_ok"] = True
        windows = None
        if left_coords and right_coords:
            windows = [
                (c1, max(1, left_coords[0][0] - 50), left_coords[-1][1] + 50),
                (c2, max(1, right_coords[0][0] - 50), right_coords[-1][1] + 50),
            ]
        hits = self.realign_fusion_transcript(fusion_seq, breakpoint_windows=windows)
        if not hits:
            flags["is_valid"] = False
            flags["reasons"].append("Fusion transcript does not realign cleanly to genome")
        else:
            meta["realignment_hits"] = len(hits)
            try:
                meta["best_hit_mapq"] = max(h.mapq for h in hits)
            except Exception:
                meta["best_hit_mapq"] = 0

    def canonical_locus_name(self, gene_name: Optional[str]) -> Optional[str]:
        if not gene_name:
            return gene_name
        return ANTISENSE_SUFFIX_RE.sub("", gene_name)

    def _build_symbol_biotype_index(self) -> None:
        if self._symbol_biotype_cache:
            return
        try:
            genes = (
                self.interval_index.get_all_genes()
                if self.interval_index is not None
                else list(self.db.features_of_type("gene"))
            )
            for gene in genes:
                attrs = getattr(gene, "attributes", {}) or {}
                gene_name = attrs.get("gene_name", [None])[0]
                if not gene_name:
                    continue
                biotype = (
                    attrs.get("gene_type", [None])[0]
                    or attrs.get("gene_biotype", [None])[0]
                    or attrs.get("transcript_biotype", [None])[0]
                )
                if biotype:
                    self._symbol_biotype_cache[gene_name] = biotype
        except Exception as exc:
            logger.warning("Failed to build symbol->biotype index: %s", exc)

    def _get_gene_biotype_by_symbol_or_coords(
        self, gene_symbol: Optional[str], chrom: Optional[str] = None, pos: Optional[int] = None
    ) -> Optional[str]:
        if not gene_symbol:
            return None
        if not self._symbol_biotype_cache:
            self._build_symbol_biotype_index()
        if gene_symbol in self._symbol_biotype_cache:
            return self._symbol_biotype_cache[gene_symbol]
        if chrom is not None and pos is not None:
            try:
                genes = (
                    self.interval_index.get_genes_at(chrom, pos, window=2000)
                    if self.interval_index is not None
                    else list(self.db.region(region=(chrom, max(1, pos - 2000), pos + 2000), featuretype="gene"))
                )
                for g in genes:
                    attrs = getattr(g, "attributes", {}) or {}
                    g_name = attrs.get("gene_name", [None])[0]
                    if g_name and g_name.upper() == gene_symbol.upper():
                        biotype = (
                            attrs.get("gene_type", [None])[0]
                            or attrs.get("gene_biotype", [None])[0]
                            or attrs.get("transcript_biotype", [None])[0]
                        )
                        if biotype:
                            self._symbol_biotype_cache[gene_symbol] = biotype
                            return biotype
            except Exception:
                pass
        return None

    def get_gene_biotype(
        self, gene_name: Optional[str], chrom: Optional[str] = None, pos: Optional[int] = None
    ) -> Optional[str]:
        if not gene_name:
            return None
        normalized = self.normalize_gene_label(gene_name)
        return self._get_gene_biotype_by_symbol_or_coords(normalized, chrom=chrom, pos=pos)

    def _clear_read_level_data(self) -> None:
        self.fusion_assigned_pairs.clear()
        self.fusion_read_scores.clear()

    def _update_metadata_flags(self, meta: dict, flags: dict, confidence: Optional[float] = None) -> None:
        if confidence is not None:
            meta["confidence"] = confidence
        meta.update(flags)

    def validate_candidates(
        self,
        min_support: int = 2,
        window: int = 25,
        require_gene_names: bool = True,
        allow_cis_sage: bool = True,
        max_intra_chr_distance: Optional[int] = None,
        min_confidence_threshold: float = 0.30,
        low_confidence_support_cutoff: int = 2,
    ) -> None:
        self.build_metadata(min_support=min_support)
        self._clear_read_level_data()
        validator = FusionValidator(self)
        validator.filter_based_on_biotype()
        validator.filter_multicopy_artifact_pairs()
        validator.apply_frequency_filters()

        for meta in list(self.fusion_metadata.values()):
            if not meta.get("is_valid", True):
                continue
            support = meta.get("support", 0)
            flags = {"is_valid": True, "reasons": [], "class": "canonical"}
            if support < min_support:
                flags["is_valid"] = False
                flags["reasons"].append(f"Low support ({support} < {min_support})")
                self._update_metadata_flags(meta, flags)
                continue
            consensus_bp = meta.get("consensus_bp")
            if not consensus_bp:
                flags["is_valid"] = False
                flags["reasons"].append("No stable breakpoint cluster")
                self._update_metadata_flags(meta, flags)
                continue

            c1, p1, c2, p2 = consensus_bp
            g1_name, r1 = self._context_query(c1, p1)
            g2_name, r2 = self._context_query(c2, p2)
            self._apply_classification_and_filters(
                meta, flags, c1, p1, c2, p2, g1_name, r1, g2_name, r2, require_gene_names, max_intra_chr_distance
            )
            is_mito = self._is_mitochondrial_candidate(c1, c2, meta.get("left_gene"), meta.get("right_gene"))
            if flags["is_valid"] and not is_mito:
                self._attempt_reconstruction_and_realignment(meta, flags, c1, p1, c2, p2)
            elif is_mito and flags["is_valid"]:
                flags["is_valid"] = False
                flags["reasons"].append("Mitochondrial fusion candidate filtered out")

            if flags["class"] == "cis-SAGe" and not allow_cis_sage:
                flags["is_valid"] = False
                flags["reasons"].append("cis-SAGe disallowed by policy")

            conf = validator.confidence(meta, flags)
            if conf < min_confidence_threshold and support < low_confidence_support_cutoff and flags["is_valid"]:
                flags["is_valid"] = False
                flags["reasons"].append(f"Low confidence ({conf:.2f} < {min_confidence_threshold})")
            self._update_metadata_flags(meta, flags, confidence=conf)

    def report(
        self,
        output_path: str = "fusion_candidates.tsv",
        min_support: int = 2,
        include_classes: tuple = ("canonical", "cis-SAGe"),
        min_confidence: float = 0.3,
        only_valid: bool = False,
    ) -> None:
        self.validate_candidates(min_support=min_support)
        FusionValidator(self)._merge_fully_identical()
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(
                "LeftGene\tLeftBiotype\tLeftScore\tLeftChromosome\tLeftBreakpoint\t"
                "RightGene\tRightBiotype\tRightScore\tRightChromosome\tRightBreakpoint\t"
                "SupportingReads\tFusionName\tClass\tValid\tConfidence\tReasons\n"
            )
            for meta in sorted(self.fusion_metadata.values(), key=lambda x: -x.get("support", 0)):
                if meta.get("confidence", 0) < min_confidence:
                    continue
                if only_valid and not meta.get("is_valid", True):
                    continue
                if meta.get("class") not in include_classes:
                    continue
                consensus_bp = meta.get("consensus_bp")
                if not consensus_bp:
                    continue
                left_chr, left_pos, right_chr, right_pos = consensus_bp
                left_gene = meta["left_gene"]
                right_gene = meta["right_gene"]
                if left_gene == right_gene:
                    continue
                fusion_name = f"{left_gene}-{right_gene}"
                reasons = ";".join(meta.get("reasons", []))
                left_score = meta.get("left_score", 0.0)
                right_score = meta.get("right_score", 0.0)
                left_biotype = meta.get("left_biotype", "unknown")
                right_biotype = meta.get("right_biotype", "unknown")
                confidence = meta.get("confidence", 0.0) or 0.0
                f.write(
                    f"{left_gene}\t{left_biotype}\t{left_score:.2f}\t{left_chr}\t{left_pos}\t"
                    f"{right_gene}\t{right_biotype}\t{right_score:.2f}\t{right_chr}\t{right_pos}\t"
                    f"{meta.get('support', 0)}\t{fusion_name}\t{meta.get('class')}\t"
                    f"{meta.get('is_valid')}\t{confidence:.3f}\t{reasons}\n"
                )
