"""ID distributors and feature ID registries compatible with IsoQuant 4.0.0."""

from __future__ import annotations

import threading
from typing import Dict, Optional, Set, Tuple

from isoquant_next.core.utils import TranscriptNaming


class SimpleIDDistributor:
    """Sequential positive integer ID generator."""

    def __init__(self) -> None:
        self.value = 0

    def increment(self) -> int:
        self.value += 1
        return self.value


class ExcludingIdDistributor(SimpleIDDistributor):
    """Generates novel integer suffixes while skipping existing novel IDs in annotation DB or index."""

    def __init__(self, genedb=None, chr_id: Optional[str] = None) -> None:
        super().__init__()
        self.forbidden_ids: Set[int] = set()
        if not genedb or not chr_id:
            return

        # Works cleanly with both gffutils FeatureDB and isoquant_next mmap index (MagicMock-safe via `is True`)
        if getattr(genedb, "supports_fast_id_lookup", False) is True:
            self.forbidden_ids.update(genedb.forbidden_numeric_ids(chr_id))
            return

        for g in genedb.region(seqid=chr_id, start=1, featuretype="gene"):
            if g.id.startswith(TranscriptNaming.novel_gene_prefix):
                try:
                    self.forbidden_ids.add(int(g.id.split("_")[-1]))
                except (IndexError, ValueError):
                    pass

        tx_prefix_len = len(TranscriptNaming.transcript_prefix)
        for t in genedb.region(seqid=chr_id, start=1, featuretype=("transcript", "mRNA")):
            if t.id.startswith(TranscriptNaming.transcript_prefix):
                try:
                    self.forbidden_ids.add(int(t.id.split(".")[0][tx_prefix_len:]))
                except (IndexError, ValueError):
                    pass

    def increment(self) -> int:
        self.value += 1
        while self.value in self.forbidden_ids:
            self.value += 1
        return self.value


class FeatureIdStorage:
    """Caches existing or newly assigned exon/intron IDs per genomic coordinate tuple."""

    def __init__(
        self,
        id_distributor: SimpleIDDistributor,
        genedb=None,
        chr_id: Optional[str] = None,
        feature: str = "exon",
    ) -> None:
        self.id_distributor = id_distributor
        self.id_dict: Dict[Tuple[str, int, int, str], str] = {}
        self.feature_name = feature
        if not genedb or not chr_id:
            return

        if getattr(genedb, "supports_fast_id_lookup", False) is True:
            self.id_dict.update(genedb.feature_id_map(chr_id, feature))
            return

        id_attr = f"{feature}_id"
        for f in genedb.region(seqid=chr_id, start=1, featuretype=feature):
            if id_attr in f.attributes:
                try:
                    self.id_dict[(chr_id, f.start, f.end, f.strand)] = f.attributes[id_attr][0]
                except IndexError:
                    pass

    def get_id(self, chr_id: str, feature: Tuple[int, int], strand: str) -> str:
        key = (chr_id, feature[0], feature[1], strand)
        existing = self.id_dict.get(key)
        if existing is not None:
            return existing
        feature_id = f"{chr_id}.{self.id_distributor.increment()}"
        self.id_dict[key] = feature_id
        return feature_id


class AtomicIDDistributor:
    """Thread-safe sequential integer ID distributor."""

    def __init__(self) -> None:
        self.value = 0
        self._lock = threading.Lock()

    def increment(self) -> int:
        with self._lock:
            self.value += 1
            return self.value
