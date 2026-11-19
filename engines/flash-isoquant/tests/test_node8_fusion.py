"""Node 8 tests: stateless fusion detection, `.iqidx` sorted interval index, and Oracle parity."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import sys
from types import SimpleNamespace

import pysam
import pytest

from isoquant_next.annotation.iqidx import compile_gtf_to_iqidx
from isoquant_next.fusion import FusionDetector, FusionMetadata, FusionValidator, GenomicIntervalIndex

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


@dataclass
class _MockFeature:
    id: str
    chrom: str
    start: int
    end: int
    strand: str
    featuretype: str
    attributes: dict = field(default_factory=dict)


class _MockDB:
    def __init__(self, features, children_map=None) -> None:
        self._features = {f.id: f for f in features}
        self._list = list(features)
        self._children = children_map or {}

    def __getitem__(self, key):
        if key not in self._features:
            raise KeyError(key)
        return self._features[key]

    def features_of_type(self, ftype, id=None):
        for f in self._list:
            if f.featuretype == ftype and (id is None or f.id == id):
                yield f

    def children(self, parent, featuretype="exon", order_by="start"):
        pid = parent.id if hasattr(parent, "id") else parent
        out = [c for c in self._children.get(pid, []) if c.featuretype == featuretype]
        if order_by == "start":
            out.sort(key=lambda x: x.start)
        return out


def _write_test_gtf_and_fasta(tmp_path: Path):
    # Synthetic two-chromosome genome where GeneA (chr1:100-500) fuses to GeneB (chr2:100-500)
    rng_a = (
        "ACGTACGTAC" * 12
        + "ATGCGATCGTAGCTAGCTAGCTGACTGATCGATGCTAGCTA"
        + "TGCATGCATG" * 20
        + "CGTAGCTGATCGATGCTAGCTAGCATCGATCG"
        + "GGGGAAAACC" * 20
    )
    rng_b = (
        "TTTTCCCCAA" * 12
        + "CGATAGCTAGCTAGCATCGATGCTAGCTAGCTGATCGATC"
        + "ACACACGTAC" * 20
        + "GTAGCTGATCGTAGCTAGCTGATCGATCGATG"
        + "CCTTAAGGCC" * 20
    )
    seq1 = rng_a.ljust(1000, "A")[:1000]
    seq2 = rng_b.ljust(1000, "C")[:1000]
    fasta_path = tmp_path / "ref.fa"
    fasta_path.write_text(f">chr1\n{seq1}\n>chr2\n{seq2}\n", encoding="utf-8")
    pysam.faidx(str(fasta_path))

    gtf_path = tmp_path / "ann.gtf"
    gtf_path.write_text(
        "\n".join([
            'chr1\tsrc\tgene\t100\t450\t.\t+\t.\tgene_id "ENSG001"; gene_name "BCR"; gene_type "protein_coding";',
            'chr1\tsrc\ttranscript\t100\t450\t.\t+\t.\tgene_id "ENSG001"; transcript_id "ENST001"; gene_name "BCR"; gene_type "protein_coding";',
            'chr1\tsrc\texon\t100\t220\t.\t+\t.\tgene_id "ENSG001"; transcript_id "ENST001"; exon_id "EX1";',
            'chr1\tsrc\texon\t300\t450\t.\t+\t.\tgene_id "ENSG001"; transcript_id "ENST001"; exon_id "EX2";',
            'chr1\tsrc\tgene\t600\t900\t.\t-\t.\tgene_id "ENSG002"; gene_name "LINC001"; gene_type "lncRNA";',
            'chr1\tsrc\ttranscript\t600\t900\t.\t-\t.\tgene_id "ENSG002"; transcript_id "ENST002"; gene_name "LINC001"; gene_type "lncRNA";',
            'chr1\tsrc\texon\t600\t750\t.\t-\t.\tgene_id "ENSG002"; transcript_id "ENST002"; exon_id "EX3";',
            'chr2\tsrc\tgene\t100\t450\t.\t+\t.\tgene_id "ENSG003"; gene_name "ABL1"; gene_type "protein_coding";',
            'chr2\tsrc\ttranscript\t100\t450\t.\t+\t.\tgene_id "ENSG003"; transcript_id "ENST003"; gene_name "ABL1"; gene_type "protein_coding";',
            'chr2\tsrc\texon\t100\t220\t.\t+\t.\tgene_id "ENSG003"; transcript_id "ENST003"; exon_id "EX4";',
            'chr2\tsrc\texon\t300\t450\t.\t+\t.\tgene_id "ENSG003"; transcript_id "ENST003"; exon_id "EX5";',
        ]) + "\n",
        encoding="utf-8",
    )
    return fasta_path, gtf_path, seq1, seq2


def _write_synthetic_fusion_bam(bam_path: Path, seq1: str, seq2: str) -> None:
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}, {"SN": "chr2", "LN": 1000}]}
    left_seg = seq1[119:199]   # 80 bp on chr1:120-199
    right_seg = seq2[119:199]  # 80 bp on chr2:120-199
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as outf:
        for i in range(3):
            seg = pysam.AlignedSegment()
            seg.query_name = f"fusion_read_{i}"
            seg.query_sequence = left_seg + right_seg
            seg.flag = 0
            seg.reference_id = 0
            seg.reference_start = 119
            seg.mapping_quality = 60
            seg.cigarstring = "80M80S"
            seg.query_qualities = pysam.qualitystring_to_array("I" * 160)
            seg.set_tag("SA", "chr2,120,+,80M,60,0;")
            outf.write(seg)

        # Non-coding partner read (should be filtered by early non-coding check)
        nc_seg = pysam.AlignedSegment()
        nc_seg.query_name = "noncoding_read_0"
        nc_seg.query_sequence = seq1[600:680] + right_seg
        nc_seg.flag = 0
        nc_seg.reference_id = 0
        nc_seg.reference_start = 600
        nc_seg.mapping_quality = 60
        nc_seg.cigarstring = "80M80S"
        nc_seg.query_qualities = pysam.qualitystring_to_array("I" * 160)
        nc_seg.set_tag("SA", "chr2,120,+,80M,60,0;")
        outf.write(nc_seg)


def test_genomic_interval_index_parity_with_upstream_intervaltree(tmp_path: Path):
    from isoquant_lib.fusion.genomic_interval_index import GenomicIntervalIndex as OracleIndex

    g1 = _MockFeature("G1", "chr1", 100, 300, "+", "gene", {"gene_name": ["GENE_A"]})
    g2 = _MockFeature("G2", "chr1", 250, 500, "+", "gene", {"gene_name": ["GENE_B"]})
    e1 = _MockFeature("E1", "chr1", 100, 180, "+", "exon", {"gene_id": ["G1"]})
    e2 = _MockFeature("E2", "chr1", 260, 340, "+", "exon", {"gene_id": ["G2"]})
    db = _MockDB([g1, g2, e1, e2], children_map={"G1": [e1], "G2": [e2]})

    oracle_idx = OracleIndex(db)
    new_idx = GenomicIntervalIndex(db)

    for pos in (50, 100, 180, 250, 275, 300, 500, 600):
        for window in (None, 0, 20, 500):
            assert {x.id for x in new_idx.get_genes_at("chr1", pos, window)} == {
                x.id for x in oracle_idx.get_genes_at("chr1", pos, window)
            }
            assert {x.id for x in new_idx.get_exons_at("chr1", pos, window)} == {
                x.id for x in oracle_idx.get_exons_at("chr1", pos, window)
            }
    assert new_idx.get_genes_at("chrX", 100) == []
    assert [x.id for x in new_idx.find_genes_by_name("GENE_B")] == ["G2"]

    # Verify native .iqidx indexing matches
    _fa, gtf_path, _s1, _s2 = _write_test_gtf_and_fasta(tmp_path)
    iqidx_path = compile_gtf_to_iqidx(gtf_path, tmp_path / "ann.iqidx")
    from isoquant_next.annotation.iqidx import MmapAnnotationIndex

    with MmapAnnotationIndex(iqidx_path) as mmap_db:
        iq_idx = GenomicIntervalIndex(mmap_db)
        assert {g.attributes["gene_name"][0] for g in iq_idx.get_genes_at("chr1", 150)} == {"BCR"}
        assert {g.attributes["gene_name"][0] for g in iq_idx.get_genes_at("chr2", 200)} == {"ABL1"}
        assert iq_idx.get_genes_at("chr1", 500, window=0) == []
        exons_bcr = iq_idx.get_exons_of_gene(iq_idx.find_genes_by_name("BCR")[0])
        assert [(e.start, e.end) for e in exons_bcr] == [(100, 220), (300, 450)]


def test_fusion_detector_end_to_end_oracle_and_iqidx_parity(tmp_path: Path):
    import gffutils
    from isoquant_lib.fusion.fusion_detector import FusionDetector as OracleFusionDetector

    fasta_path, gtf_path, seq1, seq2 = _write_test_gtf_and_fasta(tmp_path)
    bam_path = tmp_path / "reads.bam"
    _write_synthetic_fusion_bam(bam_path, seq1, seq2)

    db_path = tmp_path / "ann.db"
    gffutils.create_db(str(gtf_path), dbfn=str(db_path), force=True, disable_infer_genes=True, disable_infer_transcripts=True)
    iqidx_path = compile_gtf_to_iqidx(gtf_path, tmp_path / "ann.iqidx")

    oracle_report = tmp_path / "oracle.fusions.tsv"
    new_report = tmp_path / "new.fusions.tsv"

    oracle_det = OracleFusionDetector(str(bam_path), str(db_path), reference_fasta=str(fasta_path))
    oracle_det.detect_fusions()
    oracle_det.report(output_path=str(oracle_report), min_support=2)

    with FusionDetector(str(bam_path), str(iqidx_path), reference_fasta=str(fasta_path)) as new_det:
        new_det.detect_fusions()
        new_det.report(output_path=str(new_report), min_support=2)
        assert "ABL1--BCR" in new_det.fusion_metadata
        meta = new_det.fusion_metadata["ABL1--BCR"]
        assert meta["is_valid"] is True
        assert meta["support"] == oracle_det.fusion_metadata["ABL1--BCR"]["support"]
        assert meta["reconstruction_ok"] is True

    oracle_text = oracle_report.read_text(encoding="utf-8")
    new_text = new_report.read_text(encoding="utf-8")
    assert new_text == oracle_text


def test_stateless_instance_caches_and_zero_forbidden_imports(tmp_path: Path):
    # Two independent detectors must not pollute one another's lookup caches
    g_a = _MockFeature("GA", "chr1", 100, 300, "+", "gene", {"gene_name": ["GENEA"], "gene_type": ["protein_coding"]})
    e_a = _MockFeature("EA", "chr1", 100, 300, "+", "exon", {"gene_id": ["GA"]})
    db_a = _MockDB([g_a, e_a], children_map={"GA": [e_a]})

    g_b = _MockFeature("GB", "chr1", 100, 300, "+", "gene", {"gene_name": ["GENEB"], "gene_type": ["protein_coding"]})
    e_b = _MockFeature("EB", "chr1", 100, 300, "+", "exon", {"gene_id": ["GB"]})
    db_b = _MockDB([g_b, e_b], children_map={"GB": [e_b]})

    det_a = FusionDetector("a.bam", db_a, reference_fasta=None)
    det_b = FusionDetector("b.bam", db_b, reference_fasta=None)

    assert det_a.assign_fusion_gene_cached("chr1", 150)[0] == "GENEA"
    assert det_b.assign_fusion_gene_cached("chr1", 150)[0] == "GENEB"
    assert det_a._context_query("chr1", 150) == ("GENEA", "exonic")
    assert det_b._context_query("chr1", 150) == ("GENEB", "exonic")

    # Verify targeted Cython SW fallback realignment when minimap2 is disabled
    fa, gtf_path, seq1, seq2 = _write_test_gtf_and_fasta(tmp_path)
    iqidx = compile_gtf_to_iqidx(gtf_path, tmp_path / "sw.iqidx")
    with FusionDetector("x.bam", str(iqidx), reference_fasta=str(fa)) as det_sw:
        det_sw.aligner.minimap2_bin = None  # force Cython sw_align_bytes fallback path
        fusion_seq, lc, rc = det_sw.reconstruct_fusion_transcript("chr1", 150, "chr2", 150, exon_padding=100)
        assert fusion_seq is not None
        hits = det_sw.realign_fusion_transcript(
            fusion_seq,
            min_match_len=30,
            min_mapq=10,
            breakpoint_windows=[("chr1", lc[0][0], lc[-1][1]), ("chr2", rc[0][0], rc[-1][1])],
        )
        assert len(hits) == 2
        assert {h.ctg for h in hits} == {"chr1", "chr2"}

    # Ensure zero forbidden dependencies in isoquant_next.fusion
    import isoquant_next.fusion.fusion_detector as fd_mod
    import isoquant_next.fusion.genomic_interval_index as gi_mod

    for forbidden in ("mappy", "intervaltree", "xgboost", "ssw", "editdistance"):
        assert forbidden not in fd_mod.__dict__
        assert forbidden not in gi_mod.__dict__
