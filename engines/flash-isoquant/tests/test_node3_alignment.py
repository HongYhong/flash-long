"""Tests for Node 3: _cy_cigar_blocks, AlignmentInfo Oracle parity, BAMStreamMerger, LocusChunker, IlluminaExonCorrector, ReadGroupers."""

from __future__ import annotations

import logging
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest

from isoquant_next.alignment import (
    AlignmentInfo,
    AlignmentTagReadGrouper,
    BAMStreamMerger,
    BarcodeSpotGrouper,
    IlluminaExonCorrector,
    LocusAlignmentSoABuffer,
    LocusChunker,
    MultiReadGrouper,
    ReadIdSplitReadGrouper,
    ReadTableGrouper,
    SharedTableData,
    create_read_grouper,
    get_grouping_pool_types,
    get_grouping_strategy_names,
    parse_grouping_spec,
    prepare_read_groups,
)
from isoquant_next.annotation.iqidx import MmapAnnotationIndex, compile_gtf_to_iqidx

UPSTREAM_ROOT = Path(__file__).resolve().parents[2]
SIMPLE_DATA_DIR = UPSTREAM_ROOT / "isoquant_tests" / "simple_data"
SHORT_TOY_BAM = UPSTREAM_ROOT / "isoquant_tests" / "short_reads_toy" / "short_toy.bam"
ENCFF400BQQ_BAM = Path(
    "/gpfs/group/data_sci_group/weiwu_group/hongyanhong/backup/longpass_benchmark250209/"
    "cDNA_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF400BQQ.sorted.k14.bam"
)
V39_IQIDX = Path(__file__).resolve().parents[1] / "tmp_bench" / "g0_probe" / "v39.iqidx"


class _SyntheticAlignment:
    def __init__(self, qname, cigartuples, seq, ref_start, ref_end, ref_name, ref_id, aligned_pairs):
        self.query_name = qname
        self.cigartuples = cigartuples
        self.seq = seq
        self.query_sequence = seq
        self.reference_start = ref_start
        self.reference_end = ref_end
        self.reference_name = ref_name
        self.reference_id = ref_id
        self.aligned_pairs = aligned_pairs
        self.is_reverse = False

    def get_aligned_pairs(self):
        return self.aligned_pairs


def test_alignment_info_oracle_synthetic_and_bam():
    import sys

    if str(UPSTREAM_ROOT) not in sys.path:
        sys.path.insert(0, str(UPSTREAM_ROOT))
    from isoquant_lib.alignment.alignment_info import AlignmentInfo as OrigAlignmentInfo

    seq1 = (
        "CTCAAGACCAAGAAGGACGACATGACCATGGCTTAAAAGAGTCTGCTCCCCACAGCCCCCTGCGAT"
        "GGATGGACGGAGGAACCAGGGTCGGACGACCTCCGATGCTAAGAGCACTCCAACTGCTGCAAACCG"
        "AAGAAGCAGGCATCGGACGCCTTTTTTTTTTAGGCGCCAAAAAAAAAAAACCAAAAAAAAAAAAA"
    )
    tuples1 = [(0, 167), (4, 30)]
    tuples2 = [(0, 167), (4, 30), (5, 5), (4, 16), (0, 30), (5, 70), (0, 201)]
    pairs = [(0, 75897), (1, 76063), (2, 75995)]

    a2 = _SyntheticAlignment("query", tuples2, seq1, 75896, 113577, "reference", 789, pairs)
    nx2 = AlignmentInfo(a2)
    or2 = OrigAlignmentInfo(a2)
    assert nx2.read_exons == or2.read_exons == [(75897, 76063), (76064, 76294)]
    assert nx2.read_blocks == or2.read_blocks == [(0, 166), (213, 443)]
    assert nx2.cigar_blocks == or2.cigar_blocks == [(0, 0), (4, 6)]
    assert nx2.get_error_count(75999, 81000, None, None, None) == (0, 0)
    assert nx2.get_error_count(75999, 81000, 0, False, None) == (0, 0)
    assert nx2.get_error_count(75999, 81000, 0, True, None) == (0, 0)
    with pytest.raises(AssertionError):
        nx2.get_seq(76000, 81000)

    a1 = _SyntheticAlignment("query", tuples1, seq1, 75896, 113577, "reference", 789, pairs)
    nx1 = AlignmentInfo(a1)
    nx1.set_aligned_pairs()
    assert len(pairs) in nx1.aligned_pairs_start_index
    assert len(nx1.aligned_pairs_start_index) == 3

    # Exact Oracle parity on every read in chr9.4M.ont.sim.polya.bam (plus mismatch checks against chr9 fasta)
    chr9_bam = SIMPLE_DATA_DIR / "chr9.4M.ont.sim.polya.bam"
    chr9_fa = SIMPLE_DATA_DIR / "chr9.4M.fa.gz"
    import gzip

    chr9_seq_parts = []
    with gzip.open(chr9_fa, "rt") as fh:
        in_chr9 = False
        for line in fh:
            if line.startswith(">"):
                in_chr9 = line[1:].strip().split()[0] == "chr9"
            elif in_chr9:
                chr9_seq_parts.append(line.strip())
    chr9_seq = "".join(chr9_seq_parts)

    checked = 0
    with pysam.AlignmentFile(str(chr9_bam), "rb") as af:
        for aln in af:
            if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                continue
            nx_a = AlignmentInfo(aln)
            or_a = OrigAlignmentInfo(aln)
            assert nx_a.read_exons == or_a.read_exons
            assert nx_a.read_blocks == or_a.read_blocks
            assert nx_a.cigar_blocks == or_a.cigar_blocks
            assert nx_a.read_start == or_a.read_start
            assert nx_a.read_end == or_a.read_end
            if len(nx_a.read_exons) > 1 and checked < 40:
                s0 = nx_a.read_exons[0][1] - 10
                e0 = nx_a.read_exons[1][0] + 10
                assert nx_a.get_error_count(s0, e0, intron_index=0, left_site=True, chr_record=chr9_seq) == or_a.get_error_count(
                    s0, e0, intron_index=0, left_site=True, chr_record=chr9_seq
                )
            checked += 1
    assert checked > 50


def test_bam_merger_locus_chunker_and_soa_buffer(tmp_path):
    chr9_bam = str(SIMPLE_DATA_DIR / "chr9.4M.ont.sim.polya.bam")
    merger = BAMStreamMerger([chr9_bam, chr9_bam], chr_id="chr9")
    coords = [aln.reference_start for _fidx, aln in merger.iter_alignments()]
    assert len(coords) > 100
    assert coords == sorted(coords)

    iqidx_p = compile_gtf_to_iqidx(SIMPLE_DATA_DIR / "chr9.4M.gtf.gz", tmp_path / "chr9.iqidx")
    with MmapAnnotationIndex(iqidx_p) as idx:
        chunker = LocusChunker("chr9", max_reads_per_locus=500)
        windows = list(chunker.iter_windows(BAMStreamMerger([chr9_bam], chr_id="chr9").iter_alignments(), iqidx=idx))

    assert len(windows) >= 5
    # Verify non-overlapping window bounds and zero read loss
    prev_end = -1
    tot_alns = 0
    for w in windows:
        assert w.start <= w.end
        assert w.start > prev_end
        prev_end = w.end
        tot_alns += len(w.alignments)
    assert tot_alns == len(coords) // 2

    # SoA compact buffer over largest chr9 locus window
    non_empty = [w for w in windows if w.alignments]
    ai_list = [AlignmentInfo(aln) for _fi, aln in non_empty[0].alignments]
    soa = LocusAlignmentSoABuffer(ai_list)
    assert len(soa) == len(ai_list)
    assert soa.exon_ptr[-1] == sum(len(a.read_exons) for a in ai_list)

    # High-density safety cap splits even overlapping reads deterministically
    capped = list(
        LocusChunker("chr9", max_reads_per_locus=25).iter_windows(
            BAMStreamMerger([chr9_bam], chr_id="chr9").iter_alignments()
        )
    )
    assert all(len(w.alignments) <= 25 for w in capped)
    assert any(w.split_for_density for w in capped)

    # Tier 1 human BAM + .iqidx LocusChunker check on chr22 if benchmark files are mounted
    if ENCFF400BQQ_BAM.exists() and V39_IQIDX.exists():
        with MmapAnnotationIndex(V39_IQIDX) as hidx:
            c22_wins = list(
                LocusChunker("chr22", max_reads_per_locus=2000).iter_windows(
                    BAMStreamMerger([str(ENCFF400BQQ_BAM)], chr_id="chr22").iter_alignments(),
                    iqidx=hidx,
                )
            )
        assert len(c22_wins) > 100
        assert sum(len(w.alignments) for w in c22_wins) > 1000


def test_illumina_exon_corrector_oracle_and_bugfixes(caplog):
    import sys

    if str(UPSTREAM_ROOT) not in sys.path:
        sys.path.insert(0, str(UPSTREAM_ROOT))
    from isoquant_lib.assignment.illumina_exon_corrector import IlluminaExonCorrector as OrigIEC

    cases_offset = [
        (
            {
                (11689307, 11851227), (11776574, 11820421), (11788720, 11851119), (11797197, 11850752),
                (11818316, 11868704), (11818325, 11868704), (11868791, 11900038), (11877909, 12147885),
                (11900116, 12003767), (11900116, 12021186), (11926162, 12147885), (12003816, 12021186),
                (12021281, 12042269), (12042327, 12045146),
            },
            [(11818191, 11818324), (11868705, 11868790), (11900039, 11900115), (12003768, 12003815),
             (12021187, 12021284), (12042270, 12042326), (12045147, 12045483)],
            [(11818191, 11818324), (11868705, 11868790), (11900039, 11900115), (12003768, 12003815),
             (12021187, 12021280), (12042270, 12042326), (12045147, 12045483)],
        ),
        (
            {
                (4797249, 4898804), (4878206, 4878677), (4878710, 4898806), (4881433, 4900490),
                (4898873, 4900490), (4900539, 4902533), (4902592, 4911329), (4902605, 4907223),
                (4907298, 4909609), (4909712, 4911178), (4909712, 5079593), (4911356, 4915185),
                (4915401, 5028533),
            },
            [(4878046, 4878209), (4898782, 4898872), (4900491, 4900538), (4902534, 4902604),
             (4907224, 4907297), (4909610, 4909711), (4911179, 4911355), (4915186, 4916833)],
            [(4878046, 4878205), (4878678, 4878709), (4898807, 4898872), (4900491, 4900538),
             (4902534, 4902604), (4907224, 4907297), (4909610, 4909711), (4911179, 4911355), (4915186, 4916833)],
        ),
    ]
    for introns, long_exons, expected in cases_offset:
        nx = IlluminaExonCorrector.from_data(introns)
        or_c = OrigIEC.from_data(introns)
        assert nx.correct_exons(long_exons) == expected
        assert nx.correct_exons(long_exons) == or_c.correct_exons(long_exons)

    # Real short_toy.bam single & double file checks
    nx_bam = IlluminaExonCorrector("chr1", 5154787, 5232327, [str(SHORT_TOY_BAM)])
    or_bam = OrigIEC("chr1", 5154787, 5232327, [str(SHORT_TOY_BAM)])
    assert nx_bam.short_introns == or_bam.short_introns
    assert nx_bam.counts == or_bam.counts
    nx_bam2 = IlluminaExonCorrector("chr1", 4561154, 4565934, [str(SHORT_TOY_BAM), str(SHORT_TOY_BAM)])
    or_bam2 = OrigIEC("chr1", 4561154, 4565934, [str(SHORT_TOY_BAM), str(SHORT_TOY_BAM)])
    assert nx_bam2.short_introns == or_bam2.short_introns
    assert nx_bam2.counts == or_bam2.counts

    # Fixed isinstance(int) + fixed logger.debug formatting (no TypeError when logging enabled)
    with caplog.at_level(logging.DEBUG, logger="IsoQuant"):
        nx_rev = IlluminaExonCorrector.from_data({(100, 120), (80, 95)})
        nx_rev.correct_exons([(50, 79), (121, 150)])


def test_read_groups_and_simple_tsv(tmp_path):
    class _MockAln:
        def __init__(self, rid="", tags=None):
            self.query_name = rid
            self._tags = tags or {}

        def has_tag(self, t):
            return t in self._tags

        def get_tag(self, t):
            return self._tags.get(t)

    ag = AlignmentTagReadGrouper("CB")
    assert ag.get_group_id(_MockAln(tags={"CB": "ACGT"})) == "ACGT"
    assert ag.get_group_id(_MockAln()) is None

    rg = ReadIdSplitReadGrouper("_")
    assert rg.get_group_id(_MockAln("read_001_groupA")) == "groupA"
    assert rg.get_group_id(_MockAln("read001")) == "NA"
    assert rg.get_group_id(_MockAln("a_b_c")) == "c"

    tsv = tmp_path / "rg.tsv"
    tsv.write_text("read_001\tcellA\tspot1\nread_002\tcellB\tspot2\n")
    sd = SharedTableData(str(tsv), 0, [1, 2], "\t")
    g0 = ReadTableGrouper(sd, 0)
    g1 = ReadTableGrouper(sd, 1)
    assert g0.get_group_id(_MockAln("read_001")) == "cellA"
    assert g1.get_group_id(_MockAln("read_002")) == "spot2"
    assert g0.get_group_id(_MockAln("missing")) == "NA"

    # Real chr9.4M.ont.sim.read_groups.tsv preparation
    rg_real = SIMPLE_DATA_DIR / "chr9.4M.ont.sim.read_groups.tsv"
    bam_real = SIMPLE_DATA_DIR / "chr9.4M.ont.sim.polya.bam"
    sample = SimpleNamespace(
        file_list=[(str(bam_real), "s0")],
        readable_names_dict={str(bam_real): "s0"},
        read_group_file=str(tmp_path / "sample_rg"),
        barcodes_split_reads=False,
    )
    args = SimpleNamespace(read_group=[f"file:{rg_real}:0:1", "tag:RG", "read_id:_"], threads=2)
    prepare_read_groups(args, sample)
    grp = create_read_grouper(args, sample, "chr9")
    assert isinstance(grp, MultiReadGrouper)
    assert get_grouping_strategy_names(args) == ["file0_col1", "tag_RG", "read_id__"]
    assert get_grouping_pool_types(args) == {0: "tsv:0:1:\t", 1: "dynamic", 2: "dynamic"}
    # Spot check first read from BAM
    with pysam.AlignmentFile(str(bam_real), "rb") as af:
        first_aln = next(iter(af))
    vals = grp.get_group_id(first_aln)
    assert vals[0] == "GR0"
