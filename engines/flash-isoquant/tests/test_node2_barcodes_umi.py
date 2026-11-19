"""Tests for Node 2: Cython SW/Myers kernels, 2-bit SharedMemory, callers, Universal MDF, and UMI filtering."""

from __future__ import annotations

import os
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from isoquant_next.barcodes import (
    Array2BitKmerIndexer,
    ArrayKmerIndexer,
    BarcodeCaller,
    BarcodeDetectionResult,
    BarcodeResult,
    CurioBarcodeDetector,
    CurioIlluminaDetector,
    DetectedElement,
    Dict2BitKmerIndexer,
    ElementType,
    ExtractionReadStats,
    ExtractionResult,
    KmerIndexer,
    LinkerBarcodeDetectionResult,
    MoleculeElement,
    MoleculeStructure,
    ReadStats,
    SharedMemoryArray2BitKmerIndexer,
    SharedMemoryIndexInfo,
    SharedMemoryStereoBarcodeDetector,
    SharedMemoryStereoSplittingBarcodeDetector,
    SplittingBarcodeDetectionResult,
    StereoBarcodeDetector,
    StereoSplittingBarcodeDetector,
    TSOBarcodeDetectionResult,
    TenXBarcodeDetectionResult,
    TenXBarcodeDetector,
    TenXSplitBarcodeDetectionResult,
    TenXSplittingBarcodeDetector,
    TenXv2BarcodeDetector,
    TenXv2SplittingBarcodeDetector,
    UMIFilter,
    UNTRUSTED_UMIS,
    UNTRUSTED_UMI_KEY,
    UniversalSingleMoleculeExtractor,
    VisiumHDBarcodeDetector,
    batch_str_to_2bit,
    batch_str_to_2bit_chunked,
    bit_to_str,
    detect_exact_positions,
    edit_distance_fast,
    find_polyt,
    find_polyt_start,
    is_untrusted_umi,
    iter_fastx_chunks,
    process_in_parallel,
    process_single_thread,
    reverese_complement,
    str_to_2bit,
)
from isoquant_next.barcodes._cy_barcode_umi import sw_align_bytes
from isoquant_next.barcodes.common import align_pattern_ssw
from isoquant_next.barcodes.table_splitter import get_chromosome_read_counts, split_read_table_by_chromosomes
from isoquant_next.barcodes.umi_filtering import create_transcript_info_dict
from isoquant_next.config.enums import IsoQuantMode

UPSTREAM_ROOT = Path(__file__).resolve().parents[2]
UNIVERSAL_DATA_DIR = UPSTREAM_ROOT / "isoquant_tests" / "universal_data"
SIMPLE_DATA_DIR = UPSTREAM_ROOT / "isoquant_tests" / "simple_data"


def _py_gotoh_sw(query: str, ref: str, match: int = 1, mismatch: int = 1, gap_open: int = 1, gap_ext: int = 1):
    """Exact Gotoh SW reference (open = gap_open + gap_ext; linear extension)."""
    qlen, rlen = len(query), len(ref)
    if qlen == 0 or rlen == 0:
        return 0, 0, 0, 0, 0
    neg = -10**9
    H = [[0] * (rlen + 1) for _ in range(qlen + 1)]
    E = [[neg] * (rlen + 1) for _ in range(qlen + 1)]
    F = [[neg] * (rlen + 1) for _ in range(qlen + 1)]
    best = 0
    best_i = best_j = 0
    for i in range(1, qlen + 1):
        qi = query[i - 1]
        for j in range(1, rlen + 1):
            E[i][j] = max(H[i][j - 1] - gap_open, E[i][j - 1] - gap_ext)
            F[i][j] = max(H[i - 1][j] - gap_open, F[i - 1][j] - gap_ext)
            sub = match if qi == ref[j - 1] else -mismatch
            h = max(0, H[i - 1][j - 1] + sub, E[i][j], F[i][j])
            H[i][j] = h
            if h > best:
                best = h
                best_i, best_j = i, j
    return best, best_i - 1, best_j - 1


def _py_levenshtein(a: str, b: str) -> int:
    if a == b:
        return 0
    prev = list(range(len(b) + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cur.append(min(cur[-1] + 1, prev[j] + 1, prev[j - 1] + (0 if ca == cb else 1)))
        prev = cur
    return prev[-1]


def test_cy_sw_and_myers_exact_vs_reference():
    pairs = [
        ("CTACACGACGCTCTTCCGATCT", "NNNNCTACACGACGGCTCTTCCGATCTNNNN"),
        ("TCCCCACTTTCACTCTCACCCCATT", "TTCCGATCTATGGCGACCTTATCAGTCCCCACTTTCACTCTCACCCCATTTTGTCTTCCTAAGAC"),
        ("ACGTACGTACGT", "ACGTACGTACGT"),
        ("AAAAAAAAAAAA", "TTTTTTTTTTTT"),
    ]
    for q, r in pairs:
        qb = np.frombuffer(q.encode("ascii"), dtype=np.uint8)
        rb = np.frombuffer(r.encode("ascii"), dtype=np.uint8)
        sc, ref_st, ref_en, q_st, q_en = sw_align_bytes(qb, rb, 1, 1, 1, 1)
        ref_sc, _, _ = _py_gotoh_sw(q, r, 1, 1, 1, 1)
        assert sc == ref_sc
        if sc > 0:
            assert 0 <= ref_st <= ref_en < len(r)
            assert 0 <= q_st <= q_en < len(q)

        # Also test exact Myers bit-parallel Levenshtein (including strings > 64 chars for DP fallback)
        assert edit_distance_fast(q, r) == _py_levenshtein(q, r)

    long_a = "ACGT" * 20
    long_b = "ACGT" * 10 + "ACTT" + "ACGT" * 9
    assert edit_distance_fast(long_a, long_b) == 1


def test_polyt_and_2bit_nep50():
    assert reverese_complement("A C G T N") == "N A C G T"
    seq = "GGGGG" + ("T" * 20) + "CCCCCCCCCC"
    st, en = find_polyt(seq)
    assert st == 5
    assert find_polyt_start(seq) == 5

    bcs = ["ACTGACTGACTGACTGACTGACTGA", "TGCATGCATGCATGCATGCATGCAT", "GGGGGGGGGGGGGGGGGGGGGGGGG"]
    enc = batch_str_to_2bit(bcs, 25)
    enc_chunked = batch_str_to_2bit_chunked(iter(bcs), 25, chunk_size=2)
    assert enc.dtype == np.uint64
    assert np.array_equal(enc, enc_chunked)
    for i, bc in enumerate(bcs):
        assert int(enc[i]) == str_to_2bit(bc)
        assert bit_to_str(int(enc[i]), 25) == bc


def test_kmer_indexers_and_shared_memory_lifecycle():
    barcodes_8 = ["ACTGACTG", "TGCATGCA", "GGGGGGGG"]
    k_idx = KmerIndexer(barcodes_8, kmer_size=4)
    a_idx = ArrayKmerIndexer(barcodes_8, kmer_size=4)
    bit_8 = batch_str_to_2bit(barcodes_8, seq_len=8)
    d2_idx = Dict2BitKmerIndexer(bit_8, kmer_size=4, seq_len=8)
    assert {r[0] for r in k_idx.get_occurrences("ACTGACTT")} == {"ACTGACTG"}
    assert {r[0] for r in a_idx.get_occurrences("ACTGACTT")} == {"ACTGACTG"}
    assert {r[0] for r in d2_idx.get_occurrences("ACTGACTT")} == {"ACTGACTG"}

    barcodes_25 = [
        "ACTGACTGACTGACTGACTGACTGA",
        "TGCATGCATGCATGCATGCATGCAT",
        "GGGGGGGGGGGGGGGGGGGGGGGGG",
    ]
    bits_25 = [str_to_2bit(b) for b in barcodes_25]
    arr2 = Array2BitKmerIndexer(bits_25, kmer_size=10, seq_len=25)
    assert arr2.total_sequences == 3
    assert len(list(arr2._get_kmer_bin_indexes(bits_25[0]))) == 16

    with SharedMemoryArray2BitKmerIndexer(bits_25, kmer_size=10, seq_len=25) as shm_idx:
        assert shm_idx.total_sequences == 3
        assert shm_idx.index_size == 48
        info = shm_idx.get_sharable_info()
        assert isinstance(info, SharedMemoryIndexInfo)
        state = info.__getstate__()
        assert len(state) == 7
        recon = SharedMemoryArray2BitKmerIndexer.from_sharable_info(info)
        try:
            q = "ACTGACTGACTGACTGACTGACTGT"
            assert {r[0] for r in recon.get_occurrences(q)} == {"ACTGACTGACTGACTGACTGACTGA"}
        finally:
            recon.close()
        shm_names = [info.barcodes_sm_name, info.index_sm_name, info.index_range_sm_name]

    # Verify owner context manager unlinked POSIX shared memory cleanly
    from multiprocessing import shared_memory as _sm

    for nm in shm_names:
        with pytest.raises(FileNotFoundError):
            _sm.SharedMemory(name=nm, create=False)


def test_tenx_visium_curio_stereo_detectors():
    # Stereo-seq real read patterns from test_barcode_detectors.py
    stereo_bcs = [
        "GAAGGAGGGATCCATAACGGCCCCA",
        "TCCCCACTTTCACTCTCACCCCATT",
        "AGTTGGCGCCCAAATGGTCAGAGGG",
        "CATCAAAGGTTCACCGAATCTCCAC",
        "GTCGCCTTCAACGTATACTGATATA",
        "ACCAATTAGCAGCCATGGGAACTCC",
        "ACATGATCCTCCGTGGCATGCGCAC",
        "TATGTACGAATTTAACCGTACCCTC",
        "GGGCACATTCCTCAGAAAGCGAAAA",
        "GGCAGGCCATTAGTTACTCGTCGAG",
        "CTCCCCTAATCGTTTGTCGGTGGAC",
    ]
    stereo = StereoBarcodeDetector(stereo_bcs)
    res = stereo.find_barcode_umi(
        "r_st1",
        "TTCCGATCTATGGCGACCTTATCAGTCCCCACTTTCACTCTCACCCCATTTTGTCTTCCTAAGACCTAGTGTAGCTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
        "CACATATAAATATGACTTTATTAACCAGTCTGTGTTTGATTTAGCACTCTAGTGGCTATCAGTGATTTTCTTTTA",
    )
    assert res.is_valid()
    assert res.barcode == "TCCCCACTTTCACTCTCACCCCATT"
    assert res.UMI == "CTAGTGTAGC"

    shared_stereo = SharedMemoryStereoBarcodeDetector(stereo_bcs)
    res_s = shared_stereo.find_barcode_umi(
        "r_st2",
        "CTACACTAGGTCTTAGGAAGACAAAATGGGGTGAGAGTGAAAGTGGGGACTGATAAGG",
    )
    assert res_s.is_valid() and res_s.barcode == "TCCCCACTTTCACTCTCACCCCATT"

    # Stereo splitting using real concatenated molecule read from test_barcode_detectors.py:184-186
    split_stereo = StereoSplittingBarcodeDetector(stereo_bcs)
    split_seq = (
        "TTCCGATCTATGGCGACCTTATCAGTATGTACGAATTTAACCGTACCCTCTTGTCTTCCTAAGACCTTGCCGTCG"
        + ("T" * 32)
        + ("ACGT" * 20)
        + ("A" * 25)
        + "CGGGCTTTCAGTCATAGGAAGACAAGTGCGCATGCCACGGAGGATCATGTCTGATAAGGTCGCCATAGATCGGAAG"
    )
    res_sp = split_stereo.find_barcode_umi("r_sp", split_seq)
    assert isinstance(res_sp, SplittingBarcodeDetectionResult)
    assert res_sp.is_valid()
    got_bcs = {p.barcode for p in res_sp.detected_patterns if p.is_valid()}
    assert "TATGTACGAATTTAACCGTACCCTC" in got_bcs
    assert "ACATGATCCTCCGTGGCATGCGCAC" in got_bcs

    # 10x v3 single & splitting + v2 splitting
    tenx_bcs = ["TCCATTGTCAACAAGG", "CAATCTAAGTACTGGT"]
    tx3 = TenXBarcodeDetector(tenx_bcs)
    r10 = tx3.find_barcode_umi(
        "r10",
        "TACACGACGCTCTTCCGATCTTCCATTGTCAACAAGGCTCATGGCCCCATTTTTTTTTTTTTTGTTTTTTTTTTTTTTTTGCTTTTTG",
    )
    assert isinstance(r10, TenXBarcodeDetectionResult)
    assert r10.barcode == "TCCATTGTCAACAAGG"
    assert r10.UMI == "CTCATGGCCCCA"

    concat_10x = (
        "CTACACGACGCTCTTCCGATCTTCCATTGTCAACAAGGCTCATGGCCCCA" + ("T" * 20) + ("GCTA" * 15)
        + "CCCATGTACTCTGCGTTGATACCACTGCTT"
        + "CTACACGACGCTCTTCCGATCTCAATCTAAGTACTGGTATCCAGGCTAGA" + ("T" * 20) + ("GCTA" * 15)
    )
    tx3_sp = TenXSplittingBarcodeDetector(tenx_bcs)
    r10_sp = tx3_sp.find_barcode_umi("r10sp", concat_10x)
    assert len(r10_sp.detected_patterns) == 2
    assert {p.get_barcode() for p in r10_sp.detected_patterns} == {"TCCATTGTCAACAAGG", "CAATCTAAGTACTGGT"}
    assert r10_sp.detected_patterns[0].tso != -1

    # Curio single & illumina
    curio = CurioBarcodeDetector((["ACGTACGT", "TGCATGCA"], ["CCCCCC", "GGGGGG"]))
    assert curio.result_type() is LinkerBarcodeDetectionResult
    cur_il = CurioIlluminaDetector(["ACGTACGTCCCCCC", "TGCATGCAGGGGGG"])
    assert isinstance(cur_il.find_barcode_umi("x", "NNNN" * 10), BarcodeDetectionResult)

    # Visium HD
    vhd = VisiumHDBarcodeDetector([["AAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCC"], ["TTTTTTTTTTTTTTT", "GGGGGGGGGGGGGGG"]])
    assert not vhd.find_barcode_umi("v", "ACGT" * 20).is_valid()


def test_universal_extraction_and_parallel_caller(tmp_path):
    # Concatenated + Duplicated + Real MDF from universal_data/
    mdf_concat = """Barcode|1:Barcode|2:UMI:PolyT:cDNA
Barcode|1\tVAR_LIST\tAAAACCCC,GGGGTTTT,ACGTACGT\t4
Barcode|2\tVAR_LIST\tAAAACCCC,GGGGTTTT,ACGTACGT\t4
UMI\tVAR_ANY\t8
"""
    ms_c = MoleculeStructure(iter(mdf_concat.splitlines()))
    ext_c = UniversalSingleMoleculeExtractor(ms_c)
    assert "barcode" in ext_c.header() and "Barcode|1" not in ext_c.header()
    det_res = {}
    ext_c.process_concatenated_elements(
        {"Barcode": [(10, 13), (14, 17)]}, det_res, "XXXXXXXXXXAAAACCCCXXXXXXXX"
    )
    assert det_res["Barcode"].seq == "AAAACCCC"

    mdf_dupl = """Barcode/1:Barcode/2:UMI:PolyT:cDNA
Barcode/1\tVAR_LIST\tAAAACCCC,GGGGTTTT,ACGTACGT
Barcode/2\tVAR_LIST\tAAAACCCC,GGGGTTTT,ACGTACGT
UMI\tVAR_ANY\t8
"""
    ext_d = UniversalSingleMoleculeExtractor(MoleculeStructure(iter(mdf_dupl.splitlines())))
    det_d1 = {}
    ext_d.process_duplicated_elements(
        {"Barcode": [(10, 17), (26, 33)]}, det_d1, "XXXXXXXXXXAAAACCCCXXXXXXXXGGGGTTTTXXXXXXXX"
    )
    assert det_d1["Barcode"].seq == ExtractionResult.NOSEQ  # Disagreeing copies resolve to NOSEQ
    det_d2 = {}
    ext_d.process_duplicated_elements(
        {"Barcode": [(10, 17), (-1, -1)]}, det_d2, "XXXXXXXXXXAAAACCCCXXXXXXXXXXXXXXXX"
    )
    assert det_d2["Barcode"].seq == "AAAACCCC"

    # Real 10x MDF with BARCODES file resolved relative to universal_data
    mdf_lines = (UNIVERSAL_DATA_DIR / "10x.mdf").read_text().splitlines()
    patched = [
        line.replace("barcodes_small.tsv", str(UNIVERSAL_DATA_DIR / "barcodes_small.tsv"))
        for line in mdf_lines
    ]
    ext_10x = UniversalSingleMoleculeExtractor(MoleculeStructure(iter(patched)))
    fq_path = UNIVERSAL_DATA_DIR / "test_reads.fq"
    chunks = list(iter_fastx_chunks(fq_path, chunk_size=2))
    assert sum(len(c) for c in chunks) > 0

    out_single = tmp_path / "single.tsv"
    out_para = tmp_path / "para.tsv"
    args_s = SimpleNamespace(
        mode=IsoQuantMode.custom_sc,
        molecule_description=str(tmp_path / "10x_abs.mdf"),
        input=[str(fq_path)],
        output=str(out_single),
        out_fasta=str(tmp_path / "single.fasta"),
        threads=1,
    )
    Path(args_s.molecule_description).write_text("\n".join(patched) + "\n")
    process_single_thread(args_s)

    args_p = SimpleNamespace(
        mode=IsoQuantMode.custom_sc,
        molecule_description=args_s.molecule_description,
        input=[str(fq_path)],
        output=str(out_para),
        out_fasta=str(tmp_path / "para.fasta"),
        threads=2,
    )
    process_in_parallel(args_p)

    s_lines = sorted(out_single.read_text().splitlines()[1:])
    p_lines = sorted(out_para.read_text().splitlines()[1:])
    assert s_lines == p_lines
    assert (tmp_path / "single.tsv.stats").is_file()
    assert "Total reads" in (tmp_path / "para.tsv.stats").read_text()

    # ReadStats vs ExtractionReadStats
    er = ExtractionResult(MoleculeStructure(iter(patched)), "r1", "+", {"PolyT": DetectedElement(10, 30, 0)})
    estats = ExtractionReadStats()
    estats.add_read(er)
    assert "Total reads:\t1" in str(estats)


def test_umi_filtering_and_table_splitter(tmp_path):
    assert is_untrusted_umi(None) and is_untrusted_umi("*") and not is_untrusted_umi("ACGTACGT")

    class _DummyType:
        def __init__(self, unique: bool):
            self._unique = unique

        def is_unique(self) -> bool:
            return self._unique

        def is_inconsistent(self) -> bool:
            return False

    class _DummyMatch:
        def __init__(self, tx, gene):
            self.assigned_transcript = tx
            self.assigned_gene = gene

    class _DummyRead:
        def __init__(self, read_id, umi, exons, unique=True):
            self.read_id = read_id
            self.umi = umi
            self.corrected_exons = exons
            self.assignment_type = _DummyType(unique)
            self.additional_attributes = {"transcript_type": "protein_coding", "polya_site": -1}
            self.isoform_matches = [_DummyMatch("TX1", "G1")]
            self.chr_id = "chr9"
            self.start = exons[0][0]
            self.end = exons[-1][1]
            self.strand = "+"
            self.barcode = "AAAA"

        def set_additional_attribute(self, k, v):
            self.additional_attributes[k] = v

    uf = UMIFilter(umi_length=10, edit_distance=2, disregard_length_diff=False)
    r1 = _DummyRead("r1", "AAAAAAAAAA", [(100, 200), (300, 400)], unique=True)
    r2 = _DummyRead("r2", "AAAAAAAAAT", [(100, 200)], unique=False)
    r3 = _DummyRead("r3", "TTTTTTTTTT", [(100, 200)], unique=True)
    r_bad = _DummyRead("r4", "*", [(100, 200)], unique=True)

    kept = uf._process_duplicates([r2, r1, r3, r_bad])
    assert sorted(r.read_id for r in kept) == ["r1", "r3"]
    assert UNTRUSTED_UMI_KEY in uf._construct_umi_dict([r1, r_bad])

    bam_path = SIMPLE_DATA_DIR / "chr9.4M.ont.sim.polya.bam"
    counts = get_chromosome_read_counts(bam_path, ["chr9"])
    assert counts["chr9"] > 0

    import pysam

    with pysam.AlignmentFile(str(bam_path), "rb") as af:
        real_ids = [a.query_name for _, a in zip(range(3), af) if a.query_name]
    tbl = tmp_path / "reads.tsv"
    tbl.write_text("read_id\tbarcode\n" + "\n".join(f"{rid}\tAAAA" for rid in real_ids) + "\n")
    out_map = split_read_table_by_chromosomes(tbl, bam_path, tmp_path / "split", ["chr9"], has_header=True)
    assert "chr9" in out_map
    assert len(Path(out_map["chr9"]).read_text().splitlines()) == len(real_ids) + 1

    # create_transcript_info_dict on real iqidx built from chr9.4M.gtf.gz
    from isoquant_next.annotation.iqidx import compile_gtf_to_iqidx

    iqidx_path = compile_gtf_to_iqidx(SIMPLE_DATA_DIR / "chr9.4M.gtf.gz", tmp_path / "chr9.iqidx")
    tx_info = create_transcript_info_dict(str(iqidx_path), ["chr9"])
    assert len(tx_info) == 54
