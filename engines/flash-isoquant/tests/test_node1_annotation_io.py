from __future__ import annotations

import io
import os
import subprocess
import sys
from pathlib import Path

import gffutils
import pytest

from isoquant_next.annotation.gene_info import GeneInfo, TranscriptModel, TranscriptModelType
from isoquant_next.annotation.gtf_validator import (
    check_gtf_duplicates_streaming,
    detect_annotation_format,
    parse_gtf_attributes,
)
from isoquant_next.annotation.iqidx import (
    MmapAnnotationIndex,
    compile_gffutils_db_to_iqidx,
    compile_gtf_to_iqidx,
)
from isoquant_next.annotation.reference import ReferenceGenome
from isoquant_next.core.concurrency import WorkerPoolExecutor
from isoquant_next.core.id_policy import ExcludingIdDistributor, FeatureIdStorage, SimpleIDDistributor
from isoquant_next.core.serialization import (
    read_bool_array,
    read_dict,
    read_int,
    read_int_neg,
    read_list,
    read_list_of_pairs,
    read_short_int,
    read_string,
    read_string_or_none,
    write_bool_array,
    write_dict,
    write_int,
    write_int_neg,
    write_list,
    write_list_of_pairs,
    write_short_int,
    write_string,
    write_string_or_none,
)
from isoquant_next.core.string_pools import StringPool, StringPoolManager
from isoquant_next.mapping.read_mapper import build_minimap2_args, index_bam_safe, stable_alignment_cache_key

REPO_ROOT = Path(__file__).resolve().parents[2]
SYNTH_DB = REPO_ROOT / "isoquant_tests" / "toy_data" / "synth.db"
SIMPLE_GTF = REPO_ROOT / "isoquant_tests" / "simple_data" / "chr9.4M.gtf.gz"
SIMPLE_BAM = REPO_ROOT / "isoquant_tests" / "simple_data" / "chr9.4M.ont.sim.polya.bam"
MAPT_FA = REPO_ROOT / "isoquant_tests" / "toy_data" / "MAPT.Mouse.reference.fasta"
MAPT_GTF = REPO_ROOT / "isoquant_tests" / "toy_data" / "MAPT.Mouse.genedb.gtf"


def _worker_read_iqidx_count(arg: tuple[str, str]) -> int:
    iqidx_path, chrom = arg
    with MmapAnnotationIndex(iqidx_path) as idx:
        return sum(len(c.gene_indices) for c in idx.iter_gene_clusters(chrom))


def test_serialization_roundtrip_and_negative_dict_fix() -> None:
    cases = [
        (write_int, read_int, 0),
        (write_int, read_int, 65536),
        (write_short_int, read_short_int, 65535),
        (write_string, read_string, ""),
        (write_string, read_string, "gene-αβγ 测试 string"),
        (write_string_or_none, read_string_or_none, None),
        (write_string_or_none, read_string_or_none, "None"),
        (write_int_neg, read_int_neg, -65536),
        (write_dict, read_dict, {"neg": -42, "pair": (-10, 65536), "text": "val"}),
        (write_dict, read_dict, {}),
    ]
    for wf, rf, val in cases:
        buf = io.BytesIO()
        wf(val, buf)
        buf.seek(0)
        assert rf(buf) == val

    buf = io.BytesIO()
    write_list([(-6, -33333, 0)], buf, lambda p, f: [write_int_neg(x, f) for x in p])
    buf.seek(0)
    assert read_list(buf, lambda f: (read_int_neg(f), read_int_neg(f), read_int_neg(f))) == [(-6, -33333, 0)]

    buf2 = io.BytesIO()
    write_list_of_pairs([("a", "b"), ("", "x")], buf2, write_string)
    buf2.seek(0)
    assert read_list_of_pairs(buf2, read_string) == [("a", "b"), ("", "x")]

    buf3 = io.BytesIO()
    write_bool_array([True, False, True, False, True], buf3)
    buf3.seek(0)
    assert read_bool_array(buf3, 5) == [True, False, True, False, True]


def test_string_pool_and_manager_roundtrip(tmp_path: Path) -> None:
    pool = StringPool()
    assert len(pool) == 0
    assert pool.add("g1") == 0
    assert pool.add("g2") == 1
    assert pool.add("g1") == 0
    assert pool.get_str(1) == "g2"
    with pytest.raises(IndexError):
        pool.get_str(99)

    mgr = StringPoolManager()
    mgr.barcode_pool.add("ACGT")
    mgr.umi_pool.add("TTTT")
    mgr.set_group_spec_pool_type(0, "dynamic")
    ids = mgr.read_group_to_ids(["cell_1"])
    assert mgr.read_group_from_ids(ids) == ["cell_1"]

    dyn_file = tmp_path / "dyn.bin"
    with open(dyn_file, "wb") as fh:
        mgr.serialize_dynamic_pools(fh)
    mgr2 = StringPoolManager()
    with open(dyn_file, "rb") as fh:
        mgr2.deserialize_dynamic_pools(fh)
    assert mgr2.read_group_from_ids(ids) == ["cell_1"]


def test_gtf_validator_fixes_gff3_and_quoted_spaces(tmp_path: Path) -> None:
    assert detect_annotation_format("a.gff3") == ("gff3", False)
    assert detect_annotation_format("a.GFF3.gz") == ("gff3", True)
    assert detect_annotation_format("a.gtf.gz") == ("gtf", True)

    attrs = parse_gtf_attributes('gene_id "G1"; gene_name "HLA A 1"; exon_number 2;')
    assert attrs["gene_id"] == "G1"
    assert attrs["gene_name"] == "HLA A 1"
    assert attrs["exon_number"] == "2"

    dup_gtf = tmp_path / "dup.gtf"
    dup_gtf.write_text(
        'chr1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id "G1"; gene_name "Foo Bar";\n'
        'chr1\tsrc\tgene\t200\t300\t.\t+\t.\tgene_id "G1"; gene_name "Foo Bar";\n',
        encoding="utf-8",
    )
    modified, cleaned = check_gtf_duplicates_streaming(dup_gtf)
    assert modified and cleaned is not None
    txt = Path(cleaned).read_text(encoding="utf-8")
    assert 'gene_id "G1.dup2";' in txt
    assert 'gene_name "Foo Bar";' in txt
    os.remove(cleaned)


def test_iqidx_matches_legacy_gene_info_oracle_on_synth_db(tmp_path: Path) -> None:
    sys.path.insert(0, str(REPO_ROOT))
    from isoquant_lib.gene_info import GeneInfo as LegacyGeneInfo

    iqidx_path = compile_gffutils_db_to_iqidx(SYNTH_DB, tmp_path / "synth.iqidx")
    db = gffutils.FeatureDB(str(SYNTH_DB), keep_order=True)
    gene_db = db["ENSMUSG00000020196.10"]
    legacy_gi = LegacyGeneInfo([gene_db], db)

    with MmapAnnotationIndex(iqidx_path) as idx:
        clusters = [c for c in idx.iter_gene_clusters("chr10") if any(
            idx.get_str(int(idx.gene_id_str[gi])) == "ENSMUSG00000020196.10" for gi in c.gene_indices
        )]
        assert len(clusters) == 1
        new_gi = GeneInfo.from_iqidx_cluster(idx, clusters[0].gene_indices)

        assert new_gi.chr_id == legacy_gi.chr_id
        assert (new_gi.start, new_gi.end) == (legacy_gi.start, legacy_gi.end)
        assert new_gi.intron_profiles.features == legacy_gi.intron_profiles.features
        assert new_gi.intron_profiles.profiles == legacy_gi.intron_profiles.profiles
        assert new_gi.exon_profiles.features == legacy_gi.exon_profiles.features
        assert new_gi.exon_profiles.profiles == legacy_gi.exon_profiles.profiles
        assert new_gi.split_exon_profiles.features == legacy_gi.split_exon_profiles.features
        assert new_gi.split_exon_profiles.profiles == legacy_gi.split_exon_profiles.profiles
        assert [(f.id, f.type) for f in new_gi.exon_property_map] == [
            (f.id, f.type) for f in legacy_gi.exon_property_map
        ]
        assert [(f.id, f.type) for f in new_gi.intron_property_map] == [
            (f.id, f.type) for f in legacy_gi.intron_property_map
        ]

        # Also validate ID distributors and StringPoolManager on iqidx directly
        mgr_idx = StringPoolManager()
        mgr_idx.build_from_iqidx(idx)
        mgr_db = StringPoolManager()
        mgr_db.build_from_gffutils(db)
        assert list(mgr_idx.gene_pool) == list(mgr_db.gene_pool)
        assert list(mgr_idx.transcript_pool) == list(mgr_db.transcript_pool)

        excl = ExcludingIdDistributor(idx, "chr10")
        assert excl.increment() == 1
        fid_store = FeatureIdStorage(SimpleIDDistributor(), idx, "chr10", "exon")
        assert fid_store.get_id("chr10", (1000, 1100), "-") != ""


def test_compile_gtf_to_iqidx_and_mmap_multiprocess(tmp_path: Path) -> None:
    iqidx_path = compile_gtf_to_iqidx(SIMPLE_GTF, tmp_path / "chr9.iqidx")
    with MmapAnnotationIndex(iqidx_path) as idx:
        assert "chr9" in idx.seqids()
        clusters = list(idx.iter_gene_clusters("chr9"))
        total_genes = sum(len(c.gene_indices) for c in clusters)
        assert total_genes == 23
        gis = [GeneInfo.from_iqidx_cluster(idx, c.gene_indices) for c in clusters]
        total_tx = sum(len(gi.all_isoforms_exons) for gi in gis)
        assert total_tx == 54

    res = WorkerPoolExecutor(max_workers=2, mp_method="spawn").map_ordered(
        _worker_read_iqidx_count, [(str(iqidx_path), "chr9"), (str(iqidx_path), "chr9")]
    )
    assert res == [23, 23]

    # Validate GeneInfo.from_models compatibility
    tm = TranscriptModel("chr1", "+", "t1", "g1", [(100, 200), (300, 400)], TranscriptModelType.known)
    gi_tm = GeneInfo.from_models([tm])
    assert gi_tm.total_transcript_length("t1") == 202


def test_reference_and_read_mapper_utils(tmp_path: Path) -> None:
    with ReferenceGenome(MAPT_FA) as ref:
        lengths = ref.get_contig_lengths()
        assert len(lengths) >= 1
        first_chrom = next(iter(lengths))
        seq = ref.fetch_1based(first_chrom, 1, 20)
        assert len(seq) == 20

    assert stable_alignment_cache_key(MAPT_GTF, MAPT_FA, "nanopore") == stable_alignment_cache_key(
        MAPT_GTF, MAPT_FA, "nanopore"
    )
    # Cross-process PYTHONHASHSEED invariance
    sub_env = {**os.environ, "PYTHONPATH": str(REPO_ROOT / "isoquant_refactored" / "src")}
    out1 = subprocess.check_output(
        [
            sys.executable,
            "-c",
            "from isoquant_next.mapping.read_mapper import stable_alignment_cache_key; "
            f"print(stable_alignment_cache_key({str(MAPT_GTF)!r}, {str(MAPT_FA)!r}, 'ont'))",
        ],
        env={**sub_env, "PYTHONHASHSEED": "42"},
        text=True,
    ).strip()
    out2 = subprocess.check_output(
        [
            sys.executable,
            "-c",
            "from isoquant_next.mapping.read_mapper import stable_alignment_cache_key; "
            f"print(stable_alignment_cache_key({str(MAPT_GTF)!r}, {str(MAPT_FA)!r}, 'ont'))",
        ],
        env={**sub_env, "PYTHONHASHSEED": "99999"},
        text=True,
    ).strip()
    assert out1 == out2

    mm_args = build_minimap2_args("pacbio_ccs", MAPT_FA, MAPT_FA, threads=4)
    assert "splice:hq" in mm_args and "--MD" in mm_args

    # Test BAM indexing helper on copied simple BAM
    bam_copy = tmp_path / "test.bam"
    bam_copy.write_bytes(SIMPLE_BAM.read_bytes())
    idx_file = index_bam_safe(bam_copy, threads=1)
    assert idx_file.exists() and idx_file.name.endswith(".bai")
