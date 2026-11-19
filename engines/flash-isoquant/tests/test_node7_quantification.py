"""Node 7 (T19) tests: quantification counters, weights, Polars/SciPy COO MTX, combine_counts, and scripts."""

from __future__ import annotations

import gzip
import os
import sys
from argparse import Namespace
from pathlib import Path
from types import SimpleNamespace

import gffutils
import numpy as np
import pysam
import pytest
from scipy.io import mmread

from isoquant_next.annotation.gene_info import GeneInfo
from isoquant_next.annotation.iqidx import MmapAnnotationIndex, compile_gffutils_db_to_iqidx
from isoquant_next.assignment.profiles import Int8ProfileView
from isoquant_next.assignment.types import (
    IsoformMatch,
    MatchClassification,
    MatchEvent,
    MatchEventSubtype,
    ReadAssignment,
    ReadAssignmentType,
)
from isoquant_next.core.file_utils import merge_counts
from isoquant_next.core.stats import combine_counts
from isoquant_next.core.string_pools import StringPoolManager
from isoquant_next.core.utils import junctions_from_blocks
from isoquant_next.quantification import (
    COUNTING_STRATEGIES,
    CompositeCounter,
    ExonCounter,
    ExonSpliceSiteCounter,
    IntronCounter,
    IntronRetentionCounter,
    JointExonCounter,
    ReadWeightCounter,
    convert_profile_to_matrix,
    convert_profile_to_mtx,
    convert_to_matrix,
    convert_to_mtx,
    create_gene_counter,
    create_transcript_counter,
    load_gene_name_map,
)
from isoquant_next.scripts import (
    convert_barcode,
    convert_exon_splice_site_to_group_lists,
    convert_to_allinfo,
    convert_to_read_assignments,
    generate_spot_ids,
    inject_tags,
    load_barcoded_reads,
    load_coordinate_file,
    load_mapping,
    load_spot_map,
    transform_reads,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
TESTS_DATA = REPO_ROOT / "isoquant_tests"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


class TestCountingStrategiesAndAssignedCounters:
    @pytest.mark.parametrize("strategy", COUNTING_STRATEGIES)
    def test_weight_counter_and_gene_tx_counters_match_oracle(self, strategy, tmp_path):
        from isoquant_lib.assignment.isoform_assignment import (
            IsoformMatch as LibIM,
            MatchClassification as LibMC,
            MatchEvent as LibME,
            MatchEventSubtype as LibMES,
            ReadAssignment as LibRA,
            ReadAssignmentType as LibRAT,
        )
        from isoquant_lib.quantification.long_read_counter import (
            ReadWeightCounter as LibRWC,
            create_gene_counter as lib_gene_c,
            create_transcript_counter as lib_tx_c,
        )
        from isoquant_lib.utils.string_pools import StringPoolManager as LibSP

        rw_nxt = ReadWeightCounter(strategy)
        rw_lib = LibRWC(strategy)
        for n in (0, 1, 2, 3):
            assert rw_nxt.process_ambiguous(n) == rw_lib.process_ambiguous(n)
        for at_n, at_l in [
            (ReadAssignmentType.inconsistent, LibRAT.inconsistent),
            (ReadAssignmentType.inconsistent_ambiguous, LibRAT.inconsistent_ambiguous),
            (ReadAssignmentType.inconsistent_non_intronic, LibRAT.inconsistent_non_intronic),
        ]:
            for n in (1, 2, 3):
                assert rw_nxt.process_inconsistent(at_n, n) == rw_lib.process_inconsistent(at_l, n)

        sp_l = LibSP()
        sp_n = StringPoolManager()
        gi_stub = SimpleNamespace(all_isoforms_introns={"T1": [(150, 199)], "T2": [(150, 199)], "T_MONO": []})

        def _make_pair(rid, rat_name, matches_spec, exons=((100, 149), (200, 250))):
            ms_l = []
            ms_n = []
            for cls_name, g, t, evs in matches_spec:
                ms_l.append(
                    LibIM(
                        LibMC[cls_name],
                        sp_l,
                        assigned_gene=g,
                        assigned_transcript=t,
                        match_subclassification=[LibME(LibMES[e]) for e in evs],
                    )
                )
                ms_n.append(
                    IsoformMatch(
                        MatchClassification[cls_name],
                        sp_n,
                        assigned_gene=g,
                        assigned_transcript=t,
                        match_subclassification=[MatchEvent(MatchEventSubtype[e]) for e in evs],
                    )
                )
            rl = LibRA(rid, LibRAT[rat_name], sp_l, match=ms_l[0])
            rn = ReadAssignment(rid, ReadAssignmentType[rat_name], sp_n, match=ms_n[0])
            for m in ms_l[1:]:
                rl.add_match(m)
            for m in ms_n[1:]:
                rn.add_match(m)
            rl.gene_info = gi_stub
            rn.gene_info = gi_stub
            rl.corrected_exons = list(exons)
            rn.corrected_exons = list(exons)
            rl.exons = list(exons)
            rn.exons = list(exons)
            rl.gene_assignment_type = (
                LibRAT.unique if len({m.assigned_gene for m in ms_l if m.assigned_gene}) == 1 else LibRAT.ambiguous
            )
            rn.gene_assignment_type = (
                ReadAssignmentType.unique
                if len({m.assigned_gene for m in ms_n if m.assigned_gene}) == 1
                else ReadAssignmentType.ambiguous
            )
            return rl, rn

        reads = [
            _make_pair("r1", "unique", [("full_splice_match", "G1", "T1", ["fsm"])]),
            _make_pair("r2", "unique_minor_difference", [("incomplete_splice_match", "G1", "T1", ["ism_internal"])]),
            _make_pair("r3", "unique", [("full_splice_match", "G2", "T_MONO", ["mono_exon_match"])], exons=((100, 300),)),
            _make_pair("r4", "ambiguous", [("full_splice_match", "G1", "T1", ["fsm"]), ("full_splice_match", "G1", "T2", ["fsm"])]),
            _make_pair("r5", "inconsistent_non_intronic", [("novel_in_catalog", "G1", "T1", ["intron_retention"])]),
            _make_pair("r6", "inconsistent", [("intergenic", "G3", None, [])]),
            _make_pair("r7", "noninformative", [("intergenic", None, None, [])]),
        ]

        d_lib = tmp_path / "lib"
        d_nxt = tmp_path / "nxt"
        d_lib.mkdir()
        d_nxt.mkdir()
        gl = lib_gene_c(str(d_lib / "gene"), strategy, complete_feature_list=["G1", "G2", "G3"])
        tl = lib_tx_c(str(d_lib / "tx"), strategy, complete_feature_list=["T1", "T2", "T_MONO"])
        gn = create_gene_counter(str(d_nxt / "gene"), strategy, complete_feature_list=["G1", "G2", "G3"])
        tn = create_transcript_counter(str(d_nxt / "tx"), strategy, complete_feature_list=["T1", "T2", "T_MONO"])

        for rl, rn in reads:
            gl.add_read_info(rl)
            tl.add_read_info(rl)
            gn.add_read_info(rn)
            tn.add_read_info(rn)

        comp_lib = [gl, tl]
        comp_nxt = CompositeCounter([gn, tn])
        for c in comp_lib:
            c.dump()
        comp_nxt.dump()

        args = Namespace(normalization_method="usable_reads", counts_format=["default"])
        for c in comp_lib:
            c.finalize(args)
        comp_nxt.finalize(args)

        for fname in ("gene_counts.tsv", "gene_tpm.tsv", "tx_counts.tsv", "tx_tpm.tsv"):
            assert (d_nxt / fname).read_text() == (d_lib / fname).read_text(), f"mismatch on {strategy}/{fname}"


class TestExonIntronJointSpliceSiteOracle:
    def test_exon_and_intron_dedup_across_iqidx_builds(self, tmp_path):
        synth_db_path = TESTS_DATA / "toy_data" / "synth.db"
        iqidx_path = compile_gffutils_db_to_iqidx(synth_db_path, tmp_path / "synth.iqidx")
        with MmapAnnotationIndex(iqidx_path) as iqidx:
            gi1 = GeneInfo.from_iqidx_cluster(iqidx, [0], delta=3)
            gi2 = GeneInfo.from_iqidx_cluster(iqidx, [0], delta=3)

        exon_c = ExonCounter(str(tmp_path / "exon"))
        intron_c = IntronCounter(str(tmp_path / "intron"))
        inc = Int8ProfileView(np.ones(len(gi1.exon_property_map), dtype=np.int8))
        exc = Int8ProfileView(np.full(len(gi2.exon_property_map), -1, dtype=np.int8))
        for _ in range(3):
            exon_c.add_read_info_from_profile(inc, "-", gi1.exon_property_map, group_id=0)
        for _ in range(2):
            exon_c.add_read_info_from_profile(exc, "-", gi2.exon_property_map, group_id=0)
        exon_c.dump()

        inc_i = Int8ProfileView(np.ones(len(gi1.intron_property_map), dtype=np.int8))
        intron_c.add_read_info_from_profile(inc_i, "-", gi1.intron_property_map, group_id=0)
        intron_c.dump()

        exon_rows = Path(exon_c.output_counts_file_name).read_text().strip().splitlines()[1:]
        assert len(exon_rows) == len(gi1.exon_property_map)
        for line in exon_rows:
            cols = line.split("\t")
            assert cols[-2:] == ["3", "2"]

        ir_lines = Path(intron_c.output_counts_file_name).read_text().strip().splitlines()[1:]
        assert len(ir_lines) == len(gi1.intron_property_map)

    def test_joint_exon_and_splice_site_counters_match_oracle_exact(self, tmp_path):
        from isoquant_lib.assignment.isoform_assignment import (
            IsoformMatch as LibIM,
            MatchClassification as LibMC,
            MatchEvent as LibME,
            MatchEventSubtype as LibMES,
            ReadAssignment as LibRA,
            ReadAssignmentType as LibRAT,
        )
        from isoquant_lib.gene_info import GeneInfo as LibGeneInfo
        from isoquant_lib.quantification.long_read_counter import (
            ExonSpliceSiteCounter as LibESSC,
            IntronRetentionCounter as LibIRC,
            JointExonCounter as LibJEC,
        )
        from isoquant_lib.utils.string_pools import StringPoolManager as LibSP

        synth_db_path = TESTS_DATA / "toy_data" / "synth.db"
        db = gffutils.FeatureDB(str(synth_db_path), keep_order=True)
        gene_id = "ENSMUSG00000020196.10"
        iqidx_path = compile_gffutils_db_to_iqidx(synth_db_path, tmp_path / "synth_joint.iqidx")
        with MmapAnnotationIndex(iqidx_path) as iqidx:
            gi_n = GeneInfo.from_iqidx_cluster(iqidx, [0], delta=3)
        gi_l = LibGeneInfo([db[gene_id]], db, delta=3)

        sp_l = LibSP()
        sp_n = StringPoolManager()
        pool_l = sp_l.get_read_group_pool(0)
        pool_n = sp_n.get_read_group_pool(0)
        g0_l, g1_l = pool_l.add("c1"), pool_l.add("c2")
        g0_n, g1_n = pool_n.add("c1"), pool_n.add("c2")
        assert (g0_l, g1_l) == (g0_n, g1_n)

        tx_id = next(iter(gi_l.all_isoforms_exons))
        tx_exons = list(gi_l.all_isoforms_exons[tx_id])
        n_ex = len(gi_l.exon_profiles.features)
        n_in = len(gi_l.intron_profiles.features)

        jc_l = LibJEC(str(tmp_path / "lib_j"), sp_l, group_index=0)
        jc_n = JointExonCounter(str(tmp_path / "nxt_j"), sp_n, group_index=0)
        ess_l = LibESSC(str(tmp_path / "lib_e"), sp_l, group_index=0, emit_read_ids=True)
        ess_n = ExonSpliceSiteCounter(str(tmp_path / "nxt_e"), sp_n, group_index=0, emit_read_ids=True)
        ir_l = LibIRC(str(tmp_path / "lib_ir"), sp_l, group_index=0)
        ir_n = IntronRetentionCounter(str(tmp_path / "nxt_ir"), sp_n, group_index=0)

        for idx, (gid, blocks, evs) in enumerate(
            [
                (g0_n, tx_exons, ["fsm", "terminal_site_match_left"]),
                (g1_n, [tx_exons[0], tx_exons[-1]], ["intron_retention"]),
            ]
        ):
            ml = LibIM(LibMC.full_splice_match, sp_l, assigned_gene=gene_id, assigned_transcript=tx_id, match_subclassification=[LibME(LibMES[e], isoform_region=(0, 0)) for e in evs])
            rl = LibRA(f"r{idx}", LibRAT.unique, sp_l, match=ml)
            rl.gene_info = gi_l
            rl.strand = gi_l.gene_strands[gene_id]
            rl.exons = blocks
            rl.corrected_exons = blocks
            rl.polyA_found = True
            rl.read_group_ids = [gid]
            rl.exon_gene_profile = [1 if i == 0 else -1 for i in range(n_ex)]
            rl.intron_gene_profile = [1] * n_in

            mn = IsoformMatch(MatchClassification.full_splice_match, sp_n, assigned_gene=gene_id, assigned_transcript=tx_id, match_subclassification=[MatchEvent(MatchEventSubtype[e], isoform_region=(0, 0)) for e in evs])
            rn = ReadAssignment(f"r{idx}", ReadAssignmentType.unique, sp_n, match=mn)
            rn.gene_info = gi_n
            rn.strand = gi_n.gene_strands[gene_id]
            rn.exons = blocks
            rn.corrected_exons = blocks
            rn.polyA_found = True
            rn.read_group_ids = [gid]
            rn.exon_gene_profile = Int8ProfileView(np.asarray(rl.exon_gene_profile, dtype=np.int8))
            rn.intron_gene_profile = Int8ProfileView(np.asarray(rl.intron_gene_profile, dtype=np.int8))

            for cl, cn in ((jc_l, jc_n), (ess_l, ess_n), (ir_l, ir_n)):
                cl.add_read_info(rl)
                cn.add_read_info(rn)

        for c in (jc_l, jc_n, ess_l, ess_n, ir_l, ir_n):
            c.dump()

        assert Path(jc_n.output_counts_file_name).read_text() == Path(jc_l.output_counts_file_name).read_text()
        assert Path(ess_n.output_counts_file_name).read_text() == Path(ess_l.output_counts_file_name).read_text()
        assert Path(ir_n.output_counts_file_name).read_text() == Path(ir_l.output_counts_file_name).read_text()


class TestPolarsSciPyMTXAndCombineOracle:
    def test_grouped_matrix_mtx_and_iqidx_gene_map_match_oracle(self, tmp_path):
        from isoquant_lib.quantification.convert_grouped_counts import (
            convert_profile_to_matrix as lib_p2m,
            convert_profile_to_mtx as lib_p2mtx,
            convert_to_matrix as lib_2m,
            convert_to_mtx as lib_2mtx,
            load_gene_name_map as lib_map,
        )

        # Compare load_gene_name_map between .iqidx and .db
        synth_db = TESTS_DATA / "toy_data" / "synth.db"
        iqidx_path = compile_gffutils_db_to_iqidx(synth_db, tmp_path / "synth_nm.iqidx")
        gn_l, tn_l, t2g_l = lib_map(str(synth_db))
        gn_n, tn_n, t2g_n = load_gene_name_map(iqidx_path)
        assert gn_n == gn_l and tn_n == tn_l and t2g_n == t2g_l

        # Linear gene counts -> full matrix + gzipped MTX (raw + TPM)
        lin_gene = tmp_path / "lin_gene.tsv"
        lin_gene.write_text(
            "feature_id\tgroup_id\tcount\n"
            "G1\tcellA\t5.00\n"
            "G1\tcellB\t2.00\n"
            "G2\tcellB\t8.00\n"
        )
        usable = {"cellA": 100.0, "cellB": 200.0}

        m_lib = str(tmp_path / "lib_mat")
        m_nxt = str(tmp_path / "nxt_mat")
        assert lib_2m(str(lin_gene), m_lib, gn_l, convert_to_tpm=True, usable_reads_per_group=usable) == 2
        assert convert_to_matrix(str(lin_gene), m_nxt, gn_n, convert_to_tpm=True, usable_reads_per_group=usable) == 2
        assert Path(m_nxt + ".tsv").read_text() == Path(m_lib + ".tsv").read_text()

        lib_2mtx(str(lin_gene), str(tmp_path / "lib_mtx"), gn_l, gzipped=True, convert_to_tpm=True, usable_reads_per_group=usable)
        convert_to_mtx(str(lin_gene), str(tmp_path / "nxt_mtx"), gn_n, gzipped=True, convert_to_tpm=True, usable_reads_per_group=usable)
        assert (tmp_path / "nxt_mtx.barcodes.tsv").read_text() == (tmp_path / "lib_mtx.barcodes.tsv").read_text()
        assert (tmp_path / "nxt_mtx.features.tsv").read_text() == (tmp_path / "lib_mtx.features.tsv").read_text()
        np.testing.assert_allclose(
            mmread(str(tmp_path / "nxt_mtx.matrix.mtx.gz")).toarray(),
            mmread(str(tmp_path / "lib_mtx.matrix.mtx.gz")).toarray(),
            rtol=1e-7,
            atol=1e-9,
        )

        # Profile linear counts -> full matrix + include/exclude MTX
        lin_prof = tmp_path / "lin_prof.tsv"
        lin_prof.write_text(
            "chr\tstart\tend\tstrand\tflags\tgene_ids\tgroup_id\tinclude_counts\texclude_counts\n"
            "chr1\t10\t20\t+\tE\tG1\tcellA\t3\t1\n"
            "chr1\t10\t20\t+\tE\tG1\tcellB\t0\t4\n"
            "chr1\t50\t80\t+\tE\tG1\tcellA\t2\t0\n"
        )
        lib_p2m(str(lin_prof), str(tmp_path / "lib_pmat"), feature_type="exon")
        convert_profile_to_matrix(str(lin_prof), str(tmp_path / "nxt_pmat"), feature_type="exon")
        assert (tmp_path / "nxt_pmat.tsv").read_text() == (tmp_path / "lib_pmat.tsv").read_text()

        lib_p2mtx(str(lin_prof), str(tmp_path / "lib_pmtx"))
        convert_profile_to_mtx(str(lin_prof), str(tmp_path / "nxt_pmtx"))
        assert (tmp_path / "nxt_pmtx.barcodes.tsv").read_text() == (tmp_path / "lib_pmtx.barcodes.tsv").read_text()
        assert (tmp_path / "nxt_pmtx.features.tsv").read_text() == (tmp_path / "lib_pmtx.features.tsv").read_text()
        np.testing.assert_array_equal(
            mmread(str(tmp_path / "nxt_pmtx.include.matrix.mtx")).toarray(),
            mmread(str(tmp_path / "lib_pmtx.include.matrix.mtx")).toarray(),
        )
        np.testing.assert_array_equal(
            mmread(str(tmp_path / "nxt_pmtx.exclude.matrix.mtx")).toarray(),
            mmread(str(tmp_path / "lib_pmtx.exclude.matrix.mtx")).toarray(),
        )

    def test_combine_counts_and_merge_counts(self, tmp_path):
        from isoquant_lib.utils.stats import combine_counts as lib_combine

        s1_dir = tmp_path / "s1"
        s2_dir = tmp_path / "s2"
        s1_dir.mkdir()
        s2_dir.mkdir()
        (s1_dir / "s1.gene_counts.tsv").write_text("feature_id\tcount\nG1\t10.00\nG2\t2.00\n__ambiguous\t0\n__no_feature\t1\n__not_aligned\t0\n")
        (s2_dir / "s2.gene_counts.tsv").write_text("feature_id\tcount\nG1\t4.00\nG3\t7.00\n__ambiguous\t0\n__no_feature\t2\n__not_aligned\t0\n")
        (s1_dir / "s1.gene_tpm.tsv").write_text("feature_id\tTPM\nG1\t800000.0\nG2\t200000.0\n")
        (s2_dir / "s2.gene_tpm.tsv").write_text("feature_id\tTPM\nG1\t300000.0\nG3\t700000.0\n")
        (s1_dir / "s1.tx_counts.tsv").write_text("feature_id\tcount\nT1\t10.00\n__ambiguous\t0\n__no_feature\t0\n__not_aligned\t0\n")
        (s2_dir / "s2.tx_counts.tsv").write_text("feature_id\tcount\nT1\t4.00\n__ambiguous\t0\n__no_feature\t0\n__not_aligned\t0\n")
        (s1_dir / "s1.tx_tpm.tsv").write_text("feature_id\tTPM\nT1\t1000000.0\n")
        (s2_dir / "s2.tx_tpm.tsv").write_text("feature_id\tTPM\nT1\t1000000.0\n")

        input_data = SimpleNamespace(
            samples=[
                SimpleNamespace(prefix="S1", out_gene_counts_tsv=str(s1_dir / "s1.gene"), out_transcript_counts_tsv=str(s1_dir / "s1.tx")),
                SimpleNamespace(prefix="S2", out_gene_counts_tsv=str(s2_dir / "s2.gene"), out_transcript_counts_tsv=str(s2_dir / "s2.tx")),
            ]
        )
        out_lib = tmp_path / "out_lib"
        out_nxt = tmp_path / "out_nxt"
        out_lib.mkdir()
        out_nxt.mkdir()
        lib_combine(input_data, str(out_lib))
        combine_counts(input_data, str(out_nxt))

        for name in ("combined_gene_counts.tsv", "combined_gene_tpm.tsv", "combined_transcript_counts.tsv", "combined_transcript_tpm.tsv"):
            l_rows = sorted((out_lib / name).read_text().splitlines()[1:])
            n_rows = sorted((out_nxt / name).read_text().splitlines()[1:])
            assert l_rows == n_rows

        # merge_counts cleanly closes handles and sums stats
        mc = create_gene_counter(str(tmp_path / "SAMPLE.gene"), "unique_only")
        c1 = tmp_path / "SAMPLE_chr1.gene_counts.tsv"
        c2 = tmp_path / "SAMPLE_chr2.gene_counts.tsv"
        c1.write_text("feature_id\tcount\nG1\t3.00\n")
        c2.write_text("feature_id\tcount\nG2\t5.00\n")
        (tmp_path / "SAMPLE_chr1.gene_counts.tsv.stats").write_text("__ambiguous\t1\n__no_feature\t2\n")
        (tmp_path / "SAMPLE_chr2.gene_counts.tsv.stats").write_text("__ambiguous\t4\n__no_feature\t1\n")
        (tmp_path / "SAMPLE_chr1.gene_counts.tsv.usable").write_text("NA\t10\n")
        (tmp_path / "SAMPLE_chr2.gene_counts.tsv.usable").write_text("NA\t20\n")
        merge_counts(mc, "SAMPLE", ["chr1", "chr2"], unaligned_reads=7)
        merged_text = Path(mc.output_counts_file_name).read_text()
        assert "G1\t3.00\n" in merged_text and "G2\t5.00\n" in merged_text
        assert "__ambiguous\t5\n" in merged_text and "__not_aligned\t7\n" in merged_text
        assert mc.reads_for_tpm[0] == 30


class TestNode7AuxiliaryScripts:
    def test_scripts_match_upstream_byte_for_byte(self, tmp_path):
        from isoquant_lib.scripts.convert_read_info import (
            convert_to_allinfo as lib_allinfo,
            convert_to_read_assignments as lib_ra,
        )
        from isoquant_lib.scripts.exon_splice_site_to_group_lists import convert as lib_ess_conv
        from isoquant_lib.scripts.prepare_visium_spot_ids import generate_spot_ids as lib_gen_spots

        # 1. convert_read_info
        ri = tmp_path / "sample.read_info.tsv.gz"
        with gzip.open(ri, "wt") as fh:
            fh.write("# header\n")
            fh.write("read_id\tchr\tstrand\tgene_id\tgene_assignment_type\tisoform_id\tisoform_assignment_type\tassignment_events\tclassification\texons\tpolyA\tCAGE\tcanonical\tbarcode\tumi\tcell_type\tgroups\tadditional\n")
            fh.write("R1\tchr1\t+\tG1\tunique\tT1\tunique\tfsm,terminal_site_match_left,correct_polya_site\tfull_splice_match\t100-200,300-400\tTrue\tFalse\tTrue\tACGT\tTTTT\tcellA\tgrp1\ttranscript_type=protein_coding;\n")

        lib_ra(str(ri), str(tmp_path / "lib.ra.tsv"))
        convert_to_read_assignments(str(ri), str(tmp_path / "nxt.ra.tsv"))
        assert (tmp_path / "nxt.ra.tsv").read_text() == (tmp_path / "lib.ra.tsv").read_text()

        lib_allinfo(str(ri), str(tmp_path / "lib.allinfo"))
        convert_to_allinfo(str(ri), str(tmp_path / "nxt.allinfo"))
        assert (tmp_path / "nxt.allinfo").read_text() == (tmp_path / "lib.allinfo").read_text()

        # 2. exon_splice_site_to_group_lists
        ess_in = tmp_path / "ess_in.tsv"
        ess_in.write_text("region_gene_candidate\tn_full\tn_left\tn_right\tgroup_id\tread_ids_full\tread_ids_left\tread_ids_right\nF1\t2\t0\t0\tc1\tr1;r2\tNA\tNA\n")
        lib_ess_conv(str(ess_in), str(tmp_path / "lib_ess.tsv"))
        convert_exon_splice_site_to_group_lists(str(ess_in), str(tmp_path / "nxt_ess.tsv"))
        assert (tmp_path / "nxt_ess.tsv").read_text() == (tmp_path / "lib_ess.tsv").read_text()

        # 3. prepare_visium_spot_ids + barcodes_seq_to_spot_id
        p1 = tmp_path / "p1.tsv"
        p2 = tmp_path / "p2.tsv"
        p1.write_text("AAAA\t1\t2\n")
        p2.write_text("CCCC\t3\t4\n")
        lib_gen_spots(str(p1), str(p2), str(tmp_path / "lib_spots.tsv"), ["s_002um_", "s_008um_"], [""], "_", "|")
        generate_spot_ids(str(p1), str(p2), str(tmp_path / "nxt_spots.tsv"), ["s_002um_", "s_008um_"], [""], "_", "|")
        assert (tmp_path / "nxt_spots.tsv").read_text() == (tmp_path / "lib_spots.tsv").read_text()

        y_map = tmp_path / "ymap.tsv"
        x_map = tmp_path / "xmap.tsv"
        smap = tmp_path / "smap.tsv"
        y_map.write_text("AA\t10\n")
        x_map.write_text("CC\t20\n")
        smap.write_text("s_002um_20_10_1\ts_008um_5_2_1\ts_016um_2_1_1\n")
        reads_tsv = tmp_path / "bc_reads.tsv"
        reads_tsv.write_text("read_id\tbarcode\tumi\nr1\tAA|CC\tTTTT\n")
        transform_reads(str(reads_tsv), 1, load_mapping(str(y_map)), load_mapping(str(x_map)), load_spot_map(str(smap)), str(tmp_path / "tr"))
        assert "s_002um_20_10" in (tmp_path / "tr_2um.tsv").read_text()
        assert convert_barcode("bad", {}, {}, {}) == (None, None, None)

        # 4. inject_bam_tags on chr9.4M BAM
        src_bam = TESTS_DATA / "simple_data" / "chr9.4M.ont.sim.polya.bam"
        with pysam.AlignmentFile(str(src_bam), "rb") as bf:
            first_qname = next(iter(bf)).query_name
        bc_file = tmp_path / "barcoded.tsv"
        bc_file.write_text(f"{first_qname}\tAAACCC\tGGGGTT\n")
        out_bam = tmp_path / "tagged.bam"
        inject_tags(str(src_bam), str(out_bam), load_barcoded_reads(str(bc_file)), suffix="-1")
        with pysam.AlignmentFile(str(out_bam), "rb") as obf:
            tagged = next(iter(obf))
            assert tagged.get_tag("CB") == "AAACCC-1"
            assert tagged.get_tag("UB") == "GGGGTT"
