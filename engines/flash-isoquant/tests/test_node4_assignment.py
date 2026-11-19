"""Node 4 (T14-T15) tests: Cython profiles, comparator, LongReadAssigner Oracle diffs, CompactMultimapTuple & I/O."""

from __future__ import annotations

import os
import sys
from functools import partial
from pathlib import Path
from types import SimpleNamespace

import gffutils
import numpy as np
import pytest

from isoquant_next.annotation.gene_info import GeneInfo
from isoquant_next.annotation.iqidx import MmapAnnotationIndex, compile_gffutils_db_to_iqidx, compile_gtf_to_iqidx
from isoquant_next.alignment.alignment_info import AlignmentInfo, PolyAInfo
from isoquant_next.alignment.bam_merger import BAMStreamMerger
from isoquant_next.assignment import (
    AmbiguityResolvingMethod,
    BasicReadAssignment,
    BasicTSVAssignmentPrinter,
    BEDPrinter,
    CombinedProfileConstructor,
    CombinedReadProfiles,
    CompactMultimapTuple,
    GeneListTmpFileAssignmentLoader,
    IOSupport,
    Int8ProfileView,
    IsoformDiff,
    IsoformMatch,
    JunctionComparator,
    LongReadAssigner,
    MatchClassification,
    MatchEvent,
    MatchEventSubtype,
    MultimapResolver,
    MultimapResolvingStrategy,
    NonOverlappingFeaturesProfileConstructor,
    NormalTmpFileAssignmentLoader,
    OverlappingFeaturesProfileConstructor,
    QuickTmpFileAssignmentLoader,
    ReadAssignment,
    ReadAssignmentCompositePrinter,
    ReadAssignmentType,
    ReadInfoPrinter,
    SqantiTSVPrinter,
    TmpFileAssignmentPrinter,
    VoidPrinter,
    resolve_compact_bytes,
    resolve_compact_multimappers,
)
from isoquant_next.assignment._cy_profile_assign import count_iqidx_exon_overlaps_cy
from isoquant_next.core.string_pools import StringPoolManager
from isoquant_next.core.utils import equal_ranges, junctions_from_blocks, overlaps_at_least

REPO_ROOT = Path(__file__).resolve().parents[2]
TESTS_DATA = REPO_ROOT / "isoquant_tests"

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


def _make_params(delta: int = 3, resolve_ambiguous=AmbiguityResolvingMethod.monoexon_only, **overrides) -> SimpleNamespace:
    defaults = dict(
        delta=delta,
        minor_exon_extension=15,
        major_exon_extension=100,
        min_abs_exon_overlap=10,
        min_rel_exon_overlap=0.2,
        max_suspicious_intron_abs_len=0,
        max_suspicious_intron_rel_len=0.0,
        max_fake_terminal_exon_len=0,
        micro_intron_length=0,
        max_intron_abs_diff=5,
        max_intron_rel_diff=0.1,
        apa_delta=50,
        minimal_exon_overlap=5,
        minimal_intron_absence_overlap=5,
        max_intron_shift=10,
        max_missed_exon_len=10,
        minimal_polya_overlap=5,
        upstream_region_len=20,
        resolve_ambiguous=resolve_ambiguous,
        correct_minor_errors=True,
        needs_reference=False,
        check_canonical=True,
        cage=None,
        count_exons=True,
    )
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


def _norm_assignment_signature(ra):
    matches = tuple(
        sorted(
            (
                m.assigned_gene,
                m.assigned_transcript,
                m.match_classification.name if m.match_classification else None,
                tuple(sorted(e.event_type.name for e in m.match_subclassifications)),
            )
            for m in ra.isoform_matches
        )
    )
    return (ra.assignment_type.name, ra.gene_assignment_type.name, matches)


class TestCythonProfilesOracle:
    @pytest.mark.parametrize(
        "known_features, gene_region, read_exons, polya_pos, polyt_pos, delta, expected_gene, expected_read",
        [
            ([(50, 60), (90, 110), (140, 180)], (0, 200), [(0, 49), (61, 89), (111, 139), (181, 200)], -1, -1, 3, [1, 1, 1], [1, 1, 1]),
            ([(50, 60), (90, 110), (140, 180)], (0, 200), [(0, 50), (60, 90), (110, 140), (181, 200)], -1, -1, 3, [1, 1, 1], [1, 1, 1]),
            ([(50, 60), (90, 110), (140, 180)], (0, 200), [(0, 50), (71, 90), (110, 140), (181, 200)], -1, -1, 3, [-1, 1, 1], [-1, 1, 1]),
            ([(50, 60), (90, 110), (140, 180)], (0, 200), [(0, 50), (71, 80), (100, 130), (181, 200)], -1, -1, 3, [-1, -1, -1], [-1, -1, -1]),
            ([(50, 60), (90, 110), (140, 180)], (0, 200), [(0, 30), (71, 80), (100, 105), (186, 200)], -1, -1, 3, [-1, -1, -1], [-1, -1, -1]),
            ([(50, 60), (90, 110), (140, 180)], (0, 200), [(0, 30), (71, 80), (100, 200)], -1, -1, 3, [-1, -1, -1], [-1, -1]),
            ([(50, 60), (90, 110), (140, 180)], (0, 200), [(1000, 1030), (1071, 1080), (1100, 1200)], -1, -1, 3, [0, 0, 0], [0, 0]),
            ([(1050, 1060), (1090, 1110), (1140, 1180)], (1000, 1200), [(0, 30), (71, 80), (100, 200)], -1, -1, 3, [0, 0, 0], [0, 0]),
            ([(50, 60), (90, 110), (140, 180)], (0, 200), [(0, 200)], -1, -1, 3, [-1, -1, -1], []),
            ([(50, 60), (90, 110), (140, 180)], (0, 200), [(80, 120)], -1, -1, 3, [0, -1, 0], []),
            ([(50, 60), (90, 110), (140, 180)], (0, 200), [(0, 50), (181, 200)], -1, -1, 3, [-1, -1, -1], [-1]),
            ([(50, 60), (90, 110), (93, 107), (140, 180)], (0, 200), [(0, 49), (61, 90), (110, 139), (181, 200)], -1, -1, 3, [1, 1, -1, 1], [1, 1, 1]),
            ([(50, 60), (90, 110), (140, 180)], (0, 200), [(0, 49), (61, 89), (111, 120)], -1, -1, 3, [1, 1, 0], [1, 1]),
            ([(50, 60), (90, 110), (140, 180)], (0, 200), [(55, 80), (111, 139), (181, 200)], -1, -1, 3, [0, -1, 1], [-1, 1]),
            ([(50, 60), (90, 110), (140, 180)], (0, 200), [(62, 89), (111, 139), (181, 200)], -1, -1, 3, [0, 1, 1], [1, 1]),
            ([(50, 60), (90, 110), (140, 180)], (0, 200), [(68, 89), (111, 139), (181, 200)], -1, 68, 3, [-2, 1, 1], [1, 1]),
            ([(50, 60), (90, 110), (140, 180)], (0, 200), [(0, 49), (61, 89), (111, 120)], 119, -1, 3, [1, 1, -2], [1, 1]),
        ],
    )
    def test_intron_profile_oracle_parity(
        self, known_features, gene_region, read_exons, polya_pos, polyt_pos, delta, expected_gene, expected_read
    ):
        from isoquant_lib.assignment.long_read_profiles import OverlappingFeaturesProfileConstructor as LibOv
        from isoquant_lib.common import equal_ranges as lib_eq, overlaps_at_least as lib_ov

        lib_c = LibOv(known_features, gene_region, comparator=partial(lib_eq, delta=delta), absence_condition=partial(lib_ov, delta=10), delta=delta)
        nxt_c = OverlappingFeaturesProfileConstructor(
            known_features, gene_region, comparator=partial(equal_ranges, delta=delta), absence_condition=partial(overlaps_at_least, delta=10), delta=delta
        )
        lib_p = lib_c.construct_intron_profile(read_exons, polya_pos, polyt_pos)
        nxt_p = nxt_c.construct_intron_profile(read_exons, polya_pos, polyt_pos)
        assert isinstance(nxt_p.gene_profile, Int8ProfileView)
        assert nxt_p.gene_profile == expected_gene == lib_p.gene_profile
        assert nxt_p.read_profile == expected_read == lib_p.read_profile

    @pytest.mark.parametrize(
        "known_features, gene_region, read_exons, polya_pos, polyt_pos, delta, expected_gene, expected_read",
        [
            ([(20, 60), (90, 110), (180, 200)], (20, 200), [(20, 60), (90, 110), (180, 200)], -1, -1, 3, [1, 1, 1], [1, 1, 1]),
            ([(20, 60), (90, 110), (180, 200)], (20, 200), [(21, 60), (91, 111), (178, 200)], -1, -1, 3, [1, 1, 1], [1, 1, 1]),
            ([(20, 60), (90, 110), (180, 200)], (20, 200), [(20, 60), (99, 110), (180, 200)], -1, -1, 3, [1, -1, 1], [1, -1, 1]),
            ([(20, 60), (90, 110), (180, 200)], (20, 200), [(20, 60), (90, 110)], -1, -1, 3, [1, 1, 0], [1, 1]),
            ([(20, 60), (90, 110), (180, 200)], (20, 200), [(90, 110), (180, 200)], -1, -1, 3, [0, 1, 1], [1, 1]),
            ([(20, 60), (90, 110), (180, 200)], (20, 200), [(20, 60), (90, 100)], -1, -1, 3, [1, 0, 0], [1, -1]),
        ],
    )
    def test_exon_profile_oracle_parity(
        self, known_features, gene_region, read_exons, polya_pos, polyt_pos, delta, expected_gene, expected_read
    ):
        from isoquant_lib.assignment.long_read_profiles import OverlappingFeaturesProfileConstructor as LibOv
        from isoquant_lib.common import equal_ranges as lib_eq

        lib_c = LibOv(known_features, gene_region, comparator=partial(lib_eq, delta=delta), delta=delta)
        nxt_c = OverlappingFeaturesProfileConstructor(
            known_features, gene_region, comparator=partial(equal_ranges, delta=delta), delta=delta
        )
        lib_p = lib_c.construct_exon_profile(read_exons, polya_pos, polyt_pos)
        nxt_p = nxt_c.construct_exon_profile(read_exons, polya_pos, polyt_pos)
        assert nxt_p.gene_profile == expected_gene == lib_p.gene_profile
        assert nxt_p.read_profile == expected_read == lib_p.read_profile

    @pytest.mark.parametrize(
        "known_features, read_exons, polya_pos, polyt_pos, delta, expected_gene, expected_read",
        [
            ([(20, 60), (90, 110), (180, 200)], [(20, 60), (90, 110), (180, 200)], -1, -1, 3, [1, 1, 1], [1, 1, 1]),
            ([(20, 60), (90, 110), (180, 200)], [(20, 60), (99, 110), (180, 200)], -1, -1, 3, [1, 1, 1], [1, 1, 1]),
            ([(20, 30), (31, 40), (41, 60), (90, 93), (94, 110), (180, 200)], [(20, 60), (92, 110), (180, 200)], -1, -1, 3, [1, 1, 1, 0, 1, 1], [1, 1, 1]),
            ([(20, 30), (31, 40), (41, 60), (90, 93), (94, 110), (180, 200)], [(20, 60), (90, 110)], 110, -1, 3, [1, 1, 1, 1, 1, -2], [1, 1]),
            ([(20, 30), (31, 40), (41, 60), (90, 93), (94, 110), (180, 200)], [(90, 110), (180, 200)], -1, 90, 3, [-2, -2, -2, 1, 1, 1], [1, 1]),
            ([(20, 30), (31, 40), (41, 60), (90, 93), (94, 110), (180, 200)], [(5, 15), (70, 80), (140, 150), (160, 180), (300, 310)], -1, -1, 3, [-1, -1, -1, -1, -1, -1], [0, -1, -1, 0, 0]),
            ([(20, 30), (31, 40), (41, 60), (90, 93), (94, 110), (180, 200)], [(30, 60), (70, 80), (100, 110), (160, 179)], -1, -1, 3, [0, 1, 1, -1, 1, 0], [1, -1, 1, -1]),
        ],
    )
    def test_split_exon_profile_oracle_parity(
        self, known_features, read_exons, polya_pos, polyt_pos, delta, expected_gene, expected_read
    ):
        from isoquant_lib.assignment.long_read_profiles import NonOverlappingFeaturesProfileConstructor as LibNonOv
        from isoquant_lib.common import overlaps_at_least as lib_ov

        lib_c = LibNonOv(known_features, comparator=partial(lib_ov, delta=5), delta=delta)
        nxt_c = NonOverlappingFeaturesProfileConstructor(
            known_features, comparator=partial(overlaps_at_least, delta=5), delta=delta
        )
        lib_p = lib_c.construct_profile(read_exons, polya_pos, polyt_pos)
        nxt_p = nxt_c.construct_profile(read_exons, polya_pos, polyt_pos)
        assert nxt_p.gene_profile == expected_gene == lib_p.gene_profile
        assert nxt_p.read_profile == expected_read == lib_p.read_profile

    def test_iqidx_mmap_zero_copy_exon_overlap_scan(self, tmp_path):
        idx_path = compile_gtf_to_iqidx(TESTS_DATA / "simple_data" / "chr9.4M.gtf.gz", tmp_path / "chr9.iqidx")
        with MmapAnnotationIndex(idx_path) as iqidx:
            clusters = list(iqidx.iter_gene_clusters("chr9"))
            assert len(clusters) > 0
            c0 = clusters[0]
            tx_0 = int(iqidx.gene_tx_ptr[c0.gene_indices[0]])
            tx_1 = int(iqidx.gene_tx_ptr[c0.gene_indices[-1] + 1])
            ex_0 = int(iqidx.tx_exon_ptr[tx_0])
            ex_1 = int(iqidx.tx_exon_ptr[tx_1])
            first_ex_s = int(iqidx.exon_start[ex_0])
            first_ex_e = int(iqidx.exon_end[ex_0])
            hits = count_iqidx_exon_overlaps_cy(
                iqidx.exon_start,
                iqidx.exon_end,
                ex_0,
                ex_1,
                np.array([first_ex_s], dtype=np.int32),
                np.array([first_ex_e], dtype=np.int32),
            )
            assert hits >= 1


class TestLongReadAssignerOracle:
    def test_profile_matching_and_junction_comparator_vectors(self):
        from isoquant_lib.assignment.long_read_assigner import LongReadAssigner as LibAssigner
        from isoquant_lib.assignment.long_read_profiles import MappedReadProfile as LibProf
        from isoquant_lib.gene_info import FeatureProfiles as LibFP

        # test_long_read_assigner match_profile cases
        with pytest.raises(AssertionError):
            LongReadAssigner.match_profile([], {"id1": [1], "id2": [-1]}, None)
        cases = [
            ([0, 1, 1], {"id1": [1, 1, -1], "id2": [1, -1, 1], "id3": [-1, -1, -1]}, None, ["id1"]),
            ([-1, 1, -1, 0], {"id1": [-1, 1, -1, -1], "id2": [1, 1, 1, -1]}, None, ["id1", "id2"]),
            ([0, 1, 1, -1, 0], {"id1": [-1, 1, 1, 1, -2], "id2": [-2, 1, 1, -1, -2], "id3": [1, -1, 1, -2, -2]}, None, ["id2", "id1", "id3"]),
            ([-1, 1, 1, 0], {"id1": [-1, 1, 1, -1], "id2": [1, 1, 1, -1], "id3": [-1, -1, -1, 1]}, {"id2", "id3"}, ["id2", "id3"]),
        ]
        for rgp, iso_profs, hint, expected_ids in cases:
            nxt_res = LongReadAssigner.match_profile(rgp, iso_profs, hint)
            assert [x.isoform_id for x in nxt_res[: len(expected_ids)]] == expected_ids
            with pytest.raises(AssertionError):
                LongReadAssigner.match_profile(rgp + [0], iso_profs, hint)

        # test_junction_comparator.py checks
        assert overlaps_at_least((100, 200), (150, 250), 10) is True
        assert equal_ranges((100, 200), (102, 198), delta=5) is True
        assert equal_ranges((100, 200), (102, 198), delta=1) is False
        assert junctions_from_blocks([(100, 200), (300, 400), (500, 600)]) == [(201, 299), (401, 499)]
        assert junctions_from_blocks([(100, 200)]) == []

        # CompareLib vs Next compare_junctions on multiple contradiction cases
        sp = StringPoolManager()
        gi = GeneInfo.from_region("chr1", 0, 500, delta=3)
        gi.intron_profiles.set_features([(50, 60), (80, 100), (80, 110), (80, 150), (200, 210)])
        assigner = LongReadAssigner(gi, _make_params(3, max_suspicious_intron_abs_len=10, max_suspicious_intron_rel_len=0.1), sp)

        lib_gi = SimpleNamespace(
            chr_id="chr1",
            start=0,
            end=500,
            intron_profiles=LibFP(),
            all_isoforms_exons={},
            all_isoforms_introns={},
            isoform_strands={},
        )
        lib_gi.intron_profiles.set_features([(50, 60), (80, 100), (80, 110), (80, 150), (200, 210)])
        from isoquant_lib.utils.string_pools import StringPoolManager as LibSP

        lib_assigner = LibAssigner(lib_gi, _make_params(3, max_suspicious_intron_abs_len=10, max_suspicious_intron_rel_len=0.1), LibSP())

        junc_cases = [
            ([], (20, 200), [], (20, 290)),
            ([], (20, 200), [(150, 170)], (20, 290)),
            ([], (20, 100), [(50, 170)], (20, 290)),
            ([], (150, 320), [(50, 170)], (20, 290)),
            ([(10, 50)], (0, 100), [(10, 25), (40, 49)], (0, 99)),
            ([(80, 110)], (50, 150), [(80, 90), (105, 110)], (50, 149)),
            ([(50, 60), (81, 149)], (0, 200), [(50, 100), (119, 150)], (0, 219)),
            ([(50, 60), (80, 110)], (0, 150), [(50, 60), (80, 100)], (9, 149)),
        ]
        for rj, rr, ij, ir in junc_cases:
            ev_lib = [e.event_type.name for e in lib_assigner.intron_comparator.compare_junctions(rj, rr, ij, ir)]
            ev_nxt = [e.event_type.name for e in assigner.intron_comparator.compare_junctions(rj, rr, ij, ir)]
            assert ev_nxt == ev_lib

    def test_synth_db_vs_iqidx_and_oracle_parity(self, tmp_path):
        from isoquant_lib.assignment.long_read_assigner import (
            AmbiguityResolvingMethod as LibARM,
            LongReadAssigner as LibAssigner,
        )
        from isoquant_lib.assignment.long_read_profiles import (
            CombinedProfileConstructor as LibCombinedC,
        )
        from isoquant_lib.gene_info import GeneInfo as LibGeneInfo
        from isoquant_lib.utils.string_pools import StringPoolManager as LibSP

        synth_db_path = TESTS_DATA / "toy_data" / "synth.db"
        db = gffutils.FeatureDB(str(synth_db_path), keep_order=True)
        gene_id = "ENSMUSG00000020196.10"
        params_lib = _make_params(3, resolve_ambiguous=LibARM.monoexon_only, max_missed_exon_len=20, max_fake_terminal_exon_len=20)
        params_nxt = _make_params(3, resolve_ambiguous=AmbiguityResolvingMethod.monoexon_only, max_missed_exon_len=20, max_fake_terminal_exon_len=20)

        # Build iqidx directly from synth.db so we prove zero-SQLite GeneInfo cluster equivalence too
        iqidx_path = compile_gffutils_db_to_iqidx(synth_db_path, tmp_path / "synth.iqidx")
        with MmapAnnotationIndex(iqidx_path) as iqidx:
            gi_nxt = None
            for clust in iqidx.iter_gene_clusters(db[gene_id].seqid):
                idxs = list(clust.gene_indices)
                if any(iqidx.get_str(int(iqidx.gene_id_str[i])) == gene_id for i in idxs):
                    gi_nxt = GeneInfo.from_iqidx_cluster(iqidx, idxs, delta=params_nxt.delta)
                    break
            assert gi_nxt is not None

        gi_lib = LibGeneInfo([db[gene_id]], db, delta=params_lib.delta)
        sp_lib = LibSP()
        sp_nxt = StringPoolManager()
        assigner_lib = LibAssigner(gi_lib, params_lib, sp_lib)
        assigner_nxt = LongReadAssigner(gi_nxt, params_nxt, sp_nxt)
        prof_c_lib = LibCombinedC(gi_lib, params_lib)
        prof_c_nxt = CombinedProfileConstructor(gi_nxt, params_nxt)

        test_cases = [
            ([(1000, 1100), (2000, 2100), (2300, 2400), (3000, 3300), (9500, 10000)], -1, -1),
            ([(1000, 1100), (2000, 2100), (2300, 2400), (3000, 3300)], -1, -1),
            ([(2000, 2098), (2301, 2400), (3001, 3300), (9500, 10003)], -1, -1),
            ([(2000, 2100), (2300, 2400), (3000, 3300)], -1, -1),
            ([(7000, 7500)], -1, -1),
            ([(1000, 1108), (2007, 2100), (2300, 2400), (3000, 3300), (9500, 10000)], -1, -1),
            ([(7998, 8201), (8500, 8800), (9900, 9907)], -1, -1),
            ([(1000, 1108), (2007, 2100)], -1, -1),
            ([(1000, 1200), (2000, 2100), (2300, 2400), (3000, 3300), (9500, 10000)], -1, -1),
            ([(1000, 1100), (2000, 2400), (3000, 3300), (9500, 10000)], -1, -1),
            ([(900, 1100), (2000, 2100), (2300, 2400), (3000, 3300), (9500, 10000)], -1, 900),
            ([(1000, 2100)], -1, -1),
            ([(4000, 4100)], -1, -1),
            ([(gi_lib.end + 500, gi_lib.end + 600)], -1, -1),
        ]

        from isoquant_lib.alignment.alignment_info import PolyAInfo as LibPA

        for i, (blocks, pa, pt) in enumerate(test_cases):
            cp_lib = prof_c_lib.construct_profiles(blocks, LibPA(pa, pt, -1, -1), [])
            cp_nxt = prof_c_nxt.construct_profiles(blocks, PolyAInfo(pa, pt, -1, -1), [])
            ra_lib = assigner_lib.assign_to_isoform(f"r{i}", cp_lib)
            ra_nxt = assigner_nxt.assign_to_isoform(f"r{i}", cp_nxt)
            assert _norm_assignment_signature(ra_nxt) == _norm_assignment_signature(ra_lib)

    def test_chr9_real_bam_reads_assignment_oracle_parity(self, tmp_path):
        from isoquant_lib.assignment.long_read_assigner import (
            AmbiguityResolvingMethod as LibARM,
            LongReadAssigner as LibAssigner,
        )
        from isoquant_lib.assignment.long_read_profiles import CombinedProfileConstructor as LibCombinedC
        from isoquant_lib.alignment.alignment_info import PolyAInfo as LibPA
        from isoquant_lib.gene_info import GeneInfo as LibGeneInfo
        from isoquant_lib.gtf2db import gtf2db
        from isoquant_lib.utils.string_pools import StringPoolManager as LibSP

        gtf_path = TESTS_DATA / "simple_data" / "chr9.4M.gtf.gz"
        bam_path = TESTS_DATA / "simple_data" / "chr9.4M.ont.sim.polya.bam"
        params_lib = _make_params(3, resolve_ambiguous=LibARM.monoexon_only)
        params_nxt = _make_params(3, resolve_ambiguous=AmbiguityResolvingMethod.monoexon_only)

        db_path = tmp_path / "chr9.db"
        gtf2db(str(gtf_path), str(db_path), complete_db=True)
        db = gffutils.FeatureDB(str(db_path), keep_order=True)
        iqidx_path = compile_gtf_to_iqidx(gtf_path, tmp_path / "chr9.iqidx")

        sp_lib = LibSP()
        sp_nxt = StringPoolManager()
        with MmapAnnotationIndex(iqidx_path) as iqidx:
            import pysam

            with pysam.AlignmentFile(str(bam_path), "rb") as bf:
                read_exons_list = []
                for aln in bf.fetch("chr9"):
                    if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                        continue
                    ai = AlignmentInfo(aln)
                    if ai.read_exons:
                        read_exons_list.append((aln.query_name, ai.read_exons, PolyAInfo(-1, -1, -1, -1)))
            assert len(read_exons_list) > 50

            compared = 0
            for clust in iqidx.iter_gene_clusters("chr9"):
                idxs = list(clust.gene_indices)
                gi_nxt = GeneInfo.from_iqidx_cluster(iqidx, idxs, delta=params_nxt.delta)
                gene_ids = [iqidx.get_str(int(iqidx.gene_id_str[i])) for i in idxs]
                gi_lib = LibGeneInfo([db[gid] for gid in gene_ids], db, delta=params_lib.delta)

                assigner_lib = LibAssigner(gi_lib, params_lib, sp_lib)
                assigner_nxt = LongReadAssigner(gi_nxt, params_nxt, sp_nxt)
                pc_lib = LibCombinedC(gi_lib, params_lib)
                pc_nxt = CombinedProfileConstructor(gi_nxt, params_nxt)

                for qname, exons, pa in read_exons_list:
                    if exons[-1][1] < gi_nxt.start or exons[0][0] > gi_nxt.end:
                        continue
                    cp_lib = pc_lib.construct_profiles(
                        exons, LibPA(pa.external_polya_pos, pa.external_polyt_pos, pa.internal_polya_pos, pa.internal_polyt_pos), []
                    )
                    cp_nxt = pc_nxt.construct_profiles(exons, pa, [])
                    ra_lib = assigner_lib.assign_to_isoform(qname, cp_lib)
                    ra_nxt = assigner_nxt.assign_to_isoform(qname, cp_nxt)
                    assert _norm_assignment_signature(ra_nxt) == _norm_assignment_signature(ra_lib)
                    compared += 1
            assert compared > 50


class TestMultimapAndContextIO:
    def test_compact_multimap_bytes_roundtrip_and_resolution(self):
        sp = StringPoolManager()
        sp.chromosome_pool.add("chr1")
        sp.gene_pool.add("G1")
        sp.gene_pool.add("G2")
        sp.transcript_pool.add("T1")
        sp.transcript_pool.add("T2")

        m1 = IsoformMatch(MatchClassification.full_splice_match, sp, assigned_gene="G1", assigned_transcript="T1")
        ra1 = ReadAssignment("r1", ReadAssignmentType.unique, sp, match=m1)
        ra1.chr_id = "chr1"
        ra1.exons = [(100, 200)]
        ra1.genomic_region = (90, 250)

        m2 = IsoformMatch(MatchClassification.novel_not_in_catalog, sp, assigned_gene="G2", assigned_transcript="T2")
        ra2 = ReadAssignment("r1", ReadAssignmentType.inconsistent, sp, match=m2)
        ra2.chr_id = "chr1"
        ra2.exons = [(500, 600)]
        ra2.genomic_region = (450, 700)
        ra2.multimapper = True

        c1 = BasicReadAssignment(ra1).to_compact_tuple()
        c2 = BasicReadAssignment(ra2).to_compact_tuple()
        resolved_bytes = resolve_compact_bytes([c1.to_bytes(), c2.to_bytes()], MultimapResolvingStrategy.take_best)
        r_recs = [CompactMultimapTuple.from_bytes(b) for b in resolved_bytes]
        assert r_recs[0].assignment_type == ReadAssignmentType.unique
        assert r_recs[1].assignment_type == ReadAssignmentType.suspended

        # Merge strategy with >1 informative records (verifies set-to-list fix too)
        r_merge = resolve_compact_multimappers(
            [BasicReadAssignment(ra1).to_compact_tuple(), BasicReadAssignment(ra2).to_compact_tuple()],
            MultimapResolvingStrategy.merge,
        )
        assert all(r.multimapper for r in r_merge)
        assert r_merge[0].assignment_type == ReadAssignmentType.ambiguous
        assert r_merge[1].assignment_type == ReadAssignmentType.inconsistent_ambiguous

    def test_context_printers_and_iqidx_iosupport_parity(self, tmp_path):
        from isoquant_lib.assignment.assignment_io import (
            BasicTSVAssignmentPrinter as LibTSVP,
            IOSupport as LibIOS,
            ReadInfoPrinter as LibRIP,
            SqantiTSVPrinter as LibSqP,
        )
        from isoquant_lib.gene_info import GeneInfo as LibGeneInfo
        from isoquant_lib.utils.string_pools import StringPoolManager as LibSP

        synth_db_path = TESTS_DATA / "toy_data" / "synth.db"
        db = gffutils.FeatureDB(str(synth_db_path), keep_order=True)
        gene_id = "ENSMUSG00000020196.10"
        params = _make_params(3, check_canonical=False)

        iqidx_path = compile_gffutils_db_to_iqidx(synth_db_path, tmp_path / "synth.iqidx")
        with MmapAnnotationIndex(iqidx_path) as iqidx:
            gi_nxt = GeneInfo.from_iqidx_cluster(iqidx, [0], delta=params.delta)
        gi_lib = LibGeneInfo([db[gene_id]], db, delta=params.delta)

        sp_lib = LibSP()
        sp_nxt = StringPoolManager()
        blocks = [(1000, 1100), (2000, 2100), (2300, 2400), (3000, 3300), (9500, 10000)]

        from isoquant_lib.assignment.isoform_assignment import IsoformMatch as LibIM, MatchClassification as LibMC, MatchEvent as LibME, MatchEventSubtype as LibMES, ReadAssignment as LibRA, ReadAssignmentType as LibRAT

        m_lib = LibIM(LibMC.full_splice_match, sp_lib, assigned_gene=gene_id, assigned_transcript="ENSMUST00000001712.7", match_subclassification=[LibME(LibMES.fsm)])
        ra_lib = LibRA("r_tsv", LibRAT.unique, sp_lib, match=m_lib)
        ra_lib.gene_info = gi_lib
        ra_lib.exons = blocks
        ra_lib.corrected_exons = blocks
        ra_lib.chr_id = gi_lib.chr_id
        ra_lib.strand = "+"
        ra_lib.mapped_strand = "+"
        ra_lib.additional_info["FSM_class"] = "A"

        m_nxt = IsoformMatch(MatchClassification.full_splice_match, sp_nxt, assigned_gene=gene_id, assigned_transcript="ENSMUST00000001712.7", match_subclassification=[MatchEvent(MatchEventSubtype.fsm)])
        ra_nxt = ReadAssignment("r_tsv", ReadAssignmentType.unique, sp_nxt, match=m_nxt)
        ra_nxt.gene_info = gi_nxt
        ra_nxt.exons = blocks
        ra_nxt.corrected_exons = blocks
        ra_nxt.chr_id = gi_nxt.chr_id
        ra_nxt.strand = "+"
        ra_nxt.mapped_strand = "+"
        ra_nxt.additional_info["FSM_class"] = "A"

        lib_ri = tmp_path / "lib.read_info.tsv"
        nxt_ri = tmp_path / "nxt.read_info.tsv"
        lp = LibRIP(str(lib_ri), params, LibIOS(params))
        lp.add_read_info(ra_lib)
        lp.output_file.close()
        with ReadInfoPrinter(str(nxt_ri), params, IOSupport(params)) as np_p:
            np_p.add_read_info(ra_nxt)
        assert nxt_ri.read_text() == lib_ri.read_text()

        lib_sq = tmp_path / "lib.sqanti.tsv"
        nxt_sq = tmp_path / "nxt.sqanti.tsv"
        lsq = LibSqP(str(lib_sq), params, LibIOS(params))
        lsq.add_read_info(ra_lib)
        lsq.output_file.close()
        with SqantiTSVPrinter(str(nxt_sq), params, IOSupport(params)) as nsq:
            nsq.add_read_info(ra_nxt)
        assert nxt_sq.read_text() == lib_sq.read_text()

        # Binary TmpFile roundtrip via MmapAnnotationIndex (zero SQLite!)
        bin_file = tmp_path / "assign.bin"
        with TmpFileAssignmentPrinter(str(bin_file), params) as tp:
            tp.add_gene_info(gi_nxt)
            tp.add_read_info(ra_nxt)

        with MmapAnnotationIndex(iqidx_path) as iqidx:
            with NormalTmpFileAssignmentLoader(str(bin_file), iqidx, None, sp_nxt) as loader:
                assert loader.has_next() and loader.is_gene_info()
                g_loaded = loader.get_object()
                assert g_loaded.chr_id == gi_nxt.chr_id
                assert loader.has_next() and loader.is_read_assignment()
                r_loaded = loader.get_object()
                assert r_loaded.read_id == "r_tsv"
                assert not loader.has_next()

            with QuickTmpFileAssignmentLoader(str(bin_file), sp_nxt) as qloader:
                assert qloader.get_object() is None
                bra = qloader.get_object()
                assert isinstance(bra, BasicReadAssignment)
                assert bra.read_id == "r_tsv"

    def test_assess_assignment_quality_compatibility(self, tmp_path):
        from misc.assess_assignment_quality import AssignmentData

        synth_db_path = TESTS_DATA / "toy_data" / "synth.db"
        gene_id = "ENSMUSG00000020196.10"
        tx_id = "ENSMUST00000001712.7"
        params = _make_params(3, check_canonical=False)
        iqidx_path = compile_gffutils_db_to_iqidx(synth_db_path, tmp_path / "synth_qa.iqidx")
        with MmapAnnotationIndex(iqidx_path) as iqidx:
            gi_nxt = GeneInfo.from_iqidx_cluster(iqidx, [0], delta=params.delta)

        sp_nxt = StringPoolManager()
        blocks = [(1000, 1100), (2000, 2100), (2300, 2400), (3000, 3300), (9500, 10000)]

        m_fsm = IsoformMatch(
            MatchClassification.full_splice_match,
            sp_nxt,
            assigned_gene=gene_id,
            assigned_transcript=tx_id,
            match_subclassification=[MatchEvent(MatchEventSubtype.fsm)],
        )
        ra_fsm = ReadAssignment("r_fsm", ReadAssignmentType.unique, sp_nxt, match=m_fsm)
        ra_fsm.gene_info = gi_nxt
        ra_fsm.exons = blocks
        ra_fsm.corrected_exons = blocks
        ra_fsm.chr_id = gi_nxt.chr_id
        ra_fsm.strand = "+"
        ra_fsm.mapped_strand = "+"

        m_nic = IsoformMatch(MatchClassification.novel_in_catalog, sp_nxt, assigned_gene=gene_id, assigned_transcript=None)
        ra_inc = ReadAssignment("r_inc", ReadAssignmentType.inconsistent, sp_nxt, match=m_nic)
        ra_inc.gene_info = gi_nxt
        ra_inc.exons = blocks
        ra_inc.corrected_exons = blocks
        ra_inc.chr_id = gi_nxt.chr_id
        ra_inc.strand = "+"
        ra_inc.mapped_strand = "+"

        ri_path = tmp_path / "nxt.read_info.tsv"
        with ReadInfoPrinter(str(ri_path), params, IOSupport(params)) as rip:
            rip.add_read_info(ra_fsm)
            rip.add_read_info(ra_inc)

        basic_path = tmp_path / "nxt.read_assignments.tsv"
        with BasicTSVAssignmentPrinter(str(basic_path), params, IOSupport(params)) as btp:
            btp.add_read_info(ra_fsm)
            btp.add_read_info(ra_inc)

        for tsv_path in (ri_path, basic_path):
            for preset in ("isoquant", "isoquant_legacy"):
                ad_unique = AssignmentData(str(tsv_path), preset=preset)
                assert dict(ad_unique.assigned_isoforms) == {"r_fsm": tx_id}
                ad_all = AssignmentData(str(tsv_path), preset=preset, assignment_types=AssignmentData.ALL_ASSIGNMENTS_TYPES)
                assert dict(ad_all.assigned_isoforms) == {"r_fsm": tx_id, "r_inc": "."}


def test_canonical_ambiguity_resolving_method_identity_from_config_enums(tmp_path) -> None:
    from isoquant_next.assignment.assigner import AmbiguityResolvingMethod as AssignerARM
    from isoquant_next.assignment import AmbiguityResolvingMethod as InitARM
    from isoquant_next.config.enums import AmbiguityResolvingMethod as ConfigARM

    assert AssignerARM is ConfigARM
    assert InitARM is ConfigARM

    # Ensure LongReadAssigner honors ConfigARM.monoexon_and_fsm (not an identity-mismatched duplicate Enum)
    gtf = tmp_path / "ambig.gtf"
    gtf.write_text(
        "\n".join(
            [
                'chr1\tT\tgene\t100\t900\t.\t+\t.\tgene_id "g1";',
                'chr1\tT\ttranscript\t100\t400\t.\t+\t.\tgene_id "g1"; transcript_id "tx_fsm";',
                'chr1\tT\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "tx_fsm";',
                'chr1\tT\texon\t300\t400\t.\t+\t.\tgene_id "g1"; transcript_id "tx_fsm";',
                'chr1\tT\ttranscript\t100\t900\t.\t+\t.\tgene_id "g1"; transcript_id "tx_ism";',
                'chr1\tT\texon\t100\t200\t.\t+\t.\tgene_id "g1"; transcript_id "tx_ism";',
                'chr1\tT\texon\t300\t400\t.\t+\t.\tgene_id "g1"; transcript_id "tx_ism";',
                'chr1\tT\texon\t500\t900\t.\t+\t.\tgene_id "g1"; transcript_id "tx_ism";',
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    iqidx_path = compile_gtf_to_iqidx(gtf, tmp_path / "ambig.iqidx")
    params = _make_params(3, resolve_ambiguous=ConfigARM.monoexon_and_fsm)
    with MmapAnnotationIndex(iqidx_path) as iqidx:
        gi = GeneInfo.from_iqidx_cluster(iqidx, [0], delta=params.delta)

    sp = StringPoolManager()
    assigner = LongReadAssigner(gi, params, sp)
    cons = CombinedProfileConstructor(gi, params)
    prof = cons.construct_profiles([(105, 200), (300, 395)], PolyAInfo(-1, -1, -1, -1), [])
    ra = assigner.assign_to_isoform("r_fsm_vs_ism", prof)
    assert ra.assignment_type == ReadAssignmentType.unique
    assert [m.assigned_transcript for m in ra.isoform_matches] == ["tx_fsm"]

