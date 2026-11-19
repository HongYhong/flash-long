"""Node 6 (T18) tests: ChromosomeDiscoveryState isolation, IntronGraph, alt-end NICs, splice-site correction, 5-stage E2E."""

from __future__ import annotations

import gzip
import shutil
import sys
import types
from argparse import Namespace
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pysam
import pytest

from isoquant_next.alignment.alignment_info import AlignmentInfo, PolyAInfo
from isoquant_next.annotation.gene_info import GeneInfo, TranscriptModel, TranscriptModelType
from isoquant_next.annotation.iqidx import MmapAnnotationIndex, compile_gtf_to_iqidx
from isoquant_next.annotation.reference import ReferenceGenome
from isoquant_next.assignment.assigner import AmbiguityResolvingMethod, LongReadAssigner
from isoquant_next.assignment.profiles import CombinedProfileConstructor
from isoquant_next.core.id_policy import FeatureIdStorage, SimpleIDDistributor
from isoquant_next.core.string_pools import StringPoolManager
from isoquant_next.core.utils import junctions_from_blocks
from isoquant_next.discovery import (
    ChromosomeDiscoveryState,
    GFFPrinter,
    GraphBasedModelConstructor,
    IntronGraph,
    ModelFilter,
    ModelStore,
    SpliceSiteCase,
    StrandnessReportingLevel,
    TerminalVertex,
    TranscriptEndProcessor,
    compute_most_common_case_of_deletions,
    correct_splice_site_errors,
    count_deletions_for_splice_site_locations,
    create_extended_storage,
    extract_relevant_deletions,
    extract_splice_site_locations_within_aligned_read,
    generate_updated_exon_list,
)
from isoquant_next.terminal.peaks import Peak
from isoquant_next.terminal.polya_finder import PolyAFinder

REPO_ROOT = Path(__file__).resolve().parents[2]
TESTS_DATA = REPO_ROOT / "isoquant_tests"

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


class _StubModel:
    def __init__(self, accept: bool = True) -> None:
        self._accept = accept

    def predict(self, X):
        return np.full(len(X), 1 if self._accept else 0, dtype=np.int32)


class _DummyCounter:
    def __init__(self) -> None:
        self.assignments = []
        self.confirmed = []

    def add_read_info(self, ra) -> None:
        self.assignments.append(ra)

    def add_confirmed_features(self, feats) -> None:
        self.confirmed.extend(list(feats))


def _make_discovery_args(**overrides) -> SimpleNamespace:
    defaults = dict(
        delta=3,
        minor_exon_extension=15,
        major_exon_extension=100,
        min_abs_exon_overlap=10,
        min_rel_exon_overlap=0.2,
        max_suspicious_intron_abs_len=20,
        max_suspicious_intron_rel_len=0.2,
        max_fake_terminal_exon_len=10,
        micro_intron_length=0,
        max_intron_abs_diff=5,
        max_intron_rel_diff=0.1,
        apa_delta=10,
        minimal_exon_overlap=5,
        minimal_intron_absence_overlap=5,
        max_intron_shift=10,
        max_missed_exon_len=10,
        minimal_polya_overlap=5,
        upstream_region_len=20,
        resolve_ambiguous=AmbiguityResolvingMethod.monoexon_only,
        correct_minor_errors=True,
        needs_reference=True,
        check_canonical=True,
        count_exons=True,
        cage=None,
        genedb="dummy.gtf",
        fl_data=False,
        fl_only=False,
        novel_apa=False,
        min_known_count=1,
        min_novel_count=2,
        min_novel_intron_count=2,
        min_novel_isolated_intron_abs=3,
        min_novel_count_rel=0.05,
        min_mono_count_rel=0.05,
        min_mono_exon_coverage=0.5,
        require_monointronic_polya=False,
        require_monoexonic_polya=False,
        requires_polya_for_construction=False,
        report_canonical_strategy=StrandnessReportingLevel.all,
        report_novel_unspliced=True,
        simple_alignments_mapq_cutoff=0,
        simple_models_mapq_cutoff=0,
        graph_clustering_distance=6,
        graph_clustering_ratio=0.2,
        singleton_adjacent_cov=100,
        terminal_position_abs=1,
        terminal_position_rel=0.1,
        terminal_internal_position_rel=0.05,
        no_splice_site_correction=False,
        sqanti_output=True,
        debug=False,
    )
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


def _bare_graph(apa_delta=10, polya_predictions=None, tss_predictions=None):
    g = IntronGraph.__new__(IntronGraph)
    g.params = SimpleNamespace(
        apa_delta=apa_delta,
        terminal_position_abs=1,
        terminal_position_rel=0.0,
        terminal_internal_position_rel=0.0,
    )
    g.polya_predictions = polya_predictions
    g.tss_predictions = tss_predictions
    g.outgoing_edges = defaultdict(set)
    g.incoming_edges = defaultdict(set)
    return g


class TestChromosomeStateIsolation:
    def test_explicit_chr_state_isolates_workers_and_dedups_within_chr(self):
        gi = GeneInfo.from_region("chr1", 1, 1000)
        state_a = ChromosomeDiscoveryState()
        state_b = ChromosomeDiscoveryState()
        store_a1 = ModelStore(gi, chr_state=state_a)
        store_a2 = ModelStore(gi, chr_state=state_a)
        store_b = ModelStore(gi, chr_state=state_b)

        store_a1.detected_known_isoforms.add("TX_KNOWN_1")
        assert "TX_KNOWN_1" in store_a2.detected_known_isoforms
        assert "TX_KNOWN_1" not in store_b.detected_known_isoforms
        # Default store with no explicit chr_state never mutates shared class globals
        store_default = ModelStore(gi)
        assert "TX_KNOWN_1" not in store_default.detected_known_isoforms


class TestIntronGraphAndEndProcessorOracle:
    def test_intron_graph_refinement_and_side_selection(self):
        from isoquant_lib.model_construction.intron_graph import (
            IntronGraph as LibIG,
            TerminalVertex as LibTV,
        )

        for hist, preds in [({105: 3}, [100, 500]), ({130: 3}, [100, 500]), ({98: 2, 103: 5}, [100]), ({105: 3}, None)]:
            g_nxt = _bare_graph(apa_delta=10)
            g_lib = LibIG.__new__(LibIG)
            g_lib.params = g_nxt.params
            assert g_nxt._refine_positions(hist, preds) == g_lib._refine_positions(hist, preds)

        intron = (100, 200)
        g = _bare_graph(apa_delta=10, polya_predictions=[1000], tss_predictions=[1010])
        g.intron_collector = SimpleNamespace(clustered_introns={intron: 10})
        g.cluster_polya_positions = lambda positions, i, read_end: {}
        g.cluster_terminal_positions = lambda extra, read_end, cutoff: {1005: 4}
        g._attach_side([intron], {intron: {}}, {intron: {1005: 4}}, read_end=True)
        verts = g.outgoing_edges[intron]
        assert (TerminalVertex.read_end, 1000) in verts
        assert (TerminalVertex.read_end, 1010) not in verts
        assert (TerminalVertex.read_end, 1005) not in verts

    def test_end_processor_alt_ends_and_real_l1_gbdt_integration(self, monkeypatch):
        import isoquant_next.discovery.end_processor as ep_mod

        gi = SimpleNamespace(chr_id="chr1", all_isoforms_exons={"REF1": [(100, 200), (300, 400)]})
        args = _make_discovery_args(apa_delta=10, min_novel_count=2, terminal_position_rel=0.1)
        ctx = SimpleNamespace(
            gene_info=gi,
            args=args,
            store=ModelStore(gi),
            use_tss_model=True,
            novel_apa=False,
            get_transcript_id=SimpleIDDistributor().increment,
        )
        ep = TranscriptEndProcessor(ctx)

        m = TranscriptModel("chr1", "+", "T1", "G1", [(100, 200), (300, 400), (500, 600)], TranscriptModelType.known)
        assert TranscriptEndProcessor._intron_chain_key(m) == ((201, 299), (401, 499))
        assert TranscriptEndProcessor._closest_inward([10, 20, 30], 15, True) == 20
        assert TranscriptEndProcessor._closest_inward([30, 20, 10], 25, False) == 20

        # Real L1 GBDT model (no XGBoost dependency) derives alt-end NIC for clustered reads
        source = TranscriptModel("chr1", "+", "ENST1", "G1", [(100, 200), (300, 400)], TranscriptModelType.known)
        reads = [
            SimpleNamespace(
                read_id=f"r{i}",
                strand="+",
                polyA_found=True,
                polya_info=PolyAInfo(external_polya_pos=500 + (i % 3)),
                corrected_exons=[(100, 200), (300, 500 + (i % 3))],
            )
            for i in range(30)
        ]
        # Verify real L1 GBDT derives an alt-end NIC around 500 without stubbing!
        derived_real = ep.derive_alternative_end_models(source, reads)
        assert len(derived_real) == 1
        nic_real, supporting_real = derived_real[0]
        assert 499 <= nic_real.exon_blocks[-1][1] <= 503
        assert len(supporting_real) == 30
        assert nic_real.transcript_id.endswith(".nic")

        # Duplicate filtering check
        known = TranscriptModel("chr1", "+", "REF1", "G1", [(100, 200), (300, 400)], TranscriptModelType.known)
        dup = TranscriptModel("chr1", "+", "DUP", "G1", [(105, 200), (300, 398)], TranscriptModelType.novel_in_catalog)
        keep = TranscriptModel("chr1", "+", "KEEP", "G1", [(100, 200), (300, 800)], TranscriptModelType.novel_in_catalog)
        ep.store.transcript_model_storage = [known, dup, keep]
        ep._drop_duplicate_alt_end_models()
        assert {x.transcript_id for x in ep.store.transcript_model_storage} == {"REF1", "KEEP"}


class TestSpliceSiteCorrectorOracle:
    def test_splice_site_corrector_vectors(self):
        from isoquant_lib.model_construction.transcript_splice_site_corrector import (
            extract_relevant_deletions as lib_rel_dels,
        )

        cigar = [(0, 50), (2, 2), (0, 40), (2, 4), (0, 6)]
        assert extract_relevant_deletions(cigar, 100, [(100, 200)]) == [(192, 4)] == lib_rel_dels(cigar, 100, [(100, 200)])
        assert extract_splice_site_locations_within_aligned_read(20, 40, [(1, 10), (20, 30), (40, 50)]) == [
            (20, False),
            (30, True),
            (40, False),
        ]
        assert compute_most_common_case_of_deletions({0: 10, 4: 20}, False) == 4
        assert compute_most_common_case_of_deletions({0: 10, 4: 20}, True) == -4
        assert compute_most_common_case_of_deletions({3: 20, 4: 20}, False) == -1

        reads = [SimpleNamespace(exons=[(0, 30)], boundary_deletions=[(10, 4)]) for _ in range(5)]
        mf = ModelFilter.__new__(ModelFilter)
        mf.chr_record = "C" * 10 + "GT" + "C" * 20
        mf.args = SimpleNamespace(no_splice_site_correction=False)
        mf.store = SimpleNamespace()
        model = SimpleNamespace(transcript_id="t1", transcript_type=TranscriptModelType.novel_not_in_catalog, strand="+", exon_blocks=[(0, 14), (20, 30)])
        mf.store.transcript_model_storage = [model]
        mf.store.transcript_read_ids = {"t1": reads}
        mf.correct_transcript_splice_sites()
        assert model.exon_blocks == [(0, 10), (20, 30)]

        # Known models remain untouched
        known_m = SimpleNamespace(transcript_id="k1", transcript_type=TranscriptModelType.known, strand="+", exon_blocks=[(0, 14), (20, 30)])
        mf.store.transcript_model_storage = [known_m]
        mf.store.transcript_read_ids = {"k1": reads}
        mf.correct_transcript_splice_sites()
        assert known_m.exon_blocks == [(0, 14), (20, 30)]


class TestEndToEndDiscoveryChr9:
    def test_five_stage_discovery_and_context_gff_on_chr9(self, tmp_path):
        gtf_path = TESTS_DATA / "simple_data" / "chr9.4M.gtf.gz"
        bam_path = TESTS_DATA / "simple_data" / "chr9.4M.ont.sim.polya.bam"
        fa_path = TESTS_DATA / "simple_data" / "chr9.4M.fa.gz"

        # Copy fasta or use ReferenceGenome (write clean temp copy if .fai needed — .fai already exists read-only)
        args = _make_discovery_args()
        sp = StringPoolManager()
        id_dist = SimpleIDDistributor()
        exon_ids = FeatureIdStorage(SimpleIDDistributor())
        chr_state = ChromosomeDiscoveryState()

        plain_fa = tmp_path / "chr9.4M.fa"
        with gzip.open(fa_path, "rb") as src, open(plain_fa, "wb") as dst:
            shutil.copyfileobj(src, dst)

        iqidx_path = compile_gtf_to_iqidx(gtf_path, tmp_path / "chr9.iqidx")
        with MmapAnnotationIndex(iqidx_path) as iqidx, ReferenceGenome(plain_fa) as ref:
            chr_len = ref.get_contig_lengths()["chr9"]
            chr_record = ref.fetch_1based("chr9", 1, chr_len)
            lazy_view = ref.contig_view("chr9")
            assert len(lazy_view) == chr_len
            assert lazy_view[0:25] == chr_record[0:25]
            assert lazy_view[3000250:3000280] == chr_record[3000250:3000280]
            assert lazy_view[-10:] == chr_record[-10:]
            clusters = list(iqidx.iter_gene_clusters("chr9"))
            assert clusters

            # Gather real reads from chr9.4M BAM
            pf = PolyAFinder()
            read_infos = []
            with pysam.AlignmentFile(str(bam_path), "rb") as bf:
                for aln in bf.fetch("chr9"):
                    if aln.is_unmapped or aln.is_secondary or aln.is_supplementary:
                        continue
                    ai = AlignmentInfo(aln)
                    if not ai.read_exons:
                        continue
                    pa = pf.detect_polya(aln)
                    read_infos.append((aln, ai.read_exons, pa))
            assert len(read_infos) > 50

            # Run 5-stage discovery across all clusters sharing ChromosomeDiscoveryState
            all_built_models = []
            tx_counter = _DummyCounter()
            gene_counter = _DummyCounter()
            with GFFPrinter(str(tmp_path), "chr9_disc", exon_ids, check_canonical=True, header="# test\n") as printer:
                for clust in clusters:
                    idxs = list(clust.gene_indices)
                    gi = GeneInfo.from_iqidx_cluster(iqidx, idxs, delta=args.delta)
                    assigner = LongReadAssigner(gi, args, sp)
                    pc = CombinedProfileConstructor(gi, args)
                    cluster_ras = []
                    for aln, exons, pa in read_infos:
                        if exons[-1][1] < gi.start or exons[0][0] > gi.end:
                            continue
                        cp = pc.construct_profiles(exons, pa, [])
                        ra = assigner.assign_to_isoform(aln.query_name, cp)
                        ra.chr_id = "chr9"
                        ra.strand = "-" if aln.is_reverse else "+"
                        ra.exons = exons
                        ra.corrected_exons = exons
                        ra.corrected_introns = junctions_from_blocks(exons)
                        ra.polyA_found = pa.external_polya_pos != -1 or pa.external_polyt_pos != -1
                        ra.polya_info = pa
                        ra.mapping_quality = int(aln.mapping_quality)
                        cluster_ras.append(ra)
                    if not cluster_ras:
                        continue
                    builder = GraphBasedModelConstructor(
                        gi,
                        chr_record,
                        args,
                        tx_counter,
                        gene_counter,
                        id_dist,
                        string_pools=sp,
                        chr_state=chr_state,
                    )
                    builder.process(cluster_ras)
                    printer.dump(gi, builder.transcript_model_storage)
                    all_built_models.extend(builder.transcript_model_storage)

            assert len(all_built_models) >= 1
            assert len(tx_counter.assignments) > 0
            gtf_out = (tmp_path / "chr9_disc.transcript_models.gtf").read_text()
            assert "\ttranscript\t" in gtf_out and "\texon\t" in gtf_out

            # Verify create_extended_storage directly from .iqidx without SQLite
            ext_models, ext_gi = create_extended_storage(iqidx, "chr9", chr_record, [])
            assert len(ext_models) == 54  # all 54 known transcripts on chr9.4M!
            assert ext_gi.chr_id == "chr9"


class TestIntronCollectorAndGraphConstruction:
    def test_intron_collector_clustering_and_similar_map(self):
        from isoquant_next.discovery import IntronCollector

        gi = SimpleNamespace(intron_profiles=SimpleNamespace(features=[(100, 200), (400, 500)]))
        collector = IntronCollector(gi, delta=3)

        def _ra(introns, multimapper=False):
            return SimpleNamespace(
                corrected_introns=introns,
                multimapper=multimapper,
            )

        reads = (
            [_ra([(100, 200), (400, 500)]) for _ in range(5)]
            + [_ra([(102, 199), (400, 500)]) for _ in range(2)]
            + [_ra([(700, 800)]) for _ in range(3)]
            + [_ra([(1100, 1200)])]
            + [_ra([(1500, 1600)], multimapper=True) for _ in range(4)]
        )
        all_introns = collector.collect_introns(reads)
        assert all_introns[(1500, 1600)] == 0
        assert all_introns[(100, 200)] == 5
        assert all_introns[(102, 199)] == 2

        sim_map = collector.construct_similar_intron_map(all_introns)
        assert (102, 199) in sim_map[(100, 200)]
        assert (100, 200) in sim_map[(102, 199)]

        collector.process(reads, min_count=2)
        assert collector.clustered_introns[(100, 200)] == 7
        assert collector.intron_correction_map[(102, 199)] == (100, 200)
        assert collector.clustered_introns[(700, 800)] == 3
        assert (1100, 1200) in collector.discarded_introns
        assert (1100, 1200) not in collector.clustered_introns

    def test_intron_graph_multi_isoform_simplify_pruning_and_polya(self):
        i1 = (101, 199)
        i2a = (301, 399)
        i2b = (501, 599)
        i2a_bulge = (303, 398)  # distance (2, 1) < 6 of i2a; cov=2 < 10 * 0.3=3 -> collapsed in clean_tips_and_bulges
        i_iso_noise = (1501, 1599)  # isolated novel intron with cov=2 < min_novel_isolated_intron_abs=3

        gi = SimpleNamespace(
            intron_profiles=SimpleNamespace(features=[i1, i2a, i2b]),
            all_isoforms_introns={"ISO_A": [i1, i2a], "ISO_B": [i1, i2b]},
            all_isoforms_exons={
                "ISO_A": [(50, 100), (200, 300), (400, 900)],
                "ISO_B": [(50, 100), (200, 500), (600, 800)],
            },
        )
        params = _make_discovery_args(
            delta=0,
            apa_delta=10,
            min_novel_intron_count=2,
            min_novel_isolated_intron_abs=3,
            graph_clustering_distance=6,
            graph_clustering_ratio=0.3,
            singleton_adjacent_cov=5,
            terminal_position_abs=1,
            terminal_position_rel=0.1,
            terminal_internal_position_rel=0.05,
        )

        def _ra(exons, strand="+", polya_pos=-1):
            return SimpleNamespace(
                corrected_exons=exons,
                corrected_introns=junctions_from_blocks(exons),
                multimapper=False,
                strand=strand,
                polya_info=PolyAInfo(polya_pos, -1, -1, -1),
            )

        reads = (
            [_ra([(50, 100), (200, 300), (400, 902)], polya_pos=902) for _ in range(10)]
            + [_ra([(50, 100), (200, 500), (600, 750)], polya_pos=-1) for _ in range(8)]
            + [_ra([(50, 100), (200, 302), (399, 450)], polya_pos=-1) for _ in range(2)]  # bulge i2a_bulge
            + [_ra([(1450, 1500), (1600, 1650)], polya_pos=-1) for _ in range(2)]  # isolated noise i_iso_noise
            + [_ra([(50, 100), (200, 300), (400, 420), (1800, 1850)], polya_pos=-1)]  # count=1 pruned before construct
        )

        graph = IntronGraph(params, gi, reads, polya_predictions=[900], tss_predictions=None)
        assert graph.get_outgoing(i1) == [i2a, i2b]
        assert graph.get_incoming(i2a) == [i1]
        assert graph.get_incoming(i2b) == [i1]
        assert i2a_bulge not in graph.intron_collector.clustered_introns
        assert graph.intron_collector.substitute(i2a_bulge) == i2a
        assert i_iso_noise in graph.intron_collector.discarded_introns
        assert (TerminalVertex.polya, 900) in graph.outgoing_edges[i2a]


def test_canonical_strandness_reporting_level_identity() -> None:
    from isoquant_next.config.enums import StrandnessReportingLevel as ConfigSRL
    from isoquant_next.discovery.state import StrandnessReportingLevel as StateSRL
    from isoquant_next.discovery import StrandnessReportingLevel as InitSRL

    assert StateSRL is ConfigSRL
    assert InitSRL is ConfigSRL
    assert ConfigSRL.only_canonical == StateSRL.only_canonical

