#!/usr/bin/env python3
"""Build and validate the full-cohort Slurm Job Array manifest across all datasets in backup/."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Dict, List, Optional


REPO_ROOT = Path("/gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant")
BENCH_ROOT = REPO_ROOT / "isoquant_refactored"
DATA_ROOT = Path("/gpfs/group/data_sci_group/weiwu_group/hongyanhong/backup/longpass_benchmark250209")
GENOME_DIR = DATA_ROOT / "genome"

V39_GTF = GENOME_DIR / "lrgasp_gencode_v39_annotation_sirvs.human.gtf"
V46_GTF = GENOME_DIR / "gencode.v46.annotation.gtf"
VM28_GTF = GENOME_DIR / "lrgasp_gencode_vM28_sirvs.mouse.gtf"
VM36_GTF = DATA_ROOT / "cDNA_pacbio_fulpass/mouse/F121_cell/genome/gencode.vM36.annotation.gtf"
SEQUIN_GTF_GZ = GENOME_DIR / "lrgasp_sirvs4.gtf.gz"

GRCH38_FA = GENOME_DIR / "GRCh38.p14.genome.fa"
GRCM39_FA = GENOME_DIR / "GRCm39.genome.fa"
MANATEE_FA = DATA_ROOT / "cDNA_pacbio_fulpass/Trichechus_manatus/genome/GCA_030013775.1_ASM3001377v1_genomic.fa"
SEQUIN_FA = DATA_ROOT / "sequin_benchmark/genome/lrgasp_sirv4.fa"
CAGE_BED = GENOME_DIR / "CAGEPAS/CAGE_WTC11.all_reps_v2.bed"


def _make_task(
    *,
    slug: str,
    gate: str,
    protocol: str,
    species: str,
    cell_line: str,
    sample_id: str,
    mode: str = "pair",
    data_type: str,
    fl_data: bool = True,
    stranded: Optional[str] = None,
    threads: int = 16,
    repeats: int = 1,
    bam: Optional[Path] = None,
    fastq: Optional[Path] = None,
    gtf: Optional[Path] = None,
    ref: Optional[Path] = None,
    illumina_bam: Optional[Path] = None,
    cage_bed: Optional[Path] = None,
    yaml_entries: Optional[List[dict]] = None,
    eval_sirv: bool = False,
    complete_genedb: bool = True,
    notes: str = "",
) -> dict:
    # Validate inputs on disk immediately
    for p in (bam, fastq, gtf, ref, illumina_bam, cage_bed):
        if p is not None and not p.exists():
            raise FileNotFoundError(f"[{slug}] Missing required input path: {p}")
    if yaml_entries:
        for ent in yaml_entries:
            for b in ent["bams"]:
                if not Path(b).exists():
                    raise FileNotFoundError(f"[{slug}] Missing YAML BAM path: {b}")

    needs_bai_link = False
    if bam is not None:
        has_bai = Path(f"{bam}.bai").exists() or bam.with_suffix(".bai").exists()
        needs_bai_link = not has_bai

    size_gib = 0.0
    if bam is not None:
        size_gib = round(bam.stat().st_size / (1024**3), 3)
    elif fastq is not None:
        size_gib = round(fastq.stat().st_size / (1024**3), 3)
    elif yaml_entries:
        size_gib = round(
            sum(Path(b).stat().st_size for ent in yaml_entries for b in ent["bams"]) / (1024**3), 3
        )

    queue_tier = "heavy" if threads >= 32 else "standard"
    return {
        "slug": slug,
        "gate": gate,
        "protocol": protocol,
        "species": species,
        "cell_line": cell_line,
        "sample_id": sample_id,
        "mode": mode,
        "data_type": data_type,
        "fl_data": fl_data,
        "stranded": stranded,
        "threads": threads,
        "repeats": repeats,
        "queue_tier": queue_tier,
        "size_gib": size_gib,
        "needs_bai_link": needs_bai_link,
        "bam": str(bam) if bam else None,
        "fastq": str(fastq) if fastq else None,
        "gtf": str(gtf) if gtf else None,
        "ref": str(ref) if ref else None,
        "illumina_bam": str(illumina_bam) if illumina_bam else None,
        "cage_bed": str(cage_bed) if cage_bed else None,
        "yaml_entries": yaml_entries,
        "eval_sirv": eval_sirv,
        "complete_genedb": complete_genedb,
        "notes": notes,
    }


def _find_single_bam(dir_path: Path, sample_id: str) -> Path:
    candidates = [
        p
        for p in sorted(dir_path.glob(f"{sample_id}*.bam"))
        if "nosupplementary" not in p.name and "reads_discarded" not in p.name
    ]
    if len(candidates) != 1:
        raise FileNotFoundError(f"Expected 1 BAM matching {sample_id} in {dir_path}, found {candidates}")
    return candidates[0]


def build_manifest() -> List[dict]:
    tasks: List[dict] = []

    # ------------------------------------------------------------------
    # Gate G0/G1: Thread scaling & warm cache tests (preserve/include existing slugs)
    # ------------------------------------------------------------------
    bqq = DATA_ROOT / "cDNA_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF400BQQ.sorted.k14.bam"
    tasks.append(
        _make_task(
            slug="g1_encff400bqq_t1",
            gate="G0_G1",
            protocol="cDNA_PacBio",
            species="Human",
            cell_line="H1_ES_cell",
            sample_id="ENCFF400BQQ",
            data_type="pacbio_ccs",
            threads=1,
            bam=bqq,
            gtf=V39_GTF,
            ref=GRCH38_FA,
            notes="1-thread scaling baseline",
        )
    )
    tasks.append(
        _make_task(
            slug="g1_encff400bqq_t8_abba3",
            gate="G0_G1",
            protocol="cDNA_PacBio",
            species="Human",
            cell_line="H1_ES_cell",
            sample_id="ENCFF400BQQ",
            data_type="pacbio_ccs",
            threads=8,
            repeats=3,
            bam=bqq,
            gtf=V39_GTF,
            ref=GRCH38_FA,
            notes="3-repeat A-B-B-A cold GTF",
        )
    )
    tasks.append(
        _make_task(
            slug="g1_encff400bqq_t8_warm_abba3",
            gate="G0_G1",
            protocol="cDNA_PacBio",
            species="Human",
            cell_line="H1_ES_cell",
            sample_id="ENCFF400BQQ",
            mode="warm_idx",
            data_type="pacbio_ccs",
            threads=8,
            repeats=3,
            bam=bqq,
            gtf=V39_GTF,
            ref=GRCH38_FA,
            notes="3-repeat A-B-B-A warm precompiled .db vs .iqidx",
        )
    )

    # ------------------------------------------------------------------
    # Gate G1: cDNA_pacbio_fulpass (all 14 replicates)
    # ------------------------------------------------------------------
    cdna_pb_groups = [
        ("Human", "H1_DE_cell", DATA_ROOT / "cDNA_pacbio_fulpass/human/H1_DE_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF235QXW", "ENCFF712CBL", "ENCFF731HST"]),
        ("Human", "H1_ES_cell", DATA_ROOT / "cDNA_pacbio_fulpass/human/H1_ES_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF400BQQ", "ENCFF436GKZ", "ENCFF853OFP"]),
        ("Human", "WTC11_cell", DATA_ROOT / "cDNA_pacbio_fulpass/human/WTC11_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF245IPA", "ENCFF370NFS", "ENCFF563QZR"]),
        ("Mouse", "F121_cell", DATA_ROOT / "cDNA_pacbio_fulpass/mouse/F121_cell/alignments", VM28_GTF, GRCM39_FA, ["ENCFF313VYZ", "ENCFF667VXS", "ENCFF874VSI"]),
        ("Manatee", "Trichechus_manatus", DATA_ROOT / "cDNA_pacbio_fulpass/Trichechus_manatus/alignments", None, MANATEE_FA, ["ENCFF589LDN", "ENCFF835YMN"]),
    ]
    existing_slugs_map = {
        ("cDNA_PacBio", "ENCFF400BQQ"): "g1_encff400bqq_t16",
        ("cDNA_PacBio", "ENCFF436GKZ"): "g2_encff436gkz_t16",
        ("capped_PacBio", "ENCFF684YOO"): "g1_encff684yoo_t16",
        ("capped_PacBio", "ENCFF339FMQ"): "g1_encff339fmq_t16",
        ("directRNA", "ENCFF846RID"): "g2_encff846rid_drna_t16",
        ("directRNA", "ENCFF207KVO"): "g2_encff207kvo_drna_t16",
        ("RBMX_KD", "con_rep1"): "g2_con_rep1_v46_t16",
        ("RBMX_KD", "shRBMX_rep2"): "g3_shrbmx_rep2_v46_t32",
        ("capped_ONT", "ENCFF089AMF"): "g3_encff089amf_t32",
        ("R2C2_ONT", "ENCFF598YQO"): "g3_encff598yqo_mouse_t16",
        ("SIRV_simulation", "simulated1"): "g3_simulated1_t16",
    }

    for sp, cl, adir, gtf, ref, sids in cdna_pb_groups:
        for sid in sids:
            bam = _find_single_bam(adir, sid)
            slug = existing_slugs_map.get(("cDNA_PacBio", sid), f"g1_cdna_pb_{sp.lower()}_{cl.lower()}_{sid.lower()}_t16")
            tasks.append(
                _make_task(
                    slug=slug,
                    gate="G1",
                    protocol="cDNA_PacBio",
                    species=sp,
                    cell_line=cl,
                    sample_id=sid,
                    data_type="pacbio_ccs",
                    threads=16,
                    bam=bam,
                    gtf=gtf,
                    ref=ref,
                    notes="De novo no-genedb" if gtf is None else "",
                )
            )

    # Extra reference-version contrasts for G1
    tasks.append(
        _make_task(
            slug="g1_encff370nfs_wtc11_v46_t16",
            gate="G1",
            protocol="cDNA_PacBio_v46",
            species="Human",
            cell_line="WTC11_cell",
            sample_id="ENCFF370NFS",
            data_type="pacbio_ccs",
            threads=16,
            bam=_find_single_bam(DATA_ROOT / "cDNA_pacbio_fulpass/human/WTC11_cell/alignments", "ENCFF370NFS"),
            gtf=V46_GTF,
            ref=GRCH38_FA,
            notes="GENCODE v46 reference contrast",
        )
    )
    tasks.append(
        _make_task(
            slug="g1_encff313vyz_mouse_vm36_t16",
            gate="G1",
            protocol="cDNA_PacBio_vM36",
            species="Mouse",
            cell_line="F121_cell",
            sample_id="ENCFF313VYZ",
            data_type="pacbio_ccs",
            threads=16,
            bam=_find_single_bam(DATA_ROOT / "cDNA_pacbio_fulpass/mouse/F121_cell/alignments", "ENCFF313VYZ"),
            gtf=VM36_GTF,
            ref=GRCM39_FA,
            notes="GENCODE vM36 reference contrast",
        )
    )

    # ------------------------------------------------------------------
    # Gate G1: capped_pacbio_fulpass (all 12 replicates)
    # ------------------------------------------------------------------
    cap_pb_groups = [
        ("Human", "H1_DE_cell", DATA_ROOT / "capped_pacbio_fulpass/human/H1_DE_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF054KCY", "ENCFF296KQK", "ENCFF761BFK"]),
        ("Human", "H1_ES_cell", DATA_ROOT / "capped_pacbio_fulpass/human/H1_ES_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF143ICB", "ENCFF339FMQ", "ENCFF684YOO"]),
        ("Human", "WTC11_cell", DATA_ROOT / "capped_pacbio_fulpass/human/WTC11_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF003QZT", "ENCFF105WIJ", "ENCFF212HLP"]),
        ("Mouse", "F121_cell", DATA_ROOT / "capped_pacbio_fulpass/mouse/F121_cell/alignments", VM28_GTF, GRCM39_FA, ["ENCFF310IPO", "ENCFF535DQA", "ENCFF654JHQ"]),
    ]
    for sp, cl, adir, gtf, ref, sids in cap_pb_groups:
        for sid in sids:
            bam = _find_single_bam(adir, sid)
            slug = existing_slugs_map.get(("capped_PacBio", sid), f"g1_cap_pb_{sp.lower()}_{cl.lower()}_{sid.lower()}_t16")
            tasks.append(
                _make_task(
                    slug=slug,
                    gate="G1",
                    protocol="capped_PacBio",
                    species=sp,
                    cell_line=cl,
                    sample_id=sid,
                    data_type="pacbio_ccs",
                    threads=16,
                    bam=bam,
                    gtf=gtf,
                    ref=ref,
                )
            )

    # Multi-sample YAML merges in G1
    tasks.append(
        _make_task(
            slug="g1_multisample_yaml_t16",
            gate="G1",
            protocol="MultiSample_YAML",
            species="Human",
            cell_line="H1_ES_cell",
            sample_id="ENCFF684YOO+ENCFF339FMQ",
            mode="yaml",
            data_type="pacbio_ccs",
            threads=16,
            gtf=V39_GTF,
            ref=GRCH38_FA,
            yaml_entries=[
                {"name": "ENCFF684YOO", "bams": [str(DATA_ROOT / "capped_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF684YOO.bam")]},
                {"name": "ENCFF339FMQ", "bams": [str(DATA_ROOT / "capped_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF339FMQ.bam")]},
            ],
            notes="2-rep H1_ES YAML combine_counts",
        )
    )
    tasks.append(
        _make_task(
            slug="g1_cross_cell_yaml_t16",
            gate="G1",
            protocol="MultiSample_YAML",
            species="Human",
            cell_line="CrossCell_H1DE_H1ES_WTC11",
            sample_id="ENCFF731HST+ENCFF400BQQ+ENCFF370NFS",
            mode="yaml",
            data_type="pacbio_ccs",
            threads=16,
            gtf=V39_GTF,
            ref=GRCH38_FA,
            yaml_entries=[
                {"name": "H1DE_ENCFF731HST", "bams": [str(DATA_ROOT / "cDNA_pacbio_fulpass/human/H1_DE_cell/alignments/ENCFF731HST.sorted.k14.bam")]},
                {"name": "H1ES_ENCFF400BQQ", "bams": [str(DATA_ROOT / "cDNA_pacbio_fulpass/human/H1_ES_cell/alignments/ENCFF400BQQ.sorted.k14.bam")]},
                {"name": "WTC11_ENCFF370NFS", "bams": [str(DATA_ROOT / "cDNA_pacbio_fulpass/human/WTC11_cell/alignments/ENCFF370NFS.sorted.k14.bam")]},
            ],
            notes="Cross-cell-line PacBio YAML merge",
        )
    )

    # ------------------------------------------------------------------
    # Gate G2: directRNA_fulpass (all 12 replicates)
    # ------------------------------------------------------------------
    drna_groups = [
        ("Human", "H1_DE_cell", DATA_ROOT / "directRNA_fulpass/human/H1_DE_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF008HMY", "ENCFF937WUA", "ENCFF965CGM"]),
        ("Human", "H1_ES_cell", DATA_ROOT / "directRNA_fulpass/human/H1_ES_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF207KVO", "ENCFF373LMV", "ENCFF846RID"]),
        ("Human", "WTC11_cell", DATA_ROOT / "directRNA_fulpass/human/WTC11_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF155CFF", "ENCFF600LIU", "ENCFF771DIX"]),
        ("Mouse", "F121_cell", DATA_ROOT / "directRNA_fulpass/mouse/F121_cell/alignments", VM28_GTF, GRCM39_FA, ["ENCFF349BIN", "ENCFF412NKJ", "ENCFF765AEC"]),
    ]
    for sp, cl, adir, gtf, ref, sids in drna_groups:
        for sid in sids:
            bam = _find_single_bam(adir, sid)
            slug = existing_slugs_map.get(("directRNA", sid), f"g2_drna_{sp.lower()}_{cl.lower()}_{sid.lower()}_t16")
            tasks.append(
                _make_task(
                    slug=slug,
                    gate="G2",
                    protocol="directRNA",
                    species=sp,
                    cell_line=cl,
                    sample_id=sid,
                    data_type="nanopore",
                    stranded="forward",
                    fl_data=True,
                    threads=16,
                    bam=bam,
                    gtf=gtf,
                    ref=ref,
                )
            )

    # ------------------------------------------------------------------
    # Gate G2: cDNA_ONT_fulpass (all 15 replicates)
    # ------------------------------------------------------------------
    cdna_ont_groups = [
        ("Human", "H1_DE_cell", DATA_ROOT / "cDNA_ONT_fulpass/human/H1_DE_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF180OMN", "ENCFF277PKW", "ENCFF754BZS"]),
        ("Human", "H1_ES_cell", DATA_ROOT / "cDNA_ONT_fulpass/human/H1_ES_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF326IBB", "ENCFF714YOZ", "ENCFF955CIA"]),
        ("Human", "WTC11_cell", DATA_ROOT / "cDNA_ONT_fulpass/human/WTC11_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF023EXJ", "ENCFF263YFG", "ENCFF961HLO"]),
        ("Mouse", "F121_cell", DATA_ROOT / "cDNA_ONT_fulpass/mouse/F121_cell/alignments", VM28_GTF, GRCM39_FA, ["ENCFF232YSU", "ENCFF288PBL", "ENCFF683TBO"]),
        ("Manatee", "Trichechus_manatus", DATA_ROOT / "cDNA_ONT_fulpass/trichechus_manatus/Trichechus_manatus/alignments", None, MANATEE_FA, ["ENCFF090GGQ", "ENCFF620ASK", "ENCFF629AWL"]),
    ]
    for sp, cl, adir, gtf, ref, sids in cdna_ont_groups:
        for sid in sids:
            bam = _find_single_bam(adir, sid)
            sz_gb = bam.stat().st_size / (1024**3)
            threads = 32 if sz_gb >= 4.5 else 16
            tasks.append(
                _make_task(
                    slug=f"g2_cdna_ont_{sp.lower()}_{cl.lower()}_{sid.lower()}_t{threads}",
                    gate="G2",
                    protocol="cDNA_ONT",
                    species=sp,
                    cell_line=cl,
                    sample_id=sid,
                    data_type="nanopore",
                    fl_data=True,
                    threads=threads,
                    bam=bam,
                    gtf=gtf,
                    ref=ref,
                    notes="De novo no-genedb" if gtf is None else "",
                )
            )

    # ------------------------------------------------------------------
    # Gate G2 / G3: RBMX_KD (all 6 replicates + CAGE + YAML)
    # ------------------------------------------------------------------
    rbmx_dir = DATA_ROOT / "RBMX_KD/alignments"
    for sid in ["ama_rep1", "ama_rep2", "con_rep1", "con_rep2", "shRBMX_rep1", "shRBMX_rep2"]:
        bam = rbmx_dir / f"{sid}.fastq.sorted.nosecondary.bam"
        sz_gb = bam.stat().st_size / (1024**3)
        threads = 32 if sz_gb >= 4.5 else 16
        slug = existing_slugs_map.get(("RBMX_KD", sid), f"g2_rbmx_{sid.lower()}_v46_t{threads}")
        gate = "G3" if sid == "shRBMX_rep2" else "G2"
        tasks.append(
            _make_task(
                slug=slug,
                gate=gate,
                protocol="RBMX_KD",
                species="Human",
                cell_line="WTC11_RBMX_KD",
                sample_id=sid,
                data_type="nanopore",
                fl_data=False,
                threads=threads,
                bam=bam,
                gtf=V46_GTF,
                ref=GRCH38_FA,
            )
        )

    # RBMX_KD CAGE-assisted isoquant-next runs
    tasks.append(
        _make_task(
            slug="g2_con_rep1_cage_t16",
            gate="G2",
            protocol="RBMX_KD_CAGE",
            species="Human",
            cell_line="WTC11_RBMX_KD",
            sample_id="con_rep1",
            mode="cage",
            data_type="nanopore",
            fl_data=True,
            threads=16,
            bam=rbmx_dir / "con_rep1.fastq.sorted.nosecondary.bam",
            gtf=V46_GTF,
            ref=GRCH38_FA,
            cage_bed=CAGE_BED,
            notes="isoquant-next CAGE-assisted TSS discovery",
        )
    )
    for sid in ["ama_rep1", "shRBMX_rep1"]:
        tasks.append(
            _make_task(
                slug=f"g2_{sid.lower()}_cage_t16",
                gate="G2",
                protocol="RBMX_KD_CAGE",
                species="Human",
                cell_line="WTC11_RBMX_KD",
                sample_id=sid,
                mode="cage",
                data_type="nanopore",
                fl_data=True,
                threads=16,
                bam=rbmx_dir / f"{sid}.fastq.sorted.nosecondary.bam",
                gtf=V46_GTF,
                ref=GRCH38_FA,
                cage_bed=CAGE_BED,
                notes="isoquant-next CAGE-assisted TSS discovery",
            )
        )

    tasks.append(
        _make_task(
            slug="g2_rbmx_yaml_t32",
            gate="G2",
            protocol="MultiSample_YAML",
            species="Human",
            cell_line="WTC11_RBMX_KD",
            sample_id="con_rep1+ama_rep1+shRBMX_rep1",
            mode="yaml",
            data_type="nanopore",
            fl_data=False,
            threads=32,
            gtf=V46_GTF,
            ref=GRCH38_FA,
            yaml_entries=[
                {"name": "con_rep1", "bams": [str(rbmx_dir / "con_rep1.fastq.sorted.nosecondary.bam")]},
                {"name": "ama_rep1", "bams": [str(rbmx_dir / "ama_rep1.fastq.sorted.nosecondary.bam")]},
                {"name": "shRBMX_rep1", "bams": [str(rbmx_dir / "shRBMX_rep1.fastq.sorted.nosecondary.bam")]},
            ],
            notes="Multi-condition RBMX_KD YAML + gencode.v46",
        )
    )

    # ------------------------------------------------------------------
    # Gate G3: capped_ONT_fulpass (all 12 replicates)
    # ------------------------------------------------------------------
    cap_ont_groups = [
        ("Human", "H1_DE_cell", DATA_ROOT / "capped_ONT_fulpass/human/H1_DE_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF430MWW", "ENCFF626JLU", "ENCFF933LQH"]),
        ("Human", "H1_ES_cell", DATA_ROOT / "capped_ONT_fulpass/human/H1_ES_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF085EVG", "ENCFF089AMF", "ENCFF661WLV"]),
        ("Human", "WTC11_cell", DATA_ROOT / "capped_ONT_fulpass/human/WTC11_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF104BNW", "ENCFF654SNK", "ENCFF934KDM"]),
        ("Mouse", "F121_cell", DATA_ROOT / "capped_ONT_fulpass/mouse/F121_cell/alignments", VM28_GTF, GRCM39_FA, ["ENCFF056EOI", "ENCFF275RMO", "ENCFF356OJC"]),
    ]
    for sp, cl, adir, gtf, ref, sids in cap_ont_groups:
        for sid in sids:
            bam = _find_single_bam(adir, sid)
            sz_gb = bam.stat().st_size / (1024**3)
            threads = 32 if sz_gb >= 4.5 else 16
            slug = existing_slugs_map.get(("capped_ONT", sid), f"g3_cap_ont_{sp.lower()}_{cl.lower()}_{sid.lower()}_t{threads}")
            tasks.append(
                _make_task(
                    slug=slug,
                    gate="G3",
                    protocol="capped_ONT",
                    species=sp,
                    cell_line=cl,
                    sample_id=sid,
                    data_type="nanopore",
                    fl_data=True,
                    threads=threads,
                    bam=bam,
                    gtf=gtf,
                    ref=ref,
                )
            )

    # ------------------------------------------------------------------
    # Gate G3: R2C2_ONT_fulpass (all 24 primary + 12 seqkitrmdup replicates)
    # ------------------------------------------------------------------
    r2c2_groups = [
        ("R2C2_ONT", "Human", "H1_DE_cell", DATA_ROOT / "R2C2_ONT_fulpass/human/H1_DE_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF083OEF", "ENCFF498ZHI", "ENCFF605OXP", "ENCFF887KHY", "ENCFF930WML", "ENCFF988MFC"]),
        ("R2C2_ONT", "Human", "H1_ES_cell", DATA_ROOT / "R2C2_ONT_fulpass/human/H1_ES_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF095SQG", "ENCFF243MUU", "ENCFF337GCH", "ENCFF367JPM", "ENCFF387KPE", "ENCFF531FUR"]),
        ("R2C2_ONT", "Human", "WTC11_cell", DATA_ROOT / "R2C2_ONT_fulpass/human/WTC11_cell/alignments", V39_GTF, GRCH38_FA, ["ENCFF063ASB", "ENCFF178BYM", "ENCFF357UQD", "ENCFF408XXR", "ENCFF542VPN", "ENCFF679LUJ"]),
        ("R2C2_ONT", "Mouse", "F121_cell", DATA_ROOT / "R2C2_ONT_fulpass/mouse/F121_cell/alignments", VM28_GTF, GRCM39_FA, ["ENCFF055REA", "ENCFF335WMC", "ENCFF598YQO", "ENCFF717LLT", "ENCFF742GCO", "ENCFF797PJV"]),
        ("R2C2_ONT_rmdup", "Human", "H1_DE_cell_seqkitrmdup", DATA_ROOT / "R2C2_ONT_fulpass/human/H1_DE_cell_seqkitrmdup/alignments", V39_GTF, GRCH38_FA, ["ENCFF083OEF", "ENCFF498ZHI", "ENCFF605OXP", "ENCFF887KHY", "ENCFF930WML", "ENCFF988MFC"]),
        ("R2C2_ONT_rmdup", "Human", "H1_ES_cell_seqkitrmdup", DATA_ROOT / "R2C2_ONT_fulpass/human/H1_ES_cell_seqkitrmdup/alignments", V39_GTF, GRCH38_FA, ["ENCFF095SQG", "ENCFF243MUU", "ENCFF337GCH", "ENCFF367JPM", "ENCFF387KPE", "ENCFF531FUR"]),
    ]
    for proto, sp, cl, adir, gtf, ref, sids in r2c2_groups:
        for sid in sids:
            bam = _find_single_bam(adir, sid)
            sz_gb = bam.stat().st_size / (1024**3)
            threads = 32 if sz_gb >= 4.5 else 16
            tag = "r2c2_rmdup" if proto == "R2C2_ONT_rmdup" else "r2c2"
            slug = existing_slugs_map.get((proto, sid), f"g3_{tag}_{sp.lower()}_{cl.lower()}_{sid.lower()}_t{threads}")
            tasks.append(
                _make_task(
                    slug=slug,
                    gate="G3",
                    protocol=proto,
                    species=sp,
                    cell_line=cl,
                    sample_id=sid,
                    data_type="nanopore",
                    fl_data=True,
                    threads=threads,
                    bam=bam,
                    gtf=gtf,
                    ref=ref,
                )
            )

    # ------------------------------------------------------------------
    # Gate G3: SIRV_simulation (all 3 replicates) + sequin_benchmark
    # ------------------------------------------------------------------
    for i in (1, 2, 3):
        sid = f"simulated{i}"
        slug = f"g3_{sid}_t16"
        tasks.append(
            _make_task(
                slug=slug,
                gate="G3",
                protocol="SIRV_simulation",
                species="Human_SIRV",
                cell_line="NanoSim",
                sample_id=sid,
                data_type="nanopore",
                fl_data=False,
                threads=16,
                bam=DATA_ROOT / f"SIRV_simulation/alignments/{sid}_aligned_reads.sorted.bam",
                gtf=V46_GTF,
                ref=GRCH38_FA,
                eval_sirv=True,
                notes="NanoSim ground-truth Recall/Precision/F1",
            )
        )

    tasks.append(
        _make_task(
            slug="g3_sequin_benchmark_t16",
            gate="G3",
            protocol="Sequin_SIRV4",
            species="Sequin_SIRV4",
            cell_line="SpikeIn_Combined",
            sample_id="ont_combined_no_supp",
            data_type="nanopore",
            fl_data=False,
            threads=16,
            bam=DATA_ROOT / "sequin_benchmark/combined/ont_combined_no_supp.bam",
            gtf=SEQUIN_GTF_GZ,
            ref=SEQUIN_FA,
            complete_genedb=False,
            notes="Spike-in Sequin+ERCC+SIRV4 benchmark (exon-only GTF)",
        )
    )

    # ------------------------------------------------------------------
    # Gate G4: FASTQ de novo + Illumina hybrid + Resume/Visualize verification
    # ------------------------------------------------------------------
    tasks.append(
        _make_task(
            slug="g4_encff008hmy_fastq_t16",
            gate="G4",
            protocol="FASTQ_directRNA",
            species="Human",
            cell_line="H1_DE_cell",
            sample_id="ENCFF008HMY",
            mode="fastq",
            data_type="nanopore",
            stranded="forward",
            fl_data=True,
            threads=16,
            fastq=DATA_ROOT / "ENCFF008HMY.fastq",
            gtf=V39_GTF,
            ref=GRCH38_FA,
            notes="FASTQ from-scratch minimap2 alignment",
        )
    )
    simple_dir = REPO_ROOT / "isoquant_tests/simple_data"
    tasks.append(
        _make_task(
            slug="g4_chr9_illumina_t8",
            gate="G4",
            protocol="Illumina_Hybrid",
            species="Human_chr9",
            cell_line="Sim_chr9_4M",
            sample_id="ONT_Simulated_chr9",
            mode="illumina_pair",
            data_type="nanopore",
            fl_data=False,
            threads=8,
            bam=simple_dir / "chr9.4M.ont.sim.polya.bam",
            illumina_bam=simple_dir / "chr9.4M.Illumina.bam",
            gtf=simple_dir / "chr9.4M.gtf.gz",
            ref=simple_dir / "chr9.4M.fa.gz",
            notes="Hybrid short-read Illumina splice correction pair",
        )
    )
    tasks.append(
        _make_task(
            slug="g4_resume_visualize_t16",
            gate="G4",
            protocol="Resume_Visualize",
            species="Human",
            cell_line="H1_ES_cell",
            sample_id="ENCFF400BQQ",
            mode="resume_vis",
            data_type="pacbio_ccs",
            fl_data=True,
            threads=16,
            bam=bqq,
            gtf=V39_GTF,
            ref=GRCH38_FA,
            notes="End-to-end --resume state parity + isoquant-next-visualize plots",
        )
    )

    # Assert unique slugs
    seen: Dict[str, dict] = {}
    for t in tasks:
        if t["slug"] in seen:
            raise ValueError(f"Duplicate task slug: {t['slug']}")
        seen[t["slug"]] = t
    return tasks


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate full-cohort benchmark manifest.")
    parser.add_argument("--tier", choices=["all", "standard", "heavy", "pending_standard", "pending_heavy"], default="all")
    parser.add_argument("--out-json", type=Path, default=BENCH_ROOT / "all_samples_manifest.json")
    parser.add_argument("--out-tsv", type=Path, default=BENCH_ROOT / "all_samples_manifest.tsv")
    args = parser.parse_args()

    all_tasks = build_manifest()
    bench_dir = BENCH_ROOT / "benchmarks"

    def is_done(slug: str) -> bool:
        b = bench_dir / slug / "benchmark_summary.json"
        if not b.is_file():
            return False
        try:
            d = json.loads(b.read_text(encoding="utf-8"))
        except Exception:
            return False
        s = d.get("summary", {})
        if not s:
            return False
        return all(v.get("failures", 1) == 0 for v in s.values())

    for idx, t in enumerate(all_tasks):
        t["global_index"] = idx
        t["already_done"] = is_done(t["slug"])

    if args.tier == "all":
        selected = all_tasks
    elif args.tier == "standard":
        selected = [t for t in all_tasks if t["queue_tier"] == "standard"]
    elif args.tier == "heavy":
        selected = [t for t in all_tasks if t["queue_tier"] == "heavy"]
    elif args.tier == "pending_standard":
        selected = [t for t in all_tasks if t["queue_tier"] == "standard" and not t["already_done"]]
    elif args.tier == "pending_heavy":
        selected = [t for t in all_tasks if t["queue_tier"] == "heavy" and not t["already_done"]]
    else:
        selected = all_tasks

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(selected, indent=2) + "\n", encoding="utf-8")

    cols = ["array_idx", "slug", "gate", "protocol", "species", "cell_line", "sample_id", "mode", "threads", "repeats", "queue_tier", "size_gib", "needs_bai_link", "already_done"]
    with open(args.out_tsv, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for i, t in enumerate(selected):
            w.writerow([i] + [t[c] for c in cols[1:]])

    by_proto: Dict[str, int] = {}
    for t in selected:
        by_proto[t["protocol"]] = by_proto.get(t["protocol"], 0) + 1
    print(f"Manifest written: {args.out_json} ({len(selected)} tasks, done={sum(t['already_done'] for t in selected)}, pending={sum(not t['already_done'] for t in selected)})")
    for k, v in sorted(by_proto.items()):
        print(f"  {k:20s}: {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
