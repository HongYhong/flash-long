"""Context-managed GTF model printer and zero-SQLite extended storage builder."""

from __future__ import annotations

import logging
import os
from collections import defaultdict, namedtuple
from typing import Optional

from isoquant_next.annotation.gene_info import GeneInfo, TranscriptModel
from isoquant_next.annotation.iqidx import MmapAnnotationIndex
from isoquant_next.core.utils import max_range

logger = logging.getLogger("IsoQuant")


def validate_exons(novel_exons) -> bool:
    return novel_exons == sorted(novel_exons) and all(0 < x[0] <= x[1] for x in novel_exons)


class VoidTranscriptPrinter:
    def __enter__(self) -> "VoidTranscriptPrinter":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def dump(self, gene_info, transcript_model_storage=None) -> None:
        pass

    def close(self) -> None:
        pass


class GFFPrinter:
    """Context-managed GTF printer for discovered and reference transcript models."""

    def __init__(
        self,
        outf_prefix: str,
        sample_name: str,
        exon_id_storage,
        gtf_suffix: str = ".transcript_models.gtf",
        check_canonical: bool = False,
        header: str = "",
    ) -> None:
        os.makedirs(outf_prefix, exist_ok=True)
        self.model_fname = os.path.join(outf_prefix, sample_name + gtf_suffix)
        self.exon_id_storage = exon_id_storage
        self.printed_gene_ids: set[str] = set()
        self._closed = False
        self.out_gff = open(self.model_fname, "w", encoding="utf-8")
        if header:
            self.out_gff.write(f"# {sample_name} IsoQuant generated GTF\n{header}")
            self.out_gff.flush()
        self.check_canonical = check_canonical

    def __enter__(self) -> "GFFPrinter":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        self.close()

    def close(self) -> None:
        if not self._closed:
            self._closed = True
            f = getattr(self, "out_gff", None)
            if f is not None and not f.closed:
                f.close()

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def dump(self, gene_info, transcript_model_storage) -> None:
        if not transcript_model_storage:
            return
        gene_to_model = defaultdict(list)
        gene_regions = gene_info.get_gene_regions() if not gene_info.empty() else {}
        GFFGeneInfo = namedtuple("GFFGeneInfo", ("chr_id", "strand", "gene_region"))
        gene_info_dict = {}

        for i, model in enumerate(transcript_model_storage):
            if not validate_exons(model.exon_blocks):
                logger.warning(
                    "Transcript model %s has incorrect coordinates and will be ignored: %s",
                    model.transcript_id,
                    str(model.exon_blocks),
                )
                continue
            gene_id = model.gene_id
            gene_to_model[gene_id].append(i)
            tx_region = (model.exon_blocks[0][0], model.exon_blocks[-1][1])
            if gene_id not in gene_info_dict:
                assert model.chr_id == gene_info.chr_id
                g_range = max_range(gene_regions[gene_id], tx_region) if gene_id in gene_regions else tx_region
                gene_info_dict[gene_id] = GFFGeneInfo(model.chr_id, model.strand, g_range)
            else:
                rec = gene_info_dict[gene_id]
                assert model.chr_id == rec.chr_id
                gene_info_dict[gene_id] = GFFGeneInfo(model.chr_id, rec.strand, max_range(rec.gene_region, tx_region))

        for gene_id, coords in sorted(((g, gene_info_dict[g].gene_region) for g in gene_info_dict), key=lambda x: x[1]):
            if gene_id not in self.printed_gene_ids:
                g_add = gene_info.feature_attributes[gene_id] if gene_info and gene_id in gene_info.feature_attributes else ""
                src = gene_info.sources[gene_id] if gene_info and gene_id in gene_info.sources else "IsoQuant"
                g_rec = gene_info_dict[gene_id]
                self.out_gff.write(
                    f'{g_rec.chr_id}\t{src}\tgene\t{coords[0]}\t{coords[1]}\t.\t{g_rec.strand}\t.\t'
                    f'gene_id "{gene_id}"; transcripts "{len(gene_to_model[gene_id])}"; {g_add}\n'
                )
                self.printed_gene_ids.add(gene_id)

            for m_idx in gene_to_model[gene_id]:
                model = transcript_model_storage[m_idx]
                if not model.check_additional("exons"):
                    model.add_additional_attribute("exons", str(len(model.exon_blocks)))
                t_add = (
                    " " + gene_info.feature_attributes[model.transcript_id]
                    if gene_info and model.transcript_id in gene_info.feature_attributes
                    else ""
                )
                self.out_gff.write(
                    f'{model.chr_id}\t{model.source}\ttranscript\t{model.exon_blocks[0][0]}\t{model.exon_blocks[-1][1]}'
                    f'\t.\t{model.strand}\t.\tgene_id "{model.gene_id}"; transcript_id "{model.transcript_id}"; '
                    f"{model.additional_attributes_str()}{t_add}\n"
                )
                prefix = f"{model.chr_id}\t{model.source}\t"
                suffix = f'.\t{model.strand}\t.\tgene_id "{model.gene_id}"; transcript_id "{model.transcript_id}";'
                exons_to_print = list(model.other_features) + [(e[0], e[1], "exon") for e in model.exon_blocks]
                exons_to_print = sorted(exons_to_print, reverse=(model.strand == "-"))
                for i, e in enumerate(exons_to_print):
                    exon_str_id = self.exon_id_storage.get_id(model.chr_id, e, model.strand)
                    exon_key = f"{model.transcript_id}_{e[0]}_{e[1]}_{model.strand}"
                    e_add = (
                        " " + gene_info.feature_attributes[model.transcript_id]
                        if gene_info and exon_key in gene_info.feature_attributes
                        else ""
                    )
                    self.out_gff.write(
                        f"{prefix}{e[2]}\t{e[0]}\t{e[1]}\t{suffix}"
                        f' exon_number "{i + 1}"; exon_id "{exon_str_id}"; {e_add}\n'
                    )
        self.out_gff.flush()


def create_extended_storage(genedb, chr_id: str, chr_record, novel_model_storage):
    """Build combined known + discovered TranscriptModel storage from MmapAnnotationIndex or gffutils DB."""
    all_models = []
    if genedb is None:
        all_models.extend(novel_model_storage)
        gi = GeneInfo.from_region(chr_id, 1, len(chr_record))
        if chr_record is not None:
            gi.set_reference_sequence(1, len(chr_record), chr_record)
        return all_models, gi

    if isinstance(genedb, MmapAnnotationIndex):
        g0, g1 = genedb.chrom_gene_slice(chr_id)
        if g0 >= g1:
            all_models.extend(novel_model_storage)
            gi = GeneInfo.from_region(chr_id, 1, len(chr_record))
            if chr_record is not None:
                gi.set_reference_sequence(1, len(chr_record), chr_record)
            return all_models, gi
        gi = GeneInfo.from_iqidx_cluster(genedb, list(range(g0, g1)), prepare_profiles=False)
    else:
        gene_list = list(genedb.region(seqid=chr_id, start=1, featuretype="gene"))
        if not gene_list:
            all_models.extend(novel_model_storage)
            gi = GeneInfo.from_region(chr_id, 1, len(chr_record))
            if chr_record is not None:
                gi.set_reference_sequence(1, len(chr_record), chr_record)
            return all_models, gi
        gi = GeneInfo(gene_list, genedb, prepare_profiles=False)

    if chr_record is not None:
        gi.set_reference_sequence(1, len(chr_record), chr_record)
    for isoform_id in gi.all_isoforms_exons.keys():
        all_models.append(TranscriptModel.from_reference_transcript(gi, isoform_id))
    all_models.extend(novel_model_storage)
    return all_models, gi
