"""Deterministic external aligner runner using blake2b fingerprints, ExitStack pipes, and CSI fallback."""

from __future__ import annotations

import hashlib
import logging
import os
import shlex
import subprocess
from contextlib import ExitStack
from pathlib import Path
from typing import Iterable, List, Optional

import pysam

from isoquant_next.config.enums import ASSEMBLY, DATATYPE_TO_ALIGNER, NANOPORE_DATA, PACBIO_CCS_DATA
from isoquant_next.core.error_codes import IsoQuantError, IsoQuantExitCode
from isoquant_next.core.utils import get_path_to_program

logger = logging.getLogger("IsoQuant")
_BAI_CONTIG_LIMIT = (1 << 29) - 1  # ~536.87 Mb BAM BAI ceiling


def stable_alignment_cache_key(
    reads_path: str | Path,
    ref_path: str | Path,
    data_type: str,
    extra_opts: str = "",
) -> str:
    """Order-independent, PYTHONHASHSEED-immune blake2b key replacing Python's salted hash()."""
    h = hashlib.blake2b(digest_size=12)
    h.update(str(Path(reads_path).resolve()).encode("utf-8"))
    h.update(b"|")
    h.update(str(Path(ref_path).resolve()).encode("utf-8"))
    h.update(b"|")
    h.update(data_type.encode("utf-8"))
    h.update(b"|")
    h.update(extra_opts.encode("utf-8"))
    return h.hexdigest()


def index_bam_safe(bam_path: str | Path, threads: int = 1) -> Path:
    """Create `.bai`, automatically falling back to `.csi` when any contig exceeds BAI's 536Mb limit."""
    p = Path(bam_path)
    needs_csi = False
    with pysam.AlignmentFile(str(p), "rb", check_sq=False) as af:
        if any(length >= _BAI_CONTIG_LIMIT for length in af.lengths):
            needs_csi = True

    if needs_csi:
        csi_path = Path(str(p) + ".csi")
        if not csi_path.exists():
            pysam.index("-c", "-@", str(max(1, threads)), str(p))
        return csi_path

    bai_path = Path(str(p) + ".bai")
    if not bai_path.exists():
        try:
            pysam.index("-@", str(max(1, threads)), str(p))
        except pysam.SamtoolsError:
            pysam.index("-c", "-@", str(max(1, threads)), str(p))
            return Path(str(p) + ".csi")
    return bai_path


def build_minimap2_args(
    data_type: str,
    ref_path: str | Path,
    reads_path: str | Path,
    threads: int = 8,
    junc_bed: Optional[str | Path] = None,
    secondary: bool = True,
    stranded: str = "none",
    extra_options: Optional[List[str]] = None,
) -> List[str]:
    preset = {
        PACBIO_CCS_DATA: "splice:hq",
        ASSEMBLY: "splice:hq",
        NANOPORE_DATA: "splice",
    }.get(data_type, "splice")
    kmer = "14" if data_type == NANOPORE_DATA else "15"
    cmd = [
        "-a",
        "-x",
        preset,
        "-k",
        kmer,
        "-w",
        "5",
        "-t",
        str(max(1, threads)),
        "--MD",
        "-Y",
        f"--secondary={'yes' if secondary else 'no'}",
    ]
    if stranded == "forward":
        cmd.append("-uf")
    if junc_bed is not None:
        cmd.extend(["--junc-bed", str(junc_bed)])
    if extra_options:
        cmd.extend(extra_options)
    cmd.extend([str(ref_path), str(reads_path)])
    return cmd


class ReadMapper:
    """Safe subprocess orchestrator for minimap2 + samtools sort without SIGPIPE or fd leaks."""

    def __init__(
        self,
        data_type: str,
        threads: int = 8,
        minimap2_path: Optional[str] = None,
        samtools_path: Optional[str] = None,
    ) -> None:
        self.data_type = data_type
        self.threads = max(1, int(threads))
        self.minimap2_bin = minimap2_path or get_path_to_program(DATATYPE_TO_ALIGNER.get(data_type, "minimap2"))
        self.samtools_bin = samtools_path or get_path_to_program("samtools")

    def align_to_sorted_bam(
        self,
        reads_path: str | Path,
        ref_path: str | Path,
        out_dir: str | Path,
        junc_bed: Optional[str | Path] = None,
        sample_label: Optional[str] = None,
        *,
        stranded: str = "none",
        mapping_options: Optional[str] = None,
    ) -> Path:
        if not self.minimap2_bin:
            raise IsoQuantError(IsoQuantExitCode.ALIGNER_NOT_FOUND, "minimap2 executable not found")
        if not self.samtools_bin:
            raise IsoQuantError(IsoQuantExitCode.SAMTOOLS_FAILED, "samtools executable not found")

        out_base = Path(out_dir)
        out_base.mkdir(parents=True, exist_ok=True)
        fp = stable_alignment_cache_key(reads_path, ref_path, self.data_type, str(junc_bed or ""))
        stem = sample_label or Path(reads_path).stem
        out_bam = out_base / f"{stem}.{fp}.sorted.bam"
        tmp_bam = out_bam.with_suffix(".tmp.bam")
        log_path = out_base / f"{stem}.{fp}.align.log"

        extra = shlex.split(mapping_options) if mapping_options else None
        mm_cmd = [self.minimap2_bin] + build_minimap2_args(
            self.data_type,
            ref_path,
            reads_path,
            threads=self.threads,
            junc_bed=junc_bed,
            stranded=stranded,
            extra_options=extra,
        )
        sort_cmd = [self.samtools_bin, "sort", "-@", str(min(self.threads, 8)), "-o", str(tmp_bam), "-"]

        with ExitStack() as stack:
            log_fh = stack.enter_context(open(log_path, "w", encoding="utf-8"))
            p_mm = stack.enter_context(
                subprocess.Popen(mm_cmd, stdout=subprocess.PIPE, stderr=log_fh)
            )
            assert p_mm.stdout is not None
            p_sort = stack.enter_context(
                subprocess.Popen(sort_cmd, stdin=p_mm.stdout, stdout=log_fh, stderr=log_fh)
            )
            p_mm.stdout.close()  # Allow SIGPIPE delivery if samtools exits prematurely

            sort_rc = p_sort.wait()
            mm_rc = p_mm.wait()

        if sort_rc != 0:
            raise IsoQuantError(IsoQuantExitCode.SAMTOOLS_FAILED, f"samtools sort failed ({sort_rc})")
        if mm_rc != 0:
            raise IsoQuantError(IsoQuantExitCode.ALIGNMENT_FAILED, f"minimap2 failed ({mm_rc})")

        os.replace(tmp_bam, out_bam)
        index_bam_safe(out_bam, threads=self.threads)
        return out_bam
