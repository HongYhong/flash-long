"""Output directory configuration loader preferring `run_manifest.json` with read-only legacy `.params` fallback."""

from __future__ import annotations

from dataclasses import dataclass, field
import glob
import json
import os
from pathlib import Path
import pickle
from typing import Any, Dict, List, Optional


@dataclass
class SamplePaths:
    """Resolved artifact paths for a single IsoQuant output sample directory."""

    prefix: str
    sample_dir: Path
    transcript_models_gtf: Optional[Path] = None
    extended_annotation_gtf: Optional[Path] = None
    transcript_counts_tsv: Optional[Path] = None
    discovered_transcript_counts_tsv: Optional[Path] = None
    gene_counts_tsv: Optional[Path] = None
    transcript_tpm_tsv: Optional[Path] = None
    gene_tpm_tsv: Optional[Path] = None


def _first_existing(candidates: List[Path]) -> Optional[Path]:
    for p in candidates:
        if p.is_file():
            return p
    return None


@dataclass
class OutputConfig:
    """Inspects an IsoQuant / isoquant-next output directory without requiring live pipeline state."""

    output_directory: str
    use_counts: bool = False
    ref_only: Optional[bool] = None
    gtf: Optional[str] = None
    manifest: Dict[str, Any] = field(default_factory=dict)
    cli_options: Dict[str, Any] = field(default_factory=dict)
    data_type: str = "nanopore"
    samples: List[SamplePaths] = field(default_factory=list)
    reference_gtf: Optional[Path] = None
    iqidx_path: Optional[Path] = None
    source_kind: str = "directory_scan"

    def __post_init__(self) -> None:
        self.load()

    def load(self) -> None:
        root = Path(os.path.abspath(self.output_directory))
        manifest_path = root / "run_manifest.json"
        cli_path = root / "run_cli_options.json"
        legacy_params = root / ".params"

        prefix_hint: Optional[str] = None
        if manifest_path.is_file():
            self.manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            self.source_kind = "run_manifest"
            cfg = self.manifest.get("config", {})
            self.data_type = str(cfg.get("data_type", "nanopore"))
            rt = cfg.get("runtime", {})
            prefix_hint = rt.get("prefix")
            genedb = cfg.get("genedb")
            if cli_path.is_file():
                self.cli_options = json.loads(cli_path.read_text(encoding="utf-8"))
                prefix_hint = self.cli_options.get("prefix") or prefix_hint
                genedb = self.cli_options.get("genedb") or genedb
            self._resolve_annotation_paths(root, genedb)
        elif legacy_params.is_file():
            self.source_kind = "legacy_params"
            with open(legacy_params, "rb") as fh:
                params_obj = pickle.load(fh)
            self.data_type = str(getattr(params_obj, "data_type", "nanopore"))
            prefix_hint = getattr(params_obj, "prefix", None)
            genedb = getattr(params_obj, "original_annotation", None) or getattr(params_obj, "genedb", None)
            self._resolve_annotation_paths(root, genedb)
        else:
            self._resolve_annotation_paths(root, None)

        if self.gtf:
            gtf_p = Path(os.path.abspath(self.gtf))
            if gtf_p.suffix.lower() == ".iqidx":
                self.iqidx_path = gtf_p
            else:
                self.reference_gtf = gtf_p

        self.samples = self._discover_samples(root, prefix_hint)
        if self.ref_only is None:
            has_discovered = any(s.transcript_models_gtf is not None for s in self.samples)
            self.ref_only = not has_discovered

    def _resolve_annotation_paths(self, root: Path, genedb: Optional[str]) -> None:
        if genedb:
            gp = Path(genedb)
            if gp.is_file():
                if gp.suffix.lower() == ".iqidx":
                    self.iqidx_path = gp
                elif not gp.name.lower().endswith(".db"):
                    self.reference_gtf = gp
        if self.iqidx_path is None:
            found_iqidx = sorted(root.glob("*.iqidx"))
            if found_iqidx:
                self.iqidx_path = found_iqidx[0]
        if self.reference_gtf is None:
            found_gtf = [
                p for p in sorted(root.glob("*.gtf")) + sorted(root.glob("*.gtf.gz"))
                if ".cleaned." not in p.name
            ]
            if found_gtf:
                self.reference_gtf = found_gtf[0]

    def _discover_samples(self, root: Path, prefix_hint: Optional[str]) -> List[SamplePaths]:
        candidates: List[Path] = []
        for entry in sorted(root.iterdir()):
            if not entry.is_dir() or entry.name.startswith("."):
                continue
            if list(entry.glob("*.gtf")) or list(entry.glob("*.tsv*")):
                candidates.append(entry)
        if not candidates:
            candidates = [root]

        resolved: List[SamplePaths] = []
        for sdir in candidates:
            sample_prefix = sdir.name if sdir != root else (prefix_hint or "OUT")
            # Auto-detect file prefix if directory holds files under another prefix
            if not (sdir / f"{sample_prefix}.gene_counts.tsv").exists() and not (
                sdir / f"{sample_prefix}.transcript_models.gtf"
            ).exists():
                models = sorted(sdir.glob("*.transcript_models.gtf"))
                counts = sorted(sdir.glob("*.gene_counts.tsv"))
                probe = models[0] if models else (counts[0] if counts else None)
                if probe is not None:
                    sample_prefix = probe.name.split(".", 1)[0]

            sp = SamplePaths(
                prefix=sample_prefix,
                sample_dir=sdir,
                transcript_models_gtf=_first_existing([sdir / f"{sample_prefix}.transcript_models.gtf"]),
                extended_annotation_gtf=_first_existing([sdir / f"{sample_prefix}.extended_annotation.gtf"]),
                transcript_counts_tsv=_first_existing(
                    [
                        sdir / f"{sample_prefix}.transcript_counts.tsv",
                        sdir / f"{sample_prefix}.transcript_grouped_counts.tsv",
                    ]
                ),
                discovered_transcript_counts_tsv=_first_existing(
                    [
                        sdir / f"{sample_prefix}.discovered_transcript_counts.tsv",
                        sdir / f"{sample_prefix}.discovered_transcript_grouped_counts.tsv",
                    ]
                ),
                gene_counts_tsv=_first_existing(
                    [
                        sdir / f"{sample_prefix}.gene_counts.tsv",
                        sdir / f"{sample_prefix}.gene_grouped_counts.tsv",
                    ]
                ),
                transcript_tpm_tsv=_first_existing([sdir / f"{sample_prefix}.transcript_tpm.tsv"]),
                gene_tpm_tsv=_first_existing([sdir / f"{sample_prefix}.gene_tpm.tsv"]),
            )
            resolved.append(sp)
        return resolved

    def list_combined_count_tables(self) -> List[Path]:
        root = Path(os.path.abspath(self.output_directory))
        return sorted(Path(p) for p in glob.glob(str(root / "combined_*.tsv")))
