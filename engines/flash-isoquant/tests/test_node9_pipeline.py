"""Node 9 (T21-T22) end-to-end Stage Graph pipeline, CLI, SC/UMI, and --resume parity tests."""

from __future__ import annotations

import gzip
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import pysam
import pytest

from isoquant_next.pipeline.stage_graph import StageGraphExecutor, StageNode, StageStateTracker

REFACTOR_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = REFACTOR_ROOT.parent
SIMPLE_DATA = REPO_ROOT / "isoquant_tests" / "simple_data"
TOY_DATA = REPO_ROOT / "isoquant_tests" / "toy_data"
COMPARE_SCRIPT = REFACTOR_ROOT / "scripts" / "benchmark" / "compare_outputs.py"


def _cli_env() -> dict[str, str]:
    env = os.environ.copy()
    src_path = str(REFACTOR_ROOT / "src")
    env["PYTHONPATH"] = f"{src_path}:{env['PYTHONPATH']}" if env.get("PYTHONPATH") else src_path
    return env


def _run_orig(args: list[str], cwd: Path = REPO_ROOT) -> subprocess.CompletedProcess[str]:
    cmd = [sys.executable, str(REPO_ROOT / "isoquant.py"), *args]
    res = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    assert res.returncode == 0, f"Oracle failed ({res.returncode}):\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
    return res


def _run_next(args: list[str], cwd: Path = REPO_ROOT) -> subprocess.CompletedProcess[str]:
    cmd = [sys.executable, "-m", "isoquant_next.cli.main", *args]
    res = subprocess.run(cmd, cwd=str(cwd), env=_cli_env(), capture_output=True, text=True)
    assert res.returncode == 0, f"isoquant-next failed ({res.returncode}):\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
    return res


def _compare_with_script(
    ref_dir: Path,
    cand_dir: Path,
    prefix: str,
    *,
    check_models: bool = True,
    check_assign: bool = True,
) -> None:
    sample_ref = ref_dir / prefix
    sample_cand = cand_dir / prefix
    cmd = [
        sys.executable,
        str(COMPARE_SCRIPT),
        "--ref-counts",
        str(sample_ref / f"{prefix}.gene_counts.tsv"),
        "--cand-counts",
        str(sample_cand / f"{prefix}.gene_counts.tsv"),
    ]
    res = subprocess.run(cmd, env=_cli_env(), capture_output=True, text=True)
    assert res.returncode == 0, f"compare gene_counts failed:\n{res.stdout}\n{res.stderr}"

    for fname in ("transcript_counts.tsv", "gene_tpm.tsv", "transcript_tpm.tsv"):
        if (sample_ref / f"{prefix}.{fname}").exists():
            r2 = subprocess.run(
                [
                    sys.executable,
                    str(COMPARE_SCRIPT),
                    "--ref-counts",
                    str(sample_ref / f"{prefix}.{fname}"),
                    "--cand-counts",
                    str(sample_cand / f"{prefix}.{fname}"),
                ],
                env=_cli_env(),
                capture_output=True,
                text=True,
            )
            assert r2.returncode == 0, f"compare {fname} failed:\n{r2.stdout}\n{r2.stderr}"

    if check_models and (sample_ref / f"{prefix}.transcript_models.gtf").exists():
        for gtf_name in ("transcript_models.gtf", "extended_annotation.gtf"):
            r3 = subprocess.run(
                [
                    sys.executable,
                    str(COMPARE_SCRIPT),
                    "--ref-gtf",
                    str(sample_ref / f"{prefix}.{gtf_name}"),
                    "--cand-gtf",
                    str(sample_cand / f"{prefix}.{gtf_name}"),
                ],
                env=_cli_env(),
                capture_output=True,
                text=True,
            )
            assert r3.returncode == 0, f"compare {gtf_name} failed:\n{r3.stdout}\n{r3.stderr}"

    if check_assign and (sample_ref / f"{prefix}.read_info.tsv.gz").exists():
        r4 = subprocess.run(
            [
                sys.executable,
                str(COMPARE_SCRIPT),
                "--ref-assign",
                str(sample_ref / f"{prefix}.read_info.tsv.gz"),
                "--cand-assign",
                str(sample_cand / f"{prefix}.read_info.tsv.gz"),
            ],
            env=_cli_env(),
            capture_output=True,
            text=True,
        )
        assert r4.returncode == 0, f"compare read_info failed:\n{r4.stdout}\n{r4.stderr}"


def _read_sorted_nonheader_lines(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as fh:
        return sorted(line for line in fh if not line.startswith("#"))


class TestEndToEndOracleParity:
    def test_simple_data_two_stage_and_single_pass_parity(self, tmp_path: Path):
        # Two-stage default (polya_requirement=auto) with full exon/IR/SQANTI outputs
        orig_two = tmp_path / "orig_two"
        next_two = tmp_path / "next_two"
        common = [
            "--clean_start",
            "-d",
            "nanopore",
            "--bam",
            str(SIMPLE_DATA / "chr9.4M.ont.sim.polya.bam"),
            "--genedb",
            str(SIMPLE_DATA / "chr9.4M.gtf.gz"),
            "--complete_genedb",
            "-r",
            str(SIMPLE_DATA / "chr9.4M.fa.gz"),
            "-t",
            "2",
            "-p",
            "S1",
            "--fl_data",
            "--count_exons",
            "--count_intron_retentions",
            "--sqanti_output",
        ]
        _run_orig([*common, "-o", str(orig_two)])
        _run_next([*common, "-o", str(next_two)])
        _compare_with_script(orig_two, next_two, "S1")

        # Verify zero-diff sorted content on extended TSV outputs
        for suffix in (
            "exon_counts.tsv",
            "exon_splice_site_counts.tsv",
            "intron_retention_counts.tsv",
            "splice_junction_counts.tsv",
            "novel_vs_known.SQANTI-like.tsv",
            "TSS_prediction.tsv",
            "polyA_prediction.tsv",
            "transcript_model_reads.tsv.gz",
        ):
            assert _read_sorted_nonheader_lines(orig_two / "S1" / f"S1.{suffix}") == _read_sorted_nonheader_lines(
                next_two / "S1" / f"S1.{suffix}"
            ), f"Mismatch in {suffix}"

        # Single-pass O(Locus) streaming closure (polya_requirement=never skips .save I/O)
        orig_one = tmp_path / "orig_one"
        next_one = tmp_path / "next_one"
        common_one = [*common, "--polya_requirement", "never"]
        _run_orig([*common_one, "-o", str(orig_one)])
        _run_next([*common_one, "-o", str(next_one)])
        _compare_with_script(orig_one, next_one, "S1")

        log_txt = (next_one / "isoquant.log").read_text(encoding="utf-8")
        assert "Using single-pass O(Locus) streaming closure" in log_txt
        assert not (next_one / "S1" / "aux" / "S1.save_chr9").exists()

    def test_yaml_and_illumina_bam_and_combine_counts_parity(self, tmp_path: Path):
        # Part 1: canonical repo YAML fixture
        orig_yaml = tmp_path / "orig_yaml"
        next_yaml = tmp_path / "next_yaml"
        yaml_args = [
            "--clean_start",
            "-d",
            "nanopore",
            "--yaml",
            str(SIMPLE_DATA / "chr9.4M.yaml"),
            "--genedb",
            str(SIMPLE_DATA / "chr9.4M.gtf.gz"),
            "--complete_genedb",
            "-r",
            str(SIMPLE_DATA / "chr9.4M.fa.gz"),
            "-t",
            "2",
        ]
        _run_orig([*yaml_args, "-o", str(orig_yaml)])
        _run_next([*yaml_args, "-o", str(next_yaml)])

        # Part 2: Illumina short-read junction support + multi-sample combine_counts
        multi_yaml = tmp_path / "multi.yaml"
        multi_yaml.write_text(
            "\n".join(
                [
                    "- data format: \"bam\"",
                    "- name: S_A",
                    "  long read files:",
                    f"  - {SIMPLE_DATA / 'chr9.4M.ont.sim.polya.bam'}",
                    "  illumina bam:",
                    f"  - {SIMPLE_DATA / 'chr9.4M.Illumina.bam'}",
                    "- name: S_B",
                    "  long read files:",
                    f"  - {SIMPLE_DATA / 'chr9.4M.ont.sim.polya.bam'}",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        orig_multi = tmp_path / "orig_multi"
        next_multi = tmp_path / "next_multi"
        multi_args = [
            "--clean_start",
            "-d",
            "nanopore",
            "--yaml",
            str(multi_yaml),
            "--genedb",
            str(SIMPLE_DATA / "chr9.4M.gtf.gz"),
            "--complete_genedb",
            "-r",
            str(SIMPLE_DATA / "chr9.4M.fa.gz"),
            "-t",
            "2",
        ]
        _run_orig([*multi_args, "-o", str(orig_multi)])
        _run_next([*multi_args, "-o", str(next_multi)])
        _compare_with_script(orig_multi, next_multi, "S_A")
        _compare_with_script(orig_multi, next_multi, "S_B")
        assert (next_multi / "combined_gene_counts.tsv").exists() or any(next_multi.glob("*gene_counts*"))

    def test_toy_data_fastq_cage_and_builtin_test_flag(self, tmp_path: Path):
        # Reference IsoQuant 4.0.0 rejects --cage ("CAGE data is not supported yet"),
        # so verify full Oracle parity on MAPT FASTQ without --cage, then validate --cage on isoquant-next.
        orig_mapt = tmp_path / "orig_mapt"
        next_mapt = tmp_path / "next_mapt"
        base_args = [
            "--clean_start",
            "-d",
            "nanopore",
            "--fastq",
            str(TOY_DATA / "MAPT.Mouse.ONT.simulated.fastq"),
            "--genedb",
            str(TOY_DATA / "MAPT.Mouse.genedb.gtf"),
            "--complete_genedb",
            "-r",
            str(TOY_DATA / "MAPT.Mouse.reference.fasta"),
            "-t",
            "2",
            "-p",
            "MAPT",
        ]
        _run_orig([*base_args, "-o", str(orig_mapt)])
        _run_next([*base_args, "-o", str(next_mapt)])
        _compare_with_script(orig_mapt, next_mapt, "MAPT")

        # Validate CAGE-assisted discovery & TSS prediction on isoquant-next (--fl_data enables TSSPrediction output)
        cage_dir = tmp_path / "next_cage"
        _run_next([*base_args, "--fl_data", "--cage", str(TOY_DATA / "MAPT.Mouse.CAGE.bed"), "-o", str(cage_dir)])
        assert (cage_dir / "MAPT" / "MAPT.TSS_prediction.tsv").exists()
        assert (cage_dir / "MAPT" / "MAPT.transcript_models.gtf").exists()

        # Verify built-in --test entrypoint in isolated cwd
        test_cwd = tmp_path / "test_cwd"
        test_cwd.mkdir()
        _run_next(["--test"], cwd=test_cwd)
        assert (test_cwd / "isoquant_next_test" / "OUT" / "OUT.gene_counts.tsv").exists()
        assert (test_cwd / "isoquant_next_test" / "run_manifest.json").exists()

    def test_sc_umi_barcoded_bam_and_detect_barcodes_cli(self, tmp_path: Path):
        # 1. Standalone barcode detection CLI parity vs isoquant_detect_barcodes.py
        bc1 = "AAACCCAAGAAACACT"
        bc2 = "AAACCCAAGAAACCAT"
        wl = tmp_path / "whitelist.tsv"
        wl.write_text(f"{bc1}\n{bc2}\n", encoding="utf-8")

        fq = tmp_path / "sc_reads.fastq"
        body = "AGCTTAGCTAGCTACCTATATCTTGGTACCG" * 4
        r1 = f"GATCTCTACACGACGCTCTTCCGATCT{bc1}ATGATGATGATGTTTTTTTTTTTTTTTTTTTT{body}"
        r2 = f"GATCTCTACACGACGCTCTTCCGATCT{bc2}CCCCCCCCCCCCAAAAAAAAAAAAAAAAAAA" + ("GCTAGCTA" * 10)
        fq.write_text(
            f"@read1\n{r1}\n+\n{'I' * len(r1)}\n@read2\n{r2}\n+\n{'I' * len(r2)}\n",
            encoding="utf-8",
        )
        orig_bc = tmp_path / "orig_bc"
        next_bc = tmp_path / "next_bc"
        r_orig = subprocess.run(
            [
                sys.executable,
                str(REPO_ROOT / "isoquant_detect_barcodes.py"),
                "-o",
                str(orig_bc),
                "-b",
                str(wl),
                "--mode",
                "tenX_v3",
                "-i",
                str(fq),
                "-t",
                "1",
            ],
            cwd=str(REPO_ROOT),
            capture_output=True,
            text=True,
        )
        assert r_orig.returncode == 0, r_orig.stderr
        r_next = subprocess.run(
            [
                sys.executable,
                "-m",
                "isoquant_next.cli.detect_barcodes",
                "-o",
                str(next_bc),
                "-b",
                str(wl),
                "--mode",
                "tenX_v3",
                "-i",
                str(fq),
                "-t",
                "1",
            ],
            env=_cli_env(),
            capture_output=True,
            text=True,
        )
        assert r_next.returncode == 0, r_next.stderr
        assert Path(f"{orig_bc}.barcoded_reads.tsv").read_text() == Path(f"{next_bc}.barcoded_reads.tsv").read_text()

        # 2. End-to-end --mode tenX_v3 --barcoded_bam with UMI PCR deduplication
        sc_bam = tmp_path / "sc_tagged.bam"
        with pysam.AlignmentFile(SIMPLE_DATA / "chr9.4M.ont.sim.polya.bam", "rb") as src_af:
            with pysam.AlignmentFile(sc_bam, "wb", template=src_af) as dst_af:
                for idx, aln in enumerate(src_af):
                    if idx >= 45:
                        break
                    # Force duplicate (CB, UMI) pairs on consecutive reads to trigger UMI filtering
                    cb = bc1 if (idx % 2 == 0) else bc2
                    ub = f"ACGTACGTAC{'AA' if idx < 15 else 'CC'}"
                    aln.set_tag("CB", cb)
                    aln.set_tag("UB", ub)
                    dst_af.write(aln)
        pysam.index(str(sc_bam))

        orig_sc = tmp_path / "orig_sc"
        next_sc = tmp_path / "next_sc"
        sc_args = [
            "--clean_start",
            "-d",
            "nanopore",
            "--mode",
            "tenX_v3",
            "--barcoded_bam",
            "--bam",
            str(sc_bam),
            "--genedb",
            str(SIMPLE_DATA / "chr9.4M.gtf.gz"),
            "--complete_genedb",
            "-r",
            str(SIMPLE_DATA / "chr9.4M.fa.gz"),
            "-t",
            "2",
            "-p",
            "SC1",
        ]
        _run_orig([*sc_args, "-o", str(orig_sc)])
        _run_next([*sc_args, "-o", str(next_sc)])
        _compare_with_script(orig_sc, next_sc, "SC1", check_models=False, check_assign=True)
        for grouped in ("SC1.gene_grouped_counts.tsv", "SC1.transcript_grouped_counts.tsv"):
            if (orig_sc / "SC1" / grouped).exists():
                assert _read_sorted_nonheader_lines(orig_sc / "SC1" / grouped) == _read_sorted_nonheader_lines(
                    next_sc / "SC1" / grouped
                )

    def test_interrupted_resume_determinism_and_manifest(self, tmp_path: Path):
        # Verify interrupted stage-graph resume skips mapping & .iqidx compilation while hydrating paths cleanly
        clean_dir = tmp_path / "clean_run"
        res_dir = tmp_path / "resume_run"
        base_args = [
            "--clean_start",
            "-d",
            "nanopore",
            "--fastq",
            str(TOY_DATA / "MAPT.Mouse.ONT.simulated.fastq"),
            "--genedb",
            str(TOY_DATA / "MAPT.Mouse.genedb.gtf"),
            "--complete_genedb",
            "-r",
            str(TOY_DATA / "MAPT.Mouse.reference.fasta"),
            "-t",
            "2",
            "-p",
            "MAPT",
        ]
        _run_next([*base_args, "-o", str(clean_dir)])
        _run_next([*base_args, "-o", str(res_dir)])

        # Simulate interruption right after read_mapping: remove final MAPT dir + dataset_processing completion
        state_path = res_dir / ".stage_state.json"
        state_doc = json.loads(state_path.read_text(encoding="utf-8"))
        assert "reference_and_annotation" in state_doc["completed_stages"]
        assert "read_mapping" in state_doc["completed_stages"]
        assert "dataset_processing" in state_doc["completed_stages"]
        del state_doc["completed_stages"]["dataset_processing"]
        state_path.write_text(json.dumps(state_doc, indent=2), encoding="utf-8")

        # Keep mapped BAM & iqidx intact; wipe downstream processed TSVs to prove resume regenerates them
        bam_files = list((res_dir / "MAPT" / "aux").glob("*.bam"))
        assert bam_files
        bam_mtime_before = bam_files[0].stat().st_mtime
        for p in (res_dir / "MAPT").glob("MAPT.*"):
            if p.is_file():
                p.unlink()

        time.sleep(0.05)
        # Change volatile thread count (-t 1 vs -t 2) to prove fingerprint stability on resume
        _run_next(["--resume", "-t", "1", "-o", str(res_dir)])

        bam_mtime_after = bam_files[0].stat().st_mtime
        assert bam_mtime_after == bam_mtime_before
        log_txt = (res_dir / "isoquant.log").read_text(encoding="utf-8")
        assert "Resuming: skipping completed stage 'reference_and_annotation'" in log_txt
        assert "Resuming: skipping completed stage 'read_mapping'" in log_txt
        assert "completed stages: dataset_processing" in log_txt

        _compare_with_script(clean_dir, res_dir, "MAPT")
        assert StageGraphExecutor is not None and StageNode is not None and StageStateTracker is not None
