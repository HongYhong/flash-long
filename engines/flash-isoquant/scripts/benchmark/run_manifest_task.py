#!/usr/bin/env python3
"""Execute one task entry from an all-samples benchmark manifest inside Slurm scratch."""

from __future__ import annotations

import argparse
import gzip
import json
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List


REPO_ROOT = Path("/gpfs/group/data_sci_group/weiwu_group/hongyanhong/projects/long_bench/IsoQuant")
BENCH_ROOT = REPO_ROOT / "isoquant_refactored"
CONDA_SH = "/gpfs/group/data_sci_group/weiwu_group/hongyanhong/miniconda3/miniconda3/etc/profile.d/conda.sh"


def _copy_or_gunzip(src: Path, dst_dir: Path, decompress: bool = True) -> Path:
    dst_dir.mkdir(parents=True, exist_ok=True)
    if src.name.endswith(".gz") and decompress:
        out = dst_dir / src.name[:-3]
        if not out.exists():
            with gzip.open(src, "rb") as f_in, open(out, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out, length=1024 * 1024)
        return out
    out = dst_dir / src.name
    if not out.exists():
        os.symlink(src, out)
    return out


def _prepare_bam(src_bam: Path, scratch: Path, link_name: str, threads: int) -> Path:
    has_bai = Path(f"{src_bam}.bai").exists() or src_bam.with_suffix(".bai").exists()
    if has_bai:
        return src_bam
    link_dir = scratch / "links"
    link_dir.mkdir(parents=True, exist_ok=True)
    dst_bam = link_dir / f"{link_name}.bam"
    if not dst_bam.exists():
        os.symlink(src_bam, dst_bam)
    if not Path(f"{dst_bam}.bai").exists():
        subprocess.run(["samtools", "index", "-@", str(min(8, threads)), str(dst_bam)], check=True)
    return dst_bam


def _prepare_ref_and_gtf(task: dict, scratch: Path) -> tuple[Path, Path | None]:
    in_dir = scratch / "inputs"
    ref = Path(task["ref"])
    if task["mode"] == "illumina_pair":
        # Copy FA into scratch so pyfaidx writes .fai in scratch instead of isoquant_tests/
        in_dir.mkdir(parents=True, exist_ok=True)
        dst_ref = in_dir / ref.name[:-3] if ref.name.endswith(".gz") else in_dir / ref.name
        if ref.name.endswith(".gz"):
            _copy_or_gunzip(ref, in_dir, decompress=True)
        elif not dst_ref.exists():
            shutil.copy2(ref, dst_ref)
        ref = dst_ref

    gtf = Path(task["gtf"]) if task.get("gtf") else None
    if gtf is not None and gtf.name.endswith(".gz"):
        gtf = _copy_or_gunzip(gtf, in_dir, decompress=True)
    return ref, gtf


def _common_cli(task: dict, ref: Path, gtf: Path | None, scratch: Path, prefix: str | None) -> str:
    parts: List[str] = ["--clean_start", "-d", task["data_type"]]
    if task.get("stranded"):
        parts.extend(["--stranded", task["stranded"]])
    if task.get("fl_data"):
        parts.append("--fl_data")

    mode = task["mode"]
    if mode == "fastq":
        parts.extend(["--fastq", str(task["fastq"])])
    elif mode == "yaml":
        yaml_path = scratch / f"{task['slug']}.yaml"
        lines = ['- data format: "bam"']
        for ent in task["yaml_entries"]:
            lines.append(f"- name: {ent['name']}")
            lines.append("  long read files:")
            for b in ent["bams"]:
                pb = _prepare_bam(Path(b), scratch, f"{ent['name']}_{Path(b).stem}", int(task["threads"]))
                lines.append(f"  - {pb}")
        yaml_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        parts.extend(["--yaml", str(yaml_path)])
    else:
        bam = _prepare_bam(Path(task["bam"]), scratch, task["sample_id"], int(task["threads"]))
        parts.extend(["--bam", str(bam)])

    if task.get("illumina_bam"):
        ill_bam = _prepare_bam(Path(task["illumina_bam"]), scratch, f"{task['sample_id']}_illumina", int(task["threads"]))
        parts.extend(["--illumina_bam", str(ill_bam)])

    if gtf is not None:
        parts.extend(["--genedb", str(gtf)])
        if task.get("complete_genedb", True):
            parts.append("--complete_genedb")
    if ref is not None:
        parts.extend(["-r", str(ref)])
    parts.extend(["-t", str(task["threads"])])
    if prefix:
        parts.extend(["-p", prefix])
    return " ".join(shlex.quote(x) for x in parts)


def _wrap_conda(body: str) -> str:
    return f"source {shlex.quote(CONDA_SH)} && conda activate sqanti3 && cd {shlex.quote(str(REPO_ROOT))} && {body}"


def _run_bench(log_dir: Path, repeats: int, baseline_cmd: str | None, candidate_cmd: str) -> None:
    cmd = [
        sys.executable,
        str(BENCH_ROOT / "scripts/benchmark/run_tiered_benchmark.py"),
        "--repeats",
        str(repeats),
        "--interval",
        "1.0",
        "--log-dir",
        str(log_dir),
        "--candidate-cmd",
        candidate_cmd,
    ]
    if baseline_cmd is not None:
        cmd.extend(["--baseline-cmd", baseline_cmd])
    subprocess.run(cmd, check=True)


def _run_summarize(ref_dir: Path, cand_dir: Path, prefix: str, out_json: Path) -> None:
    subprocess.run(
        [
            sys.executable,
            str(BENCH_ROOT / "scripts/benchmark/summarize_pair.py"),
            "--ref-dir",
            str(ref_dir),
            "--cand-dir",
            str(cand_dir),
            "--prefix",
            prefix,
            "--json-out",
            str(out_json),
        ],
        check=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description="Run a single manifest benchmark task.")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--index", type=int, required=True)
    parser.add_argument("--scratch", type=Path, required=True)
    parser.add_argument("--skip-existing", action="store_true", default=False)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    if args.index < 0 or args.index >= len(manifest):
        raise IndexError(f"Task index {args.index} out of range [0, {len(manifest)})")
    task = manifest[args.index]
    slug = task["slug"]
    log_dir = BENCH_ROOT / "benchmarks" / slug
    summary_path = log_dir / "benchmark_summary.json"

    if args.skip_existing and summary_path.is_file():
        data = json.loads(summary_path.read_text(encoding="utf-8"))
        s = data.get("summary", {})
        if s and all(v.get("failures", 1) == 0 for v in s.values()):
            print(f"[SKIP] {slug} already completed with 0 failures.")
            return 0

    log_dir.mkdir(parents=True, exist_ok=True)
    args.scratch.mkdir(parents=True, exist_ok=True)
    (log_dir / "task_meta.json").write_text(json.dumps(task, indent=2) + "\n", encoding="utf-8")

    ref, gtf = _prepare_ref_and_gtf(task, args.scratch)
    mode = task["mode"]
    repeats = int(task.get("repeats", 1))

    if mode == "cage":
        cage_flag = f"--cage {shlex.quote(str(task['cage_bed']))}"
        common = _common_cli(task, ref, gtf, args.scratch, task["sample_id"])
        cand_cmd = _wrap_conda(
            f"rm -rf {shlex.quote(str(args.scratch / 'next'))} && isoquant-next {common} {cage_flag} -o {shlex.quote(str(args.scratch / 'next'))}"
        )
        _run_bench(log_dir, repeats, None, cand_cmd)
        sdir = args.scratch / "next" / task["sample_id"]
        tss = sdir / f"{task['sample_id']}.TSS_prediction.tsv"
        gtf_out = sdir / f"{task['sample_id']}.transcript_models.gtf"
        gc = sdir / f"{task['sample_id']}.gene_counts.tsv"
        stats = {
            "tss_prediction_lines": sum(1 for _ in open(tss, encoding="utf-8")) if tss.exists() else 0,
            "transcript_models_count": sum(1 for ln in open(gtf_out, encoding="utf-8") if "\ttranscript\t" in ln) if gtf_out.exists() else 0,
            "gene_counts_lines": sum(1 for _ in open(gc, encoding="utf-8")) if gc.exists() else 0,
        }
        (log_dir / "cage_summary.json").write_text(json.dumps(stats, indent=2) + "\n", encoding="utf-8")
    elif mode == "warm_idx":
        db_path = args.scratch / "pre_ref.db"
        iq_path = args.scratch / "pre_ref.iqidx"
        subprocess.run(
            ["bash", "-lc", _wrap_conda(
                f"python -c \"from isoquant_lib.gtf2db import gtf2db; gtf2db({repr(str(gtf))}, {repr(str(db_path))}, True, False)\" && "
                f"python -c \"from isoquant_next.annotation.iqidx import compile_gtf_to_iqidx; compile_gtf_to_iqidx({repr(str(gtf))}, {repr(str(iq_path))})\""
            )],
            check=True,
        )
        bam = _prepare_bam(Path(task["bam"]), args.scratch, task["sample_id"], int(task["threads"]))
        base_flags = f"-d {shlex.quote(task['data_type'])} --fl_data --bam {shlex.quote(str(bam))} --complete_genedb -r {shlex.quote(str(ref))} -t {task['threads']} -p {shlex.quote(task['sample_id'])}"
        orig_cmd = _wrap_conda(f"rm -rf {shlex.quote(str(args.scratch / 'orig'))} && python isoquant.py {base_flags} --genedb {shlex.quote(str(db_path))} -o {shlex.quote(str(args.scratch / 'orig'))}")
        next_cmd = _wrap_conda(f"rm -rf {shlex.quote(str(args.scratch / 'next'))} && isoquant-next {base_flags} --genedb {shlex.quote(str(iq_path))} -o {shlex.quote(str(args.scratch / 'next'))}")
        _run_bench(log_dir, repeats, orig_cmd, next_cmd)
        _run_summarize(args.scratch / "orig", args.scratch / "next", task["sample_id"], log_dir / "accuracy_summary.json")
    elif mode == "resume_vis":
        clean_dir = args.scratch / "clean"
        res_dir = args.scratch / "resumed"
        common = _common_cli(task, ref, gtf, args.scratch, task["sample_id"])
        clean_cmd = _wrap_conda(f"rm -rf {shlex.quote(str(clean_dir))} && isoquant-next {common} -o {shlex.quote(str(clean_dir))}")
        prep_cmd = _wrap_conda(
            f"rm -rf {shlex.quote(str(res_dir))} && "
            f"isoquant-next {common} -o {shlex.quote(str(res_dir))} && "
            f"python -c \"import json, pathlib; p=pathlib.Path('{res_dir}/.stage_state.json'); d=json.loads(p.read_text()); "
            f"cs=d.get('completed_stages',{{}}); d['completed_stages']={{k:v for k,v in cs.items() if k=='reference_and_annotation'}}; p.write_text(json.dumps(d, indent=2))\" && "
            f"rm -rf {shlex.quote(str(res_dir / task['sample_id']))}"
        )
        subprocess.run(["bash", "-lc", prep_cmd], check=True)
        # Drop --clean_start for resume run so StageStateTracker reuses prepared .iqidx
        resume_common = common.replace("--clean_start ", "")
        resume_cmd = _wrap_conda(f"isoquant-next {resume_common} --resume -o {shlex.quote(str(res_dir))}")
        _run_bench(log_dir, 1, clean_cmd, resume_cmd)
        # Also verify the resumed StageStateTracker skipped reference_and_annotation
        next_log = log_dir / "candidate.rep0.stderr.log"
        res_iq_log = res_dir / "isoquant.log"
        log_blob = ""
        if next_log.exists():
            log_blob += next_log.read_text(encoding="utf-8", errors="replace")
        if res_iq_log.exists():
            log_blob += "\n" + res_iq_log.read_text(encoding="utf-8", errors="replace")
        skipped_ref = "skipping completed stage 'reference_and_annotation'" in log_blob
        _run_summarize(clean_dir, res_dir, task["sample_id"], log_dir / "accuracy_summary.json")

        # Run isoquant-next-visualize end-to-end (PNG + PDF)
        vis_png = log_dir / "vis_png"
        vis_pdf = log_dir / "vis_pdf"
        subprocess.run(
            ["bash", "-lc", _wrap_conda(
                f"isoquant-next-visualize -o {shlex.quote(str(clean_dir))} --visual_output {shlex.quote(str(vis_png))} --counts --find_genes --num_genes 6 --format png && "
                f"isoquant-next-visualize -o {shlex.quote(str(clean_dir))} --visual_output {shlex.quote(str(vis_pdf))} --counts --find_genes --num_genes 6 --format pdf"
            )],
            check=True,
        )
        pngs = sorted(p.name for p in vis_png.glob("*.png")) if vis_png.exists() else []
        pdfs = sorted(p.name for p in vis_pdf.glob("*.pdf")) if vis_pdf.exists() else []
        (log_dir / "visualize_summary.json").write_text(
            json.dumps(
                {
                    "resume_skipped_reference_and_annotation": skipped_ref,
                    "png_count": len(pngs),
                    "pdf_count": len(pdfs),
                    "png_files": pngs,
                    "pdf_files": pdfs,
                },
                indent=2,
            )
            + "\n",
            encoding="utf-8",
        )
    else:
        prefix = None if mode == "yaml" else task["sample_id"]
        common = _common_cli(task, ref, gtf, args.scratch, prefix)
        orig_cmd = _wrap_conda(f"rm -rf {shlex.quote(str(args.scratch / 'orig'))} && python isoquant.py {common} -o {shlex.quote(str(args.scratch / 'orig'))}")
        next_cmd = _wrap_conda(f"rm -rf {shlex.quote(str(args.scratch / 'next'))} && isoquant-next {common} -o {shlex.quote(str(args.scratch / 'next'))}")
        _run_bench(log_dir, repeats, orig_cmd, next_cmd)

        if mode == "yaml":
            for ent in task["yaml_entries"]:
                _run_summarize(
                    args.scratch / "orig",
                    args.scratch / "next",
                    ent["name"],
                    log_dir / f"accuracy_{ent['name']}.json",
                )
        else:
            _run_summarize(args.scratch / "orig", args.scratch / "next", task["sample_id"], log_dir / "accuracy_summary.json")
            if task.get("eval_sirv"):
                subprocess.run(
                    [
                        sys.executable,
                        str(BENCH_ROOT / "scripts/benchmark/eval_simulated_sirv.py"),
                        "--ref-dir",
                        str(args.scratch / "orig"),
                        "--cand-dir",
                        str(args.scratch / "next"),
                        "--prefix",
                        task["sample_id"],
                        "--json-out",
                        str(log_dir / f"eval_{task['sample_id']}.json"),
                    ],
                    check=True,
                )

    for src, dst_name in (
        (args.scratch / "orig/isoquant.log", "orig.isoquant.log"),
        (args.scratch / "next/isoquant.log", "next.isoquant.log"),
        (args.scratch / "next/run_manifest.json", "next.run_manifest.json"),
    ):
        if src.exists():
            shutil.copy2(src, log_dir / dst_name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
