#!/usr/bin/env python3
"""Normalize AdapterRemoval outputs into MSTFlow fragment units and QC."""

from __future__ import annotations

import argparse
import csv
import importlib.util
from contextlib import ExitStack
from pathlib import Path
from types import SimpleNamespace


def load_fragment_tools(path: Path):
    spec = importlib.util.spec_from_file_location("mstflow_fragment_tools", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import fragment helper: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def rewrite_fastqs(fragment_tools, candidates, output: Path, gzip_level: int):
    """Concatenate existing FASTQs into one validated deterministic gzip."""
    with fragment_tools.open_text(str(output), "w", gzip_level) as destination:
        for candidate in candidates:
            if not candidate.is_file():
                continue
            for record in fragment_tools.fastq_records(str(candidate)):
                fragment_tools.write_record(destination, record)


def collect_stats(
    fragment_tools,
    output: Path,
    *,
    single=(),
    r1=None,
    r2=None,
    kmer_length=35,
):
    fragment_tools.collect_stats(
        SimpleNamespace(
            single=[str(path) for path in single],
            r1=str(r1) if r1 else None,
            r2=str(r2) if r2 else None,
            output=str(output),
            prefix="",
            minimum_usable_span=0,
            kmer_length=kmer_length,
        )
    )


def read_metrics(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return {row["metric"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}


def write_rows(path: Path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument(
        "--mode",
        required=True,
        choices=("adapterremoval_se", "collapsed_only", "all_fragments"),
    )
    parser.add_argument("--kmer-length", type=int, required=True)
    parser.add_argument(
        "--gzip-level",
        type=int,
        default=6,
        choices=range(0, 10),
        metavar="{0..9}",
    )
    parser.add_argument("--fragment-helper", type=Path, required=True)
    args = parser.parse_args()

    if args.kmer_length <= 0:
        parser.error("--kmer-length must be positive")
    if not args.raw_dir.is_dir():
        parser.error(f"--raw-dir is not a directory: {args.raw_dir}")

    fragment_tools = load_fragment_tools(args.fragment_helper)
    sample_id = args.sample_id
    prefix = args.raw_dir / f"{sample_id}.ar"
    qc = Path(f"{sample_id}.adapterremoval.fragment_qc.tsv")

    if args.mode == "adapterremoval_se":
        output = Path(f"{sample_id}.prehost.fq.gz")
        rewrite_fastqs(
            fragment_tools,
            [Path(f"{prefix}.truncated.gz"), Path(f"{prefix}.singleton.truncated.gz")],
            output,
            args.gzip_level,
        )
        collect_stats(fragment_tools, qc, single=[output], kmer_length=args.kmer_length)
        retained = int(read_metrics(qc)["retained_fragments"])
    elif args.mode == "collapsed_only":
        output = Path(f"{sample_id}.prehost.fq.gz")
        rewrite_fastqs(
            fragment_tools,
            [Path(f"{prefix}.collapsed.gz"), Path(f"{prefix}.collapsed.truncated.gz")],
            output,
            args.gzip_level,
        )
        singleton = Path(f"{sample_id}.discarded_singleton.fq.gz")
        pair1 = Path(f"{sample_id}.discarded_R1.fq.gz")
        pair2 = Path(f"{sample_id}.discarded_R2.fq.gz")
        rewrite_fastqs(
            fragment_tools,
            [Path(f"{prefix}.singleton.truncated.gz")],
            singleton,
            args.gzip_level,
        )
        rewrite_fastqs(
            fragment_tools,
            [Path(f"{prefix}.pair1.truncated.gz")],
            pair1,
            args.gzip_level,
        )
        rewrite_fastqs(
            fragment_tools,
            [Path(f"{prefix}.pair2.truncated.gz")],
            pair2,
            args.gzip_level,
        )
        collapsed_stats = Path(f"{sample_id}.collapsed.stats.tsv")
        singleton_stats = Path(f"{sample_id}.singleton.stats.tsv")
        unmerged_stats = Path(f"{sample_id}.unmerged.stats.tsv")
        collect_stats(fragment_tools, collapsed_stats, single=[output], kmer_length=args.kmer_length)
        collect_stats(fragment_tools, singleton_stats, single=[singleton], kmer_length=args.kmer_length)
        collect_stats(
            fragment_tools,
            unmerged_stats,
            r1=pair1,
            r2=pair2,
            kmer_length=args.kmer_length,
        )
        collapsed = int(read_metrics(collapsed_stats)["retained_fragments"])
        singleton_count = int(read_metrics(singleton_stats)["retained_fragments"])
        unmerged = int(read_metrics(unmerged_stats)["retained_fragments"])
        write_rows(
            qc,
            [
                ("ancient_pe_mode", "collapsed_only"),
                ("fastq_gzip_level", args.gzip_level),
                ("retained_collapsed_fragments", collapsed),
                ("discarded_singleton_fragments", singleton_count),
                ("discarded_unmerged_fragments", unmerged),
                ("retained_fragments", collapsed),
            ],
        )
        retained = collapsed
    else:
        output_se = Path(f"{sample_id}.prehost_se.fq.gz")
        output_r1 = Path(f"{sample_id}.prehost_R1.fq.gz")
        output_r2 = Path(f"{sample_id}.prehost_R2.fq.gz")
        rewrite_fastqs(
            fragment_tools,
            [Path(f"{prefix}.collapsed.gz"), Path(f"{prefix}.collapsed.truncated.gz")],
            output_se,
            args.gzip_level,
        )
        rewrite_fastqs(
            fragment_tools,
            [Path(f"{prefix}.pair1.truncated.gz")],
            output_r1,
            args.gzip_level,
        )
        rewrite_fastqs(
            fragment_tools,
            [Path(f"{prefix}.pair2.truncated.gz")],
            output_r2,
            args.gzip_level,
        )
        singleton = Path(f"{sample_id}.discarded_singleton.fq.gz")
        rewrite_fastqs(
            fragment_tools,
            [Path(f"{prefix}.singleton.truncated.gz")],
            singleton,
            args.gzip_level,
        )
        all_stats = Path(f"{sample_id}.all_fragments.stats.tsv")
        singleton_stats = Path(f"{sample_id}.singleton.stats.tsv")
        collect_stats(
            fragment_tools,
            all_stats,
            single=[output_se],
            r1=output_r1,
            r2=output_r2,
            kmer_length=args.kmer_length,
        )
        collect_stats(fragment_tools, singleton_stats, single=[singleton], kmer_length=args.kmer_length)
        all_metrics = read_metrics(all_stats)
        singleton_count = int(read_metrics(singleton_stats)["retained_fragments"])
        rows = [
            ("ancient_pe_mode", "all_fragments"),
            ("fastq_gzip_level", args.gzip_level),
        ]
        rows.extend(all_metrics.items())
        rows.append(("discarded_singleton_fragments", singleton_count))
        write_rows(qc, rows)
        retained = int(all_metrics["retained_fragments"])

    if retained == 0:
        raise SystemExit(f"AdapterRemoval produced no retained molecular fragments for {sample_id}")


if __name__ == "__main__":
    main()
