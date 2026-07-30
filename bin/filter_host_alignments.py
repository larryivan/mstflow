#!/usr/bin/env python3
"""Filter a BWA SAM bundle and emit pair-aware dehosted MSTFlow reads."""

from __future__ import annotations

import argparse
import csv
import importlib.util
from pathlib import Path
from types import SimpleNamespace


def load_fragment_tools(path: Path):
    spec = importlib.util.spec_from_file_location("mstflow_fragment_tools", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import fragment helper: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def filter_sam(
    fragment_tools,
    sam,
    stats,
    minimum_mapq,
    gzip_level,
    *,
    output=None,
    r1=None,
    r2=None,
):
    fragment_tools.filter_host_sam(
        SimpleNamespace(
            sam=str(sam),
            output=str(output) if output else None,
            output_r1=str(r1) if r1 else None,
            output_r2=str(r2) if r2 else None,
            minimum_mapq=minimum_mapq,
            stats=str(stats),
            gzip_level=gzip_level,
        )
    )


def metrics(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        return {row["metric"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}


def write_metrics(path: Path, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--alignment-dir", type=Path, required=True)
    parser.add_argument("--minimum-mapq", type=int, required=True)
    parser.add_argument(
        "--gzip-level",
        type=int,
        default=6,
        choices=range(0, 10),
        metavar="{0..9}",
    )
    parser.add_argument("--fragment-helper", type=Path, required=True)
    args = parser.parse_args()

    if args.minimum_mapq < 0:
        parser.error("--minimum-mapq must be non-negative")

    mode_path = args.alignment_dir / "mode.txt"
    try:
        mode = mode_path.read_text(encoding="utf-8").strip()
    except OSError as exc:
        parser.error(str(exc))
    allowed = {"modern_se", "modern_pe", "ancient_se", "ancient_pe", "ancient_mixed"}
    if mode not in allowed:
        parser.error(f"invalid alignment mode in {mode_path}: {mode!r}")

    fragment_tools = load_fragment_tools(args.fragment_helper)
    sample_id = args.sample_id
    final_qc = Path(f"{sample_id}.host_removal.fragment_qc.tsv")

    if mode == "ancient_mixed":
        se_stats = Path(f"{sample_id}.host_se.stats.tsv")
        pair_stats = Path(f"{sample_id}.host_pair.stats.tsv")
        filter_sam(
            fragment_tools,
            args.alignment_dir / "se.sam.gz",
            se_stats,
            args.minimum_mapq,
            args.gzip_level,
            output=Path(f"{sample_id}.dehosted_se.fq.gz"),
        )
        filter_sam(
            fragment_tools,
            args.alignment_dir / "pair.sam.gz",
            pair_stats,
            args.minimum_mapq,
            args.gzip_level,
            r1=Path(f"{sample_id}.dehosted_R1.fq.gz"),
            r2=Path(f"{sample_id}.dehosted_R2.fq.gz"),
        )
        se_values = metrics(se_stats)
        pair_values = metrics(pair_stats)
        input_fragments = int(se_values["input_fragments"]) + int(pair_values["input_fragments"])
        retained_fragments = int(se_values["retained_fragments"]) + int(
            pair_values["retained_fragments"]
        )
        write_metrics(
            final_qc,
            [
                ("host_min_mapq", args.minimum_mapq),
                ("fastq_gzip_level", args.gzip_level),
                ("input_fragments", input_fragments),
                ("retained_fragments", retained_fragments),
                ("removed_host_fragments", input_fragments - retained_fragments),
                ("retained_single_fragments", se_values["retained_fragments"]),
                ("retained_paired_fragments", pair_values["retained_fragments"]),
            ],
        )
    elif mode in {"ancient_pe", "modern_pe"}:
        filter_sam(
            fragment_tools,
            args.alignment_dir / "pair.sam.gz",
            final_qc,
            args.minimum_mapq,
            args.gzip_level,
            r1=Path(f"{sample_id}.dehosted_R1.fq.gz"),
            r2=Path(f"{sample_id}.dehosted_R2.fq.gz"),
        )
    else:
        filter_sam(
            fragment_tools,
            args.alignment_dir / "single.sam.gz",
            final_qc,
            args.minimum_mapq,
            args.gzip_level,
            output=Path(f"{sample_id}.dehosted.fq.gz"),
        )


if __name__ == "__main__":
    main()
