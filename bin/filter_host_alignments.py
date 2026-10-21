#!/usr/bin/env python3
"""Keep the fragments of one library without a confident host alignment.

A pair is removed whole when either mate aligns. HOST_ALIGN_BWA records in mode.txt
which SAM files it wrote: single.sam.gz, pair.sam.gz, or se.sam.gz and pair.sam.gz
for an ancient library's collapsed fragments and unmerged pairs.
"""

from __future__ import annotations

import argparse
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--alignment-dir", type=Path, required=True)
    parser.add_argument("--minimum-mapq", type=int, required=True)
    parser.add_argument("--gzip-level", type=int, default=6, choices=range(0, 10), metavar="{0..9}")
    parser.add_argument("--fragment-helper", type=Path, required=True)
    args = parser.parse_args()

    if args.minimum_mapq < 0:
        parser.error("--minimum-mapq must be non-negative")

    mode_path = args.alignment_dir / "mode.txt"
    try:
        mode = mode_path.read_text(encoding="utf-8").strip()
    except OSError as exc:
        parser.error(str(exc))
    if mode not in {"modern_se", "modern_pe", "ancient_se", "ancient_mixed"}:
        parser.error(f"invalid alignment mode in {mode_path}: {mode!r}")

    tools = load_fragment_tools(args.fragment_helper)
    sample_id = args.sample_id

    def filter_sam(name, output=None, r1=None, r2=None):
        return tools.filter_host_sam(
            SimpleNamespace(
                sam=str(args.alignment_dir / name),
                output=output,
                output_r1=r1,
                output_r2=r2,
                minimum_mapq=args.minimum_mapq,
                gzip_level=args.gzip_level,
            )
        )

    pair_outputs = dict(r1=f"{sample_id}.dehosted_R1.fq.gz", r2=f"{sample_id}.dehosted_R2.fq.gz")
    if mode == "ancient_mixed":
        single = filter_sam("se.sam.gz", output=f"{sample_id}.dehosted_se.fq.gz")
        paired = filter_sam("pair.sam.gz", **pair_outputs)
        qc = {
            "host_min_mapq": args.minimum_mapq,
            **{
                key: single[key] + paired[key]
                for key in ("input_fragments", "retained_fragments", "removed_host_fragments")
            },
            "retained_single_fragments": single["retained_fragments"],
            "retained_paired_fragments": paired["retained_fragments"],
        }
    elif mode == "modern_pe":
        qc = filter_sam("pair.sam.gz", **pair_outputs)
    else:
        qc = filter_sam("single.sam.gz", output=f"{sample_id}.dehosted.fq.gz")
    tools.write_metrics(f"{sample_id}.host_removal.fragment_qc.tsv", qc)


if __name__ == "__main__":
    main()
