#!/usr/bin/env python3
"""Prepare one library's reads for Kraken2 or KrakenUniq.

One FASTQ is classified single-end. Two (R1, R2) or three (collapsed se, R1, R2)
are written as one native pair, so every fragment stays one classification record.
Fragments without a usable span of at least --minimum-usable-span are dropped.
"""

from __future__ import annotations

import argparse
import importlib.util
import re
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
    parser.add_argument("--read", action="append", type=Path, required=True)
    parser.add_argument("--minimum-usable-span", type=int, required=True)
    parser.add_argument("--kmer-length", type=int, required=True)
    parser.add_argument("--gzip-level", type=int, default=6, choices=range(0, 10), metavar="{0..9}")
    parser.add_argument("--fragment-helper", type=Path, required=True)
    args = parser.parse_args()

    if args.minimum_usable_span < args.kmer_length:
        parser.error("--minimum-usable-span must be at least --kmer-length")

    fragment_tools = load_fragment_tools(args.fragment_helper)
    sample_id = args.sample_id
    common = dict(
        minimum_usable_span=args.minimum_usable_span,
        kmer_length=args.kmer_length,
        gzip_level=args.gzip_level,
    )
    if len(args.read) == 1:
        metrics = fragment_tools.prepare_single(
            SimpleNamespace(single=[str(args.read[0])], output=f"{sample_id}.kraken_input.fq.gz", **common)
        )
    else:
        # Earlier steps name their files <sample>.<step>_{se,R1,R2}.fq.gz.
        parts = {}
        for path in args.read:
            match = re.search(r"_(se|R1|R2)\.fq\.gz$", path.name)
            if not match or match.group(1) in parts:
                parser.error(f"expected one each of *_se, *_R1 and *_R2 FASTQs, got {path.name}")
            parts[match.group(1)] = str(path)
        if "R1" not in parts or "R2" not in parts:
            parser.error("paired input needs an *_R1 and an *_R2 FASTQ")
        metrics = fragment_tools.pair_fragments(
            SimpleNamespace(
                single=[parts["se"]] if "se" in parts else [],
                r1=parts["R1"],
                r2=parts["R2"],
                output_r1=f"{sample_id}.kraken_input_R1.fq.gz",
                output_r2=f"{sample_id}.kraken_input_R2.fq.gz",
                dummy_mate_length=1,
                **common,
            )
        )
    fragment_tools.write_metrics(f"{sample_id}.kraken_input.read_length.tsv", metrics)
    if metrics["retained_fragments"] == 0:
        raise SystemExit(
            f"No {sample_id} fragments contain a usable span of at least {args.minimum_usable_span} bp"
        )


if __name__ == "__main__":
    main()
