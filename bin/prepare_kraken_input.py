#!/usr/bin/env python3
"""Prepare an explicit fragment-aware input contract for Kraken2."""

from __future__ import annotations

import argparse
import csv
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


def retained_fragments(path: Path):
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("metric") == "retained_fragments":
                return int(row["value"])
    raise ValueError(f"{path}: missing retained_fragments metric")


def ordered_pair(paths):
    """Return one explicit R1/R2 pair without relying on glob/list ordering."""
    roles = {1: [], 2: []}
    for path in paths:
        matched = [
            role
            for role in (1, 2)
            if re.search(rf"(?:^|[._-])R{role}(?:[._-]|$)", path.name)
        ]
        if len(matched) != 1:
            raise ValueError(f"paired FASTQ has no unambiguous R1/R2 role: {path}")
        roles[matched[0]].append(path)
    if len(roles[1]) != 1 or len(roles[2]) != 1:
        raise ValueError("paired input must contain exactly one R1 and one R2 FASTQ")
    return roles[1][0], roles[2][0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--age-group", choices=("modern", "ancient"), required=True)
    parser.add_argument("--read", action="append", type=Path, required=True)
    parser.add_argument("--minimum-usable-span", type=int, required=True)
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

    if len(args.read) not in (1, 2, 3):
        parser.error("one, two, or three --read values are required")
    if args.minimum_usable_span < args.kmer_length:
        parser.error("--minimum-usable-span must be at least --kmer-length")
    if len(args.read) == 3 and args.age_group != "ancient":
        parser.error("three-file mixed input is only valid for ancient samples")

    fragment_tools = load_fragment_tools(args.fragment_helper)
    sample_id = args.sample_id
    stats = Path(f"{sample_id}.kraken_input.read_length.tsv")
    mode_path = Path(f"{sample_id}.kraken_input.mode.txt")

    if len(args.read) == 3:
        single = [
            path
            for path in args.read
            if ".prehost_se." in path.name or ".dehosted_se." in path.name
        ]
        r1 = [
            path
            for path in args.read
            if ".prehost_R1." in path.name or ".dehosted_R1." in path.name
        ]
        r2 = [
            path
            for path in args.read
            if ".prehost_R2." in path.name or ".dehosted_R2." in path.name
        ]
        if len(single) != 1 or len(r1) != 1 or len(r2) != 1:
            parser.error("ancient three-file input must contain se, R1, and R2 FASTQs")
        fragment_tools.pair_fragments(
            SimpleNamespace(
                single=[str(single[0])],
                r1=str(r1[0]),
                r2=str(r2[0]),
                output_r1=f"{sample_id}.kraken_input_R1.fq.gz",
                output_r2=f"{sample_id}.kraken_input_R2.fq.gz",
                stats=str(stats),
                minimum_usable_span=args.minimum_usable_span,
                kmer_length=args.kmer_length,
                dummy_mate_length=1,
                gzip_level=args.gzip_level,
            )
        )
        mode = "paired"
    elif len(args.read) == 2:
        try:
            r1, r2 = ordered_pair(args.read)
        except ValueError as exc:
            parser.error(str(exc))
        fragment_tools.prepare_paired(
            SimpleNamespace(
                r1=str(r1),
                r2=str(r2),
                output_r1=f"{sample_id}.kraken_input_R1.fq.gz",
                output_r2=f"{sample_id}.kraken_input_R2.fq.gz",
                stats=str(stats),
                minimum_usable_span=args.minimum_usable_span,
                kmer_length=args.kmer_length,
                gzip_level=args.gzip_level,
            )
        )
        mode = "paired"
    else:
        fragment_tools.prepare_single(
            SimpleNamespace(
                single=[str(args.read[0])],
                output=f"{sample_id}.kraken_input.fq.gz",
                stats=str(stats),
                minimum_usable_span=args.minimum_usable_span,
                kmer_length=args.kmer_length,
                gzip_level=args.gzip_level,
            )
        )
        mode = "single"

    if retained_fragments(stats) == 0:
        raise SystemExit(
            f"No {sample_id} fragments contain a usable span of at least "
            f"{args.minimum_usable_span} bp"
        )
    mode_path.write_text(f"{mode}\n", encoding="utf-8")


if __name__ == "__main__":
    main()
