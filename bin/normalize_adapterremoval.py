#!/usr/bin/env python3
"""Turn AdapterRemoval's outputs into MSTFlow fragment FASTQs and QC.

Single-end libraries keep their trimmed reads. Paired-end libraries keep the
collapsed fragments (collapsed_only), or those and the unmerged pairs
(all_fragments); the reads they leave out are only counted.
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
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument(
        "--mode",
        required=True,
        choices=("adapterremoval_se", "collapsed_only", "all_fragments"),
    )
    parser.add_argument("--kmer-length", type=int, required=True)
    parser.add_argument("--gzip-level", type=int, default=6, choices=range(0, 10), metavar="{0..9}")
    parser.add_argument("--fragment-helper", type=Path, required=True)
    args = parser.parse_args()

    if args.kmer_length <= 0:
        parser.error("--kmer-length must be positive")
    if not args.raw_dir.is_dir():
        parser.error(f"--raw-dir is not a directory: {args.raw_dir}")

    tools = load_fragment_tools(args.fragment_helper)
    sample_id = args.sample_id

    def raw(*kinds):
        """AdapterRemoval's output files of these kinds that exist."""
        paths = [args.raw_dir / f"{sample_id}.ar.{kind}.gz" for kind in kinds]
        return [path for path in paths if path.is_file()]

    def rewrite(sources, output):
        """Concatenate FASTQs into one validated, reproducible gzip."""
        with tools.open_text(str(output), "w", args.gzip_level) as destination:
            for source in sources:
                for record in tools.fastq_records(str(source)):
                    tools.write_record(destination, record)
        return output

    def stats(single=(), pair=()):
        r1, r2 = map(str, pair) if pair else (None, None)
        return tools.collect_stats(
            SimpleNamespace(
                single=[str(path) for path in single],
                r1=r1,
                r2=r2,
                minimum_usable_span=0,
                kmer_length=args.kmer_length,
            )
        )

    if args.mode == "adapterremoval_se":
        output = rewrite(raw("truncated", "singleton.truncated"), Path(f"{sample_id}.prehost.fq.gz"))
        qc = stats(single=[output])
    elif args.mode == "collapsed_only":
        output = rewrite(raw("collapsed", "collapsed.truncated"), Path(f"{sample_id}.prehost.fq.gz"))
        collapsed = stats(single=[output])["retained_fragments"]
        qc = {
            "ancient_pe_mode": "collapsed_only",
            "retained_collapsed_fragments": collapsed,
            "discarded_singleton_fragments": stats(single=raw("singleton.truncated"))["retained_fragments"],
            "discarded_unmerged_fragments": stats(pair=raw("pair1.truncated", "pair2.truncated"))[
                "retained_fragments"
            ],
            "retained_fragments": collapsed,
        }
    else:
        single = rewrite(raw("collapsed", "collapsed.truncated"), Path(f"{sample_id}.prehost_se.fq.gz"))
        r1 = rewrite(raw("pair1.truncated"), Path(f"{sample_id}.prehost_R1.fq.gz"))
        r2 = rewrite(raw("pair2.truncated"), Path(f"{sample_id}.prehost_R2.fq.gz"))
        qc = {
            "ancient_pe_mode": "all_fragments",
            **stats(single=[single], pair=[r1, r2]),
            "discarded_singleton_fragments": stats(single=raw("singleton.truncated"))["retained_fragments"],
        }

    tools.write_metrics(f"{sample_id}.adapterremoval.fragment_qc.tsv", qc)
    if qc["retained_fragments"] == 0:
        raise SystemExit(f"AdapterRemoval produced no retained molecular fragments for {sample_id}")


if __name__ == "__main__":
    main()
