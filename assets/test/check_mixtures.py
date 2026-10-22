#!/usr/bin/env python3
"""Check SourceTracker2's estimates for the test data's mixtures against their known
composition.

A mixture holds a share of the read pairs of each library in it (mixtures.tsv), and
each library comes from a known environment. SourceTracker2 estimates shares of the
sink's counts, so a library's expected share is its share of read pairs times the
counts its whole library brings to the table, over the same for all of the mixture's
libraries. Fails when any estimate is further than --tolerance from that.

    python3 check_mixtures.py results/
"""
import argparse
import csv
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("results", type=Path)
    parser.add_argument("--label", default="default")
    parser.add_argument("--tolerance", type=float, default=0.10)
    args = parser.parse_args()

    attribution = args.results / "attribution" / args.label
    estimates = {}
    for row in read_tsv(attribution / "sourcetracker2" / "main.sourcetracker2.attribution.tsv"):
        estimates.setdefault(row["sink"], {})[row["source"]] = float(row["estimate"])
    counts = {row["input_sample_id"]: float(row["input_total_count"])
              for row in read_tsv(attribution / "inputs" / "attribution_sample_qc.tsv")}
    mixtures = {}
    for part in read_tsv(HERE / "mixtures.tsv"):
        mixtures.setdefault(part["mixture"], []).append(part)

    failed = False
    print(f"{'mixture':22} {'source':12} {'expected':>8} {'observed':>8}")
    for mixture, parts in mixtures.items():
        brought = {part["library"]: float(part["share"]) * counts[part["library"]] for part in parts}
        expected = {}
        for part in parts:
            expected[part["environment"]] = expected.get(part["environment"], 0.0) + brought[part["library"]] / sum(brought.values())
        for source, observed in estimates[mixture].items():
            off = abs(observed - expected.get(source, 0.0)) > args.tolerance
            failed |= off
            print(f"{mixture:22} {source:12} {expected.get(source, 0.0):8.2f} {observed:8.2f}{'  <- off' if off else ''}")
    raise SystemExit(1 if failed else 0)


if __name__ == "__main__":
    main()
