#!/usr/bin/env python3
"""Convert SourceTracker2 mixing proportions and their standard deviations into the
standard attribution table (see mstflow_attribution.py)."""
import argparse
import csv
import math
from collections import OrderedDict

from mstflow_attribution import UNKNOWN, read_map, roles, write_results


def read_table(path):
    with open(path, newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or rows[0][0] != "SampleID":
        raise SystemExit(f"{path} is not a SourceTracker2 result table")
    sources = rows[0][1:]
    duplicates = sorted({source for source in sources if sources.count(source) > 1})
    if duplicates:
        raise SystemExit(f"{path} has duplicate columns: {', '.join(duplicates)}")
    table = OrderedDict()
    for row in rows[1:]:
        values = [float(value) for value in row[1:]]
        if len(values) != len(sources) or not all(math.isfinite(value) for value in values):
            raise SystemExit(f"{path}: malformed row for {row[0]}")
        table[row[0]] = dict(zip(sources, values))
    return sources, table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--proportions", required=True)
    parser.add_argument("--stds", required=True)
    parser.add_argument("--map", required=True)
    parser.add_argument("--analysis", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    envs, sinks = roles(read_map(args.map))
    sources, means = read_table(args.proportions)
    std_sources, stds = read_table(args.stds)
    if sorted(sources) != sorted([*envs, UNKNOWN]) or std_sources != sources:
        raise SystemExit(f"SourceTracker2 reported sources {sources}, expected {[*envs, UNKNOWN]}")
    if sorted(means) != sorted(sinks) or sorted(stds) != sorted(sinks):
        raise SystemExit("SourceTracker2 results do not cover exactly the sinks in the map")

    estimates = OrderedDict(
        (sink, {source: (means[sink][source], stds[sink][source]) for source in sources}) for sink in sinks
    )
    write_results(args.output, args.analysis, "sourcetracker2", "proportion", estimates)


if __name__ == "__main__":
    main()
