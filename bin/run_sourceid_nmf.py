#!/usr/bin/env python3
"""Run SourceID-NMF (Huang et al. 2024, Bioinformatics) and write the standard
attribution table (see mstflow_attribution.py).

SourceID-NMF factorises the sink into the source profiles and one unknown profile
with a weighted non-negative matrix factorisation. Sources of one environment are
pooled first, as for the other methods. Its code is not redistributed with MSTFlow;
it is staged from the authors' repository at a pinned commit (--code-dir).
"""
import argparse
import csv
import os
import subprocess
import sys
from collections import OrderedDict
from pathlib import Path

from mstflow_attribution import UNKNOWN, pooled_sources, read_counts, read_map, roles, write_results


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--map", required=True)
    parser.add_argument("--analysis", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--code-dir", required=True, help="Directory holding SourceID-NMF.py and its modules")
    parser.add_argument("--iterations", type=int, default=2000)
    parser.add_argument("--rho", type=int, default=1)
    parser.add_argument("--weight", type=int, default=1)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    units = read_map(args.map)
    taxa, _, counts = read_counts(args.counts)
    envs, sinks = roles(units)
    pooled = pooled_sources(units, taxa, counts)
    # SourceID-NMF drops all-zero sources without dropping their labels, which would
    # shift every later column; the pooled sources are never empty, but say so loudly.
    empty = [env for env, column in pooled.items() if not any(column)]
    if empty:
        raise SystemExit("pooled sources with no counts: " + ", ".join(empty))

    # Neutral column names keep environment names with spaces or symbols out of
    # SourceID-NMF's space-separated output.
    source_ids = OrderedDict((f"S{index}", env) for index, env in enumerate(envs, start=1))
    sink_ids = OrderedDict((f"K{index}", sink) for index, sink in enumerate(sinks, start=1))
    with open("nmf_data.txt", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxon", *source_ids, *sink_ids])
        for index, taxon in enumerate(taxa):
            writer.writerow([f"t{index}",
                             *[pooled[env][index] for env in source_ids.values()],
                             *[counts[taxon][sink] for sink in sink_ids.values()]])
    with open("name.txt", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["SampleID", "Env", "SourceSink"])
        writer.writerows([[sid, env, "Source"] for sid, env in source_ids.items()])
        writer.writerows([[kid, "sink", "Sink"] for kid in sink_ids])

    command = [
        sys.executable, str(Path(args.code_dir) / "SourceID-NMF.py"),
        "-i", "nmf_data.txt", "-n", "name.txt", "-o", "nmf_proportions.txt",
        "-t", str(args.threads), "-e", str(args.iterations), "-r", str(args.rho),
        "-a", str(args.weight), "-c", str(args.tolerance),
    ]
    # The code files may be links into separate download directories, and Python
    # looks for a script's modules beside the file the link points to; name the
    # directory explicitly so they are found.
    code_dir = str(Path(args.code_dir).absolute())
    subprocess.run(command, check=True, env={**os.environ, "PYTHONPATH": code_dir})

    with open("nmf_proportions.txt", newline="") as handle:
        rows = [line.split() for line in handle if line.strip()]
    header, body = rows[0], rows[1:]
    if header != [*source_ids, "unknown"] or [row[0] for row in body] != list(sink_ids):
        raise SystemExit("SourceID-NMF returned columns or rows in an unexpected order")
    estimates = OrderedDict()
    for row in body:
        values = [float(value) for value in row[1:]]
        shares = {env: (value, None) for env, value in zip(source_ids.values(), values)}
        shares[UNKNOWN] = (values[-1], None)
        estimates[sink_ids[row[0]]] = shares
    write_results(args.output, args.analysis, "sourceid_nmf", "proportion", estimates)


if __name__ == "__main__":
    main()
