#!/usr/bin/env python3
"""Run SourcePredict (Borry et al. 2019, JOSS) and write the standard attribution
table (see mstflow_attribution.py).

SourcePredict is a classifier, not a mixture model: it embeds the sources and each
sink with a UniFrac distance on the NCBI taxonomy and asks a k-nearest-neighbour
classifier which environment the sink belongs to, after estimating how likely it is
to belong to none of them. Its estimates are therefore class probabilities
(quantity "probability"), not shares of the sink. Unlike the mixture models it
learns from every source library separately, so environments need several each.
"""
import argparse
import csv
import re
import subprocess
from collections import Counter, OrderedDict

from mstflow_attribution import UNKNOWN, read_counts, read_map, roles, write_results

RANKS = {"D": "superkingdom", "P": "phylum", "C": "class", "O": "order",
         "F": "family", "G": "genus", "S": "species"}
# SourcePredict's automatic choice of neighbours searches 10, 20 and 50 with
# five-fold cross-validation on 64% of the sources, and t-SNE's perplexity of 30 must
# stay below the number of samples, so both need a large reference panel.
MIN_SOURCES_FOR_NEIGHBOUR_SEARCH = 100
MIN_SOURCES_FOR_TSNE = 60


def taxid(feature):
    match = re.match(r"taxid_(\d+)__", feature)
    return int(match.group(1)) if match else None


def write_table(path, rows, units, counts):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["TAXID", *units])
        for feature, tid in rows:
            writer.writerow([tid, *[counts[feature][unit] for unit in units]])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--map", required=True)
    parser.add_argument("--analysis", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--rank", required=True, choices=sorted(RANKS), help="Rank of the count table")
    parser.add_argument("--neighbors", type=int, default=0, help="0 lets SourcePredict choose, if there are enough sources")
    parser.add_argument("--unknown-alpha", type=float, default=0.1)
    # GMPR, SourcePredict's own default, silently drops a sink that shares too few
    # taxa with the sources and then fails; subsampling has no such case.
    parser.add_argument("--normalization", default="SUBSAMPLE", choices=["SUBSAMPLE", "GMPR", "RLE"])
    parser.add_argument("--distance", default="weighted_unifrac", choices=["weighted_unifrac", "unweighted_unifrac"])
    parser.add_argument("--embedding", default="auto", choices=["auto", "TSNE", "UMAP", "MDS"],
                        help="auto: t-SNE for large panels, MDS for small ones")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--drop-lonely-envs", action="store_true",
                        help="Leave out environments with one source library instead of failing (cross-validation)")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()

    units = read_map(args.map)
    features, _, counts = read_counts(args.counts)
    envs, sinks = roles(units)
    sources = [unit for unit, row in units.items() if row["SourceSink"] == "source"]
    per_env = Counter(units[unit]["Env"] for unit in sources)
    lonely = sorted(env for env, n in per_env.items() if n < 2)
    if lonely and not args.drop_lonely_envs:
        raise SystemExit("SourcePredict needs at least two source libraries per environment; one only: "
                         + ", ".join(lonely))
    if lonely:
        print("SourcePredict cannot learn an environment from one library; left out: " + ", ".join(lonely))
        sources = [unit for unit in sources if units[unit]["Env"] not in lonely]
        envs = [env for env in envs if env not in lonely]
        per_env = Counter(units[unit]["Env"] for unit in sources)
        if len(per_env) < 2:
            raise SystemExit("fewer than two environments have two or more source libraries")

    rows = [(feature, taxid(feature)) for feature in features if taxid(feature) is not None]
    write_table("sources.csv", rows, sources, counts)
    write_table("sinks.csv", rows, sinks, counts)
    with open("labels.csv", "w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["sample", "labels"])
        writer.writerows([[unit, units[unit]["Env"]] for unit in sources])

    neighbors = args.neighbors
    if neighbors == 0 and len(sources) < MIN_SOURCES_FOR_NEIGHBOUR_SEARCH:
        neighbors = max(1, min(5, int(len(sources) * 0.64) - 1))
        print(f"{len(sources)} sources are too few for SourcePredict's neighbour search; using {neighbors}")
    embedding = args.embedding
    if embedding == "auto":
        embedding = "TSNE" if len(sources) + 1 >= MIN_SOURCES_FOR_TSNE else "MDS"
    kfold = max(2, min(5, min(per_env.values())))
    subprocess.run([
        "sourcepredict", "sinks.csv",
        "-s", "sources.csv", "-l", "labels.csv",
        "-r", RANKS[args.rank], "-dt", args.distance, "-n", args.normalization,
        "-me", embedding, "-kne", str(neighbors), "-k", str(kfold),
        "-a", str(args.unknown_alpha), "-se", str(args.seed), "-t", str(args.threads),
        "-o", "sourcepredict.csv",
    ], check=True)

    # One column per sink, one row per environment plus "unknown".
    with open("sourcepredict.csv", newline="") as handle:
        table = list(csv.reader(handle))
    columns = table[0][1:]
    probabilities = {row[0]: [float(value) for value in row[1:]] for row in table[1:]}
    if sorted(columns) != sorted(sinks):
        raise SystemExit("SourcePredict returned results for different sinks")
    estimates = OrderedDict()
    for sink in sinks:
        index = columns.index(sink)
        shares = {env: (probabilities.get(env, [0.0] * len(columns))[index], None) for env in envs}
        shares[UNKNOWN] = (probabilities.get("unknown", [0.0] * len(columns))[index], None)
        estimates[sink] = shares
    write_results(args.output, args.analysis, "sourcepredict", "probability", estimates)


if __name__ == "__main__":
    main()
