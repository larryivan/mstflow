#!/usr/bin/env python3
"""Make leave-one-group-out folds for testing source attribution on known sources.

Sources that share a value of the grouping column (by default the study they came
from) are held out together and turned into sinks whose environment is known. Every
method is then run on each fold exactly as on the real sinks, so its accuracy is
measured on data it has not seen, from studies it has not seen.

Holding out whole studies rather than single libraries matters: libraries of one
study share laboratory, protocol and site, and a method that has seen a study's
other libraries is judged too kindly.

Writes <fold>.counts.tsv and <fold>.map.tsv for each fold, and a summary of what each
fold holds out. Held-out units keep their true environment in ``expected_env``.
"""
import argparse
import csv
import re
from collections import OrderedDict
from pathlib import Path

from mstflow_attribution import is_na, read_counts, read_map


def fold_name(value, taken):
    base = "cv_" + (re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("_") or "group")
    name, suffix = base, 2
    while name in taken:
        name, suffix = f"{base}_{suffix}", suffix + 1
    taken.add(name)
    return name


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--map", required=True)
    parser.add_argument("--group-column", required=True,
                        help="Map column whose value groups sources; units with no value form their own group.")
    parser.add_argument("--min-sink-count", type=int, default=0,
                        help="Held-out units with fewer counts are skipped, as real sinks would be.")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    units = read_map(args.map)
    taxa, columns, counts = read_counts(args.counts)
    if set(columns) != set(units):
        raise SystemExit("the attribution count table and map list different units")
    fields = list(next(iter(units.values())).keys())
    if args.group_column not in fields:
        raise SystemExit(
            f"cross-validation groups by '{args.group_column}', which is not a samplesheet column; "
            "add it, or set attribution_cv_group to another column"
        )

    sources = [unit for unit, row in units.items() if row["SourceSink"] == "source"]
    groups = OrderedDict()
    for unit in sources:
        value = units[unit][args.group_column]
        groups.setdefault(unit if is_na(value) else value.strip(), []).append(unit)
    if len(groups) < 2:
        raise SystemExit(
            f"cross-validation needs sources from at least two '{args.group_column}' groups; found {len(groups)}"
        )

    totals = {unit: sum(counts[taxon][unit] for taxon in taxa) for unit in sources}
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    taken = set()
    summary = []
    for group, held_out in sorted(groups.items()):
        name = fold_name(group, taken)
        training = [unit for unit in sources if unit not in held_out]
        training_envs = {units[unit]["Env"] for unit in training}
        tested = [unit for unit in held_out if totals[unit] >= args.min_sink_count]
        skipped = [unit for unit in held_out if unit not in tested]
        if not tested:
            summary.append([name, group, ",".join(held_out), "", "", "skipped:held_out_units_below_min_sink_count"])
            continue

        kept = training + tested
        kept_taxa = [taxon for taxon in taxa if any(counts[taxon][unit] for unit in kept)]
        with open(out_dir / f"{name}.counts.tsv", "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["#OTU ID", *kept])
            for taxon in kept_taxa:
                writer.writerow([taxon, *[counts[taxon][unit] for unit in kept]])
        with open(out_dir / f"{name}.map.tsv", "w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=[*fields, "expected_env"],
                                    lineterminator="\n")
            writer.writeheader()
            for unit in kept:
                row = dict(units[unit])
                if unit in tested:
                    row.update({"SourceSink": "sink", "Env": "NA", "expected_env": units[unit]["Env"]})
                else:
                    row["expected_env"] = "NA"
                writer.writerow(row)

        novel = sorted({units[unit]["Env"] for unit in tested} - training_envs)
        note = "ok" if not skipped else "skipped_units_below_min_sink_count:" + ",".join(skipped)
        summary.append([name, group, ",".join(tested), ",".join(sorted(training_envs)), ",".join(novel), note])

    with open(args.summary, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["fold", "group", "held_out_units", "training_envs", "envs_absent_from_training", "note"])
        writer.writerows(summary)


if __name__ == "__main__":
    main()
