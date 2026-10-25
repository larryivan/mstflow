#!/usr/bin/env python3
"""Combine the attribution results of every method for the real sinks.

Writes
  proportions.tsv   every estimate, one row per method, sink and source
  assignments.tsv   per method and sink: the leading sources, the unknown share and the
                    call (an environment, or ABSTAIN with the reasons)
  comparison.tsv    per sink: each method's call, and whether the methods agree
  attribution.png   one panel per method, one bar per sink
"""
import argparse
import csv
from collections import OrderedDict

from mstflow_attribution import (METHOD_NAMES, RESULT_FIELDS, UNKNOWN, decide, group_results, method_key,
                                 read_map, read_results)

# Okabe-Ito: distinguishable with every common form of colour blindness.
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#000000"]
UNKNOWN_COLOUR = "#D6DCE4"


def fmt(value):
    return "NA" if value is None else f"{value:.6g}"


def colours(sources, plt):
    known = [source for source in sources if source != UNKNOWN]
    if len(known) > len(PALETTE):
        cmap = plt.get_cmap("tab20")
        palette = [plt.matplotlib.colors.to_hex(cmap(index % 20)) for index in range(len(known))]
    else:
        palette = PALETTE
    mapping = {source: palette[index] for index, source in enumerate(known)}
    mapping[UNKNOWN] = UNKNOWN_COLOUR
    return mapping


def plot(path, grouped, decisions, methods, sinks, sources, title):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    colour = colours(sources, plt)
    width = 3.0 * len(methods) + 1.6
    height = max(2.2, 0.32 * len(sinks) + 1.4)
    figure, axes = plt.subplots(1, len(methods), figsize=(width, height), sharey=True, squeeze=False)
    for axis, method in zip(axes[0], methods):
        quantity = next(entry["quantity"] for (_, m, _), entry in grouped.items() if m == method)
        for row, sink in enumerate(sinks):
            entry = grouped[("main", method, sink)]
            left = 0.0
            for source in sources:
                value = entry["estimates"].get(source, 0.0)
                axis.barh(row, value, left=left, color=colour[source], height=0.72, linewidth=0)
                left += value
            if decisions[(method, sink)]["call"] == "ABSTAIN":
                axis.text(1.02, row, "abstain", va="center", ha="left", fontsize=7, color="#6B7280")
        axis.set_title(f"{METHOD_NAMES.get(method, method)}\n({quantity})", fontsize=9)
        axis.set_xlim(0, 1)
        axis.set_xticks([0, 0.5, 1])
        axis.tick_params(labelsize=8)
        for side in ("top", "right"):
            axis.spines[side].set_visible(False)
    axes[0][0].set_yticks(range(len(sinks)))
    axes[0][0].set_yticklabels(sinks, fontsize=8)
    axes[0][0].invert_yaxis()
    handles = [Patch(color=colour[source], label=source) for source in sources]
    figure.legend(handles=handles, loc="lower center", ncol=min(len(handles), 6), frameon=False, fontsize=8,
                  bbox_to_anchor=(0.5, 0))
    figure.suptitle(title, fontsize=10)
    figure.tight_layout(rect=(0, 0.1, 0.97, 0.95))
    figure.savefig(path, dpi=200)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--map", required=True)
    parser.add_argument("--max-unknown", type=float, required=True)
    parser.add_argument("--min-top-known", type=float, required=True)
    parser.add_argument("--min-margin", type=float, required=True)
    parser.add_argument("--title", default="Source attribution")
    parser.add_argument("--no-plot", action="store_true", help="Write the tables only")
    parser.add_argument("results", nargs="+")
    args = parser.parse_args()

    units = read_map(args.map)
    sinks = [unit for unit, row in units.items() if row["SourceSink"] == "sink"]
    rows = read_results(args.results)
    grouped = group_results(rows)
    methods = sorted({method for _, method, _ in grouped}, key=method_key)
    sources = sorted({source for entry in grouped.values() for source in entry["estimates"]},
                     key=lambda source: (source == UNKNOWN, source))
    missing = [(method, sink) for method in methods for sink in sinks if ("main", method, sink) not in grouped]
    if missing:
        raise SystemExit("results are missing for: " + ", ".join(f"{m}/{s}" for m, s in missing))

    with open("proportions.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(RESULT_FIELDS)
        for row in sorted(rows, key=lambda r: (method_key(r["method"]), sinks.index(r["sink"]))):
            writer.writerow([row["analysis"], row["method"], row["sink"], row["source"],
                             fmt(row["estimate"]), fmt(row["sd"]), row["quantity"]])

    decisions = {}
    with open("assignments.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sink", "method", "quantity", "call", "top_known", "top_known_estimate", "top_known_sd",
                         "second_known", "second_known_estimate", "margin", "unknown", "unknown_sd",
                         "abstention_reasons"])
        for method in methods:
            for sink in sinks:
                entry = grouped[("main", method, sink)]
                decision = decide(entry["estimates"], args.max_unknown, args.min_top_known, args.min_margin)
                decisions[(method, sink)] = decision
                writer.writerow([
                    sink, method, entry["quantity"], decision["call"], decision["top_known"],
                    fmt(decision["top_known_estimate"]), fmt(entry["sd"].get(decision["top_known"])),
                    decision["second_known"], fmt(decision["second_known_estimate"]), fmt(decision["margin"]),
                    fmt(decision["unknown"]), fmt(entry["sd"].get(UNKNOWN)), decision["reasons"],
                ])

    with open("comparison.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sink", *methods, "methods_with_a_call", "agreement"])
        for sink in sinks:
            calls = OrderedDict((method, decisions[(method, sink)]["call"]) for method in methods)
            named = {call for call in calls.values() if call != "ABSTAIN"}
            called = sum(call != "ABSTAIN" for call in calls.values())
            if not named:
                agreement = "all_abstain"
            elif len(named) > 1:
                agreement = "disagree"
            elif called == 1 and len(methods) > 1:
                agreement = "one_method"
            else:
                agreement = "agree:" + next(iter(named))
            writer.writerow([sink, *calls.values(), called, agreement])

    if not args.no_plot:
        plot("attribution.png", grouped, decisions, methods, sinks, sources, args.title)


if __name__ == "__main__":
    main()
