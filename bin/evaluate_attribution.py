#!/usr/bin/env python3
"""Score every method on the cross-validation folds made by split_attribution_folds.py.

Each held-out source is a sink whose environment is known. It is called with the same
rule as the real sinks, then scored:

  coverage            share of held-out sources the method names an environment for
  accuracy            share of those calls that are right
  balanced_accuracy   accuracy averaged over environments, so a large environment
                      cannot hide a small one that is always wrong
  correct_overall     right calls over all held-out sources (an abstention counts as
                      not right); the single number to compare methods by

A held-out source whose environment has no other source in its fold cannot be named
correctly; abstaining is the right answer there, and it is scored separately.

Writes evaluation.tsv (one row per method), evaluation_by_env.tsv,
predictions.tsv (one row per method and held-out source) and evaluation.png.
"""
import argparse
import csv
from collections import OrderedDict, defaultdict
from pathlib import Path

from mstflow_attribution import METHOD_NAMES, decide, group_results, is_na, method_key, read_map, read_results


def ratio(numerator, denominator):
    return f"{numerator / denominator:.4g}" if denominator else "NA"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--folds", required=True, help="fold summary from split_attribution_folds.py")
    parser.add_argument("--maps", nargs="+", required=True, help="<fold>.map.tsv files")
    parser.add_argument("--results", nargs="+", required=True)
    parser.add_argument("--max-unknown", type=float, required=True)
    parser.add_argument("--min-top-known", type=float, required=True)
    parser.add_argument("--min-margin", type=float, required=True)
    parser.add_argument("--no-plot", action="store_true", help="Write the tables only")
    args = parser.parse_args()

    with open(args.folds, newline="") as handle:
        folds = {row["fold"]: row for row in csv.DictReader(handle, delimiter="\t")}
    truth = {}
    for path in args.maps:
        fold = Path(path).name[: -len(".map.tsv")]
        absent = set(filter(None, folds[fold]["envs_absent_from_training"].split(",")))
        for unit, row in read_map(path).items():
            if row["SourceSink"] == "sink" and not is_na(row.get("expected_env")):
                truth[(fold, unit)] = (row["expected_env"], row["expected_env"] in absent)

    grouped = group_results(read_results(args.results))
    methods = sorted({method for _, method, _ in grouped}, key=method_key)
    predictions = []
    for (fold, method, unit), entry in grouped.items():
        if (fold, unit) not in truth:
            raise SystemExit(f"{method} reported {unit} in {fold}, which was not held out there")
        expected, novel = truth[(fold, unit)]
        decision = decide(entry["estimates"], args.max_unknown, args.min_top_known, args.min_margin)
        called = decision["call"] != "ABSTAIN"
        correct = (not called) if novel else (decision["call"] == expected)
        predictions.append({
            "method": method, "fold": fold, "held_out": unit, "expected_env": expected,
            "env_in_training": str(not novel).lower(), "call": decision["call"],
            "top_known": decision["top_known"], "top_known_estimate": f"{decision['top_known_estimate']:.4g}",
            "unknown": f"{decision['unknown']:.4g}", "correct": str(correct).lower(),
        })
    expected_pairs = {(method, fold, unit) for method in methods for fold, unit in truth}
    reported = {(p["method"], p["fold"], p["held_out"]) for p in predictions}
    if expected_pairs - reported:
        missing = sorted(expected_pairs - reported)[:5]
        raise SystemExit("cross-validation results are missing, e.g. " + ", ".join("/".join(m) for m in missing))
    predictions.sort(key=lambda p: (method_key(p["method"]), p["fold"], p["held_out"]))

    with open("predictions.tsv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(predictions[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(predictions)

    summary, by_env = OrderedDict(), []
    for method in methods:
        rows = [p for p in predictions if p["method"] == method and p["env_in_training"] == "true"]
        novel = [p for p in predictions if p["method"] == method and p["env_in_training"] == "false"]
        called = [p for p in rows if p["call"] != "ABSTAIN"]
        right = [p for p in called if p["correct"] == "true"]
        recalls = []
        for env in sorted({p["expected_env"] for p in rows}):
            env_rows = [p for p in rows if p["expected_env"] == env]
            env_called = [p for p in env_rows if p["call"] != "ABSTAIN"]
            env_right = [p for p in env_called if p["correct"] == "true"]
            recalls.append(len(env_right) / len(env_rows))
            by_env.append([method, env, len(env_rows), len(env_called), len(env_right),
                           ratio(len(env_right), len(env_called)), ratio(len(env_right), len(env_rows))])
        summary[method] = [
            method, len(rows), len(called), ratio(len(called), len(rows)), ratio(len(right), len(called)),
            f"{sum(recalls) / len(recalls):.4g}" if recalls else "NA", ratio(len(right), len(rows)),
            len(novel), ratio(sum(p["correct"] == "true" for p in novel), len(novel)),
        ]

    with open("evaluation.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["method", "held_out", "called", "coverage", "accuracy", "balanced_accuracy",
                         "correct_overall", "held_out_unseen_env", "abstained_on_unseen_env"])
        writer.writerows(summary.values())
    with open("evaluation_by_env.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["method", "expected_env", "held_out", "called", "correct", "accuracy", "recall"])
        writer.writerows(by_env)

    if not args.no_plot:
        plot(summary, methods)


def plot(summary, methods):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axis = plt.subplots(figsize=(1.1 * len(methods) + 2.4, 2.8))
    positions = range(len(methods))
    for offset, column, label, colour in ((-0.2, 3, "coverage", "#9CA3AF"), (0.2, 6, "correct overall", "#0072B2")):
        values = [float(summary[m][column]) if summary[m][column] != "NA" else 0.0 for m in methods]
        axis.bar([p + offset for p in positions], values, width=0.38, color=colour, label=label)
    axis.set_xticks(list(positions))
    axis.set_xticklabels([METHOD_NAMES.get(m, m) for m in methods], fontsize=8)
    axis.set_ylim(0, 1)
    axis.set_ylabel("share of held-out sources", fontsize=8)
    axis.tick_params(axis="y", labelsize=8)
    axis.legend(frameon=False, fontsize=8, loc="upper right")
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    axis.set_title("Cross-validation on held-out sources", fontsize=9)
    figure.tight_layout()
    figure.savefig("evaluation.png", dpi=200)


if __name__ == "__main__":
    main()
