#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt


PALETTE = {
    "human_fecal_gastrointestinal": "#2F7DD3",
    "canid_paleofeces": "#20A46B",
    "archaeological_sediment": "#C89211",
    "terrestrial_archaeological_sediment": "#9C6B00",
    "aquatic_marine_lake_sediment": "#C89211",
    "submerged_archaeological_sediment": "#C89211",
    "extraction_library_blank": "#E0493E",
    "human_dental_calculus": "#D85C4A",
    "human_bone_tooth": "#8C63D2",
    "nonhuman_mammal_bone": "#6B7788",
    "Unknown": "#D6DCE4",
}


def read_map(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {row["#SampleID"]: row for row in rows}


def read_proportions(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        sources = [field for field in reader.fieldnames if field != "SampleID"]
    return sources, rows


def write_assignments(path, sources, rows, mapping):
    known = [source for source in sources if source != "Unknown"]
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "SampleID",
                "expected_env",
                "top_known_env",
                "top_known_prop",
                "top_all_env",
                "top_all_prop",
                "top_known_correct",
                "unknown_prop",
            ]
        )
        for row in rows:
            sample = row["SampleID"]
            expected = mapping.get(sample, {}).get("Env", "NA")
            top_known = max(known, key=lambda source: float(row.get(source, 0) or 0)) if known else "NA"
            top_all = max(sources, key=lambda source: float(row.get(source, 0) or 0)) if sources else "NA"
            writer.writerow(
                [
                    sample,
                    expected,
                    top_known,
                    f"{float(row.get(top_known, 0) or 0):.12g}",
                    top_all,
                    f"{float(row.get(top_all, 0) or 0):.12g}",
                    str(top_known == expected).lower(),
                    f"{float(row.get('Unknown', 0) or 0):.12g}",
                ]
            )


def plot(path, sources, rows, mapping, title):
    if not rows:
        Path(path).touch()
        return
    rows = sorted(rows, key=lambda row: (mapping.get(row["SampleID"], {}).get("Env", "ZZZ"), row["SampleID"]))
    samples = [row["SampleID"] for row in rows]
    height = max(4.5, 0.34 * len(samples) + 2.8)
    fig, ax = plt.subplots(figsize=(12.5, height))
    left = [0.0] * len(rows)
    for source in sources:
        values = [float(row.get(source, 0) or 0) for row in rows]
        ax.barh(samples, values, left=left, color=PALETTE.get(source, "#999999"), edgecolor="white", linewidth=0.4, label=source)
        left = [a + b for a, b in zip(left, values)]

    for index, row in enumerate(rows):
        expected = mapping.get(row["SampleID"], {}).get("Env", "NA")
        ax.scatter(-0.025, index, s=55, color=PALETTE.get(expected, "#999999"), edgecolor="black", linewidth=0.8, clip_on=False)

    ax.set_xlim(-0.05, 1.0)
    ax.set_xlabel("Mixing proportion")
    ax.set_title(title, loc="left", fontsize=15, fontweight="bold")
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax.grid(axis="x", alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=3, frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot SourceTracker2 mixing proportions.")
    parser.add_argument("--proportions", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--assignments", required=True)
    parser.add_argument("--plot", required=True)
    parser.add_argument("--title", default="SourceTracker2 source attribution")
    args = parser.parse_args()

    mapping = read_map(args.mapping)
    sources, rows = read_proportions(args.proportions)
    write_assignments(args.assignments, sources, rows, mapping)
    plot(args.plot, sources, rows, mapping, args.title)


if __name__ == "__main__":
    main()
