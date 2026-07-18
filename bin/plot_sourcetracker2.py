#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
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
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if "#SampleID" not in fieldnames:
            raise SystemExit(f"{path} is missing mapping column: #SampleID")
        if "SourceSink" not in fieldnames:
            raise SystemExit(f"{path} is missing mapping column: SourceSink")
        rows = list(reader)
    mapping = {}
    for line_no, row in enumerate(rows, start=2):
        sample = (row.get("#SampleID") or "").strip()
        if not sample:
            raise SystemExit(f"{path}:{line_no}: #SampleID is empty")
        if sample in mapping:
            raise SystemExit(f"{path}:{line_no}: duplicate #SampleID '{sample}'")
        role = (row.get("SourceSink") or "").strip().lower()
        if role not in {"source", "sink"}:
            raise SystemExit(
                f"{path}:{line_no}: SourceSink must be exactly source or sink"
            )
        row["SourceSink"] = role
        mapping[sample] = row
    return mapping


def numeric_proportion(value, context):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise SystemExit(f"{context} is not numeric: {value!r}")
    if not math.isfinite(number) or number < 0:
        raise SystemExit(f"{context} must be finite and non-negative: {value!r}")
    if number > 1.000001:
        raise SystemExit(f"{context} exceeds 1.0: {number}")
    return number


def read_result_table(path, label):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if not fieldnames or fieldnames[0] != "SampleID":
            raise SystemExit(f"{path} must have SampleID as its first column")
        duplicates = sorted(
            {field for field in fieldnames if fieldnames.count(field) > 1}
        )
        if duplicates:
            raise SystemExit(
                f"{path} contains duplicate result column(s): {', '.join(duplicates)}"
            )
        sources = fieldnames[1:]
        if not sources:
            raise SystemExit(f"{path} contains no source columns")
        rows = []
        seen = set()
        for line_no, row in enumerate(reader, start=2):
            sample = (row.get("SampleID") or "").strip()
            if not sample:
                raise SystemExit(f"{path}:{line_no}: SampleID is empty")
            if sample in seen:
                raise SystemExit(f"{path}:{line_no}: duplicate SampleID '{sample}'")
            seen.add(sample)
            parsed = {"SampleID": sample}
            for source in sources:
                parsed[source] = numeric_proportion(
                    row.get(source), f"{path}:{line_no}:{source} ({label})"
                )
            rows.append(parsed)
    return sources, rows


def read_proportions(path):
    """Backward-compatible public helper returning parsed source results."""
    return read_result_table(path, "mean")


def uncertainty_by_sample(path, sources, rows):
    if not path:
        return {}
    uncertainty_sources, uncertainty_rows = read_result_table(path, "standard deviation")
    if uncertainty_sources != sources:
        raise SystemExit(
            "uncertainty source columns do not exactly match mixing-proportion columns"
        )
    expected_samples = [row["SampleID"] for row in rows]
    observed_samples = [row["SampleID"] for row in uncertainty_rows]
    if set(expected_samples) != set(observed_samples) or len(expected_samples) != len(observed_samples):
        raise SystemExit("uncertainty samples do not exactly match mixing-proportion samples")
    return {row["SampleID"]: row for row in uncertainty_rows}


def validate_result_targets(rows, mapping, loo):
    expected_role = "source" if loo else "sink"
    expected = {
        sample
        for sample, metadata in mapping.items()
        if metadata["SourceSink"] == expected_role
    }
    observed = {row["SampleID"] for row in rows}
    missing = sorted(expected - observed)
    unexpected = sorted(observed - expected)
    if missing or unexpected:
        details = []
        if missing:
            details.append(f"missing {expected_role} result(s): {', '.join(missing)}")
        if unexpected:
            details.append(f"unexpected non-{expected_role} result(s): {', '.join(unexpected)}")
        raise SystemExit(
            "SourceTracker2 results do not exactly match expected mapping targets; "
            + "; ".join(details)
        )


def validate_proportion_sums(sources, rows, tolerance):
    unknown_columns = [source for source in sources if source.casefold() == "unknown"]
    if len(unknown_columns) != 1:
        raise SystemExit(
            "mixing proportions must contain exactly one case-insensitive Unknown column"
        )
    for row in rows:
        total = sum(row[source] for source in sources)
        if abs(1.0 - total) > tolerance:
            raise SystemExit(
                f"mixing proportions for {row['SampleID']} sum to {total:.12g}; "
                f"allowed absolute error is {tolerance:.12g}"
            )


def format_number(value):
    return "NA" if value is None else f"{value:.12g}"


def write_assignments(
    path,
    sources,
    rows,
    mapping,
    uncertainty=None,
    max_unknown_prop=0.5,
    min_top_known_prop=0.5,
    min_top_known_margin=0.05,
):
    uncertainty = uncertainty or {}
    known = [source for source in sources if source.casefold() != "unknown"]
    fields = [
        "SampleID",
        "expected_env",
        "top_known_env",
        "top_known_prop",
        "second_known_env",
        "second_known_prop",
        "top_known_margin",
        "top_all_env",
        "top_all_prop",
        "assignment",
        "abstained",
        "abstention_reason",
        "top_known_correct",
        "unknown_prop",
        "proportion_sum",
        "proportion_sum_abs_error",
        "top_known_std",
        "top_all_std",
        "unknown_std",
        "max_source_std",
        "uncertainty_available",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            sample = row["SampleID"]
            if sample not in mapping:
                raise SystemExit(f"mixing proportions contain sample absent from mapping: {sample}")
            expected = (mapping[sample].get("Env") or "NA").strip() or "NA"
            known_ranked = sorted(known, key=lambda source: (-row[source], source))
            all_ranked = sorted(sources, key=lambda source: (-row[source], source))
            top_known = known_ranked[0] if known_ranked else "NA"
            second_known = known_ranked[1] if len(known_ranked) > 1 else "NA"
            top_all = all_ranked[0] if all_ranked else "NA"
            top_known_prop = row.get(top_known) if top_known != "NA" else None
            second_known_prop = row.get(second_known) if second_known != "NA" else None
            top_known_margin = (
                top_known_prop - second_known_prop
                if top_known_prop is not None and second_known_prop is not None
                else None
            )
            top_all_prop = row.get(top_all) if top_all != "NA" else None
            unknown_name = next(
                (source for source in sources if source.casefold() == "unknown"), None
            )
            stds = uncertainty.get(sample, {})
            proportion_sum = sum(row[source] for source in sources)
            max_std = max((stds.get(source, 0.0) for source in sources), default=None) if stds else None
            unknown_prop = row.get(unknown_name) if unknown_name else 0.0
            abstention_reasons = []
            if unknown_prop >= max_unknown_prop:
                abstention_reasons.append(
                    f"unknown_prop>={max_unknown_prop:.12g}"
                )
            if top_known_prop is None or top_known_prop < min_top_known_prop:
                abstention_reasons.append(
                    f"top_known_prop<{min_top_known_prop:.12g}"
                )
            if top_known_margin is not None and top_known_margin <= 0:
                abstention_reasons.append("top_known_tie")
            elif top_known_margin is not None and top_known_margin < min_top_known_margin:
                abstention_reasons.append(
                    f"top_known_margin<{min_top_known_margin:.12g}"
                )
            abstained = bool(abstention_reasons)
            assignment = "ABSTAIN" if abstained else top_known
            writer.writerow(
                {
                    "SampleID": sample,
                    "expected_env": expected,
                    "top_known_env": top_known,
                    "top_known_prop": format_number(top_known_prop),
                    "second_known_env": second_known,
                    "second_known_prop": format_number(second_known_prop),
                    "top_known_margin": format_number(top_known_margin),
                    "top_all_env": top_all,
                    "top_all_prop": format_number(top_all_prop),
                    "assignment": assignment,
                    "abstained": str(abstained).lower(),
                    "abstention_reason": ";".join(abstention_reasons) or "NA",
                    "top_known_correct": (
                        str(top_known == expected).lower()
                        if not abstained
                        and expected.upper() not in {"", "NA", "N/A", "NONE", "NULL"}
                        else "NA"
                    ),
                    "unknown_prop": format_number(row.get(unknown_name) if unknown_name else None),
                    "proportion_sum": format_number(proportion_sum),
                    "proportion_sum_abs_error": format_number(abs(1.0 - proportion_sum)),
                    "top_known_std": format_number(stds.get(top_known) if stds and top_known != "NA" else None),
                    "top_all_std": format_number(stds.get(top_all) if stds and top_all != "NA" else None),
                    "unknown_std": format_number(stds.get(unknown_name) if stds and unknown_name else None),
                    "max_source_std": format_number(max_std),
                    "uncertainty_available": str(bool(stds)).lower(),
                }
            )


def source_colors(sources):
    cmap = plt.get_cmap("tab20")
    fallback = {
        source: matplotlib.colors.to_hex(cmap(index % cmap.N))
        for index, source in enumerate(sorted(source for source in sources if source not in PALETTE))
    }
    return {source: PALETTE.get(source, fallback.get(source, "#999999")) for source in sources}


def plot(path, sources, rows, mapping, title, uncertainty=None):
    uncertainty = uncertainty or {}
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        fig, ax = plt.subplots(figsize=(9, 4.5))
        ax.axis("off")
        ax.text(
            0.5,
            0.5,
            "No sink mixing proportions were produced.",
            ha="center",
            va="center",
            fontsize=14,
        )
        ax.set_title(title, loc="left", fontsize=15, fontweight="bold")
        fig.tight_layout()
        fig.savefig(output, dpi=220)
        plt.close(fig)
        if not output.is_file() or output.stat().st_size == 0:
            raise SystemExit(f"failed to create a non-empty plot: {output}")
        return

    rows = sorted(
        rows,
        key=lambda row: (mapping.get(row["SampleID"], {}).get("Env", "ZZZ"), row["SampleID"]),
    )
    samples = [row["SampleID"] for row in rows]
    height = max(4.5, 0.34 * len(samples) + 2.8)
    fig, ax = plt.subplots(figsize=(12.5, height))
    left = [0.0] * len(rows)
    colors = source_colors(sources)
    for source in sources:
        values = [row[source] for row in rows]
        ax.barh(
            samples,
            values,
            left=left,
            color=colors[source],
            edgecolor="white",
            linewidth=0.4,
            label=source,
        )
        if uncertainty:
            centers = [start + value / 2 for start, value in zip(left, values)]
            errors = [uncertainty.get(row["SampleID"], {}).get(source, 0.0) for row in rows]
            ax.errorbar(
                centers,
                range(len(rows)),
                xerr=errors,
                fmt="none",
                ecolor="#222222",
                elinewidth=0.55,
                capsize=1.4,
                alpha=0.65,
            )
        left = [a + b for a, b in zip(left, values)]

    for index, row in enumerate(rows):
        expected = mapping.get(row["SampleID"], {}).get("Env", "NA")
        ax.scatter(
            -0.025,
            index,
            s=55,
            color=colors.get(expected, "#999999"),
            edgecolor="black",
            linewidth=0.8,
            clip_on=False,
        )

    max_sum = max(sum(row[source] for source in sources) for row in rows)
    ax.set_xlim(-0.05, max(1.0, max_sum * 1.02))
    ax.set_xlabel("Mixing proportion (error bars: SourceTracker standard deviation)")
    ax.set_title(title, loc="left", fontsize=15, fontweight="bold")
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
    ax.grid(axis="x", alpha=0.25)
    ax.set_axisbelow(True)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=3, frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=220)
    plt.close(fig)
    if not output.is_file() or output.stat().st_size == 0:
        raise SystemExit(f"failed to create a non-empty plot: {output}")


def main():
    parser = argparse.ArgumentParser(
        description="Validate and plot SourceTracker2 mean and uncertainty outputs."
    )
    parser.add_argument("--proportions", required=True)
    parser.add_argument("--uncertainty")
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--assignments", required=True)
    parser.add_argument("--plot", required=True)
    parser.add_argument("--title", default="SourceTracker2 source attribution")
    parser.add_argument("--max-unknown-prop", type=float, default=0.5)
    parser.add_argument("--min-top-known-prop", type=float, default=0.5)
    parser.add_argument("--min-top-known-margin", type=float, default=0.05)
    parser.add_argument("--proportion-sum-tolerance", type=float, default=0.001)
    parser.add_argument("--loo", action="store_true")
    args = parser.parse_args()

    if not 0 <= args.max_unknown_prop <= 1:
        raise SystemExit("--max-unknown-prop must be between 0 and 1")
    if not 0 <= args.min_top_known_prop <= 1:
        raise SystemExit("--min-top-known-prop must be between 0 and 1")
    if not 0 <= args.min_top_known_margin <= 1:
        raise SystemExit("--min-top-known-margin must be between 0 and 1")
    if not 0 <= args.proportion_sum_tolerance <= 0.1:
        raise SystemExit("--proportion-sum-tolerance must be between 0 and 0.1")

    mapping = read_map(args.mapping)
    sources, rows = read_proportions(args.proportions)
    validate_result_targets(rows, mapping, args.loo)
    validate_proportion_sums(sources, rows, args.proportion_sum_tolerance)
    uncertainty = uncertainty_by_sample(args.uncertainty, sources, rows)
    write_assignments(
        args.assignments,
        sources,
        rows,
        mapping,
        uncertainty,
        max_unknown_prop=args.max_unknown_prop,
        min_top_known_prop=args.min_top_known_prop,
        min_top_known_margin=args.min_top_known_margin,
    )
    plot(args.plot, sources, rows, mapping, args.title, uncertainty)


if __name__ == "__main__":
    main()
