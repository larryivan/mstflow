"""Shared tables for the source attribution steps.

Every method reads the same two inputs and writes the same long result table, so the
summary, the plot and the cross-validation never need to know which method produced
a number.

attribution_counts.tsv
    ``#OTU ID`` followed by one column per analysis unit; one row per taxon.
attribution_map.tsv
    ``#SampleID``, ``SourceSink`` (source or sink) and ``Env`` (the source environment;
    ``NA`` for sinks), then any other samplesheet columns. Cross-validation folds add
    ``expected_env``, the environment a held-out source really came from.
<method>.attribution.tsv
    One row per sink and source environment, ``Unknown`` included:
    analysis, method, sink, source, estimate, sd, quantity. ``quantity`` is
    ``proportion`` for mixture models (estimates are shares of the sink) and
    ``probability`` for classifiers (estimates are class probabilities). ``sd`` is
    ``NA`` when a method reports no uncertainty.

Only the standard library is used: this runs in every method's container, including
SourceTracker2's Python 3.8.
"""

import csv
import math
from collections import OrderedDict

UNKNOWN = "Unknown"
RESULT_FIELDS = ["analysis", "method", "sink", "source", "estimate", "sd", "quantity"]
NA_VALUES = {"", "NA", "N/A", "NONE", "NULL"}
METHOD_NAMES = OrderedDict([
    ("sourcetracker2", "SourceTracker2"),
    ("feast", "FEAST"),
    ("fastst", "FastST"),
    ("sourceid_nmf", "SourceID-NMF"),
    ("sourcepredict", "SourcePredict"),
])


def method_key(method):
    """Sort methods in a fixed order, unknown ones last."""
    order = list(METHOD_NAMES)
    return (order.index(method) if method in order else len(order), method)


def is_na(value):
    return value is None or str(value).strip().upper() in NA_VALUES


def read_map(path):
    """Return the analysis units in file order, keyed by #SampleID."""
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = {"#SampleID", "SourceSink", "Env"} - set(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"{path} is missing column(s): {', '.join(sorted(missing))}")
        units = OrderedDict()
        for line_no, row in enumerate(reader, start=2):
            unit = (row["#SampleID"] or "").strip()
            role = (row["SourceSink"] or "").strip().lower()
            if not unit:
                raise SystemExit(f"{path}:{line_no}: #SampleID is empty")
            if unit in units:
                raise SystemExit(f"{path}:{line_no}: duplicate #SampleID '{unit}'")
            if role not in ("source", "sink"):
                raise SystemExit(f"{path}:{line_no}: SourceSink must be source or sink, not '{role}'")
            if role == "source" and is_na(row["Env"]):
                raise SystemExit(f"{path}:{line_no}: source '{unit}' has no Env")
            row["SourceSink"] = role
            units[unit] = row
    return units


def read_counts(path):
    """Return (taxa, units, counts) with counts[taxon][unit] as non-negative integers."""
    with open(path, newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or not rows[0] or rows[0][0] != "#OTU ID":
        raise SystemExit(f"{path} must start with a '#OTU ID' header")
    units = rows[0][1:]
    if len(set(units)) != len(units):
        raise SystemExit(f"{path} has duplicate sample columns")
    taxa, counts = [], OrderedDict()
    for line_no, row in enumerate(rows[1:], start=2):
        if len(row) != len(units) + 1:
            raise SystemExit(f"{path}:{line_no}: expected {len(units) + 1} columns, found {len(row)}")
        values = {}
        for unit, value in zip(units, row[1:]):
            number = float(value)
            if not math.isfinite(number) or number < 0 or number != int(number):
                raise SystemExit(f"{path}:{line_no}: count for {unit} is not a non-negative integer: {value}")
            values[unit] = int(number)
        taxa.append(row[0])
        counts[row[0]] = values
    return taxa, units, counts


def roles(units):
    """Split a map into source environments (sorted) and sink unit names (file order)."""
    envs = sorted({row["Env"] for row in units.values() if row["SourceSink"] == "source"})
    sinks = [unit for unit, row in units.items() if row["SourceSink"] == "sink"]
    if not envs:
        raise SystemExit("the attribution map has no sources")
    reserved = [env for env in envs if env.casefold() in ("unknown", "sampleid")]
    if reserved:
        raise SystemExit("source environments must not be called Unknown or SampleID: " + ", ".join(reserved))
    return envs, sinks


def pooled_sources(units, taxa, counts):
    """Sum the counts of all sources of each environment, as SourceTracker2 does."""
    envs, _ = roles(units)
    pooled = OrderedDict((env, [0] * len(taxa)) for env in envs)
    for unit, row in units.items():
        if row["SourceSink"] == "source":
            column = pooled[row["Env"]]
            for index, taxon in enumerate(taxa):
                column[index] += counts[taxon][unit]
    return pooled


def write_results(path, analysis, method, quantity, estimates, tolerance=1e-3):
    """Write {sink: {source: (estimate, sd or None)}} as a standard result table.

    Every sink must cover the same sources, Unknown included, and its estimates must
    sum to one within ``tolerance``; they are then rescaled to sum to exactly one.
    """
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(RESULT_FIELDS)
        expected_sources = None
        for sink, by_source in estimates.items():
            sources = sorted(by_source, key=lambda source: (source == UNKNOWN, source))
            if UNKNOWN not in by_source:
                raise SystemExit(f"{method}: sink {sink} has no {UNKNOWN} estimate")
            if expected_sources is None:
                expected_sources = sources
            elif sources != expected_sources:
                raise SystemExit(f"{method}: sink {sink} covers different sources from the others")
            total = 0.0
            for source in sources:
                estimate = by_source[source][0]
                if not math.isfinite(estimate) or estimate < -1e-9:
                    raise SystemExit(f"{method}: {sink}/{source} estimate is not a finite share: {estimate}")
                total += max(estimate, 0.0)
            if abs(total - 1.0) > tolerance:
                raise SystemExit(f"{method}: estimates for {sink} sum to {total:.6g}, not 1")
            for source in sources:
                estimate, sd = by_source[source]
                writer.writerow([
                    analysis, method, sink, source,
                    format_number(max(estimate, 0.0) / total),
                    "NA" if sd is None or not math.isfinite(sd) else format_number(sd / total),
                    quantity,
                ])


def decide(estimates, max_unknown, min_top_known, min_margin):
    """Name the source of one sink, or abstain.

    ``estimates`` maps each source, Unknown included, to its estimate. A sink is
    assigned to its largest known source only when the unknown share stays below
    ``max_unknown``, that source reaches ``min_top_known``, and it leads the runner-up
    by at least ``min_margin``; otherwise the result is ABSTAIN, with every reason.
    """
    known = sorted((source for source in estimates if source != UNKNOWN),
                   key=lambda source: (-estimates[source], source))
    top = known[0] if known else None
    second = known[1] if len(known) > 1 else None
    top_value = estimates[top] if top else 0.0
    margin = top_value - estimates[second] if second else None
    unknown = estimates.get(UNKNOWN, 0.0)
    reasons = []
    if unknown >= max_unknown:
        reasons.append(f"unknown>={max_unknown:g}")
    if top_value < min_top_known:
        reasons.append(f"top_known<{min_top_known:g}")
    if margin is not None and margin <= 0:
        reasons.append("tie")
    elif margin is not None and margin < min_margin:
        reasons.append(f"margin<{min_margin:g}")
    return {
        "top_known": top or "NA",
        "top_known_estimate": top_value,
        "second_known": second or "NA",
        "second_known_estimate": estimates[second] if second else None,
        "margin": margin,
        "unknown": unknown,
        "call": "ABSTAIN" if reasons else top,
        "reasons": ";".join(reasons) or "NA",
    }


def group_results(rows):
    """{(analysis, method, sink): {"quantity": q, "estimates": {...}, "sd": {...}}} in input order."""
    grouped = OrderedDict()
    for row in rows:
        key = (row["analysis"], row["method"], row["sink"])
        entry = grouped.setdefault(key, {"quantity": row["quantity"], "estimates": OrderedDict(), "sd": {}})
        entry["estimates"][row["source"]] = row["estimate"]
        entry["sd"][row["source"]] = row["sd"]
    return grouped


def read_results(paths):
    rows = []
    for path in paths:
        with open(path, newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames != RESULT_FIELDS:
                raise SystemExit(f"{path} is not an attribution result table")
            for row in reader:
                row["estimate"] = float(row["estimate"])
                row["sd"] = None if is_na(row["sd"]) else float(row["sd"])
                rows.append(row)
    return rows


def format_number(value):
    return f"{value:.6g}"
