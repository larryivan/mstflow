#!/usr/bin/env python3
"""Select, aggregate and depth-check the analysis units for source attribution.

Writes the count table and map that every attribution method reads (see
mstflow_attribution.py), and records for every input sample whether it was used and,
if not, why.
"""
import argparse
import csv
from collections import OrderedDict


NA_VALUES = {"", "NA", "N/A", "NONE", "NULL"}
# Runs merged into one analysis unit must agree on these.
CONSISTENT_AGGREGATION_FIELDS = ("source_sink", "model_source_env")
QC_FIELDS = [
    "input_sample_id",
    "analysis_sample_id",
    "unit_id",
    "source_sink",
    "model_source_env",
    "kraken_classified_fraction",
    "bracken_total_reads",
    "taxon_filter_retained_reads",
    "taxon_filter_retained_fraction",
    "input_total_count",
    "analysis_total_count",
    "required_min_count",
    "aggregation_n",
    "aggregation_members",
    "included",
    "exclusion_reason",
]


def norm(value):
    return "" if value is None else str(value).strip()


def is_na(value):
    return norm(value).upper() in NA_VALUES


def first_present(row, column):
    value = norm(row.get(column))
    return "" if is_na(value) else value


def read_count(value, context):
    """A read count. The tables read here were written by earlier MSTFlow steps as
    plain integers, so anything else means a damaged input."""
    try:
        return int(value)
    except (TypeError, ValueError):
        raise SystemExit(f"{context} is not a whole read count: {value!r}")


def read_metadata(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if "sample_id" not in fieldnames:
            raise SystemExit(f"{path} is missing metadata column: sample_id")
        rows = list(reader)

    by_sample = OrderedDict()
    for line_no, row in enumerate(rows, start=2):
        sample_id = norm(row.get("sample_id"))
        if not sample_id:
            raise SystemExit(f"{path}:{line_no}: sample_id is empty")
        if sample_id in by_sample:
            raise SystemExit(f"{path}:{line_no}: duplicate sample_id '{sample_id}'")
        by_sample[sample_id] = row
    if not by_sample:
        raise SystemExit(f"metadata contains no samples: {path}")
    return fieldnames, by_sample


def read_counts(path):
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise SystemExit(f"count table is empty: {path}")
        if not header or header[0] != "#OTU ID":
            raise SystemExit(
                f"{path} must be a tab-delimited feature table with first column '#OTU ID'"
            )
        samples = [norm(value) for value in header[1:]]
        if any(not sample for sample in samples):
            raise SystemExit(f"{path}: sample column names must not be empty")
        duplicates = sorted({sample for sample in samples if samples.count(sample) > 1})
        if duplicates:
            raise SystemExit(f"{path}: duplicate sample column(s): {', '.join(duplicates)}")

        features = []
        matrix = {}
        for line_no, row in enumerate(reader, start=2):
            if not row:
                continue
            if len(row) != len(header):
                raise SystemExit(
                    f"{path}:{line_no}: expected {len(header)} columns, found {len(row)}"
                )
            feature = norm(row[0])
            if not feature:
                raise SystemExit(f"{path}:{line_no}: feature ID is empty")
            if feature in matrix:
                raise SystemExit(f"{path}:{line_no}: duplicate feature ID '{feature}'")
            features.append(feature)
            matrix[feature] = {
                sample: read_count(value, f"{path}:{line_no}: feature '{feature}', sample '{sample}'")
                for sample, value in zip(samples, row[1:])
            }
    return samples, features, matrix


def read_feature_sample_qc(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        required = {
            "sample_id",
            "kraken_total_reads",
            "kraken_classified_reads",
            "bracken_total_reads",
            "taxon_filter_retained_reads",
        }
        missing = sorted(required.difference(fields))
        if missing:
            raise SystemExit(f"{path} is missing feature-QC column(s): {', '.join(missing)}")
        rows = list(reader)
    by_sample = {}
    for line_no, row in enumerate(rows, start=2):
        sample = norm(row.get("sample_id"))
        if not sample:
            raise SystemExit(f"{path}:{line_no}: sample_id is empty")
        if sample in by_sample:
            raise SystemExit(f"{path}:{line_no}: duplicate sample_id '{sample}'")
        by_sample[sample] = row
    return by_sample


def complete_sum(rows, field):
    """The sum of an integer QC field over the runs of a unit, or None if any run lacks it."""
    values = [norm(row.get(field)) for row in rows]
    if any(is_na(value) for value in values):
        return None
    return sum(read_count(value, f"feature QC field {field} for sample '{row['sample_id']}'")
               for row, value in zip(rows, values))


def ratio(numerator, denominator):
    return numerator / denominator if numerator is not None and denominator else None


def aggregate_feature_qc(members, feature_qc):
    rows = [feature_qc[sample] for sample in members]
    kraken_total = complete_sum(rows, "kraken_total_reads")
    kraken_classified = complete_sum(rows, "kraken_classified_reads")
    bracken_total = complete_sum(rows, "bracken_total_reads")
    retained_reads = complete_sum(rows, "taxon_filter_retained_reads")
    return {
        "kraken_classified_fraction": ratio(kraken_classified, kraken_total),
        "bracken_total_reads": bracken_total,
        "taxon_filter_retained_reads": retained_reads,
        "taxon_filter_retained_fraction": ratio(retained_reads, bracken_total),
    }


def quality_exclusion_reason(qc, taxon_scope, min_classified, min_retained_reads, min_retained_fraction):
    classified = qc["kraken_classified_fraction"]
    bracken_reads = qc["bracken_total_reads"]
    retained_reads = qc["taxon_filter_retained_reads"]
    retained_fraction = qc["taxon_filter_retained_fraction"]
    if classified is None:
        return "missing_kraken_classified_fraction"
    if bracken_reads is None:
        return "missing_bracken_total_reads"
    if taxon_scope == "filtered" and retained_reads is None:
        return "missing_taxon_filter_retained_reads"
    if taxon_scope == "filtered" and retained_fraction is None:
        return "missing_taxon_filter_retained_fraction"
    if classified < min_classified:
        return f"below_min_kraken_classified_fraction:{min_classified:.12g}"
    if taxon_scope == "all":
        if bracken_reads < min_retained_reads:
            return f"below_min_bracken_total_reads:{min_retained_reads:.12g}"
    else:
        if retained_reads < min_retained_reads:
            return f"below_min_taxon_filter_retained_reads:{min_retained_reads:.12g}"
        if retained_fraction < min_retained_fraction:
            return f"below_min_taxon_filter_retained_fraction:{min_retained_fraction:.12g}"
    return ""


def require_one_to_one(metadata_samples, count_samples):
    metadata_set = set(metadata_samples)
    count_set = set(count_samples)
    missing_metadata = sorted(count_set.difference(metadata_set))
    missing_counts = sorted(metadata_set.difference(count_set))
    problems = []
    if missing_metadata:
        problems.append(
            "count table sample(s) absent from metadata: " + ", ".join(missing_metadata)
        )
    if missing_counts:
        problems.append(
            "metadata sample(s) absent from count table: " + ", ".join(missing_counts)
        )
    if problems:
        raise SystemExit("; ".join(problems))


def check_group_consistency(group_id, members, metadata):
    for field in CONSISTENT_AGGREGATION_FIELDS:
        values = {first_present(metadata[sample], field) or "<MISSING>" for sample in members}
        if len(values) > 1:
            raise SystemExit(
                f"aggregation group '{group_id}' has conflicting {field}: {', '.join(sorted(values))}"
            )


def collapse_values(rows, field):
    values = []
    for row in rows:
        value = first_present(row, field)
        if value and value not in values:
            values.append(value)
    return ";".join(values) or "NA"


def format_number(value, count=False):
    """NA for an undefined value, such as a fraction of zero reads; counts exactly."""
    if value is None:
        return "NA"
    return str(value) if count else f"{value:.12g}"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--feature-sample-qc", required=True)
    parser.add_argument("--output-counts", required=True)
    parser.add_argument("--output-map", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--sample-qc", required=True)
    parser.add_argument("--aggregate-by", choices=["sample_id", "unit_id"], default="sample_id")
    parser.add_argument("--min-source-count", type=int, default=0,
                        help="Exclude source units with fewer counts, e.g. below the SourceTracker2 rarefaction depth.")
    parser.add_argument("--min-sink-count", type=int, default=0,
                        help="Exclude sink units with fewer counts.")
    parser.add_argument("--taxon-scope", choices=["all", "filtered"], default="filtered")
    parser.add_argument("--min-kraken-classified-fraction", type=float, default=0.001)
    parser.add_argument("--min-taxon-filter-retained-reads", type=float, default=1000.0)
    parser.add_argument("--min-taxon-filter-retained-fraction", type=float, default=0.001)
    args = parser.parse_args()

    metadata_fields, metadata = read_metadata(args.metadata)
    samples, features, matrix = read_counts(args.counts)
    require_one_to_one(metadata, samples)
    feature_qc = read_feature_sample_qc(args.feature_sample_qc)
    require_one_to_one(feature_qc, samples)
    if args.aggregate_by not in metadata_fields:
        raise SystemExit(f"--aggregate-by column is absent from metadata: {args.aggregate_by}")

    input_totals = {sample: sum(matrix[feature][sample] for feature in features) for sample in samples}
    group_of = {sample: norm(metadata[sample][args.aggregate_by]) for sample in samples}
    groups = OrderedDict()
    for sample in samples:
        groups.setdefault(group_of[sample], []).append(sample)

    reasons = {}
    group_metadata = {}
    group_feature_qc = {}
    final_groups = []
    for group_id, members in groups.items():
        # Validate the complete biological unit before applying any exclusion.
        # Otherwise one conflicting or low-quality run could be removed first
        # and leave a scientifically mislabelled partial unit behind.
        check_group_consistency(group_id, members, metadata)
        rows = [metadata[sample] for sample in members]
        group_metadata[group_id] = {field: collapse_values(rows, field) for field in metadata_fields}
        group_metadata[group_id].update({"sample_id": group_id, args.aggregate_by: group_id})
        group_feature_qc[group_id] = aggregate_feature_qc(members, feature_qc)
        total = sum(input_totals[sample] for sample in members)
        role = group_metadata[group_id]["source_sink"]
        required_depth = args.min_source_count if role == "source" else args.min_sink_count
        if role not in ("source", "sink"):
            group_reason = "source_sink_not_selected"
        else:
            group_reason = quality_exclusion_reason(
                group_feature_qc[group_id],
                args.taxon_scope,
                args.min_kraken_classified_fraction,
                args.min_taxon_filter_retained_reads,
                args.min_taxon_filter_retained_fraction,
            )
        if not group_reason and total <= 0:
            group_reason = "zero_total_count"
        elif not group_reason and total < required_depth:
            group_reason = f"below_min_{role}_count:{required_depth}"
        if group_reason:
            for sample in members:
                reasons[sample] = group_reason
        else:
            final_groups.append(group_id)

    if not final_groups:
        raise SystemExit("no samples remain for source attribution after role, exclusion, and depth QC")
    source_groups = [
        group for group in final_groups if group_metadata[group]["source_sink"] == "source"
    ]
    sink_groups = [
        group for group in final_groups if group_metadata[group]["source_sink"] == "sink"
    ]
    if not source_groups:
        raise SystemExit("source attribution needs at least one source after exclusions and depth QC")
    if not sink_groups:
        raise SystemExit("source attribution needs at least one sink after exclusions and depth QC")
    reserved_source_envs = [
        group_metadata[group]["model_source_env"]
        for group in source_groups
        if group_metadata[group]["model_source_env"].casefold() in {"unknown", "sampleid"}
    ]
    if reserved_source_envs:
        raise SystemExit(
            "source environment names must not use the reserved result columns "
            "Unknown or SampleID; found: " + ", ".join(sorted(set(reserved_source_envs)))
        )

    aggregated = {}
    for feature in features:
        aggregated[feature] = {
            group: sum(matrix[feature][sample] for sample in groups[group])
            for group in final_groups
        }

    retained_features = [
        feature
        for feature in features
        if any(aggregated[feature][group] > 0 for group in final_groups)
    ]
    if not retained_features:
        raise SystemExit("no nonzero features remain for source attribution")
    analysis_group_totals = {
        group: sum(aggregated[feature][group] for feature in retained_features)
        for group in final_groups
    }

    with open(args.output_counts, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["#OTU ID", *final_groups])
        for feature in retained_features:
            writer.writerow([feature, *[aggregated[feature][group] for group in final_groups]])

    reserved_mapping_fields = {
        "#SampleID",
        "SourceSink",
        "Env",
        "aggregation_n",
        "aggregation_members",
        "sample_id",
        "source_sink",
        "model_source_env",
    }
    mapping_extras = [
        field for field in metadata_fields if field not in reserved_mapping_fields
    ]
    map_fields = [
        "#SampleID",
        "SourceSink",
        "Env",
        "aggregation_n",
        "aggregation_members",
        *mapping_extras,
    ]
    with open(args.output_map, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=map_fields, lineterminator="\n")
        writer.writeheader()
        for group in final_groups:
            row = group_metadata[group]
            output = {
                "#SampleID": group,
                "SourceSink": row["source_sink"],
                "Env": row["model_source_env"],
                "aggregation_n": len(groups[group]),
                "aggregation_members": ",".join(groups[group]),
            }
            output.update({field: row[field] for field in mapping_extras})
            writer.writerow(output)

    with open(args.sample_qc, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=QC_FIELDS, lineterminator="\n")
        writer.writeheader()
        for sample in samples:
            row, group = metadata[sample], group_of[sample]
            qc, members = group_feature_qc[group], groups[group]
            role = first_present(row, "source_sink")
            writer.writerow({
                "input_sample_id": sample,
                "analysis_sample_id": group,
                "unit_id": first_present(row, "unit_id") or sample,
                "source_sink": role or "NA",
                "model_source_env": first_present(row, "model_source_env") or "NA",
                "kraken_classified_fraction": format_number(qc["kraken_classified_fraction"]),
                "bracken_total_reads": format_number(qc["bracken_total_reads"], count=True),
                "taxon_filter_retained_reads": format_number(qc["taxon_filter_retained_reads"], count=True),
                "taxon_filter_retained_fraction": format_number(qc["taxon_filter_retained_fraction"]),
                "input_total_count": input_totals[sample],
                "analysis_total_count": format_number(analysis_group_totals.get(group), count=True),
                "required_min_count": {"source": args.min_source_count, "sink": args.min_sink_count}.get(role, 0),
                "aggregation_n": len(members),
                "aggregation_members": ",".join(members),
                "included": str(group in final_groups).lower(),
                "exclusion_reason": reasons.get(sample, "NA"),
            })

    with open(args.summary, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["input_samples", len(samples)])
        writer.writerow(["candidate_analysis_units", len(groups)])
        writer.writerow(["kept_analysis_units", len(final_groups)])
        writer.writerow(["kept_input_samples", sum(len(groups[group]) for group in final_groups)])
        writer.writerow(["excluded_input_samples", len(reasons)])
        writer.writerow(["source_analysis_units", len(source_groups)])
        writer.writerow(["sink_analysis_units", len(sink_groups)])
        writer.writerow(["kept_analysis_unit_ids", ",".join(final_groups)])
        excluded = [sample for sample in samples if sample in reasons]
        writer.writerow(["excluded_input_sample_ids", ",".join(excluded) or "NONE"])


if __name__ == "__main__":
    main()
