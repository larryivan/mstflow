#!/usr/bin/env python3
import argparse
import csv
import fnmatch
import hashlib
import math
import random
import re
import sys
from collections import OrderedDict


NA_VALUES = {"", "NA", "N/A", "NONE", "NULL"}
CONSISTENT_AGGREGATION_FIELDS = (
    "source_sink",
    "model_source_env",
    "category",
    "source_class",
    "panel_role",
    "panel_layer",
    "validation_group",
    "project_group",
    "host_group",
)
QC_FIELDS = [
    "analysis_label",
    "input_sample_id",
    "analysis_sample_id",
    "unit_id",
    "source_sink",
    "model_source_env",
    "source_class",
    "panel_role",
    "panel_layer",
    "category",
    "validation_group",
    "project_group",
    "host_group",
    "taxon_scope",
    "kraken_classified_fraction",
    "bracken_total_reads",
    "taxon_filter_retained_reads",
    "taxon_filter_retained_fraction",
    "input_total_count",
    "analysis_total_count",
    "required_rarefaction_depth",
    "aggregation_n",
    "aggregation_members",
    "included",
    "exclusion_reason",
]


def norm(value):
    return "" if value is None else str(value).strip()


def is_na(value):
    return norm(value).upper() in NA_VALUES


def first_present(row, *columns):
    for column in columns:
        value = norm(row.get(column))
        if not is_na(value):
            return value
    return ""


def metadata_env(row):
    return first_present(row, "model_source_env", "category")


def metadata_source_class(row):
    return first_present(row, "source_class", "class", "category")


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


def finite_nonnegative(value, context):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise SystemExit(f"{context} is not numeric: {value!r}")
    if not math.isfinite(number) or number < 0:
        raise SystemExit(f"{context} must be finite and non-negative: {value!r}")
    return number


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
                sample: finite_nonnegative(value, f"{path}:{line_no}:{sample}")
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
            "kraken_classified_fraction",
            "bracken_total_reads",
            "taxon_filter_retained_reads",
            "taxon_filter_retained_fraction",
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


def optional_qc_number(row, field):
    value = norm(row.get(field))
    if is_na(value):
        return None
    return finite_nonnegative(value, f"feature QC field {field}")


def complete_sum(rows, field):
    values = [optional_qc_number(row, field) for row in rows]
    return None if any(value is None for value in values) else sum(values)


def aggregate_feature_qc(members, feature_qc):
    rows = [feature_qc[sample] for sample in members]
    kraken_total = complete_sum(rows, "kraken_total_reads")
    kraken_classified = complete_sum(rows, "kraken_classified_reads")
    bracken_total = complete_sum(rows, "bracken_total_reads")
    retained_reads = complete_sum(rows, "taxon_filter_retained_reads")
    classified_fraction = (
        kraken_classified / kraken_total
        if kraken_classified is not None and kraken_total is not None and kraken_total > 0
        else None
    )
    retained_fraction = (
        retained_reads / bracken_total
        if retained_reads is not None and bracken_total is not None and bracken_total > 0
        else None
    )
    return {
        "kraken_total_reads": kraken_total,
        "kraken_classified_reads": kraken_classified,
        "kraken_classified_fraction": classified_fraction,
        "bracken_total_reads": bracken_total,
        "taxon_filter_retained_reads": retained_reads,
        "taxon_filter_retained_fraction": retained_fraction,
    }


def quality_exclusion_reason(
    row,
    taxon_scope,
    min_classified,
    min_retained_reads,
    min_retained_fraction,
):
    classified = optional_qc_number(row, "kraken_classified_fraction")
    bracken_reads = optional_qc_number(row, "bracken_total_reads")
    retained_reads = optional_qc_number(row, "taxon_filter_retained_reads")
    retained_fraction = optional_qc_number(row, "taxon_filter_retained_fraction")
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


def parse_patterns(values):
    patterns = []
    for value in values or []:
        for item in re.split(r"[,;]", value):
            item = item.strip()
            if item:
                patterns.append(item)
    return patterns


def matches_pattern(value, patterns):
    folded = norm(value).casefold()
    return any(fnmatch.fnmatchcase(folded, pattern.casefold()) for pattern in patterns)


def panel_role_excluded(value, excluded_patterns):
    role = norm(value).casefold()
    tokens = {token for token in re.split(r"[^a-z0-9]+", role) if token}
    for pattern in excluded_patterns:
        folded = pattern.casefold()
        if fnmatch.fnmatchcase(role, folded) or folded in tokens:
            return True
    return False


def sample_total(sample, features, matrix):
    return sum(matrix[feature][sample] for feature in features)


def rarefy_counts(values, depth, seed, group_id):
    """Subsample one analysis unit without expanding its count vector."""
    integer_values = [math.ceil(value) for value in values]
    total = sum(integer_values)
    if depth == 0:
        return integer_values
    if total < depth:
        raise SystemExit(
            f"analysis unit '{group_id}' has {total} counts, below LOO rarefaction depth {depth}"
        )
    digest = hashlib.sha256(f"{seed}\0{group_id}".encode("utf-8")).digest()
    group_seed = int.from_bytes(digest[:8], "big")
    selected = sorted(random.Random(group_seed).sample(range(total), depth))
    rarefied = [0] * len(integer_values)
    feature_index = 0
    boundary = integer_values[0] if integer_values else 0
    for position in selected:
        while position >= boundary:
            feature_index += 1
            boundary += integer_values[feature_index]
        rarefied[feature_index] += 1
    return rarefied


def initial_exclusion_reason(
    row,
    include_panel_roles,
    exclude_panel_roles,
    exclude_source_envs,
    exclude_source_classes,
    loo,
):
    panel_role = first_present(row, "panel_role")
    if panel_role_excluded(panel_role, exclude_panel_roles):
        return f"panel_role_excluded:{panel_role or 'NA'}"
    if include_panel_roles and not matches_pattern(panel_role, include_panel_roles):
        return f"panel_role_not_included:{panel_role or 'NA'}"

    source_sink = norm(row.get("source_sink")).lower()
    if is_na(source_sink) or source_sink not in {"source", "sink"}:
        return "source_sink_not_selected"
    if loo and source_sink != "source":
        return "loo_sources_only"
    source_env = metadata_env(row)
    source_class = metadata_source_class(row)
    if source_sink == "source" and is_na(source_env):
        return "missing_source_environment"
    if source_sink == "source" and matches_pattern(source_env, exclude_source_envs):
        return f"source_environment_excluded:{source_env}"
    if source_sink == "source" and matches_pattern(source_class, exclude_source_classes):
        return f"source_class_excluded:{source_class}"
    return ""


def canonical_group_value(row, field):
    if field == "model_source_env":
        return metadata_env(row)
    if field == "source_class":
        return metadata_source_class(row)
    return first_present(row, field)


def check_group_consistency(group_id, members, metadata):
    for field in CONSISTENT_AGGREGATION_FIELDS:
        values = {
            canonical_group_value(metadata[sample], field)
            if not is_na(canonical_group_value(metadata[sample], field))
            else "<MISSING>"
            for sample in members
        }
        if len(values) > 1:
            raise SystemExit(
                f"aggregation group '{group_id}' has conflicting {field}: {', '.join(sorted(values))}"
            )


def collapse_values(rows, field):
    values = []
    for row in rows:
        value = norm(row.get(field))
        if not is_na(value) and value not in values:
            values.append(value)
    return "NA" if not values else values[0] if len(values) == 1 else ";".join(values)


def aggregate_metadata(group_id, members, metadata, metadata_fields, aggregate_by):
    rows = [metadata[sample] for sample in members]
    result = {field: collapse_values(rows, field) for field in metadata_fields}
    result["sample_id"] = group_id
    if aggregate_by in result:
        result[aggregate_by] = group_id
    result["source_sink"] = norm(rows[0].get("source_sink")).lower()
    result["model_source_env"] = metadata_env(rows[0]) or "NA"
    result["source_class"] = metadata_source_class(rows[0]) or "NA"
    result["unit_id"] = (
        group_id if aggregate_by == "unit_id" else first_present(rows[0], "unit_id") or group_id
    )
    return result


def format_number(value):
    return "NA" if value is None else f"{value:.12g}"


def make_qc_row(
    analysis_label,
    sample,
    analysis_sample,
    row,
    input_total,
    analysis_total,
    required_depth,
    members,
    included,
    reason,
    feature_qc,
    taxon_scope,
):
    return {
        "analysis_label": analysis_label,
        "input_sample_id": sample,
        "analysis_sample_id": analysis_sample or "NA",
        "unit_id": first_present(row, "unit_id") or sample,
        "source_sink": first_present(row, "source_sink") or "NA",
        "model_source_env": metadata_env(row) or "NA",
        "source_class": metadata_source_class(row) or "NA",
        "panel_role": first_present(row, "panel_role") or "NA",
        "panel_layer": first_present(row, "panel_layer", "layer") or "NA",
        "category": first_present(row, "category") or "NA",
        "validation_group": first_present(row, "validation_group") or "NA",
        "project_group": first_present(row, "project_group") or "NA",
        "host_group": first_present(row, "host_group") or "NA",
        "taxon_scope": taxon_scope,
        "kraken_classified_fraction": feature_qc.get("kraken_classified_fraction", "NA"),
        "bracken_total_reads": feature_qc.get("bracken_total_reads", "NA"),
        "taxon_filter_retained_reads": feature_qc.get("taxon_filter_retained_reads", "NA"),
        "taxon_filter_retained_fraction": feature_qc.get("taxon_filter_retained_fraction", "NA"),
        "input_total_count": format_number(input_total),
        "analysis_total_count": format_number(analysis_total),
        "required_rarefaction_depth": required_depth,
        "aggregation_n": len(members) if members else 0,
        "aggregation_members": ",".join(members) if members else "NA",
        "included": str(bool(included)).lower(),
        "exclusion_reason": reason or "NA",
    }


def write_qc_rows(path, rows):
    if not path:
        return
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=QC_FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description="Select, aggregate, depth-check, and document SourceTracker2 inputs."
    )
    parser.add_argument("--counts", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--feature-sample-qc", required=True)
    parser.add_argument("--output-counts", required=True)
    parser.add_argument("--output-map", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--exclusions")
    parser.add_argument("--sample-qc")
    parser.add_argument("--aggregate-by", default="sample_id")
    parser.add_argument("--source-rarefaction-depth", type=int, default=0)
    parser.add_argument("--sink-rarefaction-depth", type=int, default=0)
    parser.add_argument("--include-panel-role", action="append", default=[])
    parser.add_argument("--exclude-panel-role", action="append", default=[])
    parser.add_argument("--exclude-source-env", "--exclude-source-envs", action="append", default=[])
    parser.add_argument("--exclude-source-class", "--exclude-source-classes", action="append", default=[])
    parser.add_argument("--analysis-label", default="default")
    parser.add_argument("--random-seed", type=int, default=1)
    parser.add_argument("--taxon-scope", choices=["all", "filtered"], default="filtered")
    parser.add_argument("--min-kraken-classified-fraction", type=float, default=0.001)
    parser.add_argument("--min-taxon-filter-retained-reads", type=float, default=1000.0)
    parser.add_argument("--min-taxon-filter-retained-fraction", type=float, default=0.001)
    parser.add_argument(
        "--loo",
        action="store_true",
        help="Prepare source-only inputs for SourceTracker2 leave-one-out validation.",
    )
    args = parser.parse_args()

    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", args.analysis_label):
        raise SystemExit(
            "--analysis-label must be filesystem-safe: letters/numbers followed by letters, "
            "numbers, dot, underscore, or hyphen"
        )
    if args.source_rarefaction_depth < 0 or args.sink_rarefaction_depth < 0:
        raise SystemExit("rarefaction depths must be non-negative")
    if not 0 <= args.random_seed <= 2_147_483_647:
        raise SystemExit("--random-seed must be between 0 and 2147483647")
    if args.loo and args.source_rarefaction_depth != 0:
        raise SystemExit(
            "SourceTracker2 2.0.1 LOO requires --source-rarefaction-depth 0: "
            "the upstream CLI subtracts unrarefied held-out counts from already-rarefied "
            "source counts"
        )
    if not 0 <= args.min_kraken_classified_fraction <= 1:
        raise SystemExit("--min-kraken-classified-fraction must be between 0 and 1")
    if args.min_taxon_filter_retained_reads < 0:
        raise SystemExit("--min-taxon-filter-retained-reads must be non-negative")
    if not 0 <= args.min_taxon_filter_retained_fraction <= 1:
        raise SystemExit("--min-taxon-filter-retained-fraction must be between 0 and 1")

    metadata_fields, metadata = read_metadata(args.metadata)
    samples, features, matrix = read_counts(args.counts)
    require_one_to_one(metadata, samples)
    feature_qc = read_feature_sample_qc(args.feature_sample_qc)
    require_one_to_one(feature_qc, samples)
    if args.aggregate_by != "sample_id" and args.aggregate_by not in metadata_fields:
        raise SystemExit(f"--aggregate-by column is absent from metadata: {args.aggregate_by}")

    include_panel_roles = parse_patterns(args.include_panel_role)
    exclude_panel_roles = sorted(parse_patterns(args.exclude_panel_role))
    exclude_source_envs = parse_patterns(args.exclude_source_env)
    exclude_source_classes = parse_patterns(args.exclude_source_class)
    input_totals = {sample: sample_total(sample, features, matrix) for sample in samples}

    reasons = {}
    candidate_group = {}
    groups = OrderedDict()
    for sample in samples:
        row = metadata[sample]
        group_id = sample if args.aggregate_by == "sample_id" else first_present(row, args.aggregate_by)
        if not group_id:
            reasons[sample] = f"missing_aggregation_key:{args.aggregate_by}"
            continue
        candidate_group[sample] = group_id
        groups.setdefault(group_id, []).append(sample)

    group_metadata = {}
    group_feature_qc = {}
    group_totals = {}
    final_groups = []
    for group_id, members in groups.items():
        # Validate the complete biological unit before applying any exclusion.
        # Otherwise one conflicting or low-quality run could be removed first
        # and leave a scientifically mislabelled partial unit behind.
        check_group_consistency(group_id, members, metadata)
        group_metadata[group_id] = aggregate_metadata(
            group_id, members, metadata, metadata_fields, args.aggregate_by
        )
        group_feature_qc[group_id] = aggregate_feature_qc(members, feature_qc)
        total = sum(input_totals[sample] for sample in members)
        group_totals[group_id] = total
        role = group_metadata[group_id]["source_sink"]
        required_depth = (
            args.sink_rarefaction_depth
            if args.loo and role == "source"
            else args.source_rarefaction_depth if role == "source" else args.sink_rarefaction_depth
        )
        group_reason = initial_exclusion_reason(
            group_metadata[group_id],
            include_panel_roles,
            exclude_panel_roles,
            exclude_source_envs,
            exclude_source_classes,
            args.loo,
        )
        if not group_reason:
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
            depth_label = "loo_analysis_unit" if args.loo and role == "source" else role
            group_reason = f"below_{depth_label}_rarefaction_depth:{required_depth}"
        if group_reason:
            for sample in members:
                reasons[sample] = group_reason
        else:
            final_groups.append(group_id)

    if not final_groups:
        raise SystemExit("no SourceTracker2 samples remain after role, exclusion, and depth QC")
    source_groups = [
        group for group in final_groups if group_metadata[group]["source_sink"] == "source"
    ]
    sink_groups = [
        group for group in final_groups if group_metadata[group]["source_sink"] == "sink"
    ]
    if not source_groups:
        raise SystemExit("SourceTracker2 requires at least one source after exclusions and depth QC")
    if not args.loo and not sink_groups:
        raise SystemExit("SourceTracker2 requires at least one sink after exclusions and depth QC")
    reserved_source_envs = [
        group_metadata[group]["model_source_env"]
        for group in source_groups
        if group_metadata[group]["model_source_env"].casefold() in {"unknown", "sampleid"}
    ]
    if reserved_source_envs:
        raise SystemExit(
            "SourceTracker2 source environment names must not use reserved result columns "
            "Unknown or SampleID; found: " + ", ".join(sorted(set(reserved_source_envs)))
        )
    if args.loo:
        source_groups_by_env = OrderedDict()
        for group in source_groups:
            source_groups_by_env.setdefault(group_metadata[group]["model_source_env"], []).append(
                group
            )
        singleton_envs = [
            env for env, env_groups in source_groups_by_env.items() if len(env_groups) < 2
        ]
        if singleton_envs:
            raise SystemExit(
                "SourceTracker2 LOO requires at least two independent analysis units in every "
                "source environment; insufficient: " + ", ".join(singleton_envs)
            )

    aggregated = {}
    for feature in features:
        aggregated[feature] = {
            group: sum(matrix[feature][sample] for sample in groups[group])
            for group in final_groups
        }

    # SourceTracker2 2.0.1 does not rarefy held-out source samples in --loo
    # mode. Rarefy every analysis unit here before it pools source
    # environments, so held-out sinks have equal depth and subtraction is
    # mathematically valid. Stable per-group seeds prevent panel ordering from
    # changing an existing unit's subsample.
    if args.loo and args.sink_rarefaction_depth > 0:
        for group in final_groups:
            rarefied = rarefy_counts(
                [aggregated[feature][group] for feature in features],
                args.sink_rarefaction_depth,
                args.random_seed,
                group,
            )
            for feature, value in zip(features, rarefied):
                aggregated[feature][group] = value

    retained_features = [
        feature
        for feature in features
        if any(aggregated[feature][group] > 0 for group in final_groups)
    ]
    if not retained_features:
        raise SystemExit("no nonzero features remain for SourceTracker2")
    analysis_group_totals = {
        group: sum(aggregated[feature][group] for feature in retained_features)
        for group in final_groups
    }

    with open(args.output_counts, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["#OTU ID", *final_groups])
        for feature in retained_features:
            writer.writerow(
                [feature, *[format_number(aggregated[feature][group]) for group in final_groups]]
            )

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
            output.update({field: row.get(field, "NA") for field in mapping_extras})
            writer.writerow(output)

    qc_rows = []
    for sample in samples:
        group = candidate_group.get(sample)
        members = groups.get(group, []) if group else []
        role = norm(metadata[sample].get("source_sink")).lower()
        required_depth = (
            args.sink_rarefaction_depth
            if args.loo and role == "source"
            else args.source_rarefaction_depth if role == "source"
            else args.sink_rarefaction_depth if role == "sink" else 0
        )
        included = bool(group in final_groups)
        qc_rows.append(
            make_qc_row(
                args.analysis_label,
                sample,
                group,
                metadata[sample],
                input_totals[sample],
                analysis_group_totals.get(group),
                required_depth,
                members,
                included,
                reasons.get(sample, ""),
                group_feature_qc.get(group, feature_qc[sample]),
                args.taxon_scope,
            )
        )
    write_qc_rows(args.sample_qc, qc_rows)
    write_qc_rows(args.exclusions, [row for row in qc_rows if row["included"] == "false"])

    with open(args.summary, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["analysis_label", args.analysis_label])
        writer.writerow(["loo", str(args.loo).lower()])
        writer.writerow(["sample_completeness", "one_to_one"])
        writer.writerow(["aggregate_by", args.aggregate_by])
        writer.writerow(["taxon_scope", args.taxon_scope])
        writer.writerow(["input_samples", len(samples)])
        writer.writerow(["candidate_analysis_units", len(groups)])
        writer.writerow(["kept_analysis_units", len(final_groups)])
        writer.writerow(["kept_input_samples", sum(len(groups[group]) for group in final_groups)])
        writer.writerow(["excluded_input_samples", sum(sample in reasons for sample in samples)])
        writer.writerow(["source_analysis_units", len(source_groups)])
        writer.writerow(["sink_analysis_units", len(sink_groups)])
        writer.writerow(["source_rarefaction_depth", args.source_rarefaction_depth])
        writer.writerow(["sink_rarefaction_depth", args.sink_rarefaction_depth])
        writer.writerow(
            [
                "loo_analysis_unit_rarefaction_depth",
                args.sink_rarefaction_depth if args.loo else "NA",
            ]
        )
        writer.writerow(["random_seed", args.random_seed])
        writer.writerow(["min_kraken_classified_fraction", args.min_kraken_classified_fraction])
        writer.writerow(["min_taxon_filter_retained_reads", args.min_taxon_filter_retained_reads])
        writer.writerow(["min_taxon_filter_retained_fraction", args.min_taxon_filter_retained_fraction])
        writer.writerow(["include_panel_roles", ",".join(include_panel_roles) or "ALL_SAFE_ROLES"])
        writer.writerow(["exclude_panel_roles", ",".join(exclude_panel_roles)])
        writer.writerow(["exclude_source_envs", ",".join(exclude_source_envs) or "NONE"])
        writer.writerow(["exclude_source_classes", ",".join(exclude_source_classes) or "NONE"])
        writer.writerow(["kept_analysis_unit_ids", ",".join(final_groups)])
        excluded = [sample for sample in samples if sample in reasons]
        writer.writerow(["excluded_input_sample_ids", ",".join(excluded) or "NONE"])


if __name__ == "__main__":
    main()
