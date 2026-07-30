#!/usr/bin/env python3
import argparse
import csv
import os
import re
import sys


REQUIRED_COLUMNS = ["sample_id", "fastq_1", "fastq_2", "layout", "age_group"]
ALLOWED_LAYOUTS = {"SE", "PE"}
ALLOWED_AGE_GROUPS = {"modern", "ancient"}
ALLOWED_PREPROCESS_MODES = {
    "fastp_pe",
    "fastp_se",
    "adapterremoval_pe_collapse",
    "adapterremoval_se",
}
ALLOWED_SOURCE_SINK = {"source", "sink"}
ALLOWED_HOST_ACTIONS = {"remove", "skip"}
SAFE_SAMPLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
STANDARD_METADATA_COLUMNS = [
    "preprocess_mode",
    "bracken_read_len",
    "host_action",
    "host_ref",
    "source_sink",
    "model_source_env",
    "unit_id",
    "panel_role",
    "panel_layer",
    "category",
    "validation_group",
    "project_group",
    "host_group",
    "source_major",
    "source_class",
    "sample_label",
    "project",
    "experiment",
    "run",
]
OUTPUT_COLUMNS = REQUIRED_COLUMNS + STANDARD_METADATA_COLUMNS


def norm(value):
    return "" if value is None else str(value).strip()


def is_na(value):
    return norm(value).upper() in {"", "NA", "N/A", "NONE", "NULL"}


def fail(errors):
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    sys.exit(1)


def host_reference_errors(path):
    if is_na(path):
        return ["default host FASTA is not configured"]
    errors = []
    if not os.path.isfile(path):
        return [f"host FASTA does not exist: {path}"]
    for suffix in (".amb", ".ann", ".bwt", ".pac", ".sa"):
        sidecar = f"{path}{suffix}"
        if not os.path.isfile(sidecar):
            errors.append(f"BWA index sidecar does not exist: {sidecar}")
    return errors


def choose_preprocess_mode(layout, age_group, requested):
    if requested and requested != "auto":
        return requested
    if age_group == "modern" and layout == "PE":
        return "fastp_pe"
    if age_group == "modern" and layout == "SE":
        return "fastp_se"
    if age_group == "ancient" and layout == "PE":
        return "adapterremoval_pe_collapse"
    if age_group == "ancient" and layout == "SE":
        return "adapterremoval_se"
    raise ValueError(f"cannot infer preprocess mode for layout={layout}, age_group={age_group}")


def first_present(row, *columns):
    for column in columns:
        value = norm(row.get(column))
        if not is_na(value):
            return value
    return ""


def output_columns(input_columns):
    """Return canonical columns followed by every user-provided extra column."""
    return OUTPUT_COLUMNS + [column for column in input_columns if column not in OUTPUT_COLUMNS]


def normalized_output_row(
    row,
    sample_id,
    fastq_1,
    fastq_2,
    layout,
    age_group,
    preprocess_mode,
    host_action,
    host_group,
    host_ref,
    input_columns,
):
    source_sink = norm(row.get("source_sink")).lower()
    if is_na(source_sink):
        source_sink = ""
    category = first_present(row, "category", "source_class", "class")
    model_source_env = first_present(row, "model_source_env")
    # ``category`` is the Reference Panel's canonical environment label.  It is
    # a safe alias only after a row has been explicitly selected as a source or
    # sink; it never opts an unrelated row into SourceTracker.
    if not model_source_env and source_sink in ALLOWED_SOURCE_SINK:
        model_source_env = category

    # Preserve the complete input schema.  This prevents panel design fields
    # from disappearing between metadata validation and SourceTracker QC.
    output = {column: norm(row.get(column)) for column in input_columns}
    requested_bracken_length = norm(row.get("bracken_read_len"))
    if is_na(requested_bracken_length) or requested_bracken_length.lower() in {
        "auto",
        "per_sample",
    }:
        requested_bracken_length = "auto"
    elif requested_bracken_length.isdigit() and int(requested_bracken_length) > 0:
        requested_bracken_length = str(int(requested_bracken_length))
    output.update({
        "sample_id": sample_id,
        "fastq_1": fastq_1,
        "fastq_2": "NA" if is_na(fastq_2) else fastq_2,
        "layout": layout,
        "age_group": age_group,
        "preprocess_mode": preprocess_mode,
        "bracken_read_len": requested_bracken_length,
        "host_action": host_action,
        "host_ref": (
            "NA" if host_action == "skip" else ("default" if is_na(host_ref) else host_ref)
        ),
        "source_sink": "NA" if is_na(source_sink) else source_sink,
        "model_source_env": "NA" if is_na(model_source_env) else model_source_env,
        "unit_id": first_present(row, "unit_id") or sample_id,
        "panel_role": first_present(row, "panel_role") or "NA",
        "panel_layer": first_present(row, "panel_layer", "layer") or "NA",
        "category": category or "NA",
        "validation_group": first_present(row, "validation_group") or "NA",
        "project_group": first_present(row, "project_group") or "NA",
        "host_group": host_group,
        "source_major": first_present(row, "source_major", "major") or "NA",
        "source_class": first_present(row, "source_class", "class", "category") or "NA",
        "sample_label": norm(row.get("sample_label")) or sample_id,
        "project": norm(row.get("project")) or "NA",
        "experiment": norm(row.get("experiment")) or "NA",
        "run": norm(row.get("run")) or sample_id,
    })
    return output


def main():
    parser = argparse.ArgumentParser(description="Validate and normalize mstflow metadata.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--default-host-group", default="human")
    parser.add_argument(
        "--default-host-fasta",
        help=(
            "Default host FASTA used when host_ref=default. Supplying this option "
            "enables early FASTA and classic-BWA sidecar validation."
        ),
    )
    parser.add_argument(
        "--require-explicit-host-for-nondefault",
        action="store_true",
        help=(
            "Reject rows whose host_group differs from --default-host-group when "
            "host_ref is default/NA. This prevents non-human samples from silently "
            "using the default human reference."
        ),
    )
    parser.add_argument(
        "--require-explicit-host-action",
        action="store_true",
        help=(
            "Reject rows that omit host_action. Use this for mixed host-aware runs "
            "so every row explicitly declares remove or skip."
        ),
    )
    parser.add_argument(
        "--skip-filesystem-checks",
        action="store_true",
        help=(
            "Validate metadata semantics while leaving FASTQ/host paths opaque. "
            "Use only when the workflow controller independently stages and "
            "checks every referenced file."
        ),
    )
    args = parser.parse_args()

    default_host_group = norm(args.default_host_group).lower()
    if not default_host_group:
        fail(["--default-host-group must not be empty"])

    errors = []
    seen = set()
    rows_out = []
    checked_host_references = {}

    with open(args.input, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail([f"metadata is empty: {args.input}"])

        duplicate_columns = sorted(
            {column for column in reader.fieldnames if reader.fieldnames.count(column) > 1}
        )
        if duplicate_columns:
            fail([f"metadata contains duplicate column name(s): {', '.join(duplicate_columns)}"])

        missing = [column for column in REQUIRED_COLUMNS if column not in reader.fieldnames]
        if missing:
            fail([f"metadata is missing required column(s): {', '.join(missing)}"])

        for line_no, row in enumerate(reader, start=2):
            sample_id = norm(row.get("sample_id"))
            if is_na(sample_id):
                errors.append(
                    f"line {line_no}: sample_id is empty or uses a reserved missing-value token"
                )
                continue
            if not SAFE_SAMPLE_ID.fullmatch(sample_id):
                errors.append(
                    f"line {line_no}: sample_id '{sample_id}' is not filesystem-safe; "
                    "use only letters, digits, '.', '_' and '-' and start with a letter or digit"
                )
            if sample_id in seen:
                errors.append(f"line {line_no}: duplicate sample_id '{sample_id}'")
            seen.add(sample_id)

            layout = norm(row.get("layout")).upper()
            if layout not in ALLOWED_LAYOUTS:
                errors.append(f"line {line_no}: layout must be SE or PE")

            age_group = norm(row.get("age_group")).lower()
            if age_group not in ALLOWED_AGE_GROUPS:
                errors.append(f"line {line_no}: age_group must be modern or ancient")

            fastq_1 = norm(row.get("fastq_1"))
            fastq_2 = norm(row.get("fastq_2"))
            if is_na(fastq_1):
                errors.append(f"line {line_no}: fastq_1 is required")
            elif not args.skip_filesystem_checks and not os.path.isfile(fastq_1):
                errors.append(f"line {line_no}: fastq_1 is not a regular file: {fastq_1}")

            if layout == "PE":
                if is_na(fastq_2):
                    errors.append(f"line {line_no}: PE samples require fastq_2")
                elif not args.skip_filesystem_checks and not os.path.isfile(fastq_2):
                    errors.append(f"line {line_no}: fastq_2 is not a regular file: {fastq_2}")
            elif not is_na(fastq_2):
                errors.append(f"line {line_no}: SE samples must set fastq_2 to NA")

            requested_mode = norm(row.get("preprocess_mode")).lower()
            if is_na(requested_mode):
                requested_mode = "auto"

            try:
                preprocess_mode = choose_preprocess_mode(layout, age_group, requested_mode)
            except ValueError as exc:
                errors.append(f"line {line_no}: {exc}")
                preprocess_mode = requested_mode

            if preprocess_mode not in ALLOWED_PREPROCESS_MODES:
                errors.append(
                    f"line {line_no}: preprocess_mode must be auto or one of "
                    f"{', '.join(sorted(ALLOWED_PREPROCESS_MODES))}"
                )

            if preprocess_mode.startswith("fastp") and age_group != "modern":
                errors.append(f"line {line_no}: fastp preprocessing is only allowed for modern samples")
            if preprocess_mode.startswith("adapterremoval") and age_group != "ancient":
                errors.append(f"line {line_no}: AdapterRemoval preprocessing is only allowed for ancient samples")
            if preprocess_mode.endswith("_pe") or preprocess_mode == "adapterremoval_pe_collapse":
                if layout != "PE":
                    errors.append(f"line {line_no}: {preprocess_mode} requires PE layout")
            if preprocess_mode.endswith("_se") and layout != "SE":
                errors.append(f"line {line_no}: {preprocess_mode} requires SE layout")

            requested_bracken_length = norm(row.get("bracken_read_len"))
            if not is_na(requested_bracken_length) and requested_bracken_length.lower() not in {
                "auto",
                "per_sample",
            }:
                if not requested_bracken_length.isdigit() or int(requested_bracken_length) <= 0:
                    errors.append(
                        f"line {line_no}: bracken_read_len must be auto or a positive integer"
                    )

            requested_host_action = norm(row.get("host_action")).lower()
            if args.require_explicit_host_action and is_na(requested_host_action):
                errors.append(
                    f"line {line_no}: host_action must be explicitly set to remove or skip"
                )
            host_action = "remove" if is_na(requested_host_action) else requested_host_action
            if host_action not in ALLOWED_HOST_ACTIONS:
                errors.append(
                    f"line {line_no}: host_action must be remove or skip"
                )
                # Continue validation deterministically without treating an
                # invalid value as an accidental skip.
                host_action = "remove"

            host_ref = norm(row.get("host_ref")) or "default"
            requested_host_group = norm(row.get("host_group")).lower()
            host_group = (
                "none" if host_action == "skip" else default_host_group
            ) if is_na(requested_host_group) else requested_host_group
            if host_action == "skip":
                if not is_na(host_ref):
                    errors.append(
                        f"line {line_no}: host_action=skip requires host_ref=NA"
                    )
                host_ref = "NA"
            else:
                if (
                    not args.skip_filesystem_checks
                    and not is_na(host_ref)
                    and host_ref != "default"
                    and not os.path.exists(host_ref)
                ):
                    errors.append(
                        f"line {line_no}: host_ref must be default, NA, or an existing FASTA path: {host_ref}"
                    )
                if args.default_host_fasta is not None and not args.skip_filesystem_checks:
                    effective_host_ref = (
                        norm(args.default_host_fasta)
                        if is_na(host_ref) or host_ref == "default"
                        else host_ref
                    )
                    if effective_host_ref not in checked_host_references:
                        checked_host_references[effective_host_ref] = host_reference_errors(
                            effective_host_ref
                        )
                        errors.extend(checked_host_references[effective_host_ref])
                if (
                    args.require_explicit_host_for_nondefault
                    and host_group != default_host_group
                    and (is_na(host_ref) or host_ref == "default")
                ):
                    errors.append(
                        f"line {line_no}: host_group '{host_group}' requires an explicit host_ref; "
                        f"the default reference is for host_group '{default_host_group}'"
                    )

            source_sink = norm(row.get("source_sink")).lower()
            if is_na(source_sink):
                source_sink = ""
            category = first_present(row, "category", "source_class", "class")
            model_source_env = first_present(row, "model_source_env")
            if not model_source_env and source_sink in ALLOWED_SOURCE_SINK:
                model_source_env = category
            if source_sink and source_sink not in ALLOWED_SOURCE_SINK:
                errors.append(f"line {line_no}: source_sink must be source or sink")
            if source_sink == "source" and is_na(model_source_env):
                errors.append(
                    f"line {line_no}: source rows require model_source_env or category"
                )

            unit_id = first_present(row, "unit_id") or sample_id
            if not SAFE_SAMPLE_ID.fullmatch(unit_id):
                errors.append(
                    f"line {line_no}: unit_id '{unit_id}' is not filesystem-safe; "
                    "use only letters, digits, '.', '_' and '-' and start with a letter or digit"
                )

            rows_out.append(
                normalized_output_row(
                    row,
                    sample_id,
                    fastq_1,
                    fastq_2,
                    layout,
                    age_group,
                    preprocess_mode,
                    host_action,
                    host_group,
                    host_ref,
                    reader.fieldnames,
                )
            )

    if not rows_out:
        errors.append(f"metadata contains no sample rows: {args.input}")

    if errors:
        fail(errors)

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=output_columns(reader.fieldnames),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows_out)


if __name__ == "__main__":
    main()
