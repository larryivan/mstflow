#!/usr/bin/env python3
"""Check a samplesheet and write the normalized copy the workflow reads.

Layout is inferred from fastq_2 when it is not given, the trimming mode from age
group and layout, and host removal defaults to the default host reference. Every
problem is reported, one per line, before the run fails. Columns the workflow does
not know are kept, so they reach the attribution map.
"""
import argparse
import csv
import re
import sys


REQUIRED_COLUMNS = ["sample_id", "fastq_1", "age_group"]
OUTPUT_COLUMNS = [
    "sample_id",
    "fastq_1",
    "fastq_2",
    "layout",
    "age_group",
    "preprocess_mode",
    "host_action",
    "host_ref",
    "host_group",
    "source_sink",
    "model_source_env",
    "unit_id",
    "project",
]
# Modern libraries are trimmed by fastp; ancient ones by AdapterRemoval, which also
# collapses overlapping read pairs into the molecules they came from.
PREPROCESS_MODES = {
    ("modern", "PE"): "fastp_pe",
    ("modern", "SE"): "fastp_se",
    ("ancient", "PE"): "adapterremoval_pe_collapse",
    ("ancient", "SE"): "adapterremoval_se",
}
SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")
SAFE_ID_RULE = "use only letters, digits, '.', '_' and '-' and start with a letter or digit"


def norm(value):
    return "" if value is None else str(value).strip()


def is_na(value):
    return norm(value).upper() in {"", "NA", "N/A", "NONE", "NULL"}


def fail(errors):
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--default-host-group", default="human")
    parser.add_argument(
        "--require-explicit-host-for-nondefault",
        action="store_true",
        help="Reject rows whose host_group differs from --default-host-group but use the default host "
             "reference, so non-human samples are not silently filtered against the human genome.",
    )
    parser.add_argument(
        "--require-explicit-host-action",
        action="store_true",
        help="Reject rows that leave host_action empty, so every row declares remove or skip.",
    )
    args = parser.parse_args()

    default_host_group = norm(args.default_host_group).lower()
    if not default_host_group:
        fail(["--default-host-group must not be empty"])

    with open(args.input, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames
        if columns is None:
            fail([f"metadata is empty: {args.input}"])
        duplicate_columns = sorted({column for column in columns if columns.count(column) > 1})
        if duplicate_columns:
            fail([f"metadata contains duplicate column name(s): {', '.join(duplicate_columns)}"])
        missing = [column for column in REQUIRED_COLUMNS if column not in columns]
        if missing:
            fail([f"metadata is missing required column(s): {', '.join(missing)}"])
        rows = list(reader)

    errors = []
    seen = set()
    rows_out = []
    for line_no, row in enumerate(rows, start=2):
        where = f"line {line_no}"
        sample_id = norm(row.get("sample_id"))
        if is_na(sample_id):
            errors.append(f"{where}: sample_id is empty or uses a reserved missing-value token")
            continue
        if not SAFE_ID.fullmatch(sample_id):
            errors.append(f"{where}: sample_id '{sample_id}' is not filesystem-safe; {SAFE_ID_RULE}")
        if sample_id in seen:
            errors.append(f"{where}: duplicate sample_id '{sample_id}'")
        seen.add(sample_id)

        fastq_1 = norm(row.get("fastq_1"))
        fastq_2 = norm(row.get("fastq_2"))
        layout = norm(row.get("layout")).upper()
        if is_na(layout):
            layout = "SE" if is_na(fastq_2) else "PE"
        elif layout not in ("SE", "PE"):
            errors.append(f"{where}: layout must be SE or PE, or left empty to infer it from fastq_2")
        age_group = norm(row.get("age_group")).lower()
        if age_group not in ("modern", "ancient"):
            errors.append(f"{where}: age_group must be modern or ancient")
        if is_na(fastq_1):
            errors.append(f"{where}: fastq_1 is required")
        if layout == "PE" and is_na(fastq_2):
            errors.append(f"{where}: PE samples require fastq_2")
        elif layout == "SE" and not is_na(fastq_2):
            errors.append(f"{where}: SE samples must set fastq_2 to NA")

        host_action = norm(row.get("host_action")).lower()
        if is_na(host_action):
            if args.require_explicit_host_action:
                errors.append(f"{where}: host_action must be explicitly set to remove or skip")
            host_action = "remove"
        elif host_action not in ("remove", "skip"):
            errors.append(f"{where}: host_action must be remove or skip")
        host_ref = norm(row.get("host_ref"))
        host_group = norm(row.get("host_group")).lower()
        if host_action == "skip":
            if not is_na(host_ref):
                errors.append(f"{where}: host_action=skip conflicts with host_ref={host_ref}; leave host_ref empty")
            host_ref = "NA"
            host_group = "NA" if is_na(host_group) else host_group
        else:
            host_ref = "default" if is_na(host_ref) else host_ref
            host_group = default_host_group if is_na(host_group) else host_group
            if args.require_explicit_host_for_nondefault and host_group != default_host_group and host_ref == "default":
                errors.append(
                    f"{where}: host_group '{host_group}' requires an explicit host_ref; "
                    f"the default reference is for host_group '{default_host_group}'"
                )

        source_sink = norm(row.get("source_sink")).lower()
        source_sink = "NA" if is_na(source_sink) else source_sink
        model_source_env = norm(row.get("model_source_env"))
        model_source_env = "NA" if is_na(model_source_env) else model_source_env
        if source_sink not in ("source", "sink", "NA"):
            errors.append(f"{where}: source_sink must be source or sink")
        if source_sink == "source" and model_source_env == "NA":
            errors.append(f"{where}: source rows require model_source_env")

        unit_id = norm(row.get("unit_id"))
        unit_id = sample_id if is_na(unit_id) else unit_id
        if not SAFE_ID.fullmatch(unit_id):
            errors.append(f"{where}: unit_id '{unit_id}' is not filesystem-safe; {SAFE_ID_RULE}")

        output = {column: norm(row.get(column)) for column in columns}
        output.update({
            "sample_id": sample_id,
            "fastq_1": fastq_1,
            "fastq_2": "NA" if is_na(fastq_2) else fastq_2,
            "layout": layout,
            "age_group": age_group,
            "preprocess_mode": PREPROCESS_MODES.get((age_group, layout), "NA"),
            "host_action": host_action,
            "host_ref": host_ref,
            "host_group": host_group,
            "source_sink": source_sink,
            "model_source_env": model_source_env,
            "unit_id": unit_id,
            "project": norm(row.get("project")) or "NA",
        })
        rows_out.append(output)

    if not rows_out:
        errors.append(f"metadata contains no sample rows: {args.input}")
    if errors:
        fail(errors)

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=OUTPUT_COLUMNS + [column for column in columns if column not in OUTPUT_COLUMNS],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows_out)


if __name__ == "__main__":
    main()
