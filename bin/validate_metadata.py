#!/usr/bin/env python3
import argparse
import csv
import os
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
OUTPUT_COLUMNS = REQUIRED_COLUMNS + [
    "preprocess_mode",
    "host_ref",
    "source_sink",
    "model_source_env",
    "source_major",
    "source_class",
    "sample_label",
    "project",
    "experiment",
    "run",
]


def norm(value):
    return "" if value is None else str(value).strip()


def is_na(value):
    return norm(value).upper() in {"", "NA", "N/A", "NONE", "NULL"}


def fail(errors):
    for error in errors:
        print(f"ERROR: {error}", file=sys.stderr)
    sys.exit(1)


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
        if value:
            return value
    return ""


def normalized_output_row(
    row,
    sample_id,
    fastq_1,
    fastq_2,
    layout,
    age_group,
    preprocess_mode,
    host_ref,
):
    source_sink = norm(row.get("source_sink")).lower()
    model_source_env = norm(row.get("model_source_env"))
    return {
        "sample_id": sample_id,
        "fastq_1": fastq_1,
        "fastq_2": "NA" if is_na(fastq_2) else fastq_2,
        "layout": layout,
        "age_group": age_group,
        "preprocess_mode": preprocess_mode,
        "host_ref": "default" if is_na(host_ref) else host_ref,
        "source_sink": "NA" if is_na(source_sink) else source_sink,
        "model_source_env": "NA" if is_na(model_source_env) else model_source_env,
        "source_major": first_present(row, "source_major", "major") or "NA",
        "source_class": first_present(row, "source_class", "class") or "NA",
        "sample_label": norm(row.get("sample_label")) or sample_id,
        "project": norm(row.get("project")) or "NA",
        "experiment": norm(row.get("experiment")) or "NA",
        "run": norm(row.get("run")) or sample_id,
    }


def main():
    parser = argparse.ArgumentParser(description="Validate and normalize mstflow metadata.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    errors = []
    seen = set()
    rows_out = []

    with open(args.input, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail([f"metadata is empty: {args.input}"])

        missing = [column for column in REQUIRED_COLUMNS if column not in reader.fieldnames]
        if missing:
            fail([f"metadata is missing required column(s): {', '.join(missing)}"])

        for line_no, row in enumerate(reader, start=2):
            sample_id = norm(row.get("sample_id"))
            if not sample_id:
                errors.append(f"line {line_no}: sample_id is empty")
                continue
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
            elif not os.path.exists(fastq_1):
                errors.append(f"line {line_no}: fastq_1 does not exist: {fastq_1}")

            if layout == "PE":
                if is_na(fastq_2):
                    errors.append(f"line {line_no}: PE samples require fastq_2")
                elif not os.path.exists(fastq_2):
                    errors.append(f"line {line_no}: fastq_2 does not exist: {fastq_2}")
            elif not is_na(fastq_2):
                errors.append(f"line {line_no}: SE samples must set fastq_2 to NA")

            requested_mode = norm(row.get("preprocess_mode")).lower()
            if not requested_mode:
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

            host_ref = norm(row.get("host_ref")) or "default"
            if not is_na(host_ref) and host_ref != "default" and not os.path.exists(host_ref):
                errors.append(f"line {line_no}: host_ref must be default, NA, or an existing FASTA path: {host_ref}")

            source_sink = norm(row.get("source_sink")).lower()
            model_source_env = norm(row.get("model_source_env"))
            if source_sink and source_sink not in ALLOWED_SOURCE_SINK:
                errors.append(f"line {line_no}: source_sink must be source or sink")
            if (source_sink and is_na(model_source_env)) or (model_source_env and is_na(source_sink)):
                errors.append(f"line {line_no}: source_sink and model_source_env must be provided together")

            rows_out.append(
                normalized_output_row(
                    row,
                    sample_id,
                    fastq_1,
                    fastq_2,
                    layout,
                    age_group,
                    preprocess_mode,
                    host_ref,
                )
            )

    if errors:
        fail(errors)

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows_out)


if __name__ == "__main__":
    main()
