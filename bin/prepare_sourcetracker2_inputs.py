#!/usr/bin/env python3
import argparse
import csv
import sys


def is_na(value):
    return value is None or str(value).strip().upper() in {"", "NA", "N/A", "NONE", "NULL"}


def read_metadata(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    by_sample = {}
    errors = []
    for row in rows:
        sample_id = row["sample_id"]
        source_sink = row.get("source_sink", "").strip().lower()
        env = row.get("model_source_env", "").strip()
        if source_sink not in {"source", "sink"}:
            errors.append(f"{sample_id}: source_sink must be source or sink for SourceTracker2")
        if is_na(env):
            errors.append(f"{sample_id}: model_source_env is required for SourceTracker2")
        by_sample[sample_id] = row
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
    return by_sample


def read_counts(path):
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        rows = [row for row in reader if row]
    if not header or header[0] != "#OTU ID":
        raise SystemExit(f"{path} must be a tab-delimited feature table with first column '#OTU ID'")
    return header, rows


def main():
    parser = argparse.ArgumentParser(description="Prepare BIOM count table and mapping file for SourceTracker2.")
    parser.add_argument("--counts", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--output-counts", required=True)
    parser.add_argument("--output-map", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    metadata = read_metadata(args.metadata)
    header, rows = read_counts(args.counts)
    samples = header[1:]
    missing = [sample for sample in samples if sample not in metadata]
    if missing:
        raise SystemExit(f"count table contains sample(s) absent from metadata: {', '.join(missing)}")

    sample_totals = {sample: 0.0 for sample in samples}
    for row in rows:
        for sample, value in zip(samples, row[1:]):
            sample_totals[sample] += float(value)

    kept_samples = [sample for sample in samples if sample_totals[sample] > 0]
    dropped_samples = [sample for sample in samples if sample_totals[sample] <= 0]
    kept_indices = [samples.index(sample) for sample in kept_samples]
    if not kept_samples:
        raise SystemExit("no nonzero samples remain for SourceTracker2")
    if not any(metadata[sample]["source_sink"] == "source" for sample in kept_samples):
        raise SystemExit("SourceTracker2 requires at least one nonzero source sample")
    if not any(metadata[sample]["source_sink"] == "sink" for sample in kept_samples):
        raise SystemExit("SourceTracker2 requires at least one nonzero sink sample")

    with open(args.output_counts, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["#OTU ID", *kept_samples])
        for row in rows:
            kept_values = [row[index + 1] for index in kept_indices]
            if any(float(value) > 0 for value in kept_values):
                writer.writerow([row[0], *kept_values])

    with open(args.output_map, "w", newline="") as handle:
        fieldnames = [
            "#SampleID",
            "SourceSink",
            "Env",
            "source_major",
            "source_class",
            "sample_label",
            "project",
            "experiment",
            "run",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for sample in kept_samples:
            row = metadata[sample]
            writer.writerow(
                {
                    "#SampleID": sample,
                    "SourceSink": row["source_sink"],
                    "Env": row["model_source_env"],
                    "source_major": row.get("source_major", "NA"),
                    "source_class": row.get("source_class", "NA"),
                    "sample_label": row.get("sample_label", sample),
                    "project": row.get("project", "NA"),
                    "experiment": row.get("experiment", "NA"),
                    "run": row.get("run", sample),
                }
            )

    with open(args.summary, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["input_samples", len(samples)])
        writer.writerow(["kept_samples", len(kept_samples)])
        writer.writerow(["dropped_zero_samples", len(dropped_samples)])
        writer.writerow(["source_samples", sum(1 for sample in kept_samples if metadata[sample]["source_sink"] == "source")])
        writer.writerow(["sink_samples", sum(1 for sample in kept_samples if metadata[sample]["source_sink"] == "sink")])
        if dropped_samples:
            writer.writerow(["dropped_zero_sample_ids", ",".join(dropped_samples)])


if __name__ == "__main__":
    main()
