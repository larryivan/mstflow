#!/usr/bin/env python3
"""One row per sample: how many fragments survive each step, and whether the
profile is deep enough to interpret.

Counts are fragments (read pairs for paired-end libraries, reads for single-end),
so every column is comparable down the row. A taxon counts towards the sample's taxa
only with at least --min-taxon-reads reads: taxa with a read or two say nothing about
depth. Writes a plain TSV for users and a copy with a MultiQC header so the same table
opens the MultiQC report.
"""
import argparse
import csv
import json
import re
from pathlib import Path

PROCESSED = re.compile(r"Processed a total of ([\d,]+) reads")


def read_metric_table(path):
    with open(path, newline="") as handle:
        return {row["metric"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}


def as_int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def collect_qc(files):
    """Index the per-sample QC files by sample and step, using their file names."""
    qc = {}
    for path in map(Path, files):
        name = path.name
        for suffix, step in ((".fastp.json", "fastp"),
                             (".adapterremoval.log", "ar_log"),
                             (".adapterremoval.fragment_qc.tsv", "ar_qc"),
                             (".host_removal.fragment_qc.tsv", "host_qc")):
            if name.endswith(suffix):
                qc.setdefault(name[: -len(suffix)], {})[step] = path
    return qc


def input_and_trimmed(layout, files):
    """Fragments entering preprocessing and leaving it."""
    pairs = 2 if layout == "PE" else 1
    if "fastp" in files:
        summary = json.loads(files["fastp"].read_text())["summary"]
        return (summary["before_filtering"]["total_reads"] // pairs,
                summary["after_filtering"]["total_reads"] // pairs)
    reads_in = None
    if "ar_log" in files:
        match = PROCESSED.search(files["ar_log"].read_text(errors="replace"))
        reads_in = int(match.group(1).replace(",", "")) // pairs if match else None
    trimmed = as_int(read_metric_table(files["ar_qc"])["retained_fragments"]) if "ar_qc" in files else None
    return reads_in, trimmed


def fraction(numerator, denominator):
    if numerator is None or not denominator:
        return ""
    return f"{numerator / denominator:.3f}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--sample-qc", required=True, help="feature-table sample QC TSV")
    parser.add_argument("--counts", required=True, help="count table after the taxon filter")
    parser.add_argument("--rank", default="G")
    parser.add_argument("--min-taxon-reads", type=int, default=10)
    parser.add_argument("--min-classified", type=int, default=10000)
    parser.add_argument("--min-taxa", type=int, default=20)
    parser.add_argument("--output", default="sample_summary.tsv")
    parser.add_argument("qc_files", nargs="*")
    args = parser.parse_args()

    qc = collect_qc(args.qc_files)
    with open(args.sample_qc, newline="") as handle:
        table_qc = {row["sample_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    with open(args.counts, newline="") as handle:
        table = list(csv.reader(handle, delimiter="\t"))
    taxa_per_sample = {sample: sum(float(row[column]) >= args.min_taxon_reads for row in table[1:])
                       for column, sample in enumerate(table[0][1:], start=1)}
    rank_name = {"G": "genera", "S": "species", "F": "families"}.get(args.rank, "taxa")

    header = ["sample_id", "age_group", "layout", "input", "after_trimming", "after_host_removal",
              "host_fraction", "classified", "classified_fraction", rank_name, "flags"]
    rows = []
    with open(args.metadata, newline="") as handle:
        for meta in csv.DictReader(handle, delimiter="\t"):
            sample = meta["sample_id"]
            files = qc.get(sample, {})
            reads_in, trimmed = input_and_trimmed(meta["layout"], files)
            if "host_qc" in files:
                host = read_metric_table(files["host_qc"])
                after_host = as_int(host["retained_fragments"])
                host_removed = as_int(host["removed_host_fragments"])
            else:
                after_host, host_removed = trimmed, 0
            tq = table_qc.get(sample, {})
            classified = as_int(tq.get("kraken_classified_reads"))
            taxa = taxa_per_sample.get(sample)
            flags = []
            if classified is not None and classified < args.min_classified:
                flags.append("few_reads")
            if taxa is not None and taxa < args.min_taxa:
                flags.append(f"few_{rank_name}")
            rows.append([sample, meta["age_group"], meta["layout"], reads_in, trimmed, after_host,
                         fraction(host_removed, trimmed), classified,
                         fraction(classified, as_int(tq.get("kraken_total_reads"))), taxa,
                         ",".join(flags) or "ok"])

    def write(path, preamble=""):
        with open(path, "w", newline="") as handle:
            handle.write(preamble)
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(header)
            writer.writerows([["" if v is None else v for v in row] for row in rows])

    write(args.output)
    write(Path(args.output).with_name("sample_summary_mqc.tsv"), (
        "# id: 'mstflow_sample_summary'\n"
        "# section_name: 'MSTFlow sample summary'\n"
        "# description: 'Fragments surviving each step (read pairs for paired-end "
        "libraries). Flags mark profiles too shallow to interpret: few_reads below "
        f"{args.min_classified:,} classified, few_{rank_name} below {args.min_taxa} "
        f"with {args.min_taxon_reads} reads or more.'\n"
        "# plot_type: 'table'\n"
        "# pconfig:\n"
        "#     id: 'mstflow_sample_summary_table'\n"
        "#     namespace: 'MSTFlow'\n"
        "#     sort_rows: false\n"
    ))


if __name__ == "__main__":
    main()
