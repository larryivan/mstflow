#!/usr/bin/env python3
"""Build the count tables from each sample's Kraken2 or KrakenUniq report: the reads
the classifier placed in each taxon of one rank or below it, merged across samples,
before and after the taxon filter, with classification and lineage QC.

Each read is classified on its own, so a sample's counts are the sum of what its reads
contribute, and a sink made of two sources gets the sum of their counts.
"""
import argparse
import csv
import math
from decimal import Decimal, InvalidOperation
from pathlib import Path


PROKARYOTE_ROOTS = {"2", "2157"}
# A read count wider than a 128-bit integer is not physically obtainable; the
# bound exists so an implausible magnitude fails fast instead of being expanded.
MAX_COUNT_DIGITS = 38
ROOT_TAXIDS = {"0", "1"}
RANK_WORDS = {"D": ("superkingdom", "domain"), "P": ("phylum",), "C": ("class",), "O": ("order",),
              "F": ("family",), "G": ("genus",), "S": ("species",)}


def read_samples(path):
    """The sample IDs of the samplesheet, in its order."""
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "sample_id" not in (reader.fieldnames or []):
            raise SystemExit(f"{path} is missing metadata column: sample_id")
        rows = list(reader)

    samples = []
    for line_no, row in enumerate(rows, start=2):
        sample = (row.get("sample_id") or "").strip()
        if not sample:
            raise SystemExit(f"{path}:{line_no}: sample_id is empty")
        if sample in samples:
            raise SystemExit(f"{path}:{line_no}: duplicate sample_id '{sample}'")
        samples.append(sample)
    if not samples:
        raise SystemExit(f"metadata contains no samples: {path}")
    return samples


def id_from_path(path, *suffixes):
    """The sample ID of a per-sample file named <sample><suffix>."""
    name = Path(path).name
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    raise SystemExit(f"{path}: expected a file name ending in {' or '.join(suffixes)}")


def feature_id(taxid, name):
    return f"taxid_{taxid}__{name.strip().replace(chr(9), ' ')}"


def finite_nonnegative(value, context):
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise SystemExit(f"{context} is not numeric: {value!r}")
    if not math.isfinite(number) or number < 0:
        raise SystemExit(f"{context} must be finite and non-negative: {value!r}")
    return number


def exact_count(value, context):
    """Parse a read count as an exact integer.

    Accepts the integer-valued spellings tools actually emit (``8``, ``8.0``,
    ``8e0``) and rejects fractional, negative, NaN, and infinite values. Counts
    are kept as Python ints so summing and formatting stay exact; only relative
    abundances are floating point.
    """
    text = "" if value is None else str(value).strip()
    try:
        number = Decimal(text)
    except (InvalidOperation, ValueError):
        raise SystemExit(f"{context} is not numeric: {value!r}")
    if not number.is_finite() or number < 0:
        raise SystemExit(f"{context} must be finite and non-negative: {value!r}")
    # Bound the magnitude before int() materialises the digits. Decimal treats
    # a token such as "1e999999999" as finite, and expanding it would allocate a
    # billion-digit integer inside a non-interruptible C call; the previous
    # float-based parser collapsed it to inf and rejected it cleanly.
    if number.adjusted() >= MAX_COUNT_DIGITS:
        raise SystemExit(
            f"{context} needs more than {MAX_COUNT_DIGITS} digits "
            f"and cannot be a read count: {value!r}"
        )
    if number != number.to_integral_value():
        raise SystemExit(f"{context} must be a whole read count: {value!r}")
    return int(number)


def format_count(value):
    """Render a count exactly, with no float rounding or exponent notation."""
    return "NA" if value is None else str(int(value))


def read_kraken_report(path, rank="G"):
    """Read a Kraken2 or KrakenUniq report.

    Returns the sample, its classified and unclassified totals, its clade read counts
    at `rank`, their features, and, for KrakenUniq,
    the unique k-mers each taxon's reads cover (None for Kraken2, whose distinct
    minimizers come in a separate report). Kraken2 rows are percent, clade reads, taxon
    reads, rank code, taxid, name; KrakenUniq rows, after a header, are %, reads,
    taxReads, kmers, dup, cov, taxID, rank, taxName.
    """
    sample = id_from_path(path, ".kraken2.report.txt", ".krakenuniq.report.txt")
    totals, percents = {}, {}
    counts, features = {}, {}
    distinct = None
    header = None
    with open(path) as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[:2] == ["%", "reads"]:
                if "taxID" not in parts or "kmers" not in parts:
                    raise SystemExit(f"{path}:{line_no}: not a KrakenUniq report header")
                header, distinct = parts, {}
                continue
            if header:
                row = dict(zip(header, parts))
                percent, clade_reads = row["%"], row["reads"]
                taxid, name, rank_name = row["taxID"].strip(), row["taxName"].strip(), row["rank"].strip()
                at_rank = rank_name in RANK_WORDS[rank]
                distinct[taxid] = exact_count(row["kmers"].strip(), f"{path}:{line_no}:kmers")
            else:
                if len(parts) < 6:
                    raise SystemExit(f"{path}:{line_no}: expected at least six Kraken report columns")
                # Six columns, or eight with minimizer data; both end in rank, taxid, name.
                percent, clade_reads = parts[0], parts[1]
                rank_code, taxid, name = parts[-3].strip(), parts[-2].strip(), parts[-1].strip()
                at_rank = rank_code == rank
                if rank_code == "U":
                    taxid = "0"
            reads = exact_count(clade_reads.strip(), f"{path}:{line_no}:clade_reads")
            if taxid in ("0", "1"):
                totals[taxid] = reads
                percents[taxid] = finite_nonnegative(percent.strip(), f"{path}:{line_no}:percent")
            elif at_rank:
                fid = feature_id(taxid, name)
                counts[fid] = counts.get(fid, 0) + reads
                features[fid] = {"feature_id": fid, "taxid": taxid, "name": name, "rank": rank}
    # Kraken legitimately omits the root row for an all-unclassified sample, and may
    # omit the unclassified row for a fully classified sample.
    unclassified, classified = totals.get("0"), totals.get("1")
    if classified is None and unclassified is not None and percents["0"] >= 99.999:
        classified = 0
    if unclassified is None and classified is not None and percents["1"] >= 99.999:
        unclassified = 0
    if classified is None or unclassified is None:
        raise SystemExit(f"{path}: Kraken report must contain unclassified (taxid 0) and root (taxid 1) rows")
    metrics = {
        "kraken_total_reads": classified + unclassified,
        "kraken_classified_reads": classified,
        "kraken_unclassified_reads": unclassified,
    }
    return sample, metrics, counts, features, distinct


def read_minimizer_report(path):
    """Distinct minimizers per taxid from a Kraken2 report written with
    --report-minimizer-data: percent, clade reads, taxon reads, minimizers, distinct
    minimizers, rank, taxid, name. The minimizer columns count the whole clade."""
    sample = id_from_path(path, ".kraken2.minimizer.report.txt")
    distinct = {}
    with open(path) as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                raise SystemExit(f"{path}:{line_no}: expected the eight columns of a Kraken2 report with minimizer data")
            distinct[parts[6].strip()] = exact_count(parts[4].strip(), f"{path}:{line_no}:distinct_minimizers")
    return sample, distinct


def filter_by_distinct_kmers(counts, features, distinct_by_sample, min_distinct):
    """Zero, sample by sample, every taxon whose reads cover fewer than `min_distinct`
    distinct k-mers of its genomes.

    Reads piled on a few short stretches of a genome (a conserved gene, a low-complexity
    region, a contaminant fragment in the reference) give many reads but few distinct
    k-mers; reads from an organism that is really there are spread along it. KrakenUniq
    reports unique k-mers for this; Kraken2 reports distinct minimizers, a sample of the
    same thing. Returns the removals.
    """
    removed = []
    for fid in sorted(counts):
        taxid = features[fid]["taxid"]
        for sample, count in sorted(counts[fid].items()):
            distinct = distinct_by_sample[sample].get(taxid, 0)
            if count and distinct < min_distinct:
                removed.append([sample, fid, features[fid]["name"], format_count(count), format_count(distinct)])
                counts[fid][sample] = 0
    return removed


def require_one_to_one(expected, observed, label):
    duplicates = sorted({sample for sample in observed if observed.count(sample) > 1})
    missing = sorted(set(expected).difference(observed))
    unexpected = sorted(set(observed).difference(expected))
    problems = []
    if duplicates:
        problems.append(f"duplicate {label} for: {', '.join(duplicates)}")
    if missing:
        problems.append(f"missing {label} for: {', '.join(missing)}")
    if unexpected:
        problems.append(f"{label} absent from metadata: {', '.join(unexpected)}")
    if problems:
        raise SystemExit("; ".join(problems))


def read_nodes(nodes):
    if not nodes:
        raise SystemExit("--taxdump nodes.dmp is required when --taxon-filter prokaryote")
    parents = {}
    with open(nodes) as handle:
        for line in handle:
            parts = [part.strip() for part in line.split("|")]
            if len(parts) >= 2:
                parents[parts[0]] = parts[1]
    return parents


def lineage_has(taxid, targets, parents, cache):
    if taxid in cache:
        return cache[taxid]
    seen = set()
    current = taxid
    while current and current not in seen:
        if current in targets:
            cache[taxid] = True
            return True
        if current in ROOT_TAXIDS or current not in parents:
            break
        seen.add(current)
        current = parents[current]
    cache[taxid] = False
    return False


def select_features(features, mode, taxdump, drop_taxids):
    parents = read_nodes(taxdump) if mode == "prokaryote" else {}
    cache = {}
    reasons = {}
    for fid, feature in features.items():
        taxid = feature["taxid"]
        if taxid in drop_taxids:
            reasons[fid] = "drop_taxid"
        elif mode == "prokaryote" and not lineage_has(taxid, PROKARYOTE_ROOTS, parents, cache):
            reasons[fid] = "not_prokaryote"
        else:
            reasons[fid] = "kept"
    kept = [fid for fid, reason in reasons.items() if reason == "kept"]
    return sorted(kept), reasons


def sample_totals(samples, feature_ids, counts):
    return {
        sample: sum(counts.get(fid, {}).get(sample, 0) for fid in feature_ids)
        for sample in samples
    }


def write_matrix(path, samples, features, counts, relative=False):
    totals = sample_totals(samples, features, counts)
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["#OTU ID", *samples])
        for fid in features:
            row = [fid]
            for sample in samples:
                count = counts.get(fid, {}).get(sample, 0)
                if relative:
                    total = totals[sample]
                    row.append(f"{count / total if total else 0.0:.12g}")
                else:
                    row.append(format_count(count))
            writer.writerow(row)


def write_features(path, feature_ids, features):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=["feature_id", "taxid", "name", "rank"],
            lineterminator="\n",
        )
        writer.writeheader()
        for fid in feature_ids:
            writer.writerow(features[fid])


def write_qc(path, samples, raw_features, filtered_features, reasons, mode):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerow(["taxon_filter", mode])
        writer.writerow(["samples", len(samples)])
        writer.writerow(["raw_features", len(raw_features)])
        writer.writerow(["filtered_features", len(filtered_features)])
        writer.writerow(["dropped_features", len(raw_features) - len(filtered_features)])
        for reason in sorted(set(reasons.values())):
            writer.writerow([f"features_{reason}", sum(value == reason for value in reasons.values())])


def ratio(numerator, denominator):
    return numerator / denominator if denominator else None


def format_value(value):
    return "NA" if value is None else f"{value:.12g}"


def write_sample_qc(
    path,
    samples,
    raw_features,
    filtered_features,
    counts,
    kraken_metrics,
    kmer_filter_removed=None,
):
    raw_totals = sample_totals(samples, raw_features, counts)
    filtered_totals = sample_totals(samples, filtered_features, counts)
    fields = [
        "sample_id",
        "kraken_total_reads",
        "kraken_classified_reads",
        "kraken_unclassified_reads",
        "kraken_classified_fraction",
        "rank_total_reads",
        "rank_fraction_of_classified",
        "taxon_filter_retained_reads",
        "taxon_filter_retained_fraction",
        "raw_nonzero_features",
        "filtered_nonzero_features",
        "kmer_filter_removed_reads",
    ]
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for sample in samples:
            raw_total = raw_totals[sample]
            filtered_total = filtered_totals[sample]
            km = kraken_metrics.get(sample, {})
            classified = km.get("kraken_classified_reads")
            total = km.get("kraken_total_reads")
            writer.writerow(
                {
                    "sample_id": sample,
                    "kraken_total_reads": format_count(total),
                    "kraken_classified_reads": format_count(classified),
                    "kraken_unclassified_reads": format_count(km.get("kraken_unclassified_reads")),
                    "kraken_classified_fraction": format_value(
                        ratio(classified, total) if classified is not None and total is not None else None
                    ),
                    "rank_total_reads": format_count(raw_total),
                    "rank_fraction_of_classified": format_value(
                        ratio(raw_total, classified) if classified is not None else None
                    ),
                    "taxon_filter_retained_reads": format_count(filtered_total),
                    "taxon_filter_retained_fraction": format_value(ratio(filtered_total, raw_total)),
                    "raw_nonzero_features": sum(
                        counts.get(fid, {}).get(sample, 0) > 0 for fid in raw_features
                    ),
                    "filtered_nonzero_features": sum(
                        counts.get(fid, {}).get(sample, 0) > 0 for fid in filtered_features
                    ),
                    "kmer_filter_removed_reads": (
                        "NA" if kmer_filter_removed is None else format_count(kmer_filter_removed.get(sample, 0))
                    ),
                }
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--taxon-filter", choices=["none", "prokaryote"], default="none")
    parser.add_argument("--taxdump")
    parser.add_argument("--drop-taxids", default="0,1")
    parser.add_argument(
        "--minimizer-reports",
        nargs="+",
        default=None,
        help="Per-sample Kraken2 reports with minimizer data, for --min-distinct-kmers.",
    )
    parser.add_argument("--rank", choices=sorted(RANK_WORDS), default="G")
    parser.add_argument(
        "--min-distinct-kmers",
        type=int,
        default=0,
        help="Drop a taxon from a sample when its reads cover fewer distinct k-mers or, for Kraken2, "
             "minimizers (0: keep all).",
    )
    parser.add_argument("reports", nargs="+", help="Per-sample Kraken2 or KrakenUniq reports, one per sample.")
    args = parser.parse_args()

    samples = read_samples(args.metadata)
    counts = {}
    features = {}
    kraken_metrics = {}
    distinct_by_sample = {}

    # Sorted so the table is byte-identical regardless of the order in which
    # the workflow hands over per-sample reports.
    for path in sorted(args.reports, key=lambda value: Path(value).name):
        sample, metrics, report_counts, report_features, distinct = read_kraken_report(path, args.rank)
        if sample in kraken_metrics:
            raise SystemExit(f"duplicate Kraken report for sample {sample}")
        kraken_metrics[sample] = metrics
        features.update(report_features)
        for fid, count in report_counts.items():
            counts.setdefault(fid, {})[sample] = count
        if distinct is not None:
            distinct_by_sample[sample] = distinct
    require_one_to_one(samples, list(kraken_metrics), "Kraken report")
    if args.minimizer_reports:
        distinct_by_sample = {}
        for path in args.minimizer_reports:
            sample, distinct = read_minimizer_report(path)
            distinct_by_sample[sample] = distinct
        require_one_to_one(samples, list(distinct_by_sample), "Kraken2 minimizer report")
    if args.min_distinct_kmers > 0 and not distinct_by_sample:
        raise SystemExit("--min-distinct-kmers needs KrakenUniq reports or, for Kraken2, --minimizer-reports")

    filtering = args.min_distinct_kmers > 0
    removed = []
    if filtering:
        removed = filter_by_distinct_kmers(counts, features, distinct_by_sample, args.min_distinct_kmers)
    with open(f"{args.prefix}_kmer_filter.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "feature_id", "name", "removed_reads", "distinct_kmers"])
        writer.writerows(removed)
    removed_reads = {}
    for sample, _fid, _name, reads, _distinct in removed:
        removed_reads[sample] = removed_reads.get(sample, 0) + int(reads)

    raw_features = sorted(features)
    drop_taxids = {taxid.strip() for taxid in args.drop_taxids.split(",") if taxid.strip()}
    filtered_features, reasons = select_features(features, args.taxon_filter, args.taxdump, drop_taxids)

    write_matrix(f"{args.prefix}_counts.tsv", samples, raw_features, counts, relative=False)
    write_matrix(f"{args.prefix}_relative_abundance.tsv", samples, raw_features, counts, relative=True)
    write_features(f"{args.prefix}_features.tsv", raw_features, features)

    filtered_prefix = f"{args.prefix}_{args.taxon_filter}"
    write_matrix(f"{filtered_prefix}_counts.tsv", samples, filtered_features, counts, relative=False)
    write_matrix(f"{filtered_prefix}_relative_abundance.tsv", samples, filtered_features, counts, relative=True)
    write_features(f"{filtered_prefix}_features.tsv", filtered_features, features)
    write_qc(
        f"{filtered_prefix}_filter_qc.tsv",
        samples,
        raw_features,
        filtered_features,
        reasons,
        args.taxon_filter,
    )
    write_sample_qc(
        f"{filtered_prefix}_sample_qc.tsv",
        samples,
        raw_features,
        filtered_features,
        counts,
        kraken_metrics,
        removed_reads if filtering else None,
    )


if __name__ == "__main__":
    main()
