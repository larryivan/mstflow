#!/usr/bin/env python3
"""Validate, normalise, and summarise FASTQ fragments for MSTFlow.

A paired-end fragment is one classification unit. Mixed collapsed/unmerged
ancient fragments are emitted as native Kraken2 paired input (``--paired``): a
single-end fragment is paired with a dummy mate shorter than the Kraken2 k-mer
length, which Kraken2 2.1.6 scans into zero minimizer positions, so the dummy
cannot change the call, the confidence denominator, or any statistic below.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import re
import sys
from collections import Counter
from contextlib import ExitStack
from pathlib import Path


DNA_RUN = re.compile(r"[ACGT]+", re.IGNORECASE)
BRACKEN_DISTRIBUTION = re.compile(r"^database(\d+)mers\.kmer_distrib$")
DUMMY_MATE_MARKER = "mstflow_internal_dummy_mate=1"
DEFAULT_GZIP_LEVEL = 6
REVERSE_COMPLEMENT = str.maketrans(
    "ACGTRYKMSWBDHVNacgtrykmswbdhvn",
    "TGCAYRMKSWVHDBNtgcayrmkswvhdbn",
)


def open_text(path: str, mode: str, gzip_level: int = DEFAULT_GZIP_LEVEL):
    """Open a plain or gzipped FASTQ as text.

    ``mtime=0`` is required: ``gzip.open`` stamps ``time.time()`` into the GZIP
    header, so two identical runs would otherwise produce different bytes for
    identical decompressed content.  ``gzip_level`` only affects the compressed
    representation, never the decompressed bytes.
    """
    if str(path).endswith(".gz"):
        if not 0 <= gzip_level <= 9:
            raise ValueError(f"gzip level must be between 0 and 9: {gzip_level}")
        binary = gzip.GzipFile(path, mode + "b", compresslevel=gzip_level, mtime=0)
        return io.TextIOWrapper(binary, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def fastq_records(path: str):
    with open_text(path, "r") as handle:
        record_no = 0
        while True:
            header = handle.readline()
            if not header:
                return
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            record_no += 1
            if not sequence or not plus or not quality:
                raise ValueError(f"{path}: truncated FASTQ record {record_no}")
            header = header.rstrip("\r\n")
            sequence = sequence.rstrip("\r\n")
            plus = plus.rstrip("\r\n")
            quality = quality.rstrip("\r\n")
            if not header.startswith("@"):
                raise ValueError(f"{path}: record {record_no} header does not start with '@'")
            if not plus.startswith("+"):
                raise ValueError(f"{path}: record {record_no} plus line does not start with '+'")
            if len(sequence) != len(quality):
                raise ValueError(
                    f"{path}: record {record_no} sequence/quality lengths differ "
                    f"({len(sequence)} != {len(quality)})"
                )
            yield header, sequence, plus, quality


def canonical_read_id(header: str) -> str:
    token = header[1:].split(maxsplit=1)[0]
    return re.sub(r"(?:/|\.)(?:1|2)$", "", token)


def paired_records(r1_path: str, r2_path: str):
    r1_iter = iter(fastq_records(r1_path))
    r2_iter = iter(fastq_records(r2_path))
    pair_no = 0
    while True:
        r1 = next(r1_iter, None)
        r2 = next(r2_iter, None)
        if r1 is None and r2 is None:
            return
        pair_no += 1
        if r1 is None or r2 is None:
            raise ValueError(
                f"paired FASTQs have different record counts: {r1_path}, {r2_path} "
                f"(first mismatch at pair {pair_no})"
            )
        id1 = canonical_read_id(r1[0])
        id2 = canonical_read_id(r2[0])
        if id1 != id2:
            raise ValueError(
                f"paired FASTQ name mismatch at pair {pair_no}: {r1[0]!r} != {r2[0]!r}"
            )
        yield id1, r1, r2


def sam_groups(path: str):
    """Alignments grouped by read name, from a plain or gzipped SAM file."""
    current_name = None
    records = []
    with open_text(path, "r") as handle:
        for line_no, line in enumerate(handle, start=1):
            if line.startswith("@"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 11:
                raise ValueError(f"{path}: malformed SAM record at line {line_no}")
            try:
                flag = int(fields[1])
                mapq = int(fields[4])
            except ValueError as exc:
                raise ValueError(f"{path}: invalid SAM flag/MAPQ at line {line_no}") from exc
            record = {
                "qname": fields[0],
                "flag": flag,
                "mapq": mapq,
                "sequence": fields[9],
                "quality": fields[10],
            }
            if current_name is not None and record["qname"] != current_name:
                yield current_name, records
                records = []
            current_name = record["qname"]
            records.append(record)
        if current_name is not None:
            yield current_name, records


def sam_record_to_fastq(record: dict, mate: int | None = None):
    sequence = record["sequence"]
    quality = record["quality"]
    if sequence == "*" or quality == "*":
        raise ValueError(f"SAM record {record['qname']!r} does not contain sequence/quality")
    if record["flag"] & 0x10:
        sequence = sequence.translate(REVERSE_COMPLEMENT)[::-1]
        quality = quality[::-1]
    suffix = f"/{mate}" if mate else ""
    return f"@{record['qname']}{suffix}", sequence, "+", quality


def confident_host_alignment(record: dict, minimum_mapq: int) -> bool:
    return not (record["flag"] & 0x4) and record["mapq"] >= minimum_mapq


def filter_host_sam(args):
    paired = bool(args.output_r1 or args.output_r2)
    if paired and not (args.output_r1 and args.output_r2):
        raise ValueError("paired SAM filtering requires --output-r1 and --output-r2")
    if paired and args.output:
        raise ValueError("use either --output or paired --output-r1/--output-r2")
    if not paired and not args.output:
        raise ValueError("single-end SAM filtering requires --output")

    input_fragments = retained_fragments = removed_fragments = 0
    confident_mate1 = confident_mate2 = 0
    with ExitStack() as stack:
        def writer(path):
            return stack.enter_context(open_text(path, "w", args.gzip_level)) if path else None

        output = writer(args.output)
        output_r1 = writer(args.output_r1)
        output_r2 = writer(args.output_r2)

        for qname, alignments in sam_groups(args.sam):
            primary = [
                record
                for record in alignments
                if not (record["flag"] & 0x100) and not (record["flag"] & 0x800)
            ]
            input_fragments += 1
            if paired:
                mate1 = [record for record in primary if record["flag"] & 0x40]
                mate2 = [record for record in primary if record["flag"] & 0x80]
                if len(mate1) != 1 or len(mate2) != 1:
                    raise ValueError(
                        f"SAM fragment {qname!r} has {len(mate1)} primary R1 and "
                        f"{len(mate2)} primary R2 records"
                    )
                hit1 = confident_host_alignment(mate1[0], args.minimum_mapq)
                hit2 = confident_host_alignment(mate2[0], args.minimum_mapq)
                confident_mate1 += int(hit1)
                confident_mate2 += int(hit2)
                if hit1 or hit2:
                    removed_fragments += 1
                    continue
                write_record(output_r1, sam_record_to_fastq(mate1[0], 1))
                write_record(output_r2, sam_record_to_fastq(mate2[0], 2))
            else:
                if len(primary) != 1:
                    raise ValueError(f"SAM fragment {qname!r} has {len(primary)} primary records")
                hit = confident_host_alignment(primary[0], args.minimum_mapq)
                confident_mate1 += int(hit)
                if hit:
                    removed_fragments += 1
                    continue
                write_record(output, sam_record_to_fastq(primary[0]))
            retained_fragments += 1

    return {
        "host_min_mapq": args.minimum_mapq,
        "input_fragments": input_fragments,
        "retained_fragments": retained_fragments,
        "removed_host_fragments": removed_fragments,
        "confident_host_alignments_r1_or_se": confident_mate1,
        "confident_host_alignments_r2": confident_mate2,
    }


def write_record(handle, record):
    header, sequence, plus, quality = record
    handle.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")


def dna_runs(sequence: str) -> list[int]:
    return [len(match.group(0)) for match in DNA_RUN.finditer(sequence)]


def usable_span(sequence: str) -> int:
    return max(dna_runs(sequence), default=0)


def informative_bases(sequence: str) -> int:
    return sum(dna_runs(sequence))


def kmer_windows(sequence: str, kmer_length: int) -> int:
    return sum(max(0, run_length - kmer_length + 1) for run_length in dna_runs(sequence))


def is_dummy_mate(record) -> bool:
    """True for a synthesized sub-k-mer placeholder mate, never a sequenced read."""
    return DUMMY_MATE_MARKER in record[0].split()


def kraken_confidence_denominator(records, kmer_length: int) -> int:
    """Kraken2 2.1.6 ``total_kmers`` for one classification record.

    Verified against kraken2-2.1.6 ``src/classify.cc`` ``ClassifySequence()``:
    every minimizer position is pushed onto ``taxa`` including
    ``AMBIGUOUS_SPAN_TAXON`` positions, and only the mate-pair and
    reading-frame markers are subtracted (lines 827-842).  A sequence therefore
    contributes ``max(0, length - k + 1)`` positions whether or not it contains
    ambiguous bases, and a mate shorter than k contributes exactly zero.
    ``required_score = ceil(confidence * total_kmers)`` in ``ResolveTree()``.
    """
    return sum(max(0, len(record[1]) - kmer_length + 1) for record in records)


def dummy_mate_for(record, length: int, base: str = "N", quality: str = "!"):
    """Build a placeholder mate that shares the fragment ID of ``record``."""
    return (
        f"@{canonical_read_id(record[0])} {DUMMY_MATE_MARKER}",
        base * length,
        "+",
        quality * length,
    )


def fragments(single, r1=None, r2=None):
    """Every fragment, as the tuple of its FASTQ records: single-end files first."""
    if bool(r1) != bool(r2):
        raise ValueError("--r1 and --r2 must be provided together")
    for path in single:
        for record in fastq_records(path):
            yield (record,)
    if r1:
        for _, mate1, mate2 in paired_records(r1, r2):
            yield mate1, mate2


def median(distribution: Counter):
    total = sum(distribution.values())
    if total == 0:
        return 0
    middle_left = (total - 1) // 2
    middle_right = total // 2
    cumulative = 0
    left_value = right_value = None
    for value in sorted(distribution):
        cumulative += distribution[value]
        if left_value is None and cumulative > middle_left:
            left_value = value
        if cumulative > middle_right:
            right_value = value
            break
    middle = (left_value + right_value) / 2
    return int(middle) if middle.is_integer() else middle


def mean(distribution: Counter):
    total = sum(distribution.values())
    if total == 0:
        return 0
    return round(sum(value * count for value, count in distribution.items()) / total, 3)


class FragmentMetrics:
    """Fragment counts and length distributions shared by every subcommand."""

    def __init__(self, kmer_length: int, minimum_usable_span: int = 0):
        self.kmer_length = kmer_length
        self.minimum_usable_span = minimum_usable_span
        self.input = Counter()
        self.retained = Counter()
        self.dropped = 0
        self.informative_lengths = Counter()
        self.equivalent_lengths = Counter()
        self.bracken_read_lengths = Counter()
        self.confidence_denominators = Counter()

    def keep(self, records) -> bool:
        """Account one fragment, given as the records submitted to Kraken2 for it.

        Returns whether its longest usable span reaches the minimum. Statistics
        come from the sequenced mates only; a placeholder mate adds nothing.
        """
        mates = [record[1] for record in records if not is_dummy_mate(record)]
        kind = "paired" if len(mates) == 2 else "single"
        self.input[kind] += 1
        if max(map(usable_span, mates)) < self.minimum_usable_span:
            self.dropped += 1
            return False
        self.retained[kind] += 1
        self.informative_lengths[sum(map(informative_bases, mates))] += 1
        # A contiguous read of this length supplies the same number of k-mer windows.
        windows = sum(kmer_windows(mate, self.kmer_length) for mate in mates)
        self.equivalent_lengths[windows + self.kmer_length - 1 if windows else 0] += 1
        # Bracken's -r is the length of one sequenced read, also for paired-end
        # data, so each mate is its own observation.
        self.bracken_read_lengths.update(map(len, mates))
        self.confidence_denominators[kraken_confidence_denominator(records, self.kmer_length)] += 1
        return True

    def summary(self) -> dict:
        return {
            "input_single_fragments": self.input["single"],
            "input_paired_fragments": self.input["paired"],
            "input_fragments": sum(self.input.values()),
            "retained_single_fragments": self.retained["single"],
            "retained_paired_fragments": self.retained["paired"],
            "retained_fragments": sum(self.retained.values()),
            "dropped_below_minimum_usable_span": self.dropped,
            "minimum_usable_span": self.minimum_usable_span,
            "kraken_kmer_length": self.kmer_length,
            "median_informative_bases": median(self.informative_lengths),
            "mean_informative_bases": mean(self.informative_lengths),
            "median_kmer_equivalent_read_length": median(self.equivalent_lengths),
            "mean_kmer_equivalent_read_length": mean(self.equivalent_lengths),
            "median_bracken_read_length": median(self.bracken_read_lengths),
            "mean_bracken_read_length": mean(self.bracken_read_lengths),
            # Kraken2's confidence denominator for the records this step submits:
            # required_score = ceil(kraken2_confidence * this value).
            "median_kraken_confidence_denominator": median(self.confidence_denominators),
            "mean_kraken_confidence_denominator": mean(self.confidence_denominators),
        }


def write_metrics(path: str, metrics: dict):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "value"])
        writer.writerows(metrics.items())


def prepare_single(args):
    metrics = FragmentMetrics(args.kmer_length, args.minimum_usable_span)
    with open_text(args.output, "w", args.gzip_level) as output:
        for records in fragments(args.single):
            if metrics.keep(records):
                write_record(output, records[0])
    return metrics.summary()


def pair_fragments(args):
    """Emit mixed single-end and paired fragments as native Kraken2 paired input.

    Single-end fragments are paired with a placeholder mate shorter than the
    k-mer length so Kraken2 scans zero minimizer positions from it.  One
    fragment therefore remains exactly one Kraken2 classification record with
    the same confidence denominator as native single-end classification.
    """
    if not 0 < args.dummy_mate_length < args.kmer_length:
        raise ValueError(
            f"dummy mate length {args.dummy_mate_length} must be positive and shorter than "
            f"the k-mer length {args.kmer_length}; a mate of at least k contributes "
            "minimizer positions and would change the confidence denominator"
        )
    metrics = FragmentMetrics(args.kmer_length, args.minimum_usable_span)
    with ExitStack() as stack:
        output_r1 = stack.enter_context(open_text(args.output_r1, "w", args.gzip_level))
        output_r2 = stack.enter_context(open_text(args.output_r2, "w", args.gzip_level))
        for records in fragments(args.single, args.r1, args.r2):
            if len(records) == 1:
                records = (records[0], dummy_mate_for(records[0], args.dummy_mate_length))
            if metrics.keep(records):
                write_record(output_r1, records[0])
                write_record(output_r2, records[1])
    return metrics.summary()


def collect_stats(args):
    metrics = FragmentMetrics(args.kmer_length, args.minimum_usable_span)
    for records in fragments(args.single, args.r1, args.r2):
        metrics.keep(records)
    return metrics.summary()


def read_metric(path: str, name: str) -> float:
    with open(path, encoding="utf-8", newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        for row in rows:
            if row.get("metric") == name:
                return float(row["value"])
    raise ValueError(f"{path}: missing required metric {name!r}")


def select_bracken_length(args):
    target = read_metric(args.stats, "median_bracken_read_length")
    available = {}
    for candidate in Path(args.database).glob("database*mers.kmer_distrib"):
        match = BRACKEN_DISTRIBUTION.match(candidate.name)
        if match:
            # The selected model is an interface between two separate
            # processes. Record only its stable database-relative basename:
            # task-local absolute staging paths vary across retries and would
            # otherwise perturb checksums and downstream cache keys.
            available[int(match.group(1))] = candidate.name
    if not available:
        raise ValueError(
            f"no database{{N}}mers.kmer_distrib files found in Bracken database: {args.database}"
        )

    requested = str(args.requested).strip().lower()
    if requested == "auto":
        if target <= 0:
            raise ValueError(f"cannot select Bracken distribution from zero-length reads in {args.stats}")
        # Stable tie break: prefer the shorter distribution at equal distance.
        selected = min(available, key=lambda value: (abs(value - target), value))
        selection_mode = "nearest_per_sample"
    else:
        try:
            selected = int(requested)
        except ValueError as exc:
            raise ValueError("Bracken read length must be 'auto' or a positive integer") from exc
        if selected <= 0:
            raise ValueError("Bracken read length must be a positive integer")
        if selected not in available:
            options = ", ".join(str(value) for value in sorted(available))
            raise ValueError(
                f"requested Bracken read length {selected} is unavailable in {args.database}; "
                f"available lengths: {options}"
            )
        selection_mode = "fixed"

    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "sample_id",
                "observed_median_bracken_read_length",
                "selected_bracken_read_length",
                "selection_mode",
                "distribution_file",
            ]
        )
        writer.writerow(
            [args.sample_id, f"{target:.12g}", selected, selection_mode, available[selected]]
        )


def add_gzip_level_option(parser):
    # Compression level changes only the compressed representation. It is an
    # explicit int with no truthiness fallback so that level 0 (stored) works.
    parser.add_argument(
        "--gzip-level",
        type=int,
        default=DEFAULT_GZIP_LEVEL,
        choices=range(0, 10),
        metavar="{0..9}",
        help="gzip compression level for FASTQ outputs (default: %(default)s)",
    )


def add_common_stats_options(parser):
    parser.add_argument("--minimum-usable-span", type=int, default=35)
    parser.add_argument("--kmer-length", type=int, default=35)
    parser.add_argument("--stats", required=True)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    single = subparsers.add_parser("prepare-single", help="validate, filter, and concatenate SE FASTQs")
    single.add_argument("--single", action="append", required=True)
    single.add_argument("--output", required=True)
    add_common_stats_options(single)
    add_gzip_level_option(single)
    single.set_defaults(func=prepare_single)

    pair = subparsers.add_parser(
        "pair-fragments",
        help="emit SE and paired fragments as native Kraken2 --paired input",
    )
    pair.add_argument("--single", action="append", default=[])
    pair.add_argument("--r1")
    pair.add_argument("--r2")
    pair.add_argument("--output-r1", required=True)
    pair.add_argument("--output-r2", required=True)
    pair.add_argument(
        "--dummy-mate-length",
        type=int,
        default=1,
        help="placeholder mate length for SE fragments; must be shorter than --kmer-length",
    )
    add_common_stats_options(pair)
    add_gzip_level_option(pair)
    pair.set_defaults(func=pair_fragments)

    stats = subparsers.add_parser("stats", help="validate FASTQs and report fragment-level statistics")
    stats.add_argument("--single", action="append", default=[])
    stats.add_argument("--r1")
    stats.add_argument("--r2")
    stats.add_argument("--output", dest="stats", required=True)
    stats.add_argument("--minimum-usable-span", type=int, default=0)
    stats.add_argument("--kmer-length", type=int, default=35)
    stats.set_defaults(func=collect_stats)

    host = subparsers.add_parser(
        "filter-host-sam",
        help="retain SE fragments or complete pairs without a confident host alignment",
    )
    host.add_argument("--sam", required=True, help="query-name-grouped SAM, plain or gzipped")
    host.add_argument("--output", help="single-end FASTQ output")
    host.add_argument("--output-r1", help="paired R1 FASTQ output")
    host.add_argument("--output-r2", help="paired R2 FASTQ output")
    host.add_argument("--minimum-mapq", type=int, default=0)
    host.add_argument("--stats", required=True)
    add_gzip_level_option(host)
    host.set_defaults(func=filter_host_sam)

    bracken = subparsers.add_parser("select-bracken-length", help="select a reproducible Bracken model")
    bracken.add_argument("--stats", required=True)
    bracken.add_argument("--database", required=True)
    bracken.add_argument("--requested", default="auto")
    bracken.add_argument("--sample-id", required=True)
    bracken.add_argument("--output", required=True)
    bracken.set_defaults(func=select_bracken_length)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    if getattr(args, "minimum_usable_span", 0) < 0:
        parser.error("--minimum-usable-span must be non-negative")
    if getattr(args, "kmer_length", 1) <= 0:
        parser.error("--kmer-length must be positive")
    if getattr(args, "minimum_mapq", 0) < 0:
        parser.error("--minimum-mapq must be non-negative")
    try:
        # Every subcommand but select-bracken-length returns metrics for --stats.
        metrics = args.func(args)
        if metrics is not None:
            write_metrics(args.stats, metrics)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc


if __name__ == "__main__":
    main()
