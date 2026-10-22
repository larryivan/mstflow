#!/usr/bin/env python3
"""Authenticate the chosen taxa of one ancient library.

The library's alignments are split by taxon. For each taxon this reports how much of
its reference genomes the reads cover, and fits PyDamage's damage model to all of its
reads together: ancient DNA carries C-to-T changes at read ends, a modern contaminant
does not. A taxon is called "damaged" when PyDamage finds damage (q-value below 0.05)
with a predicted accuracy of at least 0.5, PyDamage's own filter.
"""
import argparse
import csv
import subprocess
from pathlib import Path

import numpy as np
import pysam

FIELDS = ["sample_id", "taxid", "name", "table_reads", "aligned_reads", "breadth", "c_to_t_first_base",
          "damage_pmax", "qvalue", "predicted_accuracy", "verdict"]


def read_tsv(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def damage(bam, threads):
    """PyDamage's fit to every read in `bam` together."""
    pysam.index(str(bam))
    outdir = bam.with_suffix("")
    subprocess.run(["pydamage", "--outdir", str(outdir), "analyze", "--group", "--process", str(threads), str(bam)],
                   check=True, stdout=subprocess.DEVNULL)
    with open(outdir / "pydamage_results.csv", newline="") as handle:
        return next(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--references", required=True)
    parser.add_argument("--min-reads", type=int, required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    targets = [row for row in read_tsv(args.targets) if row["sample_id"] == args.sample]
    taxon_of = {row["seqid"]: row["taxid"] for row in read_tsv(args.references)}
    wanted = {row["taxid"] for row in targets if row["has_reference"] == "true"}

    # One pass over the alignments: each read goes to its taxon's file, and marks the
    # reference bases it covers.
    bam = pysam.AlignmentFile(args.bam)
    length = dict(zip(bam.references, bam.lengths))
    split = {taxid: Path(f"taxon_{taxid}.bam") for taxid in wanted}
    writers = {taxid: pysam.AlignmentFile(str(path), "wb", template=bam) for taxid, path in split.items()}
    covered, aligned = {}, dict.fromkeys(wanted, 0)
    for read in bam:
        taxid = taxon_of.get(read.reference_name)
        if taxid not in writers:
            continue
        writers[taxid].write(read)
        aligned[taxid] += 1
        covered.setdefault(read.reference_name, np.zeros(length[read.reference_name], bool))[
            read.reference_start:read.reference_end] = True
    for writer in writers.values():
        writer.close()

    genome_length = dict.fromkeys(wanted, 0)
    for seqid, taxid in taxon_of.items():
        if taxid in genome_length and seqid in length:
            genome_length[taxid] += length[seqid]
    covered_bases = dict.fromkeys(wanted, 0)
    for seqid, bases in covered.items():
        covered_bases[taxon_of[seqid]] += int(bases.sum())

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=FIELDS, restval="NA", lineterminator="\n")
        writer.writeheader()
        for target in targets:
            taxid = target["taxid"]
            row = {key: target[key] for key in ("sample_id", "taxid", "name", "table_reads")}
            if taxid not in wanted:
                writer.writerow({**row, "verdict": "no reference genome"})
                continue
            row.update(aligned_reads=aligned[taxid],
                       breadth=f"{covered_bases[taxid] / genome_length[taxid]:.4g}")
            if aligned[taxid] < args.min_reads:
                writer.writerow({**row, "verdict": "too few reads"})
                continue
            fit = damage(split[taxid], args.threads)
            damaged = float(fit["qvalue"]) < 0.05 and float(fit["predicted_accuracy"]) >= 0.5
            writer.writerow({**row, "c_to_t_first_base": fit["CtoT-0"], "damage_pmax": fit["damage_model_pmax"],
                             "qvalue": fit["qvalue"], "predicted_accuracy": fit["predicted_accuracy"],
                             "verdict": "damaged" if damaged else "no damage"})


if __name__ == "__main__":
    main()
