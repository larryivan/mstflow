#!/usr/bin/env python3
"""Choose the taxa to authenticate in each ancient library, and the reference
sequences that belong to each.

In every library whose age_group is ancient, the taxa with the most reads in the
count table are chosen (at most --top, each with at least --min-reads). A reference
sequence belongs to a taxon when its taxid, from Kraken's `kraken:taxid|N` header
tag, descends from that taxon in the NCBI taxonomy.

Writes targets.tsv (sample_id, taxid, name, table_reads, has_reference) and
references.tsv (seqid, taxid).
"""
import argparse
import csv
import gzip
import re


def table_taxa(path):
    """{taxid: name} and {sample: {taxid: reads}} from a count table."""
    with open(path, newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    samples, names, reads = rows[0][1:], {}, {sample: {} for sample in rows[0][1:]}
    for row in rows[1:]:
        taxid, name = re.fullmatch(r"taxid_(\d+)__(.*)", row[0]).groups()
        names[taxid] = name
        for sample, value in zip(samples, row[1:]):
            reads[sample][taxid] = int(float(value))
    return names, reads


def reference_taxa(fasta, nodes, taxa):
    """{seqid: taxon} for every reference sequence that descends from a taxon in `taxa`."""
    parent = {}
    with open(nodes) as handle:
        for line in handle:
            child, up = line.split("\t|\t")[:2]
            parent[child] = up
    opener = gzip.open if fasta.endswith(".gz") else open
    found = {}
    with opener(fasta, "rt") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            seqid = line[1:].split()[0]
            match = re.search(r"kraken:taxid\|(\d+)", seqid)
            if not match:
                raise SystemExit(f"{fasta}: header without a kraken:taxid| tag: {seqid}")
            taxid, seen = match.group(1), set()
            while taxid not in taxa and taxid in parent and taxid not in seen:
                seen.add(taxid)
                taxid = parent[taxid]
            if taxid in taxa:
                found[seqid] = taxid
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--nodes", required=True)
    parser.add_argument("--top", type=int, required=True)
    parser.add_argument("--min-reads", type=int, required=True)
    args = parser.parse_args()

    names, reads = table_taxa(args.counts)
    with open(args.metadata, newline="") as handle:
        ancient = [row["sample_id"] for row in csv.DictReader(handle, delimiter="\t") if row["age_group"] == "ancient"]
    references = reference_taxa(args.fasta, args.nodes, set(names))
    with_reference = set(references.values())

    with open("targets.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "taxid", "name", "table_reads", "has_reference"])
        for sample in ancient:
            ranked = sorted(reads[sample].items(), key=lambda item: (-item[1], item[0]))
            for taxid, count in [item for item in ranked if item[1] >= args.min_reads][: args.top]:
                writer.writerow([sample, taxid, names[taxid], count, str(taxid in with_reference).lower()])
    with open("references.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["seqid", "taxid"])
        writer.writerows(sorted(references.items()))


if __name__ == "__main__":
    main()
