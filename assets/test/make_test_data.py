#!/usr/bin/env python3
"""Generate the bundled MSTFlow test dataset.

Everything here is synthetic and deterministic: random genomes stand in for eight
real genera, a random 'host' chromosome stands in for the host reference, and the
reads are cut from them with fragment lengths typical of ancient and modern
libraries. Taxon IDs and names are real NCBI ones so the outputs read naturally,
but no sequence comes from those organisms.

Run this, then build_test_db.sh, to regenerate everything in this directory.
"""
import gzip
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEED = 20261021
GENOME_LEN = 40_000
HOST_LEN = 60_000

# taxid, name, rank, parent
TAXA = [
    (1, "root", "no rank", 1),
    (131567, "cellular organisms", "no rank", 1),
    (2, "Bacteria", "superkingdom", 131567),
    (2157, "Archaea", "superkingdom", 131567),
    (2759, "Eukaryota", "superkingdom", 131567),
    (561, "Escherichia", "genus", 2), (562, "Escherichia coli", "species", 561),
    (1301, "Streptococcus", "genus", 2), (1309, "Streptococcus mutans", "species", 1301),
    (1386, "Bacillus", "genus", 2), (1423, "Bacillus subtilis", "species", 1386),
    (1678, "Bifidobacterium", "genus", 2), (1680, "Bifidobacterium adolescentis", "species", 1678),
    (286, "Pseudomonas", "genus", 2), (287, "Pseudomonas aeruginosa", "species", 286),
    (1883, "Streptomyces", "genus", 2), (1902, "Streptomyces coelicolor", "species", 1883),
    (2172, "Methanobrevibacter", "genus", 2157), (2173, "Methanobrevibacter smithii", "species", 2172),
    (4930, "Saccharomyces", "genus", 2759), (4932, "Saccharomyces cerevisiae", "species", 4930),
]
SPECIES = {562: 0.50, 1309: 0.37, 1423: 0.43, 1680: 0.60, 287: 0.66, 1902: 0.72, 2173: 0.31, 4932: 0.38}

ADAPTER1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"
ADAPTER2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAA"

# sample, age, layout, read length, fragment range, composition (source -> share)
SAMPLES = [
    ("ancient_pe", "ancient", "PE", 75, (35, 80),
     {1309: .30, 2173: .15, 1680: .10, "host": .30, "noise": .15}),
    ("ancient_se", "ancient", "SE", 75, (35, 75),
     {1423: .25, 1902: .25, 287: .10, "host": .20, "noise": .20}),
    ("modern_pe", "modern", "PE", 100, (220, 400),
     {562: .40, 1680: .20, 4932: .10, "host": .20, "noise": .10}),
    ("modern_se", "modern", "SE", 100, (120, 300),
     {287: .30, 562: .20, 1423: .20, "host": .10, "noise": .20}),
]
N_READS = 3000
COMP = str.maketrans("ACGT", "TGCA")


def genome(rng, length, gc):
    weights = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return "".join(rng.choices("ACGT", weights=weights, k=length))


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def mutate(rng, seq, rate=0.004):
    out = list(seq)
    for i in range(len(out)):
        if rng.random() < rate:
            out[i] = rng.choice([b for b in "ACGT" if b != out[i]])
    return "".join(out)


def fragment(rng, source, lo, hi, sequences):
    length = rng.randint(lo, hi)
    if source == "noise":
        return genome(rng, length, 0.5)
    seq = sequences[source]
    start = rng.randrange(0, len(seq) - length)
    frag = seq[start:start + length]
    return mutate(rng, frag if rng.random() < 0.5 else revcomp(frag))


def write_fastq(path, records):
    with gzip.GzipFile(path, "wb", mtime=0) as gz:
        for name, seq in records:
            gz.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n".encode())


def main():
    rng = random.Random(SEED)
    sequences = {taxid: genome(rng, GENOME_LEN, gc) for taxid, gc in SPECIES.items()}
    sequences["host"] = genome(rng, HOST_LEN, 0.41)

    (HERE / "genomes").mkdir(exist_ok=True)
    with open(HERE / "genomes" / "library.fna", "w") as fh:
        for taxid, seq in ((t, sequences[t]) for t in SPECIES):
            name = next(n for t, n, _, _ in TAXA if t == taxid).replace(" ", "_")
            fh.write(f">{name}|kraken:taxid|{taxid}\n")
            fh.writelines(seq[i:i + 80] + "\n" for i in range(0, len(seq), 80))

    (HERE / "host").mkdir(exist_ok=True)
    with open(HERE / "host" / "host.fa", "w") as fh:
        fh.write(">host_chr1 synthetic host chromosome\n")
        fh.writelines(sequences["host"][i:i + 80] + "\n" for i in range(0, HOST_LEN, 80))

    tax = HERE / "taxonomy"
    tax.mkdir(exist_ok=True)
    with open(tax / "nodes.dmp", "w") as nodes, open(tax / "names.dmp", "w") as names:
        for taxid, name, rank, parent in TAXA:
            nodes.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\n")
            names.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")

    (HERE / "reads").mkdir(exist_ok=True)
    rows = []
    for sample, age, layout, read_len, (lo, hi), mix in SAMPLES:
        sources = list(mix)
        weights = [mix[s] for s in sources]
        r1, r2 = [], []
        for i in range(N_READS):
            frag = fragment(rng, rng.choices(sources, weights)[0], lo, hi, sequences)
            name = f"{sample}.{i + 1}"
            if layout == "PE":
                r1.append((f"{name}/1", (frag + ADAPTER1)[:read_len]))
                r2.append((f"{name}/2", (revcomp(frag) + ADAPTER2)[:read_len]))
            else:
                r1.append((name, (frag + ADAPTER1)[:read_len]))
        if layout == "PE":
            write_fastq(HERE / "reads" / f"{sample}_R1.fastq.gz", r1)
            write_fastq(HERE / "reads" / f"{sample}_R2.fastq.gz", r2)
            rows.append((sample, f"reads/{sample}_R1.fastq.gz", f"reads/{sample}_R2.fastq.gz", age))
        else:
            write_fastq(HERE / "reads" / f"{sample}.fastq.gz", r1)
            rows.append((sample, f"reads/{sample}.fastq.gz", "", age))

    with open(HERE / "samplesheet.tsv", "w") as fh:
        fh.write("sample_id\tfastq_1\tfastq_2\tage_group\n")
        fh.writelines("\t".join(row) + "\n" for row in rows)
    print("wrote synthetic genomes, host, taxonomy, reads and samplesheet to", HERE)


if __name__ == "__main__":
    main()
