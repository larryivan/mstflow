#!/usr/bin/env python3
"""Download everything the test dataset is built from, into build/.

- reads: the first 50,000 reads (pairs) of each library in libraries.tsv, from ENA
- genomes: the reference genomes pinned in genomes.tsv, from NCBI Datasets
- taxonomy: the NCBI lineages of those genomes, cut from the current taxdump (nodes.dmp,
  names.dmp, and an empty merged.dmp, which SourcePredict's taxonomy reader expects)
- host: the human mitochondrial genome (NC_012920.1)
- mixtures: sinks of known composition, made from the reads of the libraries in
  mixtures.tsv (see that section of README.md)
- samplesheet.tsv, pointing at the reads by file name

Then run build_test_db.sh to build the Kraken2/Bracken database and the release files.
Only the Python standard library is used.
"""
import csv
import gzip
import io
import math
import random
import sys
import tarfile
import time
import urllib.request
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
BUILD = HERE / "build"
READS_PER_LIBRARY = 50_000
ADAPTERS = ("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG",
            "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT")
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def fetch(url, timeout=300):
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                return response.read()
        except Exception as exc:
            print(f"  retrying ({exc})", file=sys.stderr)
            time.sleep(3 * (attempt + 1))
    raise SystemExit(f"could not download {url}")


def ena_fastq_urls(run, layout):
    digits = run[3:]
    sub = "" if len(digits) <= 6 else "/" + digits[6:].rjust(3, "0")
    base = f"https://ftp.sra.ebi.ac.uk/vol1/fastq/{run[:6]}{sub}/{run}/{run}"
    return [f"{base}_1.fastq.gz", f"{base}_2.fastq.gz"] if layout == "PE" else [f"{base}.fastq.gz"]


def head_fastq(url, out, n):
    """Stream a gzipped FASTQ and keep its first n records, stopping once they are read."""
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=300) as response, \
                    gzip.open(response, "rt") as fin, gzip.GzipFile(out, "wb", mtime=0) as fout:
                for _ in range(n):
                    record = [fin.readline() for _ in range(4)]
                    if not record[3]:
                        break
                    if not (record[0].startswith("@") and record[2].startswith("+")):
                        raise ValueError(f"not FASTQ: {url}")
                    fout.write("".join(record).encode())
            return
        except Exception as exc:
            print(f"  retrying {url} ({exc})", file=sys.stderr)
            time.sleep(3 * (attempt + 1))
    raise SystemExit(f"could not download {url}")


def reads(libraries):
    out = BUILD / "reads"
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for lib in libraries:
        urls = ena_fastq_urls(lib["run"], lib["layout"])
        names = [f"{lib['sample_id']}_R{i}.fastq.gz" for i in (1, 2)] if len(urls) == 2 else [f"{lib['sample_id']}.fastq.gz"]
        for url, name in zip(urls, names):
            head_fastq(url, out / name, READS_PER_LIBRARY)
        files = names
        print(f"reads   {lib['sample_id']:10s} {lib['run']}", flush=True)
        rows.append([lib["sample_id"], files[0], files[1] if len(files) > 1 else "", lib["age_group"],
                     lib["host_action"], lib["source_sink"], lib["model_source_env"], lib["project"], lib["run"]])
    return rows


def fastq(path):
    with gzip.open(path, "rt") as handle:
        while True:
            record = [handle.readline().rstrip("\n") for _ in range(4)]
            if not record[3]:
                return
            yield record


def ancient_like(read1, rng):
    """An ancient-like molecule from a modern read 1: its first bases, cut to a length
    typical of ancient DNA (median 65 bp) and deaminated at both ends, sequenced as a pair
    at the original read length so that both reads run into the adapters."""
    header, sequence, _, quality = read1
    length = min(len(sequence), max(30, round(rng.lognormvariate(math.log(65), 0.35))))
    bases = list(sequence[:length])
    for i in range(length):
        p = 0.30 * 0.55 ** i
        if p < 0.001:
            break
        if bases[i] == "C" and rng.random() < p:
            bases[i] = "T"
        if bases[-1 - i] == "G" and rng.random() < p:
            bases[-1 - i] = "A"
    fragment, size = "".join(bases), len(sequence)
    name = header.split()[0].removesuffix("/1")
    return ([f"{name}/1", (fragment + ADAPTERS[0] + "A" * size)[:size], "+", quality],
            [f"{name}/2", (fragment[::-1].translate(COMPLEMENT) + ADAPTERS[1] + "A" * size)[:size], "+", quality[::-1]])


def mixtures(libraries):
    """Each mixture takes its share of read pairs from the start of each library in it."""
    out = BUILD / "reads"
    by_id = {lib["sample_id"]: lib for lib in libraries}
    parts = {}
    for part in csv.DictReader(open(HERE / "mixtures.tsv"), delimiter="\t"):
        parts.setdefault(part["mixture"], []).append(part)
    rows = []
    for mixture, components in parts.items():
        rng = random.Random(mixture)
        ancient = components[0]["ancient"] == "yes"
        pairs = []
        for part in components:
            mates = zip(fastq(out / f"{part['library']}_R1.fastq.gz"), fastq(out / f"{part['library']}_R2.fastq.gz"))
            for _, (read1, read2) in zip(range(round(float(part["share"]) * READS_PER_LIBRARY)), mates):
                pairs.append(ancient_like(read1, rng) if ancient else (read1, read2))
        names = [f"{mixture}_R{i}.fastq.gz" for i in (1, 2)]
        for mate, name in enumerate(names):
            with gzip.GzipFile(out / name, "wb", mtime=0) as handle:
                handle.write("".join("\n".join(pair[mate]) + "\n" for pair in pairs).encode())
        removes_host = any(by_id[part["library"]]["host_action"] == "remove" for part in components)
        print(f"mixture {mixture:10s} " + ", ".join(f"{p['library']} {p['share']}" for p in components), flush=True)
        rows.append([mixture, *names, "ancient" if ancient else "modern", "remove" if removes_host else "skip",
                     "sink", "NA", "NA", ""])
    return rows


def samplesheet(rows):
    with open(BUILD / "samplesheet.tsv", "w") as fh:
        fh.write("sample_id\tfastq_1\tfastq_2\tage_group\thost_action\tsource_sink\tmodel_source_env\tproject\trun\n")
        fh.writelines("\t".join(r) + "\n" for r in rows)


def genomes(pinned):
    with open(BUILD / "library.fna", "w") as library:
        for g in pinned:
            url = (f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{g['accession']}/download"
                   "?include_annotation_type=GENOME_FASTA")
            archive = zipfile.ZipFile(io.BytesIO(fetch(url)))
            for name in archive.namelist():
                if not name.endswith((".fna", ".fa", ".fasta")):
                    continue
                for line in archive.read(name).decode().splitlines():
                    if line.startswith(">"):
                        library.write(f">{line[1:].split()[0]}|kraken:taxid|{g['species_taxid']}\n")
                    elif line:
                        library.write(line + "\n")
            print(f"genome  {g['accession']:18s} {g['species']}", flush=True)


def taxonomy(pinned):
    dump = tarfile.open(fileobj=io.BytesIO(fetch("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz")))
    nodes = dump.extractfile("nodes.dmp").read().decode().splitlines(keepends=True)
    parent = {line.split("\t|\t")[0]: line.split("\t|\t")[1] for line in nodes}
    keep = set()
    for g in pinned:
        taxid = g["species_taxid"]
        while taxid not in keep:
            keep.add(taxid)
            if taxid == "1":
                break
            taxid = parent[taxid]
    out = BUILD / "taxonomy"
    out.mkdir(parents=True, exist_ok=True)
    (out / "nodes.dmp").write_text("".join(sorted((l for l in nodes if l.split("\t|\t")[0] in keep),
                                                  key=lambda l: int(l.split("\t|\t")[0]))))
    names = dump.extractfile("names.dmp").read().decode().splitlines(keepends=True)
    (out / "names.dmp").write_text("".join(l for l in names if l.split("\t|\t")[0] in keep
                                           and l.split("\t|\t")[3].startswith("scientific name")))
    (out / "merged.dmp").write_text("")
    print(f"taxonomy {len(keep)} taxa", flush=True)


def host():
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text"
    (BUILD / "host.fa").write_bytes(fetch(url))
    print("host    NC_012920.1 (human mitochondrion)", flush=True)


def main():
    BUILD.mkdir(exist_ok=True)
    libraries = list(csv.DictReader(open(HERE / "libraries.tsv"), delimiter="\t"))
    pinned = list(csv.DictReader(open(HERE / "genomes.tsv"), delimiter="\t"))
    samplesheet(reads(libraries) + mixtures(libraries))
    host()
    taxonomy(pinned)
    genomes(pinned)


if __name__ == "__main__":
    main()
