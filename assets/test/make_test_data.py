#!/usr/bin/env python3
"""Download everything the test dataset is built from, into build/.

- reads: the first 50,000 reads (pairs) of each library in libraries.tsv, from ENA
- genomes: the reference genomes pinned in genomes.tsv, from NCBI Datasets
- taxonomy: the NCBI lineages of those genomes, cut from the current taxdump (nodes.dmp,
  names.dmp, and an empty merged.dmp, which SourcePredict's taxonomy reader expects)
- host: the human mitochondrial genome (NC_012920.1)
- samplesheet.tsv, pointing at the reads by file name

Then run build_test_db.sh to build the Kraken2/Bracken database and the release files.
Only the Python standard library is used.
"""
import csv
import gzip
import io
import sys
import tarfile
import time
import urllib.request
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
BUILD = HERE / "build"
READS_PER_LIBRARY = 50_000


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
    reads(libraries)
    host()
    taxonomy(pinned)
    genomes(pinned)


if __name__ == "__main__":
    main()
