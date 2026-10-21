#!/usr/bin/env bash
# Build the test Kraken2/Bracken database and host index from build/ (made by
# make_test_data.py) with the workflow's pinned containers, then collect every file
# the test profile downloads into build/release/. Needs Docker.
set -euo pipefail
cd "$(dirname "$0")/build"

image() { sed -n "s/^ *$1 = '\(.*\)'$/\1/p" ../../../nextflow.config; }
run() { docker run --rm --platform linux/amd64 -u "$(id -u):$(id -g)" -v "$PWD:/w" -w /w "$@"; }
KRAKEN2=$(image container_kraken2); BRACKEN=$(image container_bracken); BWA=$(image container_bwa)

rm -rf kraken2_db && mkdir -p kraken2_db/taxonomy && cp taxonomy/*.dmp kraken2_db/taxonomy/
run "$KRAKEN2" kraken2-build --add-to-library library.fna --db kraken2_db --no-masking
run "$KRAKEN2" kraken2-build --build --db kraken2_db --threads 8
for len in 50 75 100 150; do       # read-length models for ancient and modern libraries
    run "$BRACKEN" bracken-build -d kraken2_db -t 8 -k 35 -l "$len"
done
run "$BWA" bwa index host.fa

rm -rf release && mkdir -p release
tar -czf release/kraken2_db.tar.gz -C kraken2_db hash.k2d opts.k2d taxo.k2d $(cd kraken2_db && ls database*mers.kmer_distrib)
cp samplesheet.tsv reads/*.fastq.gz host.fa host.fa.amb host.fa.ann host.fa.bwt host.fa.pac host.fa.sa release/
cp taxonomy/nodes.dmp taxonomy/names.dmp taxonomy/merged.dmp release/
ls -l release
