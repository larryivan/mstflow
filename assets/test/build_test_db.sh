#!/usr/bin/env bash
# Build the test Kraken2 database, host index and authentication references
# from build/ (made by make_test_data.py) with the workflow's pinned containers, then
# collect every file the test profile downloads into build/release/. Needs Docker.
set -euo pipefail
cd "$(dirname "$0")/build"

image() { sed -n "s/^ *$1 = '\(.*\)'$/\1/p" ../../../nextflow.config; }
run() { docker run --rm --platform linux/amd64 -u "$(id -u):$(id -g)" -v "$PWD:/w" -w /w "$@"; }
KRAKEN2=$(image container_kraken2); BWA=$(image container_bwa)
BOWTIE2=$(image container_bowtie2)

rm -rf kraken2_db && mkdir -p kraken2_db/taxonomy && cp taxonomy/*.dmp kraken2_db/taxonomy/
run "$KRAKEN2" kraken2-build --add-to-library library.fna --db kraken2_db --no-masking
run "$KRAKEN2" kraken2-build --build --db kraken2_db --threads 8
run "$BWA" bwa index host.fa

# Authentication references: for each genus among the three most abundant in an ancient
# library of a test run, the first genome of that genus in genomes.tsv.
printf '%s\n' Actinomyces Arachnia Bradyrhizobium Delftia Desulfobulbus Elizabethkingia Escherichia \
    Faecalibacterium Kribbella Leclercia Methanobrevibacter Micrococcus Olsenella Pseudomonas \
    Romboutsia Segatella Streptomyces > authentication_genera.txt
awk -F'\t' 'NR == FNR { want[$1] = 1; next } FNR > 1 && want[$1] && !seen[$1]++ { print $4 }' \
    authentication_genera.txt ../genomes.tsv > authentication_taxids.txt
awk 'NR == FNR { keep[$1] = 1; next } /^>/ { split($0, field, "|"); show = field[3] in keep } show' \
    authentication_taxids.txt library.fna | gzip -n > authentication.fna.gz
run "$BOWTIE2" bowtie2-build --threads 8 authentication.fna.gz authentication.fna.gz

rm -rf release && mkdir -p release
tar -czf release/kraken2_db.tar.gz -C kraken2_db hash.k2d opts.k2d taxo.k2d
cp samplesheet.tsv reads/*.fastq.gz host.fa host.fa.amb host.fa.ann host.fa.bwt host.fa.pac host.fa.sa release/
cp taxonomy/nodes.dmp taxonomy/names.dmp taxonomy/merged.dmp release/
cp authentication.fna.gz authentication.fna.gz.*.bt2 release/
ls -l release
