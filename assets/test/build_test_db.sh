#!/usr/bin/env bash
# Build the bundled test Kraken2/Bracken database and BWA host index with the
# workflow's own pinned containers. Run make_test_data.py first.
set -euo pipefail
cd "$(dirname "$0")"

image() {  # read a pinned image reference from nextflow.config
    sed -n "s/^ *$1 = '\(.*\)'$/\1/p" ../../nextflow.config
}
KRAKEN2=$(image container_kraken2)
BRACKEN=$(image container_bracken)
BWA=$(image container_bwa)
run() { docker run --rm --platform linux/amd64 -u "$(id -u):$(id -g)" -v "$PWD:/w" -w /w "$@"; }

rm -rf kraken2_db build
mkdir -p build/taxonomy kraken2_db
cp taxonomy/nodes.dmp taxonomy/names.dmp build/taxonomy/

run "$KRAKEN2" kraken2-build --add-to-library genomes/library.fna --db build --no-masking
run "$KRAKEN2" kraken2-build --build --db build --threads 2
for len in 50 100; do
    run "$BRACKEN" bracken-build -d build -t 2 -k 35 -l "$len"
done

mv build/hash.k2d build/opts.k2d build/taxo.k2d build/database*mers.kmer_distrib kraken2_db/
rm -rf build

run "$BWA" bwa index host/host.fa 2>/dev/null
echo "test database:"; ls -l kraken2_db host
