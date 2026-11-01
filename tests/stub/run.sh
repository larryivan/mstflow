#!/usr/bin/env bash
# Checks the workflow's wiring and output contracts without running any tool
# (-stub-run). For a real end-to-end run use: nextflow run . -profile test,docker
set -euo pipefail

cd "$(dirname "$0")/../.."

python3 tests/stub/make_stub_data.py

nextflow run . \
    -profile test \
    -params-file tests/stub/params.yml \
    -stub-run \
    -resume

# The same wiring with KrakenUniq in place of Kraken2.
nextflow run . \
    -profile test \
    -params-file tests/stub/params.yml \
    --classifier krakenuniq \
    --krakenuniq_db tests/stub/data/krakenuniq_stub_db \
    --krakenuniq_min_kmers 0 \
    --outdir results/stub_krakenuniq \
    -stub-run \
    -resume

# The samples again, against the panel the first run wrote instead of the sources' reads,
# and the sources alone, which are scored and written as a panel. A params file can unset
# --sources and --input, which the command line cannot.
mkdir -p results
sed -e 's|^sources:.*|sources: null|' -e 's|^outdir:.*|outdir: results/stub_panel|' \
    tests/stub/params.yml > results/stub_panel.params.yml
nextflow run . \
    -profile test \
    -params-file results/stub_panel.params.yml \
    --panel results/stub/panel \
    -stub-run \
    -resume
test -s results/stub_panel/attribution/tiny_stub/comparison.tsv
test ! -e results/stub_panel/panel

sed -e 's|^input:.*|input: null|' -e 's|^outdir:.*|outdir: results/stub_sources|' \
    tests/stub/params.yml > results/stub_sources.params.yml
nextflow run . \
    -profile test \
    -params-file results/stub_sources.params.yml \
    -stub-run \
    -resume
test -s results/stub_sources/panel/sources.tsv
test -s results/stub_sources/attribution/tiny_stub/evaluation/evaluation.tsv
test ! -e results/stub_sources/attribution/tiny_stub/comparison.tsv
