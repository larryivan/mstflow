#!/usr/bin/env bash
# Checks the workflow's wiring and output contracts without running any tool
# (-stub-run). For a real end-to-end run use: nextflow run . -profile test,docker
set -euo pipefail

cd "$(dirname "$0")/../.."

python3 tests/stub/make_stub_data.py
mkdir -p tests/stub/data/kraken2_stub_db

nextflow run . \
    -profile test \
    -params-file tests/stub/params.yml \
    -stub-run \
    -resume

# The same wiring with KrakenUniq in place of Kraken2 and Bracken.
nextflow run . \
    -profile test \
    -params-file tests/stub/params.yml \
    --classifier krakenuniq \
    --krakenuniq_db tests/stub/data/krakenuniq_stub_db \
    --krakenuniq_min_kmers 0 \
    --krakenuniq_min_reads 0 \
    --outdir results/stub_krakenuniq \
    -stub-run \
    -resume
