#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

python3 bin/create_tiny_test_data.py
mkdir -p examples/test_data/kraken2_stub_db

nextflow run . \
    -profile test \
    -params-file examples/params.tiny.stub.yml \
    -stub-run \
    -resume
