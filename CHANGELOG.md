# Changelog

## Unreleased

- The samples to trace (`--input`) and the sources (`--sources`, with `environment` and `project` columns) are two sheets; the `source_sink` and `model_source_env` columns are gone.
- `--panel` reuses the sources of an earlier run without classifying them again. Every run with `--sources` writes its sources to `results/panel/`.
- Source attribution runs whenever there are sources, with SourceTracker2 unless `--attribution_methods` says otherwise. A run with sources only scores them and writes the panel.
- The sheets are read before any task starts, as TSV or CSV. Every problem is listed at once with its line, and a summary of the samples and of the sources per environment is printed.
- Changing a label reruns source attribution only: labels no longer reach trimming, host removal or classification.
- A misspelled option stops the run.

## 0.1.0 (2026-10-21)

First public release.

- Metagenomic profiling of ancient and modern libraries in one analysis, with Kraken2 or KrakenUniq (`--classifier krakenuniq`), counted in fragments, and an optional filter on distinct minimizers.
- Source attribution with SourceTracker2, FEAST, FastST, SourceID-NMF and SourcePredict side by side, and cross-validation by study on your own sources. See docs/attribution.md.
- Authentication of the leading taxa of ancient libraries with Bowtie2 and PyDamage. See docs/authentication.md.
- `-profile docker`, `singularity`, `apptainer` and `conda`, the last running every method, FEAST and FastST included. A check at the start says what is missing and which profile to use.
- A sample summary with depth flags, a MultiQC report, and a summary printed at the end of each run.
- `--kraken2_db` accepts a `.tar.gz` or URL, such as the prebuilt indexes; the samplesheet, reads, host genome and taxdump can be URLs too.
- `-profile test` runs every step on excerpts of 18 published libraries and three mixtures of known composition (about 2.1 GB, downloaded on first use). The recipe is in `assets/test`.
- `--help`, and validation of every parameter against `nextflow_schema.json`.
- MIT licence, `CITATION.cff`, CI, and documentation under `docs/`, including `docs/validation.md`: what the test data shows about shares, calls and ancient sinks.
