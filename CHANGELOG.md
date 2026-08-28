# Changelog

## 0.1.0 (2026-10-21)

First public release.

### Added
- `-profile test` runs every step on bundled synthetic data in about a minute, with no
  downloads.
- `--help`, with parameter types, defaults and accepted values, and validation of every
  parameter against `nextflow_schema.json`.
- `summary/sample_summary.tsv`: molecules after each step, classified fraction, taxa
  detected and depth flags for every sample, also shown in the MultiQC report.
  Thresholds: `--qc_min_classified_reads`, `--qc_min_taxa`.
- A summary printed at the end of each run: where the results are and which samples
  were flagged.
- Tasks killed for memory or time are retried once with twice the resources.
- MIT licence, `CITATION.cff`, CI, and documentation under `docs/`.

### Changed
- Samplesheet: only `sample_id`, `fastq_1` and `age_group` are required. `layout` is
  inferred from `fastq_2`, which may be left empty or omitted. Relative paths are read
  relative to the samplesheet (paths relative to the launch directory still work).
- Bracken runs by default, and `--bracken_db` defaults to `--kraken2_db`.
- Parameters are validated strictly against their declared types: booleans must be
  `true` or `false`, and `--bracken_read_len` takes `auto` or a positive integer.
- The `docker` profile runs the `linux/amd64` images on any host, including ARM.

### Removed
- The `slurm`, `slurm_conda`, `slurm_apptainer` and `slurm_singularity` profiles.
  Choose a scheduler in your own config instead, passed with `-c`; see docs/usage.md.
