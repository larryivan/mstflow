<h1 align="center">MSTFlow</h1>

<p align="center">
Metagenomic profiling for ancient and modern DNA that counts molecules, not reads,<br>
with optional microbial source attribution.
</p>

<p align="center">
  <a href="https://github.com/larryivan/mstflow/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/larryivan/mstflow/actions/workflows/ci.yml/badge.svg"></a>
  <a href="https://www.nextflow.io"><img alt="Nextflow ≥ 25.10.4" src="https://img.shields.io/badge/nextflow-%E2%89%A525.10.4-23aa62"></a>
  <a href="LICENSE"><img alt="MIT licence" src="https://img.shields.io/badge/licence-MIT-blue"></a>
</p>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/images/mstflow-workflow-dark.svg">
  <img alt="Ancient and modern libraries are trimmed on their own branches, then run side by side through host removal with BWA, Kraken2 and Bracken to a genus table, with SourceTracker2 as an optional last step." src="docs/images/mstflow-workflow-light.svg" width="100%">
</picture>

## Why MSTFlow

- **Molecules, not reads.** Overlapping ancient read pairs are collapsed into the molecules
  they came from, and every library is counted the same way, so ancient and modern
  samples can be compared directly.
- **Short fragments are not lost.** Mates go to Kraken2 as native pairs, never joined
  by an `N` spacer that would silently raise the confidence threshold for short ancient
  fragments. [Why this matters](docs/design-notes.md)
- **Reproducible by default.** Every tool runs from a container pinned by digest, and
  each run records tool versions, parameters and the molecules kept at every step.

## Quick start

With [Nextflow](https://www.nextflow.io/docs/latest/install.html) ≥ 25.10.4 and Docker:

```bash
nextflow run larryivan/mstflow -profile test,docker
```

This runs every step on a small bundled dataset in about a minute once the containers
are downloaded. Nothing else to install or download.

## Your own data

A samplesheet, saved as a tab-separated file with one row per library:

| sample_id | fastq_1 | fastq_2 | age_group |
|---|---|---|---|
| bone_01 | reads/bone_01_R1.fastq.gz | reads/bone_01_R2.fastq.gz | ancient |
| gut_01 | reads/gut_01.fastq.gz | | modern |

Then point MSTFlow at a Kraken2 database ([prebuilt indexes](https://benlangmead.github.io/aws-indexes/k2)
include what Bracken needs) and a BWA-indexed host genome:

```bash
nextflow run larryivan/mstflow -profile docker \
    --input samplesheet.tsv \
    --kraken2_db k2_standard/ \
    --host_fasta chm13v2.0.fa
```

`--help` lists the options. Singularity, Apptainer and Conda work in place of Docker.

## What you get

`results/summary/sample_summary.tsv` follows each sample through the workflow, here from
the test run:

| sample_id | age | input | after host removal | classified | genera | flags |
|---|---|--:|--:|--:|--:|---|
| ancient_pe | ancient | 3,000 | 2,093 | 69.8% | 3 | ok |
| ancient_se | ancient | 3,000 | 2,407 | 68.3% | 3 | ok |
| modern_pe | modern | 3,000 | 2,404 | 88.0% | 2 | few_genera |
| modern_se | modern | 3,000 | 2,706 | 78.7% | 3 | ok |

Genus abundance tables for all samples are in `results/tables/`, and every QC metric
in `results/qc/multiqc/multiqc_report.html`.

## Documentation

[Usage](docs/usage.md) · [Parameters](docs/parameters.md) · [Output](docs/output.md) ·
[Design notes](docs/design-notes.md) · [Troubleshooting](docs/troubleshooting.md)

## Citation

Please cite MSTFlow as given in [CITATION.cff](CITATION.cff), together with the tools it
runs: fastp, AdapterRemoval, BWA, Kraken2, Bracken, FastQC, MultiQC, and SourceTracker2 if
used.

MSTFlow is released under the [MIT licence](LICENSE).
