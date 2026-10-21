<h1 align="center">MSTFlow</h1>

<p align="center">
Metagenomic profiling for ancient and modern DNA that counts molecules, not reads,<br>
and source attribution with five methods, scored on your own data.
</p>

<p align="center">
  <a href="https://github.com/larryivan/mstflow/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/larryivan/mstflow/actions/workflows/ci.yml/badge.svg"></a>
  <a href="https://www.nextflow.io"><img alt="Nextflow ≥ 25.10.4" src="https://img.shields.io/badge/nextflow-%E2%89%A525.10.4-23aa62"></a>
  <a href="LICENSE"><img alt="MIT licence" src="https://img.shields.io/badge/licence-MIT-blue"></a>
</p>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/images/mstflow-workflow-dark.svg">
  <img alt="Ancient and modern libraries are trimmed on their own lines (AdapterRemoval, fastp) and run side by side through host removal with BWA. The pair then splits into two branches, Kraken2 or KrakenUniq, and rejoins at Bracken, then the k-mer and taxon filters and one count table. An optional source attribution line fans out into five methods, SourceTracker2, FEAST, FastST, SourceID-NMF and SourcePredict, and gathers again for their comparison, cross-validated by study. A QC line runs underneath: FastQC before and after, the sample summary and MultiQC." src="docs/images/mstflow-workflow-light.svg" width="100%">
</picture>

## Why MSTFlow

- **Molecules, not reads.** Overlapping ancient read pairs are collapsed into the molecules
  they came from, and every library is counted the same way, so ancient and modern
  samples can be compared directly.
- **Short fragments are not lost.** Mates go to Kraken2 as native pairs, never joined
  by an `N` spacer that would silently raise the confidence threshold for short ancient
  fragments. [Why this matters](docs/design-notes.md)
- **Source attribution you can check.** SourceTracker2, FEAST, FastST, SourceID-NMF and
  SourcePredict run on the same table and are compared sink by sink. Cross-validation
  holds out one study's sources at a time and scores every method on them, so you can
  see which to trust on your panel. [How it works](docs/attribution.md)
- **Reproducible by default.** Every tool runs from a container pinned by digest, and
  each run records tool versions, parameters and the molecules kept at every step.

## Quick start

With [Nextflow](https://www.nextflow.io/docs/latest/install.html) ≥ 25.10.4 and Docker:

```bash
nextflow run larryivan/mstflow -profile test,docker
```

This runs every step, all five attribution methods and their cross-validation included,
on excerpts of 15 published ancient and modern libraries. The first run downloads about
1.8 GB of test data and the containers; the run itself takes about 15 minutes on a
laptop. [The test data](assets/test/README.md)

## Your own data

A samplesheet, saved as a tab-separated file with one row per library:

| sample_id | fastq_1 | fastq_2 | age_group |
|---|---|---|---|
| bone_01 | reads/bone_01_R1.fastq.gz | reads/bone_01_R2.fastq.gz | ancient |
| gut_01 | reads/gut_01.fastq.gz | | modern |

Then point MSTFlow at a Kraken2 database ([prebuilt indexes](https://benlangmead.github.io/aws-indexes/k2)
include what Bracken needs; give the `.tar.gz` or its URL as it is) and a BWA-indexed
host genome:

```bash
nextflow run larryivan/mstflow -profile docker \
    --input samplesheet.tsv \
    --kraken2_db k2_standard/ \
    --host_fasta chm13v2.0.fa
```

`--help` lists the options. Singularity and Apptainer work in place of Docker; Conda
covers everything except FEAST, FastST and SourceID-NMF.

For source attribution, add `source_sink` and `model_source_env` columns and choose the
methods with `--attribution_methods`; see [source attribution](docs/attribution.md).
KrakenUniq can replace Kraken2 with `--classifier krakenuniq`.

## What you get

`results/summary/sample_summary.tsv` follows each sample through the workflow. Some rows
from the test run:

| sample_id | age | input | after trimming | classified | genera | flags |
|---|---|--:|--:|--:|--:|---|
| oral_1 | modern | 50,000 | 42,915 | 59.2% | 27 | ok |
| calculus | ancient | 50,000 | 35,438 | 16.9% | 19 | ok |
| sediment | ancient | 50,000 | 47,583 | 4.8% | 18 | ok |
| bone | ancient | 50,000 | 45,062 | 2.4% | 4 | few_genera |

`results/attribution/default/comparison.tsv` gives each method's call for every sink,
and `evaluation/evaluation.tsv` how often each was right on held-out sources:

| method | calls made | right when it called | right overall |
|---|--:|--:|--:|
| SourceTracker2 | 86% | 100% | 86% |
| FEAST | 57% | 100% | 57% |
| FastST | 100% | 29% | 29% |
| SourceID-NMF | 0% | – | 0% |
| SourcePredict | 14% | 100% | 14% |

Genus tables for all samples are in `results/tables/`, and every QC metric in
`results/qc/multiqc/multiqc_report.html`.

## Documentation

[Usage](docs/usage.md) · [Source attribution](docs/attribution.md) · [Parameters](docs/parameters.md) ·
[Output](docs/output.md) · [Design notes](docs/design-notes.md) · [Troubleshooting](docs/troubleshooting.md)

## Citation

Please cite MSTFlow as given in [CITATION.cff](CITATION.cff), together with the tools it
runs: fastp, AdapterRemoval, BWA, Kraken2 and Bracken or KrakenUniq, FastQC, MultiQC, and
each attribution method you used ([references](docs/attribution.md#methods)).

MSTFlow is released under the [MIT licence](LICENSE). The attribution methods keep
their own licences; FEAST is for non-commercial use only.
