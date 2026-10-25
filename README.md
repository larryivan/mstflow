<h1 align="center">MSTFlow</h1>

<p align="center">
Metagenomic profiling of ancient and modern DNA in one analysis,<br>
and source attribution with five methods, scored on your own data.
</p>

<p align="center">
  <a href="https://github.com/larryivan/mstflow/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/larryivan/mstflow/actions/workflows/ci.yml/badge.svg"></a>
  <a href="https://www.nextflow.io"><img alt="Nextflow ≥ 25.10.4" src="https://img.shields.io/badge/nextflow-%E2%89%A525.10.4-23aa62"></a>
  <a href="LICENSE"><img alt="MIT licence" src="https://img.shields.io/badge/licence-MIT-blue"></a>
</p>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/images/mstflow-workflow-dark.svg">
  <img alt="Ancient and modern libraries are trimmed on their own lines (AdapterRemoval, fastp) and run side by side through host removal with BWA. The pair then splits into two branches, Kraken2 or KrakenUniq, and rejoins at the read, k-mer and taxon filters and one count table. An optional authentication spur of the ancient line rises from the table to Bowtie2, PyDamage and a verdict per taxon. An optional source attribution line fans out into five methods, SourceTracker2, FEAST, FastST, SourceID-NMF and SourcePredict, and gathers again for their comparison, cross-validated by study. A QC line runs underneath: FastQC before and after, the sample summary and MultiQC." src="docs/images/mstflow-workflow-light.svg" width="100%">
</picture>

## Why MSTFlow

- **Ancient and modern samples in one analysis.** Ancient libraries are trimmed with
  AdapterRemoval, which collapses overlapping read pairs, and modern ones with fastp;
  from then on every library is counted the same way, a read pair or a collapsed read as
  one fragment, against one database. Modern libraries cut to ancient fragment lengths
  and given ancient DNA damage are attributed to the same environment as the originals,
  oral and gut ones within five percentage points.
  [The test](docs/validation.md#ancient-and-modern-libraries-in-one-analysis)
- **Source attribution you can check.** SourceTracker2, FEAST, FastST, SourceID-NMF and
  SourcePredict run on the same table and are compared sink by sink. Cross-validation
  holds out one study's sources at a time and scores every method on them, so you can
  see which to trust on your panel, and the test data's sinks of known composition
  check their shares. [How it works](docs/attribution.md)
- **Authentication.** The leading taxa of each ancient library are aligned to reference
  genomes with Bowtie2 and tested for ancient DNA damage with PyDamage, so a modern
  contaminant can be told apart. [How it works](docs/authentication.md)
- **Reproducible by default.** Every tool runs from a container pinned by digest, and
  each run records tool versions, parameters and the fragments kept at every step.

## Quick start

With [Nextflow](https://www.nextflow.io/docs/latest/install.html) ≥ 25.10.4 and Docker:

```bash
nextflow run larryivan/mstflow -profile test,docker
```

`-profile` says how the tools are installed. If you are not sure which to use:

| What your machine has | Use |
|---|---|
| Docker, running | `-profile docker` |
| A cluster with Apptainer or Singularity | `-profile apptainer` or `-profile singularity` |
| None of these | Install [Miniforge](https://github.com/conda-forge/miniforge), then `-profile conda` |

MSTFlow checks this before it starts: if the one you chose is not installed, it lists what
it found on your machine and the `-profile` to use instead.

This runs every step, authentication and all five attribution methods with their
cross-validation included, on excerpts of 18 published ancient and modern libraries
and three sinks mixed from them in known shares.
The first run downloads about 1.9 GB of test data and the containers; the run itself
takes about 15 minutes on a laptop. [The test data](assets/test/README.md)

## Your own data

A samplesheet, saved as a tab-separated file with one row per library:

| sample_id | fastq_1 | fastq_2 | age_group |
|---|---|---|---|
| bone_01 | reads/bone_01_R1.fastq.gz | reads/bone_01_R2.fastq.gz | ancient |
| gut_01 | reads/gut_01.fastq.gz | | modern |

Then point MSTFlow at a Kraken2 database (a [prebuilt index](https://benlangmead.github.io/aws-indexes/k2)
works; give the `.tar.gz` or its URL as it is) and a BWA-indexed host genome:

```bash
nextflow run larryivan/mstflow -profile docker \
    --input samplesheet.tsv \
    --kraken2_db k2_standard/ \
    --host_fasta chm13v2.0.fa
```

`--help` lists the options. Singularity and Apptainer work in place of Docker, and
`-profile conda` runs every step without containers, FEAST and FastST included.

For source attribution, add `source_sink` and `model_source_env` columns and choose the
methods with `--attribution_methods`; see [source attribution](docs/attribution.md).
To authenticate the leading taxa, give their reference genomes with
`--authentication_genomes`; see [authentication](docs/authentication.md). KrakenUniq can
replace Kraken2 with `--classifier krakenuniq`.

## What you get

`results/summary/sample_summary.tsv` follows each sample through the workflow. Some rows
from the test run:

| sample_id | age | input | after trimming | classified | genera | flags |
|---|---|--:|--:|--:|--:|---|
| oral_1 | modern | 50,000 | 47,911 | 57.0% | 28 | ok |
| calculus | ancient | 50,000 | 35,438 | 16.9% | 19 | ok |
| sediment | ancient | 50,000 | 47,583 | 4.8% | 18 | ok |
| bone | ancient | 50,000 | 45,062 | 2.4% | 4 | few_genera |

`results/attribution/default/comparison.tsv` gives each method's call for every sink,
and `evaluation/evaluation.tsv` how often each was right on held-out sources:

| method | calls made | right when it called | right overall |
|---|--:|--:|--:|
| SourceTracker2 | 86% | 100% | 86% |
| FEAST | 100% | 100% | 100% |
| FastST | 100% | 100% | 100% |
| SourceID-NMF | 0% | – | 0% |
| SourcePredict | 0% | – | 0% |

Only SourceTracker2 also gets the shares right: on sinks mixed from known sources it
comes within 0.06 of their composition, while FEAST, FastST and SourceID-NMF miss by 0.5
or more.
[Validation](docs/validation.md)

`results/authentication/authentication.tsv` tests the leading taxa of each ancient
library for ancient DNA damage:

| sample_id | taxon | aligned reads | C→T at the first base | q-value | verdict |
|---|---|--:|--:|--:|---|
| calculus | Arachnia | 327 | 8.3% | 0.001 | damaged |
| calculus | Actinomyces | 250 | 2.1% | 0.90 | no damage |
| lab_1 | Leclercia | 1,218 | 0.6% | 0.73 | no damage |

Genus tables for all samples are in `results/tables/`, and every QC metric in
`results/qc/multiqc/multiqc_report.html`.

## Documentation

[Usage](docs/usage.md) · [Source attribution](docs/attribution.md) · [Authentication](docs/authentication.md) · [Validation](docs/validation.md) · [Parameters](docs/parameters.md) ·
[Output](docs/output.md) · [Design notes](docs/design-notes.md) · [Troubleshooting](docs/troubleshooting.md)

## Citation

Please cite MSTFlow as given in [CITATION.cff](CITATION.cff), together with the tools it
runs: fastp, AdapterRemoval, BWA, Kraken2 or KrakenUniq, FastQC, MultiQC,
Bowtie2 and PyDamage if you authenticated, and each attribution method you used
([references](docs/attribution.md#methods)).

MSTFlow is released under the [MIT licence](LICENSE). The attribution methods keep
their own licences; FEAST is for non-commercial use only.
