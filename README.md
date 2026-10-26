<h1 align="center">MSTFlow</h1>

<p align="center">
Metagenomic profiling of ancient and modern DNA in one analysis,<br>
and source attribution with five methods, plus a test of which one works on your samples.
</p>

<p align="center">
  <a href="https://github.com/larryivan/mstflow/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/larryivan/mstflow/actions/workflows/ci.yml/badge.svg"></a>
  <a href="https://www.nextflow.io"><img alt="Nextflow ≥ 25.10.4" src="https://img.shields.io/badge/nextflow-%E2%89%A525.10.4-23aa62"></a>
  <a href="LICENSE"><img alt="MIT licence" src="https://img.shields.io/badge/licence-MIT-blue"></a>
</p>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="docs/images/mstflow-workflow-dark.svg">
  <img alt="The MSTFlow workflow. Ancient and modern libraries are trimmed on separate lines, host reads are removed, and Kraken2 or KrakenUniq classifies the rest into one count table. Optional authentication (Bowtie2, PyDamage) and source attribution (five methods, cross-validated) read that table, and QC runs underneath." src="docs/images/mstflow-workflow-light.svg" width="100%">
</picture>

## What it does

- **Ancient and modern together.** Each library is trimmed for its kind, then counted and classified the same way against one database. [Tested](docs/validation.md#ancient-and-modern-libraries-in-one-analysis)
- **Source attribution you can check.** SourceTracker2, FEAST, FastST, SourceID-NMF and SourcePredict run on one table and are compared sample by sample. Cross-validation on your own sources shows which to trust. [How it works](docs/attribution.md)
- **Authentication.** The leading taxa of each ancient library are tested for ancient DNA damage, so a modern contaminant stands out. [How it works](docs/authentication.md)
- **Reproducible.** Tools run from pinned containers or conda environments, and every run records versions and parameters.

A *sink* is a sample to explain, such as an ancient dental calculus; *sources* are reference samples of known origin, such as oral, gut or soil.

## Quick start

With [Nextflow](https://www.nextflow.io/docs/latest/install.html) ≥ 25.10.4:

```bash
nextflow run larryivan/mstflow -profile test,docker
```

This runs every step on short excerpts of 18 published libraries: about 2.1 GB to download, 15 minutes on a laptop. [The test data](assets/test/README.md)

`-profile` says how the tools are installed:

| Your machine has | Use |
|---|---|
| Docker, running | `-profile docker` |
| Apptainer or Singularity (clusters) | `-profile apptainer` or `-profile singularity` |
| Neither | [Miniforge](https://github.com/conda-forge/miniforge), then `-profile conda` |

If your choice is missing, MSTFlow says so at once and lists what it found.

## Your own data

A tab-separated samplesheet, one row per library:

| sample_id | fastq_1 | fastq_2 | age_group |
|---|---|---|---|
| bone_01 | reads/bone_01_R1.fastq.gz | reads/bone_01_R2.fastq.gz | ancient |
| gut_01 | reads/gut_01.fastq.gz | | modern |

```bash
nextflow run larryivan/mstflow -profile docker \
    --input samplesheet.tsv \
    --kraken2_db k2_standard/ \
    --host_fasta chm13v2.0.fa
```

`--kraken2_db` takes a [prebuilt index](https://benlangmead.github.io/aws-indexes/k2) as downloaded (`.tar.gz` or URL), and `--host_fasta` a BWA-indexed genome. That gives genus tables. `--help` lists every option.

**To attribute sinks and score the methods,** add three columns that label the libraries of known origin as *sources*, name their environment and say which study they come from:

| sample_id | source_sink | model_source_env | project |
|---|---|---|---|
| calculus_01 | sink | | |
| oral_01 | source | oral | study_A |
| gut_01 | source | gut | study_A |
| oral_02 | source | oral | study_B |
| gut_02 | source | gut | study_B |

and switch it on with `--attribution_methods sourcetracker2,feast --attribution_cv true`. MSTFlow estimates where each sink's microbes came from, then holds out one study's sources at a time, treats them as sinks of known origin and reports how often each method names the right environment. [Details](docs/attribution.md)

To test ancient taxa for damage, give reference genomes with `--authentication_genomes`; see [usage](docs/usage.md).

## What you get

`results/summary/sample_summary.tsv` follows every sample through the workflow:

| sample_id | age | input | after trimming | classified | genera | flags |
|---|---|--:|--:|--:|--:|---|
| oral_1 | modern | 50,000 | 47,911 | 57.0% | 28 | ok |
| calculus | ancient | 50,000 | 35,438 | 16.9% | 19 | ok |
| sediment | ancient | 50,000 | 47,583 | 4.8% | 18 | ok |
| bone | ancient | 50,000 | 45,062 | 2.4% | 4 | few_genera |

Next to it: genus tables, a MultiQC report, and, when switched on, `attribution/` (each method's call per sink and its accuracy on held-out sources) and `authentication/` (a damage verdict per taxon). [Output files](docs/output.md)

## Documentation

[Usage](docs/usage.md) · [Source attribution](docs/attribution.md) · [Authentication](docs/authentication.md) · [Validation](docs/validation.md) · [Output](docs/output.md) · [Parameters](docs/parameters.md) · [Design notes](docs/design-notes.md) · [Troubleshooting](docs/troubleshooting.md)

## Citation and licence

Cite MSTFlow as in [CITATION.cff](CITATION.cff), with the tools you used: fastp, AdapterRemoval, BWA, Kraken2 or KrakenUniq, FastQC, MultiQC, Bowtie2 and PyDamage (authentication), and each attribution method ([references](docs/attribution.md#methods)).

MSTFlow is [MIT-licensed](LICENSE). The attribution methods keep their own licences; FEAST is for non-commercial use only.
