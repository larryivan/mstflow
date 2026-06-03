# mstflow

Nextflow workflow for metagenomic preprocessing, host removal, Kraken2/Bracken profiling, and optional SourceTracker2 attribution.

## Workflow

```text
metadata.tsv
-> FastQC
-> fastp or AdapterRemoval
-> host removal
-> Kraken2
-> Bracken
-> optional prokaryote filter
-> optional SourceTracker2
-> MultiQC
```

Preprocessing is selected per sample:

| sample type | tool |
| --- | --- |
| modern PE/SE | fastp |
| ancient PE | AdapterRemoval --collapse |
| ancient SE | AdapterRemoval |

## Run

```bash
nextflow run . -profile conda -params-file params.yml -resume
```

Tiny DAG test:

```bash
bash examples/run_tiny_stub.sh
```

## Metadata

Required tab-delimited columns:

| column | value |
| --- | --- |
| sample_id | unique sample ID |
| fastq_1 | R1 for PE, FASTQ for SE |
| fastq_2 | R2 for PE, `NA` for SE |
| layout | `PE` or `SE` |
| age_group | `modern` or `ancient` |

Recommended columns:

| column | value |
| --- | --- |
| preprocess_mode | usually `auto` |
| host_ref | FASTA path or `default` |
| source_sink | `source` or `sink` |
| model_source_env | SourceTracker2 source label |

## Params

```yaml
input: /path/to/metadata.tsv
outdir: results/project

remove_host: true
host_fasta: /path/to/host.fa

run_kraken2: true
kraken2_db: /path/to/kraken2_db

run_bracken: true
bracken_db: /path/to/kraken2_db
bracken_rank: G

taxon_filter: prokaryote
taxdump: /path/to/taxdump
drop_taxids: "0,1"

run_sourcetracker2: true
```

Use `taxon_filter: none` to skip lineage filtering.

## Outputs

```text
results/
├── kraken2/
├── bracken/
├── tables/
├── source_attribution/
├── qc/
└── pipeline_info/
```

Main tables:

```text
tables/bracken_G_counts.tsv
tables/bracken_G_relative_abundance.tsv
tables/bracken_G_prokaryote_counts.tsv
tables/bracken_G_prokaryote_relative_abundance.tsv
source_attribution/sourcetracker2/mixing_proportions.txt
```
