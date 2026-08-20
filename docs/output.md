# Output

Everything is written under `--outdir` (default `results/`). Start with the first three.

| Path | What it is |
|---|---|
| `summary/sample_summary.tsv` | One row per sample: molecules after trimming and host removal, host fraction, classified fraction, taxa detected, and flags |
| `tables/` | Abundance tables across all samples |
| `qc/multiqc/multiqc_report.html` | Every QC metric in one report, including the sample summary |

## `summary/sample_summary.tsv`

| Column | Meaning |
|---|---|
| `input` | Molecules in the raw files (read pairs for paired-end libraries) |
| `after_trimming` | Molecules left after adapter and quality trimming; for ancient paired-end libraries, collapsed molecules |
| `after_host_removal` | Molecules left after removing host alignments |
| `host_fraction` | Share of trimmed molecules that aligned to the host |
| `classified`, `classified_fraction` | Molecules Kraken2 assigned to a taxon |
| `genera` | Taxa detected at `--bracken_rank` (named `species`, `families` or `taxa` for other ranks), after the taxon filter if one is set |
| `flags` | `ok`, or why the profile is too shallow to interpret: `few_reads` (fewer than `--qc_min_classified_reads`) or `few_genera` (fewer than `--qc_min_taxa`) |

A flagged sample is not wrong, but its profile carries little information, and
attributing it to a source is unreliable.

## `tables/`

`G` is the rank (`--bracken_rank`) and `prokaryote` the taxon filter (`--taxon_filter`).

| File | Contents |
|---|---|
| `bracken_G_counts.tsv` | Estimated reads per genus (rows) and sample (columns), all taxa |
| `bracken_G_relative_abundance.tsv` | The same as proportions of each sample |
| `bracken_G_features.tsv` | Taxon ID, name and rank of every row |
| `bracken_G_prokaryote_*.tsv` | The same three tables restricted to bacteria and archaea |
| `bracken_G_prokaryote_filter_qc.tsv`, `*_sample_qc.tsv` | What the filter removed, and per-sample Kraken2 and Bracken totals |

## Per-step files

| Directory | Contents |
|---|---|
| `preprocess/fastp/`, `preprocess/adapterremoval/` | Trimming reports and molecule counts |
| `host_removed/` | Host-removal logs and molecule counts |
| `kraken2/` | Kraken2 reports (with and without minimizer data), read-length and confidence statistics |
| `bracken/` | Per-sample Bracken estimates, the read-length model chosen, and its status |
| `qc/fastqc_raw/`, `qc/fastqc_clean/` | FastQC before trimming and after host removal |
| `source_attribution/sourcetracker2/<label>/` | Only with `--run_sourcetracker2`: mixing proportions and their uncertainty, assignments (`sourcetracker2_assignments.tsv`), a stacked bar plot, and which samples were excluded and why |
| `pipeline_info/` | Validated samplesheet, Nextflow report, timeline, trace and DAG |

Every directory also holds a `*.versions.yml` per step, recording the exact tool versions used.

Large intermediate FASTQ files stay in the Nextflow work directory unless
`--publish_intermediates true` is set. Delete the work directory once you no longer
need `-resume`.
