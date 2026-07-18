# mstflow

MSTFlow is a Nextflow DSL2 workflow for fragment-aware metagenomic preprocessing, host removal, Kraken2/Bracken profiling, and optional SourceTracker2 source attribution. It supports modern and ancient, paired- and single-end libraries while keeping the biological counting unit explicit.

## Workflow

```text
metadata.tsv
  -> metadata and reference preflight
  -> raw FastQC
  -> fastp (modern) or AdapterRemoval (ancient)
  -> pair-aware host removal
  -> cleaned FastQC
  -> Kraken2
  -> per-sample Bracken model selection
  -> count/relative-abundance tables
  -> optional SourceTracker2
  -> MultiQC and provenance reports
```

Ancient paired-end libraries default to `ancient_pe_mode: collapsed_only`, matching the molecule-level protocol in `steps.md`. The excluded singleton and unmerged counts are reported. The optional `all_fragments` mode retains all complete molecular fragments (collapsed reads and complete unmerged pairs), preserves pairing through host removal, and converts each surviving pair into one classification record separated by a non-k-mer-forming `N` gap. Orphan singletons are quantified but excluded, matching modern paired-end handling.

## Requirements

- Nextflow 25.10.4 or newer
- Bash
- Mamba/Conda for the pinned environments
- a BWA-indexed host FASTA when `remove_host: true`
- a Kraken2 database when `run_kraken2: true`
- matching Bracken distributions when `run_bracken: true`

The primary tools are version-pinned in [`envs/mstflow.yml`](envs/mstflow.yml). SourceTracker uses its separate compatibility environment in [`envs/sourcetracker2.yml`](envs/sourcetracker2.yml).

## Run

Run locally with the pinned Conda environments:

```bash
nextflow run . -profile conda -params-file /absolute/path/to/params.yml -resume
```

The portable `slurm` profile submits processes through Slurm without assuming an account, partition, module system, work path, or storage layout:

```bash
nextflow run . -profile slurm -params-file /absolute/path/to/params.yml -resume
```

Site policies belong in an external configuration supplied with Nextflow's global `-c` option; they are intentionally not hard-coded in this repository. See [`docs/server_setup.md`](docs/server_setup.md) for portable server setup, reference validation, resource tuning, and provenance guidance.

## Tests

Run the self-contained stub DAG test without installing scientific tools:

```bash
bash examples/run_tiny_stub.sh
```

Stub mode validates workflow wiring and output contracts; it does not validate real FastQC, BWA, Kraken2, Bracken, or SourceTracker2 behavior.

## Metadata

The input must resolve to exactly one tab-delimited file. Required columns are:

| column | accepted value |
| --- | --- |
| `sample_id` | unique, filesystem-safe sample identifier |
| `fastq_1` | absolute R1/SE FASTQ path |
| `fastq_2` | absolute R2 path for PE; `NA` for SE |
| `layout` | `PE` or `SE` |
| `age_group` | `modern` or `ancient` |

Optional analysis columns include:

| column | purpose |
| --- | --- |
| `preprocess_mode` | normally `auto` |
| `host_action` | `remove` or `skip`; omitted rows retain the historical `remove` behavior |
| `host_ref` | sample-specific indexed FASTA, or `default` |
| `host_group` | biological host group; non-human groups require explicit `host_ref` by default |
| `source_sink` | exactly `source` or `sink` for SourceTracker2 |
| `model_source_env` | source environment label; sinks use `NA` |
| `unit_id` | biological unit used to detect duplicate/run relationships |
| `panel_role`, `panel_layer`, `category` | preserved reference-panel design fields |

SourceTracker runs require complete, explicit source/sink assignments. Controls, holdouts, NO_BONE experiments, and production reference panels should be separate manifests rather than silently mixing roles.

## Core parameters

```yaml
input: /absolute/path/to/metadata.tsv
outdir: /absolute/path/to/results

remove_host: true
host_fasta: /absolute/path/to/chm13.fa
host_remove_max_forks: 4
host_min_mapq: 0
default_host_group: human
require_explicit_host_for_nondefault: true
require_explicit_host_action: false

run_kraken2: true
kraken2_db: /absolute/path/to/k2_pluspf_20260226
kraken2_kmer_len: 35
kraken2_min_read_len: 35
kraken2_confidence: 0.05
kraken2_minimum_hit_groups: 2
kraken2_use_names: true
kraken2_report_minimizer_data: true
kraken2_save_read_classifications: false
kraken2_max_forks: 1
kraken2_cpus: 16
kraken2_memory: 160 GB
kraken2_time: 36h

run_bracken: true
bracken_db: /absolute/path/to/k2_pluspf_20260226
bracken_rank: G
bracken_read_len: auto

taxon_filter: prokaryote
taxdump: /absolute/path/to/taxdump

ancient_pe_mode: collapsed_only
publish_intermediates: false

run_sourcetracker2: false
sourcetracker2_aggregate_by: unit_id
sourcetracker2_taxon_scope: filtered
sourcetracker2_analysis_label: production
sourcetracker2_exclude_panel_roles: "control,holdout"
sourcetracker2_min_kraken_classified_fraction: 0.001
sourcetracker2_min_taxon_filter_retained_reads: 1000
sourcetracker2_min_taxon_filter_retained_fraction: 0.001
sourcetracker2_max_unknown_prop: 0.5
sourcetracker2_min_top_known_prop: 0.5
sourcetracker2_min_top_known_margin: 0.05
sourcetracker2_loo: false
sourcetracker2_seed: 1
```

`bracken_read_len: auto` selects the nearest available `database<N>mers.kmer_distrib` per sample from the median observed single-read length and records the choice. For paired-end data, each mate contributes one read-length observation; the two mate lengths are never added together. A positive integer forces a model and fails if that distribution is absent.

If a valid Kraken report yields no estimate at the requested Bracken rank and `bracken_threshold`, MSTFlow retains the sample column with zero estimated counts, emits a zero-root Bracken report, and records `zero_estimates/no_estimates_at_requested_rank_threshold` in its per-sample status file. This means only that Bracken produced no estimate under the declared settings; it does not establish biological absence. The original Kraken report remains unchanged in the Kraken outputs. MSTFlow does not silently lower the threshold or invent an abundance estimate, and every other Bracken failure remains fatal.

`kraken2_min_read_len` must be at least `kraken2_kmer_len`. The default Kraken confidence and minimum-hit-group settings are explicit so sources and sinks cannot silently use different classifier defaults. Per-read classifications are disabled by default; reports and downstream Bracken inputs are retained.

`host_min_mapq: 0` removes an entire fragment when either mate has any primary host alignment. A threshold such as 20 is suitable only as a separately declared sensitivity analysis. Minimizer data and taxon names are enabled in the Kraken QC report; MSTFlow also emits a Bracken-compatible six-column report.

`require_explicit_host_for_nondefault: true` prevents a non-human `host_group` from silently using the default human reference. Each such sample must provide its own indexed `host_ref`.

`host_action: remove` sends the sample through BWA host subtraction. `host_action: skip` bypasses BWA and passes the preprocessed reads directly to clean-read QC and classification; skip rows must use `host_ref: NA`. Manifests that omit `host_action` remain backward compatible and default to `remove`. For a new mixed host-aware run, set `require_explicit_host_action: true` so every row must declare its policy.

Kraken memory depends on the selected database. The CPU, memory, time, and concurrency defaults are explicit and should be overridden after measuring the intended database and executor. Host removal is capped at four concurrent tasks by default to control repeated index I/O.

SourceTracker defaults to biological-unit aggregation, filtered taxa, and explicit QC. Complete units are validated and filtered as a whole, so a conflicting or poor lane cannot leave a partial unit. `control` and `holdout` panel roles are excluded before modeling; clear or replace `sourcetracker2_exclude_panel_roles` only in a separately labelled validation run when one of those roles is intentionally needed. Assignments abstain when Unknown is at least 0.5, the strongest known source is below 0.5, the top sources tie, or their margin is below 0.05. SourceTracker2 2.0.1 is run in one process with a recorded `sourcetracker2_seed` for reproducibility. Set `sourcetracker2_loo: true` in a separately labelled, source-only validation run, set `sourcetracker2_source_rarefaction_depth: 0`, and retain at least two independent analysis units per source environment. In LOO, MSTFlow deterministically rarefies each unit to `sourcetracker2_sink_rarefaction_depth` before SourceTracker, then disables both buggy upstream rarefaction paths.

`publish_intermediates: false` keeps large FASTQ and per-read intermediates in the resumable Nextflow work directory instead of duplicating them into the results tree.

## Outputs

```text
results/
├── bracken/
├── kraken2/
├── qc/
├── source_attribution/
├── tables/
└── pipeline_info/
```

Important deliverables include Bracken counts and relative-abundance tables, AdapterRemoval/host/Kraken fragment QC, Kraken minimizer and read-length summaries, per-sample Bracken read-length selections, SourceTracker mixing proportions when enabled, MultiQC, and Nextflow trace/report/timeline files.

Run outputs, work directories, local params, downloaded data, databases, and scheduler logs are excluded by `.gitignore`. Commit source/configuration changes only; archive run provenance with the scientific analysis rather than in the code repository.
