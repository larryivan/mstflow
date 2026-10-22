# Parameters

Every option, generated from `nextflow_schema.json`. `nextflow run larryivan/mstflow --help`
prints the same list; `--help <name>` explains one option. Options marked *advanced* are
hidden from `--help` unless `--showHidden` is given.

## Input and output

| Option | Description | Type | Default |
|---|---|---|---|
| `--input` | Samplesheet (TSV) listing every library; a local path or a URL. Required columns: `sample_id`, `fastq_1`, `age_group` (`ancient` or `modern`). Add `fastq_2` for paired-end libraries. Relative paths are read relative to the samplesheet. See docs/usage.md. | string | required |
| `--outdir` | Directory for results. | string | `results` |
| `--publish_intermediates` | Also copy large intermediate FASTQ files to the results directory. | boolean | `false` |
| `--multiqc_title` | Title shown on the MultiQC report. | string | `mstflow` |

## Host removal

| Option | Description | Type | Default |
|---|---|---|---|
| `--remove_host` | Remove reads that align to the host genome. | boolean | `true` |
| `--host_fasta` | Host FASTA with its BWA index (.amb .ann .bwt .pac .sa) next to it; a local path or a URL. Used for every sample whose `host_ref` is empty or `default`. Build the index once with `bwa index host.fa`. | string |  |
| `--host_min_mapq` | Minimum mapping quality for a host alignment to remove a fragment. 0 removes a fragment when either mate has any primary host alignment. | integer | `0` |
| `--default_host_group` | Host group that the default host FASTA represents. | string | `human` |
| `--require_explicit_host_for_nondefault` | Require an explicit `host_ref` for samples whose `host_group` differs from the default. | boolean | `true` |
| `--require_explicit_host_action` | Require every samplesheet row to declare `host_action` (remove or skip). | boolean | `false` |

<details><summary>Advanced</summary>

| Option | Description | Type | Default |
|---|---|---|---|
| `--host_remove_max_forks` | Maximum concurrent host-alignment tasks. | integer | `4` |

</details>

## Read preprocessing

| Option | Description | Type | Default |
|---|---|---|---|
| `--min_length_modern` | Minimum read length after trimming, modern libraries. | integer | `50` |
| `--min_length_ancient` | Minimum fragment length after trimming, ancient libraries. | integer | `35` |
| `--ancient_pe_mode` | How ancient paired-end libraries are counted. `collapsed_only` keeps merged molecules only. `all_fragments` also keeps complete unmerged pairs, each submitted to Kraken2 as one paired record. See docs/design-notes.md. | `collapsed_only` / `all_fragments` | `collapsed_only` |

<details><summary>Advanced</summary>

| Option | Description | Type | Default |
|---|---|---|---|
| `--fastp_extra_args` | Extra arguments passed to fastp. | string |  |
| `--adapterremoval_extra_args` | Extra arguments passed to AdapterRemoval. | string |  |
| `--fastq_gzip_level` | Gzip level (0-9) for FASTQ written by the workflow's own helpers. | integer | `6` |

</details>

## Classifier

| Option | Description | Type | Default |
|---|---|---|---|
| `--classifier` | Which classifier to use. `kraken2` classifies with Kraken2. `krakenuniq` uses KrakenUniq instead, as aMeta does, and removes taxa whose reads cover too few unique k-mers. Bracken then re-estimates abundances after either. Use one classifier and database for sources and sinks alike. | `kraken2` / `krakenuniq` | `kraken2` |
| `--krakenuniq_db` | KrakenUniq database: a directory with database.kdb, database.idx and taxDB (and Bracken's database<N>mers.kmer_distrib files, unless --run_bracken false), or a .tar.gz of one; local or a URL. | string |  |
| `--krakenuniq_preload_size` | Load the KrakenUniq database in chunks of this size, e.g. 16G; unset loads it whole. | string |  |
| `--krakenuniq_min_kmers` | Remove a taxon from a sample when its reads cover fewer unique k-mers. | integer | `1000` |
| `--krakenuniq_min_reads` | Remove a taxon from a sample when it has fewer reads. | integer | `200` |

## Kraken2 classification

| Option | Description | Type | Default |
|---|---|---|---|
| `--kraken2_db` | Kraken2 database: a directory (hash.k2d, opts.k2d, taxo.k2d) or a .tar.gz of one, local or a URL. Required with the default classifier. Prebuilt indexes from https://benlangmead.github.io/aws-indexes/k2 can be given as the downloaded .tar.gz or its URL; the workflow unpacks it. Bracken needs the `database<N>mers.kmer_distrib` files, which prebuilt indexes include. | string |  |
| `--kraken2_confidence` | Kraken2 confidence threshold. | number | `0.05` |
| `--kraken2_minimum_hit_groups` | Kraken2 minimum hit groups. | integer | `2` |
| `--kraken2_kmer_len` | k-mer length the database was built with. | integer | `35` |
| `--kraken2_min_read_len` | Shortest read submitted to Kraken2; must be at least `kraken2_kmer_len`. | integer | `35` |
| `--kraken2_min_distinct_minimizers` | Remove a taxon from a sample when its reads cover fewer distinct minimizers (0: keep all). False positives often rest on many reads from a few short stretches of a genome. Distinct minimizers measure how much of the genome the reads cover, as unique k-mers do in KrakenUniq. Removals are listed in `tables/bracken_<rank>_kmer_filter.tsv`. See docs/usage.md. | integer | `0` |
| `--kraken2_save_read_classifications` | Keep per-read Kraken2 classifications (large). | boolean | `false` |
| `--kraken2_cpus` | CPUs per Kraken2 task. | integer | `16` |
| `--kraken2_memory` | Memory per Kraken2 task; must hold the whole database. | string | `160 GB` |
| `--kraken2_time` | Time limit per Kraken2 task. | string | `36h` |
| `--kraken2_max_forks` | Maximum concurrent Kraken2 tasks. | integer | `1` |

<details><summary>Advanced</summary>

| Option | Description | Type | Default |
|---|---|---|---|
| `--run_kraken2` | Classify reads (with Kraken2 or KrakenUniq); false stops after preprocessing. | boolean | `true` |
| `--kraken2_use_names` | Include taxon names in Kraken2 output. | boolean | `true` |
| `--kraken2_report_minimizer_data` | Include minimizer counts in the Kraken2 report. | boolean | `true` |
| `--kraken2_extra_args` | Extra arguments passed to Kraken2. | string |  |

</details>

## Bracken abundance

| Option | Description | Type | Default |
|---|---|---|---|
| `--run_bracken` | Re-estimate abundances at --bracken_rank with Bracken; without it the tables hold the classifier's own clade counts. Bracken moves reads the classifier left at higher ranks down to --bracken_rank. It needs the database's database<N>mers.kmer_distrib files, which prebuilt Kraken2 indexes ship; for KrakenUniq, build them with `bracken-build -y krakenuniq -k 31`. | boolean | `true` |
| `--bracken_rank` | Taxonomic rank of the output tables. | `D` / `P` / `C` / `O` / `F` / `G` / `S` | `G` |
| `--bracken_read_len` | Bracken read-length model: `auto` picks the nearest available one per sample, or give a length. | string or integer | `auto` |
| `--bracken_threshold` | Minimum reads for a taxon to be re-estimated. | integer | `10` |

<details><summary>Advanced</summary>

| Option | Description | Type | Default |
|---|---|---|---|
| `--bracken_extra_args` | Extra arguments passed to Bracken. | string |  |

</details>

## Taxon filter

| Option | Description | Type | Default |
|---|---|---|---|
| `--taxon_filter` | Restrict the filtered tables to bacteria and archaea. | `none` / `prokaryote` | `none` |
| `--taxdump` | NCBI taxdump: a directory with nodes.dmp (and names.dmp, for SourcePredict), or the nodes.dmp file itself, with the others beside it; local or a URL. Needed for `--taxon_filter prokaryote`, SourcePredict and authentication. Use the snapshot the classifier database was built from. | string |  |
| `--drop_taxids` | Comma-separated taxon IDs dropped from every table. | string or integer | `0,1` |

## Authentication (optional)

| Option | Description | Type | Default |
|---|---|---|---|
| `--authentication_genomes` | Reference genomes of the taxa to authenticate: a FASTA, plain or gzipped, whose headers carry `kraken:taxid\|<taxid>`; local or a URL. A Bowtie2 index beside it (`<fasta>.1.bt2` and the rest) is used when present, and built otherwise. Needs `--taxdump`. See docs/authentication.md. | string |  |
| `--authentication_top_taxa` | Authenticate this many of the most abundant taxa in each ancient library. | integer | `10` |
| `--authentication_min_reads` | A taxon needs this many reads in the count table to be chosen, and this many aligned reads for its damage to be tested. | integer | `100` |

## Source attribution (optional)

| Option | Description | Type | Default |
|---|---|---|---|
| `--attribution_methods` | Methods to run, comma-separated; empty runs none. `sourcetracker2`, `feast`, `fastst` (experimental), `sourceid_nmf` and `sourcepredict`. All read the same table, so their results can be compared; see docs/attribution.md. | string |  |
| `--attribution_label` | Name of this analysis; results go to `attribution/<label>/`. | string | `default` |
| `--attribution_aggregate_by` | Unit each method models: each library, or all runs sharing a `unit_id`. | `sample_id` / `unit_id` | `unit_id` |
| `--attribution_taxon_scope` | Taxa given to the methods: after the taxon filter, or all. | `filtered` / `all` | `filtered` |
| `--attribution_cv` | Also score every method by cross-validation on the sources. Sources sharing a value of `attribution_cv_group` are held out together, run as sinks of known origin, and scored. Results go to `attribution/<label>/evaluation/`. | boolean | `false` |
| `--attribution_cv_group` | Samplesheet column whose values are held out together in cross-validation. Holding out whole studies tests a method on studies it has not seen; libraries without a value are held out one by one. | string | `project` |
| `--attribution_max_unknown` | Abstain when the unknown share reaches this value. | number | `0.5` |
| `--attribution_min_top_known` | Abstain when the leading known source is below this value. | number | `0.5` |
| `--attribution_min_margin` | Abstain when the two leading known sources are closer than this. | number | `0.05` |
| `--attribution_seed` | Random seed for every method. | integer | `1` |

<details><summary>Advanced</summary>

| Option | Description | Type | Default |
|---|---|---|---|
| `--attribution_min_classified_fraction` | Leave out libraries with a smaller Kraken2 classified fraction. | number | `0.001` |
| `--attribution_min_retained_reads` | Leave out libraries with fewer reads after the taxon filter. | integer | `1000` |
| `--attribution_min_retained_fraction` | Leave out libraries keeping a smaller fraction of reads after the taxon filter. | number | `0.001` |

</details>

## SourceTracker2

| Option | Description | Type | Default |
|---|---|---|---|
| `--sourcetracker2_source_rarefaction_depth` | Rarefy pooled sources to this depth (0: no rarefaction). | integer | `1000` |
| `--sourcetracker2_sink_rarefaction_depth` | Rarefy sinks to this depth (0: no rarefaction). | integer | `1000` |
| `--sourcetracker2_alpha2` | Prior on the unknown source. SourceTracker2's default. With few taxa per sink the unknown source can absorb sinks that clearly match a source; 0.1, the default of the original SourceTracker, guards against that. See docs/attribution.md. | number | `0.001` |

<details><summary>Advanced</summary>

| Option | Description | Type | Default |
|---|---|---|---|
| `--sourcetracker2_alpha1` | Prior on the known sources. | number | `0.001` |
| `--sourcetracker2_beta` | Prior on the mixing proportions. | number | `10` |
| `--sourcetracker2_restarts` | Independent Gibbs chains per sink. | integer | `10` |
| `--sourcetracker2_draws_per_restart` | Draws taken from each chain. | integer | `1` |
| `--sourcetracker2_burnin` | Passes before the first draw. | integer | `100` |
| `--sourcetracker2_delay` | Passes between draws. | integer | `10` |

</details>

## FEAST, SourceID-NMF and SourcePredict

| Option | Description | Type | Default |
|---|---|---|---|
| `--feast_em_iterations` | FEAST: expectation-maximisation iterations. | integer | `1000` |
| `--feast_coverage` | FEAST: rarefaction depth (0: the smallest sink or pooled source). | integer | `0` |
| `--sourceid_nmf_code` | SourceID-NMF: where its code is fetched from, a URL or a directory. SourceID-NMF has no licence, so MSTFlow does not redistribute it. By default its files are fetched from the authors' repository at a fixed commit; on machines without internet access, clone the repository and give its directory. | string | `https://raw.githubusercontent.com/ZiyiHuang0708/SourceID-NMF/4fd1cdff0b2781cdae8bf1f5c71aa6b72779578a` |
| `--sourceid_nmf_iterations` | SourceID-NMF: maximum iterations. | integer | `2000` |
| `--sourcepredict_neighbors` | SourcePredict: neighbours for its classifier (0: chosen automatically). | integer | `0` |

<details><summary>Advanced</summary>

| Option | Description | Type | Default |
|---|---|---|---|
| `--sourceid_nmf_rho` | SourceID-NMF: penalty parameter. | integer | `1` |
| `--sourceid_nmf_weight` | SourceID-NMF: weighting matrix coefficient. | integer | `1` |
| `--sourceid_nmf_tolerance` | SourceID-NMF: convergence threshold. | number | `1e-06` |
| `--sourcepredict_unknown_alpha` | SourcePredict: share of the sink used to model an unknown source. | number | `0.1` |
| `--sourcepredict_normalization` | SourcePredict: count normalisation. | `SUBSAMPLE` / `GMPR` / `RLE` | `SUBSAMPLE` |
| `--sourcepredict_distance` | SourcePredict: distance between samples. | `weighted_unifrac` / `unweighted_unifrac` | `weighted_unifrac` |
| `--sourcepredict_embedding` | SourcePredict: embedding (auto: t-SNE for large panels, MDS for small ones). | `auto` / `TSNE` / `UMAP` / `MDS` | `auto` |

</details>

## Sample summary flags

| Option | Description | Type | Default |
|---|---|---|---|
| `--qc_min_classified_reads` | Flag samples with fewer classified molecules than this. | integer | `10000` |
| `--qc_min_taxa` | Flag samples with fewer taxa (at bracken_rank, after the taxon filter) than this. Shallow profiles carry little information, and downstream attribution of them is unreliable. | integer | `20` |

## Containers

<details><summary>Advanced</summary>

| Option | Description | Type | Default |
|---|---|---|---|
| `--container_fastqc` | Container image for FastQC. | string |  |
| `--container_fastp` | Container image for fastp. | string |  |
| `--container_adapterremoval` | Container image for AdapterRemoval. | string |  |
| `--container_bwa` | Container image for BWA. | string |  |
| `--container_kraken2` | Container image for Kraken2. | string |  |
| `--container_bracken` | Container image for Bracken. | string |  |
| `--container_multiqc` | Container image for MultiQC. | string |  |
| `--container_python_tools` | Container image for the workflow's Python helpers. | string |  |
| `--container_sourcetracker2` | Container image for SourceTracker2. | string |  |
| `--container_sourcepredict` | Container image for SourcePredict. | string |  |
| `--container_attribution` | Container image for FEAST, FastST and SourceID-NMF. | string |  |
| `--container_krakenuniq` | Container image for KrakenUniq. | string |  |
| `--container_bowtie2` | Container image for Bowtie2 and samtools. | string |  |
| `--container_pydamage` | Container image for PyDamage. | string |  |

</details>
