# Changelog

## 0.1.0 (2026-10-21)

First public release.

### Added
- Source attribution with five methods on the same count table: SourceTracker2, FEAST,
  FastST (experimental), SourceID-NMF and SourcePredict (`--attribution_methods`). Their
  estimates, calls and agreement are written side by side, with one figure
  (`attribution/<label>/`). See docs/attribution.md.
- `-profile conda` runs every step without containers, FEAST and FastST included. Their R
  packages are built once at the start of the run from the commits in the attribution image
  (`--feast_code`, `--fastst_code`, `--npreg_code`), in an environment whose R and OpenBLAS
  are the image's; on the test data their estimates are identical to the image's. The test
  profile runs under Conda in CI as well as under Docker.
- Cross-validation of every method on the sources (`--attribution_cv`): the sources of
  one study at a time are held out as sinks of known origin, and each method's coverage
  and accuracy on them are reported.
- KrakenUniq as an alternative to Kraken2 (`--classifier krakenuniq`), with aMeta's
  unique k-mer threshold.
- Authentication (`--authentication_genomes`): the most abundant taxa of each ancient
  library are aligned to reference genomes with Bowtie2 and tested for ancient DNA
  damage with PyDamage, with the breadth of genome they cover. See
  docs/authentication.md.
- `--kraken2_min_distinct_minimizers` removes taxa whose reads cover too few distinct
  minimizers, sample by sample, and lists what it removed.
- `-profile test` runs every step, all attribution methods and their cross-validation
  included, on excerpts of 18 published ancient and modern libraries and three sinks
  mixed from them in known shares, with a 274-genome Kraken2 database,
  downloaded from a release asset (about 1.9 GB). The recipe is in `assets/test`, and
  `assets/test/check_mixtures.py` checks the mixtures' estimates.
- `--kraken2_db` accepts a `.tar.gz` archive, local or a URL, such as the prebuilt
  indexes, and unpacks it. The samplesheet, reads, host genome and taxdump can also be
  URLs.
- `--help`, with parameter types, defaults and accepted values, and validation of every
  parameter against `nextflow_schema.json`.
- `summary/sample_summary.tsv`: fragments after each step, classified fraction, taxa
  detected and depth flags for every sample, also shown in the MultiQC report.
  Thresholds: `--qc_min_classified_reads`, `--qc_min_taxa`.
- A summary printed at the end of each run: where the results are and which samples
  were flagged.
- Tasks killed for memory or time are retried once with twice the resources.
- A check at the start that the software environment chosen with `-profile` is installed
  (Docker running, Apptainer, Singularity, conda or mamba). If not, the run stops before
  anything is downloaded, listing what the machine has and the `-profile` to use, instead
  of failing in the first task; `--check_software false` skips it. The conda profile uses
  conda when mamba is not installed.
- MIT licence, `CITATION.cff`, CI, and documentation under `docs/`, with
  `docs/validation.md`: what the test data shows about the shares, the calls, and ancient
  sinks against modern sources, with figures.

### Changed
- Samplesheet: only `sample_id`, `fastq_1` and `age_group` are required. `layout` is
  inferred from `fastq_2`, which may be left empty or omitted, and libraries with
  `host_action` set to `skip` need no `host_ref`. Relative paths are read relative to
  the samplesheet (paths relative to the launch directory still work).
- The count table holds the reads the classifier placed in each taxon at `--rank` or
  below it; Bracken is no longer run. Its redistribution of reads placed above the rank
  follows the composition of each sample, so in a sink mixed from two sources it moved
  one source's reads into the other's genera. Every read counts, with no per-sample
  minimum of reads per taxon, which would break that sum as well; the depth flag counts
  only taxa with `--qc_min_taxon_reads` (10) reads. Tables are named `kraken2_<rank>_*`
  or `krakenuniq_<rank>_*`.
- Parameters are validated strictly against their declared types: booleans must be
  `true` or `false`.
  On/off options given on the command line, such as `--remove_host false`, now take
  effect; Nextflow passes them as text, which had read as on.
- The `docker` profile runs the `linux/amd64` images on any host, including ARM.
- One minimum read length, `--min_length` (35 bp), for ancient and modern libraries
  alike, in place of `--min_length_ancient` and `--min_length_modern` (50 bp).
- `--sourcetracker2_alpha2` defaults to 0.1, the default of the original SourceTracker,
  instead of the 0.001 of SourceTracker2 2.0.1, with which the unknown source can
  absorb sinks that clearly match a source.
- The workflow is split into preprocessing, profiling and attribution subworkflows.
- SourceTracker2 options are now attribution options (`--attribution_*`) shared by
  every method, and its priors and sampling settings are parameters
  (`--sourcetracker2_alpha2` and others) instead of free-form extra arguments. Results
  move from `source_attribution/sourcetracker2/` to `attribution/`.

### Fixed
- Ancient libraries sequenced on instruments that write quality 42 ('K'), such as the
  HiSeq 4000, HiSeq X and NovaSeq, no longer stop AdapterRemoval.

### Removed
- The `slurm`, `slurm_conda`, `slurm_apptainer` and `slurm_singularity` profiles.
  Choose a scheduler in your own config instead, passed with `-c`; see docs/usage.md.
- `--run_sourcetracker2`, `--sourcetracker2_extra_args` and `--sourcetracker2_loo`;
  cross-validation by study replaces SourceTracker2's leave-one-out.
- Bracken, with `--run_bracken`, `--bracken_*` and `--krakenuniq_min_reads`; and
  `--fragment_pair_gap`.
- The `panel_role` column with its `--attribution_*_panel_roles` options, and
  `--attribution_exclude_source_envs`/`_classes`: leave `source_sink` empty to keep a
  library out of attribution. `category` no longer stands in for `model_source_env`.
- The samplesheet's `preprocess_mode` and `bracken_read_len` columns: the trimming
  follows from `age_group` and `layout`.
- `attribution_excluded_samples.tsv`: `attribution_sample_qc.tsv` says which samples
  were left out and why.
