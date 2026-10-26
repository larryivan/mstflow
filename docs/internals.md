# Internals

For maintainers: how MSTFlow drives Kraken2, and what some parameters guarantee. Running the workflow needs none of it.

## `all_fragments`

`ancient_pe_mode: all_fragments` keeps all complete molecular fragments, collapsed reads and complete unmerged pairs. It preserves pairing through host removal and submits each surviving fragment as one native Kraken2 `--paired` record. Orphan singletons are quantified but excluded, as for modern paired-end libraries.

## Why native paired input

Kraken2 2.1.6 counts every minimizer position in `total_kmers`, ambiguous positions from an `N` run included, and subtracts only its own mate-pair and reading-frame markers (`src/classify.cc`, `ClassifySequence`). Since `required_score = ceil(confidence * total_kmers)`, joining two 40 bp mates with a 100 `N` spacer raises the denominator from 12 native paired positions to 146, and the 0.05 confidence requirement from 1 hit to 8. Short ancient fragments that native `--paired` classifies would be lost.

So MSTFlow submits both mates natively and gives each collapsed single-end fragment a placeholder mate shorter than `kraken2_kmer_len`. The placeholder adds exactly zero minimizer positions, which reproduces native single-end classification, and it is left out of the informative-length and k-mer-equivalent-length statistics. The effective denominator is recorded per sample as `median_kraken_confidence_denominator` in `*.kraken2.read_length.tsv`.

## `--check-names`

It is not a Kraken2 2.1.6 option. Kraken2's Perl wrapper parses long options with `Getopt::Long` and never inspects the result, so an unknown long option prints `Unknown option: <name>` and classification still succeeds with exit status 0. MSTFlow therefore

- rejects `--check-names` during parameter validation,
- checks paired record counts, normalised mate identifiers and mate order in `bin/mstflow_fragment_tools.py` before Kraken2 starts, and
- fails the task if Kraken2 reports that it ignored any option.

## Parameter semantics

- `kraken2_min_read_len` must be at least `kraken2_kmer_len`.
- The Kraken2 confidence and minimum-hit-group settings are explicit, so sources and sinks cannot silently use different classifier defaults. Per-read classifications are off by default; the reports are kept. Minimizer data and taxon names are on in the QC report, and MSTFlow also writes the standard six-column report.
- `host_min_mapq: 0` removes a whole fragment when either mate has any primary host alignment. A threshold such as 20 suits only a separately declared sensitivity analysis.
- `require_explicit_host_for_nondefault: true` stops a non-human `host_group` from silently using the default human reference; each such sample must give its own indexed `host_ref`.
- `host_action: remove` sends a sample through BWA host subtraction. `skip` bypasses BWA and passes the preprocessed reads to clean-read QC and classification, and skip rows must use `host_ref: NA`. Without a `host_action` column every row is `remove`; `require_explicit_host_action: true` makes every row declare it.
- Kraken memory depends on the database. The CPU, memory, time and concurrency defaults are explicit and should be overridden after measuring your database and executor. Host removal runs at most four tasks at once by default, to limit repeated index I/O.
- `fastq_gzip_level` (0-9, default 6) is the compression level of the FASTQ written by `bin/mstflow_fragment_tools.py`, matching AdapterRemoval and `gzip(1)` instead of Python's 9. It changes only the compressed bytes, not content, fragment membership or counts. Lowering it towards 1 costs roughly 10% more scratch space and makes host-removal and Kraken input writes much faster; 0 stores without compressing. The gzip header timestamp is fixed, so rewriting a file gives identical bytes.
- Attribution works on biological units (`attribution_aggregate_by: unit_id`): the lanes of a unit are checked and filtered together, so one poor lane cannot leave a partial unit.
- `publish_intermediates: false` keeps large FASTQ and per-read intermediates in the resumable work directory instead of copying them into the results.
