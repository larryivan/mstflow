# Design notes

Why MSTFlow does some things the way it does. None of this is needed to run the
workflow; it is here for anyone who wants to check the reasoning or change it.

## Counting fragments

Every count in MSTFlow is a sequenced fragment: a read pair from a paired-end library is
one fragment, and so are a collapsed read and a single read from a single-end library.
PCR duplicates are not removed, so a molecule copied several times before sequencing is
counted several times.

Ancient paired-end libraries default to `ancient_pe_mode: collapsed_only`: only read pairs that overlap and collapse into one fragment are kept. The excluded singleton and unmerged counts are reported. The optional `all_fragments` mode retains all complete molecular fragments (collapsed reads and complete unmerged pairs), preserves pairing through host removal, and submits each surviving fragment as one native Kraken2 `--paired` record. Orphan singletons are quantified but excluded, matching modern paired-end handling.

## Ancient and modern libraries in one analysis

Source panels for ancient samples are mostly modern, so an ancient sink has to be
attributed against modern sources. The two kinds of library differ in three technical
ways: ancient fragments are short (median 72 bp in the test data, against reads of 76
to 151 bp in its modern libraries), they carry deamination at their ends, and they are
trimmed by AdapterRemoval rather than fastp. Both are trimmed to the same minimum length
(`--min_length`, 35 bp), counted in fragments and classified against one database.

To check that the rest does not bias attribution, each of the seven modern source
libraries in the test data was made into two ancient-like libraries. Every read 1 was
cut to a length drawn from the test data's ancient fragments; in the second library its
ends were also given C-to-T and G-to-A changes, with probability 0.30 × 0.55^i at i bases
from the end. Both were written as read pairs that run into the adapters, so they took
the ancient branch of the workflow, and each was attributed with SourceTracker2 with its
own study's sources left out, as in cross-validation.

| Library | Profile change (Bray-Curtis) | Two libraries of one environment | Own environment: modern / cut / cut and damaged |
|---|--:|--:|---|
| oral_1 | 0.02 | 0.61 | 0.97 / 0.96 / 0.97 |
| oral_2 | 0.03 | 0.68 | 0.99 / 0.99 / 0.99 |
| gut_1 | 0.04 | 0.86 | 0.94 / 0.96 / 0.96 |
| gut_2 | 0.04 | 0.68 | 0.62 / 0.57 / 0.57 |
| gut_3 | 0.04 | 0.65 | 0.87 / 0.87 / 0.83 |
| soil_1 | 0.11 | 0.24 | 0.94 / 0.89 / 0.89 |
| soil_2 | 0.17 | 0.24 | 0.88 / too few counts |

The test data has a single soil study, so each soil library was attributed with the
other as the only soil source.

For oral and gut libraries, cutting to ancient lengths changed the genus profile more
than ten times less than two libraries of one environment differ, and the estimate by
at most five percentage points; damage added little. Short reads are classified less
often (57% of oral_1 fell to 43%), and soil, of which only 3 to 6% is classified with
this database, lost more of its profile: soil_1 cut short still went to soil, while
soil_2 kept too few counts to be attributed at all. What remains between an ancient
sink and a modern source is the biology that attribution is meant to measure.

## Why not Bracken

Bracken re-estimates abundances by moving the reads a classifier could only place at a
family or higher down to the genera below, in proportion to the genera found in the
same sample. That suits a profile of one sample, but a sink is a mixture, and the reads
of one source get moved into the genera of another. In the test data's gut and soil
mixture (10% and 90% of read pairs), Bracken gave Bradyrhizobium, a soil genus, 275
reads where the soil's share of its own count was 705, and Escherichia, a gut genus, 397
where 45 were expected: the soil's reads placed only to Proteobacteria went to Escherichia once
the faeces brought it in. SourceTracker2 then put 12% of the mixture in soil, against
the 24% its reads predicted.

A classifier places each read on its own, so a mixture's counts are the sum of its
sources' counts. With the classifier's own counts, every mixture in the test data came
within 0.06 of its known composition, and the soil one within 0.01; in cross-validation
FEAST rose from 57% to 100% correct and FastST from 29% to 100%. The cost is that reads
placed only above the rank are not counted, which leaves poorly covered samples, such as
soils in a small database, with fewer counts.

For the same reason the table keeps every read, with no minimum per taxon and sample.
Bracken drops taxa below 10 reads; applied to the classifier's counts, that minimum
alone brought FEAST back to 57% and FastST to 29%, and moved the ancient-like mixture
0.11 from its composition. A genus with 15 reads in a source and 7 in a mixture of it
is kept in one and dropped from the other, so the counts no longer add up, and FEAST and
FastST lean on rare genera. Rare genera do make a shallow sample look rich, so the depth
flag in the sample summary counts only taxa with `--qc_min_taxon_reads` (10) reads.

## Why `all_fragments` uses native paired input

Kraken2 2.1.6 counts every minimizer position in `total_kmers`, including the ambiguous positions produced by an `N` run, and subtracts only its own mate-pair and reading-frame markers (`src/classify.cc`, `ClassifySequence`). Because `required_score = ceil(confidence * total_kmers)`, concatenating two 40 bp mates with a 100 `N` spacer raises the denominator from 12 native paired positions to 146, and the confidence requirement at 0.05 from 1 hit to 8. Short ancient fragments that native `--paired` classifies are therefore lost. MSTFlow submits both mates natively and gives each collapsed single-end fragment a placeholder mate shorter than `kraken2_kmer_len`, which contributes exactly zero minimizer positions and so reproduces native single-end classification. Placeholder mates are excluded from the informative length and k-mer-equivalent length statistics.

`kraken2_confidence`'s effective denominator is recorded per sample as `median_kraken_confidence_denominator` in `*.kraken2.read_length.tsv`.

`--check-names` is not a Kraken2 2.1.6 option. Its Perl wrapper parses long options with `Getopt::Long` and never inspects the parse result, so an unrecognised long option prints `Unknown option: <name>` and classification still succeeds with exit status 0. MSTFlow therefore rejects `--check-names` during parameter validation, verifies paired record counts, normalised mate identifiers, and mate ordering in `bin/mstflow_fragment_tools.py` before Kraken2 starts, and fails the task if Kraken2 reports that it ignored any option.

## Parameter semantics

`kraken2_min_read_len` must be at least `kraken2_kmer_len`. The default Kraken confidence and minimum-hit-group settings are explicit so sources and sinks cannot silently use different classifier defaults. Per-read classifications are disabled by default; the reports are kept.

`host_min_mapq: 0` removes an entire fragment when either mate has any primary host alignment. A threshold such as 20 is suitable only as a separately declared sensitivity analysis. Minimizer data and taxon names are enabled in the Kraken QC report; MSTFlow also emits the standard six-column report.

`require_explicit_host_for_nondefault: true` prevents a non-human `host_group` from silently using the default human reference. Each such sample must provide its own indexed `host_ref`.

`host_action: remove` sends the sample through BWA host subtraction. `host_action: skip` bypasses BWA and passes the preprocessed reads directly to clean-read QC and classification; skip rows must use `host_ref: NA`. A samplesheet without a `host_action` column treats every row as `remove`; set `require_explicit_host_action: true` to make every row declare it.

Kraken memory depends on the selected database. The CPU, memory, time, and concurrency defaults are explicit and should be overridden after measuring the intended database and executor. Host removal is capped at four concurrent tasks by default to control repeated index I/O.

`fastq_gzip_level` (0-9, default 6) sets the compression level for FASTQ written by `bin/mstflow_fragment_tools.py`, matching the AdapterRemoval and `gzip(1)` defaults instead of Python's level 9. It changes only the compressed representation: decompressed FASTQ content, fragment membership, and counts are unaffected. Lower it towards 1 to trade roughly 10% more scratch space for markedly faster host-removal and Kraken input writes; level 0 stores without compressing. The helper fixes the gzip header timestamp, so rewriting the same file gives identical bytes.

Source attribution works on biological units (`attribution_aggregate_by: unit_id`): the lanes of a unit are checked and filtered together, so one poor lane cannot leave a partial unit. A sink is assigned to its leading source only when Unknown is below `attribution_max_unknown`, the leading source reaches `attribution_min_top_known`, and it beats the runner-up by `attribution_min_margin`; otherwise the result is reported as uncertain. [attribution.md](attribution.md) describes each method and the cross-validation.

`publish_intermediates: false` keeps large FASTQ and per-read intermediates in the resumable Nextflow work directory instead of duplicating them into the results tree.
