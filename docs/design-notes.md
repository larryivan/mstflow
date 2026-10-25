# Design notes

Why MSTFlow counts the way it does. None of this is needed to run it. How it drives Kraken2 is in [internals](internals.md).

## Counting fragments

Every count is a sequenced fragment: a read pair, a merged read and a single read each count once. PCR duplicates are not removed.

Ancient paired-end libraries keep only the pairs that overlap and merge (`ancient_pe_mode: collapsed_only`, the default), and the excluded counts are reported. `all_fragments` also keeps complete unmerged pairs.

## Ancient and modern libraries in one analysis

Source panels for ancient samples are mostly modern, so an ancient sink has to be attributed against modern sources. Both kinds of library are trimmed to the same minimum length (`--min_length`, 35 bp), counted in fragments and classified against one database. Modern libraries cut to ancient lengths and damaged in silico are attributed as the originals are ([the test](validation.md#ancient-and-modern-libraries-in-one-analysis)).

## Why not Bracken

A sink is a mixture. Bracken moves the reads a classifier could only place at family level or above into the genera below, in proportion to what the same sample holds, so one source's reads end up in another's genera. A classifier places each read on its own, so a mixture's counts are the sum of its sources' counts.

For the same reason there is no minimum of reads per taxon and sample: a genus kept in a source can fall below it in a mixture. The cost is that reads placed only above the rank are not counted, which leaves poorly covered samples with fewer counts. [Validation](validation.md#choices-the-tests-decided) measures both choices.
