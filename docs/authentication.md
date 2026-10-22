# Authentication

A taxon in an ancient library's profile may come from the ancient sample or may have
been added since: in the ground, in the museum, in the laboratory. Ancient DNA carries
a mark that later contamination lacks: cytosines at the ends of its molecules have
deaminated, so its reads show C-to-T changes against the reference near their 5' ends.
MSTFlow tests the most abundant taxa of each ancient library for that damage.

```bash
--authentication_genomes genomes.fna.gz --taxdump taxdump/
```

## How it works

1. **Choose the taxa.** In every library with `age_group` ancient, the
   `--authentication_top_taxa` (10) taxa with the most reads in the filtered count table
   are chosen, each with at least `--authentication_min_reads` (100).
2. **Align.** The library's reads, after host removal, are aligned end to end to the
   reference genomes with Bowtie2 (`--very-sensitive`). Only reads placed uniquely
   (MAPQ 30 or more) are kept, so a read that fits several related genomes counts for
   none of them.
3. **Test for damage.** The reads aligned to each taxon's genomes are pooled, and
   PyDamage fits its damage model to them. MSTFlow also reports the share of the
   genomes they cover.

## Reference genomes

`--authentication_genomes` is a FASTA, plain or gzipped, local or a URL, with one or
more genomes for each taxon you want tested. Every header carries the genome's taxon in
Kraken's format:

```
>NZ_LR134406.1|kraken:taxid|1750
```

so the genomes used to build a Kraken2 database can be reused as they are. A genome
belongs to a profiled taxon when its taxon descends from it in `--taxdump`, so species
genomes serve a genus-level table. A chosen taxon with no genome is reported as such.

Index the FASTA once with `bowtie2-build genomes.fna.gz genomes.fna.gz` and keep the
index beside it (`genomes.fna.gz.1.bt2` and the rest); without an index MSTFlow builds
one in each run. A few genomes per taxon are enough. With many closely related genomes,
fewer reads align uniquely to each of them.

## Output

`authentication/authentication.tsv`, one row per library and taxon:

| Column | Meaning |
|---|---|
| `table_reads` | Reads of the taxon in the count table |
| `aligned_reads` | Reads aligned uniquely to its genomes |
| `breadth` | Share of its genomes' bases covered by at least one read |
| `c_to_t_first_base` | How often a reference C is read as T at the first base of the reads |
| `damage_pmax` | PyDamage's estimate of the damage at the read ends |
| `qvalue` | PyDamage's q-value for damage against no damage |
| `predicted_accuracy` | PyDamage's estimate of how reliable its call is at this depth |
| `verdict` | `damaged`, `no damage`, `too few reads` or `no reference genome` |

A taxon is `damaged` when the q-value is below 0.05 and the predicted accuracy is at
least 0.5, PyDamage's own thresholds. `no damage` does not prove a taxon modern: with
little damage or few reads the test has little power, which the predicted accuracy
shows. A taxon that is damaged in an ancient sample but not in your blanks is the
clearest case. Reads spread evenly over a genome also argue for a real organism, while
reads piled on a few regions, with low breadth for their number, suggest a
misclassification.

`authentication/targets.tsv` lists the taxa chosen, and `authentication/bowtie2/` holds
each library's alignments and Bowtie2's log.

Treat the result as one line of evidence, together with fragment lengths, the
negative controls and what is known of the sample, not as proof on its own.
