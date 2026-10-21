# Usage

## Samplesheet

A tab-separated file with one row per sequencing library, for example:

| sample_id | fastq_1 | fastq_2 | age_group |
|---|---|---|---|
| bone_01 | reads/bone_01_R1.fastq.gz | reads/bone_01_R2.fastq.gz | ancient |
| bone_02 | reads/bone_02.fastq.gz | | ancient |
| gut_01 | reads/gut_01_R1.fastq.gz | reads/gut_01_R2.fastq.gz | modern |

The [test samplesheet](https://github.com/larryivan/mstflow/releases/download/test-data-v1/samplesheet.tsv)
is a working example, with the source attribution columns and reads named relative to it.

| Column | Required | Meaning |
|---|---|---|
| `sample_id` | yes | Unique name; letters, digits, `.`, `_` and `-` only |
| `fastq_1` | yes | Read 1, or the only read file for single-end libraries |
| `fastq_2` | no | Read 2; leave empty (or omit the column) for single-end libraries |
| `age_group` | yes | `ancient` or `modern`; decides how reads are trimmed and counted |
| `layout` | no | `PE` or `SE`; inferred from `fastq_2` when left out |
| `host_action` | no | `remove` (default) or `skip`, for libraries with no host DNA such as sediment |
| `host_ref` | no | A different indexed host FASTA for this library; empty or `default` uses `--host_fasta` |
| `host_group` | no | Name of the host for this library, such as `pig`; needs its own `host_ref` |
| `unit_id` | no | Groups several runs of one biological sample for source attribution |
| `source_sink`, `model_source_env` | no | Needed only for source attribution; see below |
| `project` | no | The study a library comes from; source attribution's cross-validation holds out one study at a time |

Paths can be absolute, URLs, or relative. Relative paths are read relative to the
samplesheet, so a sheet and its reads can be moved together; this also works when the
samplesheet itself is a URL.

Ancient libraries are adapter-trimmed and, if paired, their overlapping mates are
collapsed into single molecules with AdapterRemoval. Modern libraries are trimmed with
fastp. From then on every library is counted in molecules; see
[design notes](design-notes.md).

## References

### Kraken2 and Bracken database

Use a prebuilt index from <https://benlangmead.github.io/aws-indexes/k2>. Pass it with
`--kraken2_db` as the downloaded `.tar.gz`, its URL, or the directory you unpacked it
into; MSTFlow unpacks archives itself. Prebuilt indexes already contain the Bracken
`database<N>mers.kmer_distrib` files.

Kraken2 loads the whole database into memory, so set `--kraken2_memory` to at least the
size of `hash.k2d` plus a few GB. If a Kraken2 task is killed for memory it is retried
once with twice as much.

To build your own database, follow the Kraken2 manual, then run `bracken-build` once for
each read length you expect (for example `-l 50` for ancient and `-l 150` for modern
libraries). MSTFlow picks the nearest length for each sample.

#### Removing poorly supported taxa

A false positive often rests on many reads from a few short stretches of a genome: a
conserved gene, a low-complexity region, a contaminant fragment in a reference. Reads from
an organism that is really there are spread along its genome, so they cover many
distinct minimizers. `--kraken2_min_distinct_minimizers` removes a taxon from a sample
when its reads cover fewer than that; the removals are listed in
`tables/bracken_<rank>_kmer_filter.tsv`. The filter is off by default. There is no
threshold that suits every depth and database: start low, read what the table removes,
and, with source attribution, compare the cross-validation scores with and without it.

### KrakenUniq instead of Kraken2

```bash
--classifier krakenuniq --krakenuniq_db krakenuniq_db/
```

KrakenUniq classifies by exact k-mers and counts the unique k-mers each taxon's reads
cover; taxa below `--krakenuniq_min_kmers` (1,000) unique k-mers or
`--krakenuniq_min_reads` (200) reads are removed sample by sample, the thresholds aMeta
uses. Bracken then re-estimates abundances, as after Kraken2, so the two classifiers
give tables of the same kind. That needs Bracken's files in the KrakenUniq database:

```bash
bracken-build -d krakenuniq_db/ -y krakenuniq -k 31 -l 50
```

once per read length. The database must have been built without
`--taxids-for-genomes` and `--taxids-for-sequences`: the taxa those options add are
not in the taxonomy, and Bracken cannot place them. For a database without Bracken's
files, such as aMeta's prebuilt one, `--run_bracken false` uses KrakenUniq's own clade
counts at `--bracken_rank` instead, as aMeta does.
KrakenUniq databases can be far larger than memory: `--krakenuniq_preload_size 32G`
loads one in chunks of that size. Its index, `database.idx`, is always loaded whole; with
the default minimizer length of 15 it is 8 GB. `--kraken2_cpus`, `--kraken2_memory` and
`--kraken2_time` set KrakenUniq's resources too.

Run sources and sinks through the same classifier and database; profiles from different
ones cannot be compared.

### Host genome

Any FASTA works; index it once with BWA:

```bash
bwa index chm13v2.0.fa
```

Pass it with `--host_fasta`. Libraries with no host DNA can set `host_action` to `skip`,
and `--remove_host false` turns host removal off for the whole run.

### Taxonomy

`--taxon_filter prokaryote` keeps only bacteria and archaea in the filtered tables, and
SourcePredict computes distances on the taxonomy. Both read an NCBI taxdump from the
same snapshot as the classifier database, given with `--taxdump`: a directory with
`nodes.dmp` and `names.dmp`, or the `nodes.dmp` file with the others beside it.

## Running

```bash
nextflow run larryivan/mstflow -r 0.1.0 -profile docker \
    --input samplesheet.tsv \
    --kraken2_db k2_standard/ \
    --host_fasta chm13v2.0.fa \
    --outdir results
```

- `-profile` chooses how the tools are installed: `docker`, `singularity`, `apptainer`
  or `conda`.
- `-r 0.1.0` pins a release so the analysis can be rerun exactly.
- `-resume` continues an interrupted run without redoing finished steps.
- Many options are easier to keep in a file: `-params-file params.yml`. See
  [examples/params.example.yml](../examples/params.example.yml).
- `--help` lists every option; [parameters](parameters.md) has the full table.

### On a cluster

Nextflow submits each step as a cluster job when you tell it which scheduler to use.
Put this in a file, for example `site.config`:

```groovy
process {
    executor = 'slurm'          // or 'pbs', 'sge', 'lsf' ...
    queue    = 'your-queue'
}
```

and add `-c site.config` to the command. Accounts, queues, scratch paths and other
site settings belong in this file, not in the workflow.

If compute nodes cannot download container images, convert each `container_*` image
named in [nextflow.config](../nextflow.config) to a `.sif` file once, and point those
parameters at the files in the same config:

```groovy
params.container_kraken2 = '/shared/containers/kraken2-2.1.6.sif'
```

### Apple Silicon and other ARM machines

The containers are built for `linux/amd64`. The `docker` profile asks Docker to run them
under emulation, which works but is several times slower. Use an x86-64 Linux machine
for real datasets.

## Source attribution

Add `source_sink` (`source` or `sink`) and `model_source_env`, the environment each
source belongs to, to the samplesheet, and choose the methods:

```bash
--attribution_methods sourcetracker2,feast --attribution_cv true
```

Sinks the sources cannot explain well enough are reported as `ABSTAIN` rather than
forced into a source. [Source attribution](attribution.md) describes the five methods,
their licences, how calls are made and how cross-validation scores them.

## Reproducibility

Every run writes the validated samplesheet and the Nextflow report, timeline and trace
to `results/pipeline_info`, and each tool's version next to its output.
Together with a pinned release (`-r`) and the database snapshot, that is everything
needed to rerun an analysis exactly.
