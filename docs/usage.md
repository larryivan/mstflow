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
collapsed into single reads with AdapterRemoval. Modern libraries are trimmed with
fastp. From then on every library is counted in fragments, a read pair or a collapsed
read counting once; see
[design notes](design-notes.md).

## References

### Kraken2 database

Use a prebuilt index from <https://benlangmead.github.io/aws-indexes/k2>, or any Kraken2
database. Pass it with `--kraken2_db` as the downloaded `.tar.gz`, its URL, or the
directory you unpacked it into; MSTFlow unpacks archives itself.

Kraken2 loads the whole database into memory, so set `--kraken2_memory` to at least the
size of `hash.k2d` plus a few GB. If a Kraken2 task is killed for memory it is retried
once with twice as much.

The count table holds, for each genus (`--rank`), the reads Kraken2 placed in it or in a
species below it. Reads it could only place higher up, in a family say, are left out
rather than guessed into a genus, and every read counts, however rare its genus, so a
sink's counts are the sum of what its sources contribute. Bracken is not used, for the
reason in the [design notes](design-notes.md#why-not-bracken).

#### Removing poorly supported taxa

A false positive often rests on many reads from a few short stretches of a genome: a
conserved gene, a low-complexity region, a contaminant fragment in a reference. Reads from
an organism that is really there are spread along its genome, so they cover many
distinct minimizers. `--kraken2_min_distinct_minimizers` removes a taxon from a sample
when its reads cover fewer than that; the removals are listed in
`tables/kraken2_<rank>_kmer_filter.tsv`. The filter is off by default. There is no
threshold that suits every depth and database: start low, read what the table removes,
and, with source attribution, compare the cross-validation scores with and without it.

### KrakenUniq instead of Kraken2

```bash
--classifier krakenuniq --krakenuniq_db krakenuniq_db/
```

KrakenUniq classifies by exact k-mers and counts the unique k-mers each taxon's reads
cover; a taxon below `--krakenuniq_min_kmers` (1,000, aMeta's threshold) unique k-mers
is removed sample by sample. Its table is built the same way
as Kraken2's, from the reads it placed in each genus or below.
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
  or `conda`. Before anything is downloaded, MSTFlow checks that the one you chose is
  installed, and says which of the four are on the machine if it is not. If your compute
  nodes have it but the machine that starts the run does not, add `--check_software false`.
  With no `-profile`, the run goes on only if the tools are installed on the machine itself.
- `-profile conda` needs conda or mamba, for example from [Miniforge](https://github.com/conda-forge/miniforge);
  it uses mamba when there is one. The first run builds one environment per group of tools,
  which takes several minutes; they are kept in `$NXF_CONDA_CACHEDIR` (by default the
  work directory) and reused. FEAST and FastST are not in Bioconda, so the run also
  compiles them from source: it needs internet access to GitHub and CRAN, or local copies
  given with `--feast_code`, `--fastst_code` and `--npreg_code`.
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

## Authentication

```bash
--authentication_genomes genomes.fna.gz --taxdump taxdump/
```

aligns each ancient library to reference genomes of its most abundant taxa with Bowtie2
and tests each taxon for ancient DNA damage with PyDamage. The genomes need Kraken-style
headers (`>seqid|kraken:taxid|1750`); [authentication](authentication.md) has the
details.

## Reproducibility

Every run writes the validated samplesheet and the Nextflow report, timeline and trace
to `results/pipeline_info`, and each tool's version next to its output.
Together with a pinned release (`-r`) and the database snapshot, that is everything
needed to rerun an analysis exactly.
