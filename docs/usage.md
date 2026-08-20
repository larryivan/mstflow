# Usage

## Samplesheet

A tab-separated file with one row per sequencing library, for example:

| sample_id | fastq_1 | fastq_2 | age_group |
|---|---|---|---|
| bone_01 | reads/bone_01_R1.fastq.gz | reads/bone_01_R2.fastq.gz | ancient |
| bone_02 | reads/bone_02.fastq.gz | | ancient |
| gut_01 | reads/gut_01_R1.fastq.gz | reads/gut_01_R2.fastq.gz | modern |

The bundled [test samplesheet](../assets/test/samplesheet.tsv) is a working example.

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
| `unit_id` | no | Groups several runs of one biological sample for SourceTracker2 |
| `source_sink`, `model_source_env` | no | Needed only for SourceTracker2; see below |

Paths can be absolute, URLs, or relative. Relative paths are read relative to the
samplesheet, so a sheet and its reads can be moved together.

Ancient libraries are adapter-trimmed and, if paired, their overlapping mates are
collapsed into single molecules with AdapterRemoval. Modern libraries are trimmed with
fastp. From then on every library is counted in molecules; see
[design notes](design-notes.md).

## References

### Kraken2 and Bracken database

Download a prebuilt index from <https://benlangmead.github.io/aws-indexes/k2> and unpack
it into one directory. These indexes already contain the Bracken
`database<N>mers.kmer_distrib` files. Pass the directory with `--kraken2_db`.

Kraken2 loads the whole database into memory, so set `--kraken2_memory` to at least the
size of `hash.k2d` plus a few GB. If a Kraken2 task is killed for memory it is retried
once with twice as much.

To build your own database, follow the Kraken2 manual, then run `bracken-build` once for
each read length you expect (for example `-l 50` for ancient and `-l 150` for modern
libraries). MSTFlow picks the nearest length for each sample.

### Host genome

Any FASTA works; index it once with BWA:

```bash
bwa index chm13v2.0.fa
```

Pass it with `--host_fasta`. Libraries with no host DNA can set `host_action` to `skip`,
and `--remove_host false` turns host removal off for the whole run.

### Taxonomy for the prokaryote filter

`--taxon_filter prokaryote` keeps only bacteria and archaea in the filtered tables. It
needs an NCBI taxdump (`nodes.dmp`) from the same snapshot as the Kraken2 database,
given with `--taxdump`.

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

If compute nodes cannot download container images, convert each image listed in
[containers.lock.tsv](../containers.lock.tsv) to a `.sif` file once, and point the
matching `container_*` parameters at those files in the same config:

```groovy
params.container_kraken2 = '/shared/containers/kraken2-2.1.6.sif'
```

### Apple Silicon and other ARM machines

The containers are built for `linux/amd64`. The `docker` profile asks Docker to run them
under emulation, which works but is several times slower. Use an x86-64 Linux machine
for real datasets.

## SourceTracker2

Set `--run_sourcetracker2 true` and add two columns to the samplesheet: `source_sink`
(`source` or `sink`) and `model_source_env`, the environment each source belongs to.
Sinks whose estimate is too uncertain are reported as unassigned rather than forced into a source; the thresholds
are listed under [parameters](parameters.md#sourcetracker2-attribution-optional).

## Reproducibility

Every run writes the validated samplesheet, the Nextflow report, timeline and trace,
and the version of each tool to `results/pipeline_info` and next to each output.
Together with a pinned release (`-r`) and the database snapshot, that is everything
needed to rerun an analysis exactly.
