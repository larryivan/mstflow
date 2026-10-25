# Usage

## Run

```bash
nextflow run larryivan/mstflow -r 0.1.0 -profile docker \
    --input samplesheet.tsv \
    --kraken2_db k2_standard/ \
    --host_fasta chm13v2.0.fa \
    --outdir results
```

- `-profile` is `docker`, `singularity`, `apptainer` or `conda` ([which one?](../README.md#quick-start)). MSTFlow checks it before starting; add `--check_software false` if only your compute nodes have it.
- `-r 0.1.0` pins a release, `-resume` continues an interrupted run, and `-params-file params.yml` keeps options in a file ([example](../examples/params.example.yml)).
- `--help` lists every option; [parameters](parameters.md) has the full table.

## Samplesheet

A tab-separated file with one row per library:

| sample_id | fastq_1 | fastq_2 | age_group |
|---|---|---|---|
| bone_01 | reads/bone_01_R1.fastq.gz | reads/bone_01_R2.fastq.gz | ancient |
| bone_02 | reads/bone_02.fastq.gz | | ancient |
| gut_01 | reads/gut_01_R1.fastq.gz | reads/gut_01_R2.fastq.gz | modern |

The [test samplesheet](https://github.com/larryivan/mstflow/releases/download/test-data-v1/samplesheet.tsv) is a working example with the attribution columns.

| Column | Required | Meaning |
|---|---|---|
| `sample_id` | yes | Unique name: letters, digits, `.`, `_`, `-` |
| `fastq_1` | yes | Read 1, or the only file of a single-end library |
| `fastq_2` | no | Read 2; leave empty for single-end |
| `age_group` | yes | `ancient` or `modern`; decides how reads are trimmed and counted |
| `layout` | no | `PE` or `SE`; inferred from `fastq_2` |
| `host_action` | no | `remove` (default) or `skip`, for libraries without host DNA such as sediment |
| `host_ref` | no | Another indexed host FASTA for this library; empty or `default` uses `--host_fasta` |
| `host_group` | no | Name of that host, such as `pig`; needs its own `host_ref` |
| `unit_id` | no | Groups several runs of one biological sample |
| `source_sink`, `model_source_env` | no | For [source attribution](attribution.md) |
| `project` | no | The study a library comes from; cross-validation holds out one study at a time |

Paths can be absolute, URLs or relative to the samplesheet, so a sheet and its reads can move together.

Ancient libraries are trimmed, and overlapping mates merged, with AdapterRemoval; modern ones with fastp. From then on a read pair or a merged read counts once, as a *fragment* ([design notes](design-notes.md)).

## References

**Kraken2 database.** Use a [prebuilt index](https://benlangmead.github.io/aws-indexes/k2) or any Kraken2 database, as the downloaded `.tar.gz`, its URL or the unpacked directory (`--kraken2_db`). Kraken2 loads it into memory, so set `--kraken2_memory` to the size of `hash.k2d` plus a few GB; a task killed for memory is retried once with double.

**Host genome.** Any FASTA, indexed once with `bwa index chm13v2.0.fa` and given as `--host_fasta`. Set `host_action` to `skip` for libraries without host DNA, or `--remove_host false` for the whole run.

**Taxonomy.** `--taxon_filter prokaryote` and SourcePredict need an NCBI taxdump from the same snapshot as the database. `--taxdump` takes a directory with `nodes.dmp` and `names.dmp`, or `nodes.dmp` with the others beside it.

**Counts.** The table holds, per genus (`--rank`), the reads Kraken2 placed in it or below. Reads placed only higher up are left out, and every read counts, so a mixture's counts are the sum of its sources' ([why no Bracken](design-notes.md#why-not-bracken)).

## Source attribution and authentication

Add `source_sink` and `model_source_env` to the samplesheet and choose the methods:

```bash
--attribution_methods sourcetracker2,feast --attribution_cv true
```

[Source attribution](attribution.md) describes the five methods, how calls are made and how cross-validation scores them.

To test ancient taxa for damage, give reference genomes with Kraken-style headers (`>seqid|kraken:taxid|1750`):

```bash
--authentication_genomes genomes.fna.gz --taxdump taxdump/
```

See [authentication](authentication.md).

## More options

**Drop poorly supported taxa.** `--kraken2_min_distinct_minimizers N` removes a taxon from a sample when its reads cover fewer than N distinct minimizers. False positives often rest on many reads from a few short stretches of a genome, while a real organism's reads are spread along it. Off by default; removals are listed in `tables/kraken2_<rank>_kmer_filter.tsv`. Start low, read what it removes and, with attribution, compare cross-validation scores with and without it.

**KrakenUniq instead of Kraken2.** `--classifier krakenuniq --krakenuniq_db krakenuniq_db/`. It drops taxa below `--krakenuniq_min_kmers` unique k-mers (1,000, aMeta's threshold). A database larger than memory loads in chunks with `--krakenuniq_preload_size 32G`; its index, `database.idx` (8 GB at the default minimizer length), always loads whole. Run sources and sinks through the same classifier and database.

**Conda.** `-profile conda` needs conda or mamba, for example from [Miniforge](https://github.com/conda-forge/miniforge). The first run builds one environment per group of tools, which takes several minutes, and keeps them in `$NXF_CONDA_CACHEDIR` (by default the work directory). FEAST and FastST are compiled from source, so it needs internet access to GitHub and CRAN, or local copies (`--feast_code`, `--fastst_code`, `--npreg_code`).

**Cluster.** Put the scheduler in a file, say `site.config`, and add `-c site.config`. Accounts, queues and scratch paths belong there, not in the workflow.

```groovy
process {
    executor = 'slurm'          // or 'pbs', 'sge', 'lsf' ...
    queue    = 'your-queue'
}
```

If compute nodes cannot pull images, convert each `container_*` image in [nextflow.config](../nextflow.config) to a `.sif` file once and point the parameter at it:

```groovy
params.container_kraken2 = '/shared/containers/kraken2-2.1.6.sif'
```

**Apple Silicon.** The containers are `linux/amd64`, so Docker runs them under emulation, several times slower. Use x86-64 Linux for real datasets.

## Reproducibility

Each run writes the samplesheet, the Nextflow report, timeline and trace to `results/pipeline_info`, and each tool's version next to its output. With a pinned release (`-r`) and the database snapshot, that is what you need to rerun an analysis exactly.
