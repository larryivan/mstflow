# Portable server setup

MSTFlow keeps scientific workflow logic separate from site deployment policy. The repository provides local Docker, Singularity, Apptainer, Conda, test, and generic Slurm profiles; accounts, queues, module names, storage paths, scratch behavior, and submission limits belong in an external site configuration.

## 1. Software

Install Nextflow 25.10.4 or newer and choose one dependency backend.
For example:

```bash
nextflow -version
singularity --version
```

The process-level container identities are recorded in `containers.lock.tsv`.
They are locked for `linux/amd64`; ARM hosts need explicit amd64 emulation, and
production cluster runs should use Linux/amd64 nodes.
The main Conda environment is defined in `envs/mstflow.yml`. SourceTracker2
uses the separate compatibility environment in `envs/sourcetracker2.yml`
because version 2.0.1 requires an older Python/NumPy stack. Conda itself remains
open source; Miniforge is the recommended community conda-forge distribution.

Resolve both environments before a production launch:

```bash
mamba env create --dry-run --file envs/mstflow.yml
mamba env create --dry-run --file envs/sourcetracker2.yml
```

## 2. Execution profiles

Local Conda execution:

```bash
nextflow run . \
  -profile conda \
  -params-file /absolute/path/to/params.yml \
  -resume
```

Local container execution:

```bash
nextflow run . \
  -profile singularity \
  -params-file /absolute/path/to/params.yml \
  -resume
```

Use `docker` or `apptainer` in place of `singularity` for those engines.

Generic Slurm execution:

```bash
nextflow -c /absolute/path/to/site.config run . \
  -profile slurm_singularity \
  -params-file /absolute/path/to/params.yml \
  -resume
```

The alternatives are `slurm_apptainer` and `slurm_conda`. The historical
`slurm` profile remains an alias for `slurm_conda`.

These profiles deliberately do not choose an account, partition, work
directory, cache directory, module system, or scratch policy. Put those values
in the external configuration.

Keep site configuration outside this repository when it contains usernames, allocation names, private paths, or infrastructure assumptions.

Host alignment and pair-aware filtering intentionally run in separate
single-tool containers. The handoff is a compressed SAM file in the Nextflow
work directory. It is not published, but it can be larger than the input and
can contain host sequence. Put `workDir` on access-controlled storage with
enough temporary capacity, retain it only while `-resume` is needed, and remove
it under the site's approved data-retention procedure after outputs and
provenance have been verified.

### Pre-staged SIF images

On clusters where compute nodes cannot pull OCI images, convert each digest in
`containers.lock.tsv` to a SIF **serially**, verify and record the resulting SIF
checksum, then override the logical image parameters in the external site
configuration:

```groovy
params {
    container_fastqc = '/shared/containers/mstflow/fastqc-0.12.1.sif'
    container_fastp = '/shared/containers/mstflow/fastp-0.24.0.sif'
    container_adapterremoval = '/shared/containers/mstflow/adapterremoval-2.3.4.sif'
    container_bwa = '/shared/containers/mstflow/bwa-0.7.19.sif'
    container_kraken2 = '/shared/containers/mstflow/kraken2-2.1.6.sif'
    container_bracken = '/shared/containers/mstflow/bracken-2.9.sif'
    container_biom = '/shared/containers/mstflow/biom-format-2.1.17.sif'
    container_multiqc = '/shared/containers/mstflow/multiqc-1.29.sif'
    container_python_tools = '/shared/containers/mstflow/mstflow-python.sif'
    container_sourcetracker2 = '/shared/containers/mstflow/mstflow-sourcetracker2.sif'
}

singularity {
    autoMounts = true
    runOptions = '--cleanenv'
}
```

Add site bind paths through `singularity.runOptions` only when automatic mounts
are insufficient. Keep cache, temporary, and work directories on a filesystem
visible from every execution node:

```bash
export NXF_SINGULARITY_CACHEDIR=/shared/cache/mstflow/singularity
export SINGULARITY_CACHEDIR=/shared/cache/mstflow/singularity-ce
export SINGULARITY_TMPDIR=/shared/tmp/mstflow/singularity
```

Do not start many first-time OCI-to-SIF conversions concurrently. Some
Singularity versions can race in the shared cache. Pre-staging also makes the
executed SIF checksum part of the frozen run provenance.

## 3. Reference validation

### Host reference

When host removal is enabled, the FASTA needs classic BWA index sidecars. A Samtools index is recommended:

```bash
ref=/absolute/path/to/host.fa
for suffix in '' .amb .ann .bwt .pac .sa .fai; do
  test -s "${ref}${suffix}" || echo "MISSING ${ref}${suffix}"
done
```

Create missing indices in an appropriate compute allocation:

```bash
bwa index /absolute/path/to/host.fa
samtools faidx /absolute/path/to/host.fa
```

### Kraken2 and Bracken

The Kraken directory must contain non-empty `hash.k2d`, `opts.k2d`, and `taxo.k2d`. Bracken additionally needs at least one `database<N>mers.kmer_distrib` file:

```bash
db=/absolute/path/to/kraken2_db
for name in hash.k2d opts.k2d taxo.k2d; do
  test -s "${db}/${name}" || echo "MISSING ${db}/${name}"
done
find "$db" -maxdepth 1 -name 'database*mers.kmer_distrib' \
  -type f -size +0 -print | sort -V
```

Sources and sinks must use the same Kraken snapshot, classifier parameters, Bracken distributions, rank, and panel-building procedure. Tables from different database snapshots are not directly mergeable.

For `taxon_filter: prokaryote`, point `taxdump` to a compatible taxonomy snapshot containing `nodes.dmp`.

## 4. Resource tuning

Kraken memory requirements vary substantially by database. The defaults are conservative and explicit:

```yaml
kraken2_cpus: 16
kraken2_memory: 160 GB
kraken2_time: 36h
kraken2_max_forks: 1
host_remove_max_forks: 4
```

Override these after measuring the selected database and executor. Do not increase database or host-removal concurrency until shared-storage and memory behavior have been checked.

Executor-specific staging also belongs in the external site configuration. Large Kraken databases and host indexes should not be copied repeatedly unless the storage design makes that intentional.

## 5. Production configuration

Copy the portable example and keep the run-specific file out of Git:

```bash
cp examples/params.example.yml params.local.yml
```

Use absolute paths for metadata, references, databases, output, work, and cache locations. Run a small pilot containing modern and ancient libraries before expanding to the full panel. Inspect fragment retention, host-removal yield, selected Bracken lengths, classified fractions, and table sample identifiers.

## 6. Resume and provenance

Correct failed inputs or configuration and relaunch the same command with `-resume`. Keep the work directory until the analysis and provenance have been accepted.

For each frozen analysis, archive outside the code repository:

- the Git commit;
- the executed params file;
- an explicit environment package list;
- `containers.lock.tsv` plus the executed SIF checksums for container runs;
- the database snapshot identifiers and checksums;
- Nextflow trace, report, timeline, and DAG files; and
- emitted tool-version and scientific-QC tables.

Example:

```bash
git rev-parse HEAD > /absolute/path/to/provenance/mstflow.git-commit.txt
conda list --prefix /absolute/path/to/environment --explicit \
  > /absolute/path/to/provenance/mstflow.explicit.txt
cp params.local.yml /absolute/path/to/provenance/params.executed.yml
```

Run outputs, downloaded references, databases, work directories, caches, local params, and scheduler logs must remain outside Git.
