# Portable server setup

MSTFlow keeps scientific workflow logic separate from site deployment policy. The repository provides local, Conda, test, and generic Slurm profiles; accounts, queues, module names, storage paths, scratch behavior, and submission limits belong in an external site configuration.

## 1. Software

Install Nextflow 25.10.4 or newer and make Mamba or Conda available on `PATH`:

```bash
nextflow -version
mamba --version
```

The main environment is defined in `envs/mstflow.yml`. SourceTracker2 uses the separate compatibility environment in `envs/sourcetracker2.yml` because version 2.0.1 requires an older Python/NumPy stack.

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

Generic Slurm execution:

```bash
nextflow run . \
  -profile slurm \
  -params-file /absolute/path/to/params.yml \
  -resume
```

The generic Slurm profile deliberately does not choose an account, partition, work directory, cache directory, module system, or scratch policy. Put those values in an external configuration and compose it at launch:

```bash
nextflow -c /absolute/path/to/site.config run . \
  -profile slurm \
  -params-file /absolute/path/to/params.yml \
  -resume
```

Keep site configuration outside this repository when it contains usernames, allocation names, private paths, or infrastructure assumptions.

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
