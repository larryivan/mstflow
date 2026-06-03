# Server Setup

Install Nextflow and make `conda`/`mamba` available on `PATH`.

```bash
nextflow -version
mamba --version
```

Prepare references:

```bash
# Example only. Use your real CHM13 or host reference path.
bwa index /path/to/chm13.fa
samtools faidx /path/to/chm13.fa
```

Prepare a Kraken2 database and set the path in your params file:

```yaml
host_fasta: "/path/to/chm13.fa"
kraken2_db: "/path/to/kraken2_db"
run_kraken2: true
run_bracken: true
bracken_db: "/path/to/kraken2_db"
bracken_read_len: 100
```

If you want prokaryote-only source attribution, download an NCBI taxdump and set:

```yaml
taxon_filter: prokaryote
taxdump: "/path/to/taxdump"
```

For SourceTracker2, make sure the metadata contains `source_sink` and `model_source_env`, then enable:

```yaml
run_sourcetracker2: true
sourcetracker2_source_rarefaction_depth: 1000
sourcetracker2_sink_rarefaction_depth: 1000
```

Run:

```bash
nextflow run . -profile conda -params-file params.yml -resume
```

The SourceTracker2 process uses `envs/sourcetracker2.yml` with Python 3.8, separate from the main mstflow environment. This avoids the Python `open(..., 'U')` compatibility issue in SourceTracker2.
