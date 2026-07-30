# Reference-panel container quick start

The handoff contains:

- `download_samples.csv`: public-data URLs and checksums;
- `download_aria2c.py`: resumable downloader and verifier;
- `mstflow_input.tsv`: pipeline input; and
- `params.yml`: Kraken2/Bracken settings.

The supplied configuration runs preprocessing and Kraken2/Bracken. Host
removal and SourceTracker2 are disabled for this first pass.

## 1. Download the samples

Requirements: Python 3.9+, `aria2c`, SRA Toolkit, and `pigz`.

Put the four files in one empty run directory, then run:

```bash
cd /absolute/path/to/reference-panel-run

python3 ./download_aria2c.py ./download_samples.csv \
  --root "$PWD" \
  --jobs 16 \
  --connections-per-file 1 \
  --sra-jobs 2 \
  --sra-threads 8
```

The script resumes partial transfers, verifies sizes and MD5 checksums,
converts the SRA rows, and validates the FASTQs. A successful run ends with:

```text
PASS: all requested downloads are complete and verified
```

It is safe to rerun the same command.

## 2. Configure the database

Edit the three placeholder paths in `params.yml`:

```yaml
kraken2_db: /absolute/path/to/kraken_bracken_database
bracken_db: /absolute/path/to/kraken_bracken_database
taxdump: /absolute/path/to/nodes.dmp
```

The two database paths must identify the same snapshot. It must contain the
three Kraken2 `.k2d` files and at least one Bracken
`database<N>mers.kmer_distrib`. The database is not included in the handoff.

## 3. Run MSTFlow

Requirements: Nextflow 25.10.4+ and Singularity or Apptainer.

```bash
git clone https://github.com/larryivan/mstflow.git
PIPELINE="$(cd mstflow && pwd)"
cd /absolute/path/to/reference-panel-run
```

Local Singularity run:

```bash
nextflow run "$PIPELINE" \
  -profile singularity \
  -params-file "$PWD/params.yml" \
  --input "$PWD/mstflow_input.tsv" \
  -resume
```

Slurm run:

```bash
nextflow -c "$PWD/site.config" run "$PIPELINE" \
  -profile slurm_singularity \
  -params-file "$PWD/params.yml" \
  --input "$PWD/mstflow_input.tsv" \
  -resume
```

Keep account, partition, work directory, and other cluster-specific settings
in `site.config`. Replace `singularity` with `apptainer` in the profile name
when appropriate. For offline compute nodes and pre-staged SIF files, see
[`server_setup.md`](server_setup.md#pre-staged-sif-images).

An optional container smoke test is:

```bash
cd "$PIPELINE"
python3 bin/create_tiny_test_data.py
nextflow run . -profile singularity,test \
  -params-file examples/params.tiny.stub.yml -stub-run -resume
```

## 4. Results to return

Please return these directories under `results/reference_panel/`:

```text
bracken/
kraken2/
tables/
qc/
pipeline_info/
```

Please also retain `params.yml`, `mstflow_input.tsv`, the database snapshot
identifier, and the pipeline commit:

```bash
git -C "$PIPELINE" rev-parse HEAD > mstflow.git-commit.txt
```

Do not combine source and sink profiles generated with different databases or
classifier settings.
