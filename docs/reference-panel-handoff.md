# Reference-panel container quick start

The handoff contains:

- `download_samples.csv`: public-data manifest and checksums;
- `download_samples.aria2`: standard aria2 input file;
- `mstflow_input.tsv`: pipeline input; and
- `params.yml`: Kraken2/Bracken settings.

The supplied configuration runs preprocessing and Kraken2/Bracken. Host
removal and SourceTracker2 are disabled for this first pass.

## 1. Download the samples

Requirements: `aria2c`, SRA Toolkit, and `pigz`.

Put the four files in one empty run directory, then run:

```bash
cd /absolute/path/to/reference-panel-run

aria2c -c -j 16 -x 4 -s 4 \
  --file-allocation=none \
  --check-integrity=true \
  --auto-file-renaming=false \
  --input-file=download_samples.aria2
```

The aria2 input file contains every destination and the repository MD5 for all
direct FASTQs. The command is resumable; rerun it unchanged after an
interrupted transfer. The SRA archives are checked with `vdb-validate` below.

Eleven rows are SRA archives. Convert them using the target FASTQ paths in the
last two CSV columns:

```bash
mkdir -p .sra_conversion_tmp

tail -n +2 download_samples.csv |
while IFS=, read -r id cohort kind project run url md5 bytes sra fastq_1 fastq_2
do
  test "$kind" = sra || continue
  test -s "$fastq_1" && test -s "$fastq_2" && continue

  vdb-validate "$sra"
  tmp="$(mktemp -d ".sra_conversion_tmp/${run}.XXXXXX")"
  fasterq-dump "$sra" --split-files --skip-technical \
    --threads 8 --outdir "$tmp"
  pigz -n -p 8 "$tmp/${run}_1.fastq" "$tmp/${run}_2.fastq"

  mkdir -p "$(dirname "$fastq_1")" "$(dirname "$fastq_2")"
  mv "$tmp/${run}_1.fastq.gz" "$fastq_1"
  mv "$tmp/${run}_2.fastq.gz" "$fastq_2"
  rmdir "$tmp"
done

find data -type f -name '*.fastq.gz' -exec gzip -t {} +
```

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

Cluster run, with the scheduler set in your own `site.config`
(see [usage](usage.md#on-a-cluster)):

```bash
nextflow -c "$PWD/site.config" run "$PIPELINE" \
  -profile singularity \
  -params-file "$PWD/params.yml" \
  --input "$PWD/mstflow_input.tsv" \
  -resume
```

Keep account, partition, work directory, and other cluster-specific settings
in `site.config`. Replace `singularity` with `apptainer` when appropriate. For
offline compute nodes and pre-staged SIF files, see
[usage](usage.md#on-a-cluster).

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
