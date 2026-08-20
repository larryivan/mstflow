# Troubleshooting

**`Missing required parameter(s): input`**
Give the samplesheet with `--input`. To try the workflow without your own data, run
`nextflow run larryivan/mstflow -profile test,docker`.

**`Sample '…': fastq_1 '…' was not found`**
Relative paths are looked up next to the samplesheet first, then in the directory you
launched from. Use absolute paths if the reads live elsewhere.

**`kraken2_db is missing required file: …/hash.k2d`**
`--kraken2_db` must be the directory that contains `hash.k2d`, `opts.k2d` and
`taxo.k2d`, not a parent directory or the archive.

**`bracken_db contains no database<N>mers.kmer_distrib files`**
The database has no Bracken read-length models. Prebuilt indexes include them; for a
custom database run `bracken-build -d <db> -l <read length>` once per length.

**`Host FASTA is missing BWA index: ….bwt`**
Run `bwa index host.fa` and keep the index files next to the FASTA.

**A Kraken2 task fails with exit status 137 or 140**
It ran out of memory. Raise `--kraken2_memory` to the size of `hash.k2d` plus a few GB.
MSTFlow already retries once with twice the memory before giving up.

**Every sample is flagged `few_reads` or `few_genera`**
The libraries hold little classifiable DNA, which is common for capture libraries and
poorly preserved material. The tables are still written; interpret them with care. The
thresholds are `--qc_min_classified_reads` and `--qc_min_taxa`.

**Docker warns that the image platform does not match the host**
Expected on Apple Silicon and other ARM machines: the images are `linux/amd64` and run
under emulation. It works, only more slowly.

**Changing an option and rerunning repeats finished steps**
Add `-resume` so that only the steps affected by the change run again.

Anything else: run with `-resume` after fixing the cause, and check `.nextflow.log` in the
launch directory for the full error. Issues and questions are welcome at
<https://github.com/larryivan/mstflow/issues>.
