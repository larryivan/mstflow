# Troubleshooting

**`Missing required parameter(s): input`**
Give the samplesheet with `--input`. To try the workflow without your own data, run
`nextflow run larryivan/mstflow -profile test,docker`.

**`Sample '…': fastq_1 '…' was not found`**
Relative paths are looked up next to the samplesheet first, then in the directory you
launched from. Use absolute paths if the reads live elsewhere.

**`kraken2_db is missing required file: …/hash.k2d`**
`--kraken2_db` must be the directory that contains `hash.k2d`, `opts.k2d` and
`taxo.k2d`, or a `.tar.gz` of it; not a parent directory.

**`… does not contain a Kraken2 database`**
The archive given to `--kraken2_db` holds no `hash.k2d` at its top level or one
directory down. Check that it is a Kraken2 index and not, for example, a library.

**`The kraken2 database has no Bracken database<N>mers.kmer_distrib files`**
The database has no Bracken read-length models. Prebuilt indexes include them; for a
custom database run `bracken-build -d <db> -l <read length>` once per length, or set
`--run_bracken false` to use the classifier's own counts.

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

**SourceTracker2 assigns nearly every sink to Unknown**
With few taxa per sink, SourceTracker2's default prior lets its unknown source absorb
sinks that clearly match a source. Try `--sourcetracker2_alpha2 0.1` and compare the
cross-validation scores; see [source attribution](attribution.md#known-behaviour-of-each-method).

**`SourcePredict needs at least two source libraries per environment`**
SourcePredict learns from individual source libraries. Add sources, merge rare
environments, or leave SourcePredict out of `--attribution_methods`. In
cross-validation folds an environment left with one library is skipped instead.

**SourceID-NMF cannot download its code**
Its files are fetched from GitHub when the run starts. Without internet access, clone
<https://github.com/ZiyiHuang0708/SourceID-NMF> elsewhere, copy it over, and give the
directory with `--sourceid_nmf_code`.

**`feast and fastst run only from MSTFlow's attribution container`**
FEAST and FastST are installed from GitHub into MSTFlow's attribution image, so they
need `-profile docker`, `singularity` or `apptainer`; Conda cannot provide them.

**`cross-validation needs sources from at least two 'project' groups`**
Cross-validation holds out one value of `--attribution_cv_group` at a time. Fill the
`project` column, or name another column that separates your sources into groups.

**Changing an option and rerunning repeats finished steps**
Add `-resume` so that only the steps affected by the change run again.

Anything else: run with `-resume` after fixing the cause, and check `.nextflow.log` in the
launch directory for the full error. Issues and questions are welcome at
<https://github.com/larryivan/mstflow/issues>.
