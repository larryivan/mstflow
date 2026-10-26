# Troubleshooting

Each heading is the message you see; below it is what to do. After a fix, rerun with `-resume`. `.nextflow.log` in the launch directory has the full error, and questions are welcome at <https://github.com/larryivan/mstflow/issues>.

**`-profile … needs Docker, which is not installed` or `No software environment was selected`**
The tools have nowhere to run from. The message lists Docker, Apptainer, Singularity and Conda with what was found, and the `-profile` to use. If none is installed, install [Miniforge](https://github.com/conda-forge/miniforge) and use `-profile conda`. Docker must be running: open Docker Desktop. If only your compute nodes have the tool, add `--check_software false`.

**`Missing required parameter(s): input`**
Give the samplesheet with `--input`, or try `-profile test,docker`.

**`Sample '…': fastq_1 '…' was not found`**
Relative paths are read next to the samplesheet, then from the launch directory. Use absolute paths if the reads are elsewhere.

**`kraken2_db is missing required file: …/hash.k2d`**
`--kraken2_db` must be the directory holding `hash.k2d`, `opts.k2d` and `taxo.k2d`, or a `.tar.gz` of it, not a parent directory.

**`… does not contain a Kraken2 database`**
The archive has no `hash.k2d` at its top level or one directory down. Check that it is a Kraken2 index, not a library.

**`Host FASTA is missing BWA index: ….bwt`**
Run `bwa index host.fa` and keep the index beside the FASTA.

**A Kraken2 task fails with exit status 137 or 140**
It ran out of memory. Raise `--kraken2_memory` to the size of `hash.k2d` plus a few GB; MSTFlow already retries once with double.

**Every sample is flagged `few_reads` or `few_genera`**
There is little classifiable DNA, as in capture libraries and poorly preserved material. The tables are still written; read them with care. The thresholds are `--qc_min_classified_reads` and `--qc_min_taxa`.

**Docker warns that the image platform does not match the host**
Expected on Apple Silicon: the images are `linux/amd64` and run under emulation, only slower.

**SourceTracker2 assigns nearly every sink to Unknown**
With few taxa per sink, a small prior lets the unknown source absorb sinks that match a source. Keep `--sourcetracker2_alpha2` at 0.1 and check cross-validation ([details](attribution.md#known-behaviour-of-each-method)).

**`SourcePredict needs at least two source libraries per environment`**
Add sources, merge rare environments, or drop SourcePredict from `--attribution_methods`. In cross-validation folds, an environment left with one library is skipped instead.

**`cross-validation needs sources from at least two 'project' groups`**
Cross-validation holds out one `--attribution_cv_group` value at a time. Fill the `project` column, or name another column that groups your sources.

**`feast and fastst need MSTFlow's attribution container or a build from source`**
They are on GitHub, not in Bioconda. Use `-profile docker`, `singularity` or `apptainer` (the attribution image) or `-profile conda` (built when the run starts). With your own R and a C++ compiler, set `--r_methods_from_source true`.

**SourceID-NMF, FEAST or FastST cannot be fetched without internet access**
Clone <https://github.com/ZiyiHuang0708/SourceID-NMF> and give its directory with `--sourceid_nmf_code`. For FEAST and FastST, download the sources and give them with `--feast_code`, `--fastst_code` and `--npreg_code`.
