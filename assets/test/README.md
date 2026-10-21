# Test data

`-profile test` runs the whole workflow on short excerpts of 15 published sequencing
libraries. The files, about 1.8 GB, are attached to the
[`test-data-v1` release](https://github.com/larryivan/mstflow/releases/tag/test-data-v1)
and downloaded on the first run. This directory holds the recipe they were built from.

| File | Contents |
|---|---|
| [`libraries.tsv`](libraries.tsv) | The 15 libraries, with their ENA run and project accessions |
| [`genomes.tsv`](genomes.tsv) | The 274 reference genomes in the Kraken2/Bracken database |
| [`make_test_data.py`](make_test_data.py) | Downloads the reads, genomes, taxonomy and host genome |
| [`build_test_db.sh`](build_test_db.sh) | Builds the database and collects the release files |

## Libraries

Each library is represented by its first 50,000 reads (read pairs for paired-end runs).

| Sample | Role | Material | Run | Study |
|---|---|---|---|---|
| oral_1 | source: oral | modern dental calculus | ERR3579820 | [PRJEB34569](https://www.ebi.ac.uk/ena/browser/view/PRJEB34569) The evolution and changing ecology of the African hominid oral microbiome |
| oral_2 | source: oral | modern dental calculus | ERR3307051 | [PRJEB31185](https://www.ebi.ac.uk/ena/browser/view/PRJEB31185) Metagenome sequencing of modern dental calculus |
| oral_3 | source: oral | ancient dental calculus | ERR3579737 | [PRJEB34569](https://www.ebi.ac.uk/ena/browser/view/PRJEB34569) |
| gut_1 | source: gut | faeces | SRR5898909 | [PRJNA397112](https://www.ebi.ac.uk/ena/browser/view/PRJNA397112) Indian human gut microbiome |
| gut_2 | source: gut | faeces | SRR5763463 | [PRJNA392180](https://www.ebi.ac.uk/ena/browser/view/PRJNA392180) Seasonal cycling in the gut microbiome of the Hadza hunter-gatherers of Tanzania |
| gut_3 | source: gut | faeces | SRR6451716 | [PRJNA354235](https://www.ebi.ac.uk/ena/browser/view/PRJNA354235) Human gut metagenome and metatranscriptome |
| gut_4 | source: gut | human paleofaeces | ERR3678629 | [PRJEB33577](https://www.ebi.ac.uk/ena/browser/view/PRJEB33577) Comparison of extraction methods for recovering ancient microbial DNA from paleofeces |
| soil_1 | source: soil | topsoil | ERR2603035 | [PRJEB24121](https://www.ebi.ac.uk/ena/browser/view/PRJEB24121) Estonian forest and grassland topsoil samples |
| soil_2 | source: soil | topsoil | ERR2603043 | [PRJEB24121](https://www.ebi.ac.uk/ena/browser/view/PRJEB24121) |
| lab_1 | source: laboratory | extraction blank | ERR10114869 | [PRJEB55583](https://www.ebi.ac.uk/ena/browser/view/PRJEB55583) Reconstruction of Middle and Late Paleolithic metagenomes to access natural products |
| lab_2 | source: laboratory | extraction blank | ERR10114868 | [PRJEB55583](https://www.ebi.ac.uk/ena/browser/view/PRJEB55583) |
| calculus | sink | ancient dental calculus | ERR10114857 | [PRJEB55583](https://www.ebi.ac.uk/ena/browser/view/PRJEB55583) |
| coprolite | sink | human paleofaeces | ERR5766170 | [PRJEB44507](https://www.ebi.ac.uk/ena/browser/view/PRJEB44507) Metagenomics analysis of human paleofeces |
| sediment | sink | paleosoil | ERR13348407 | [PRJEB77295](https://www.ebi.ac.uk/ena/browser/view/PRJEB77295) Reconstruction of the vegetation of Morasko hill during the meteorite impact based on sedimentary ancient DNA from paleosoil |
| bone | sink | human bone | ERR1777568 | [PRJEB18722](https://www.ebi.ac.uk/ena/browser/view/PRJEB18722) Comparing ancient DNA preservation in petrous bone and tooth cementum |

The sources include ancient libraries as well as modern ones, as real source panels for
ancient samples do; every sink's own study is kept out of its environment's sources.
The `project` column of the samplesheet groups the libraries by study, so
cross-validation holds out one study at a time.

The data belong to the original studies; please cite them, through the publications
linked from each ENA project page, if you use these files for anything beyond testing.

## Database

The genomes were chosen so the database covers most of what each library contains.
Each library was classified with Kraken2 against NCBI core_nt, and the species that
together account for 80% of its species-level reads (at most 40 per library) were
kept, one reference genome each: 265 species in total. Nine genus representatives were
added for genera that had no species call. The `selected_for` column of `genomes.tsv`
records the run each genome was chosen for. The genomes were chosen before the final set
of libraries was settled: ERR16255222, an ancient calculus library first considered, was
replaced because its adapters carry inline barcodes, and oral_3, gut_4, the calculus and
the coprolite rely on the genomes chosen for the other libraries of their kind.

Genomes and their NCBI taxonomy were downloaded on 2026-10-21; the release carries
`nodes.dmp`, `names.dmp` and an empty `merged.dmp` cut down to the database's taxa. The database is built
with Kraken2 2.1.6 (k = 35) and Bracken read-length models for 50, 75, 100 and 150 bp.
The host reference is the human mitochondrial genome (NC_012920.1), which is enough to
exercise host removal without a 3 GB download.

## Rebuilding

```bash
python3 make_test_data.py   # downloads into build/, about 1.6 GB
bash build_test_db.sh       # needs Docker and about 4 GB of memory
```

The release files are written to `build/release/`. NCBI updates its taxonomy and genome
records, so a rebuild gives a database that is close to, but not byte-identical with,
the released one.
