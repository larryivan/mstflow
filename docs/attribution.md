# Source attribution

Source attribution estimates where the microbes in a sample came from: how much of an
ancient dental calculus is oral, how much is soil or laboratory contamination. The
samples to explain are *sinks*; the reference samples whose environments are known are
*sources*.

MSTFlow runs any of five published methods on the same count table, puts their results
side by side, and can score each of them on your own sources. The methods disagree
often, on real data more than in their papers, so treat agreement between methods, and
their accuracy on held-out sources, as part of the result.

## Samplesheet

| Column | Meaning |
|---|---|
| `source_sink` | `source` or `sink` |
| `model_source_env` | For sources, the environment, such as `oral`, `gut`, `soil`, `laboratory` |
| `project` | The study a library comes from; cross-validation holds out one study at a time |
| `unit_id` | Optional: runs of one biological sample, merged into one unit |

A library with `source_sink` left empty is profiled but left out of attribution.

Sources of one environment are pooled for the mixture models, as SourceTracker2 does,
so each environment needs at least one source that passes the depth checks.
SourcePredict learns from every source library separately and needs at least two per
environment, and many more to work well.

## Methods

```bash
--attribution_methods sourcetracker2,feast,fastst,sourceid_nmf,sourcepredict
```

| Method | Model | Reports | Uncertainty | Licence | Installed from |
|---|---|---|---|---|---|
| `sourcetracker2` | Bayesian mixture, Gibbs sampling ([Knights et al. 2011](https://doi.org/10.1038/nmeth.1650)) | shares of the sink | standard deviation over draws | BSD-3 | Bioconda |
| `feast` | mixture fitted by expectation-maximisation ([Shenhav et al. 2019](https://doi.org/10.1038/s41592-019-0431-x)) | shares of the sink | none | CC BY-NC-SA 4.0: **non-commercial use only** | GitHub, fixed commit, in the MSTFlow image |
| `fastst` | mixture fitted by generalised least squares ([Choi et al. 2025](https://doi.org/10.7717/peerj.20161)); experimental | shares of the sink | standard errors | MIT | GitHub, fixed commit, in the MSTFlow image |
| `sourceid_nmf` | weighted non-negative matrix factorisation ([Huang et al. 2024](https://doi.org/10.1093/bioinformatics/btae227)) | shares of the sink | none | none published | fetched from the authors' repository at run time |
| `sourcepredict` | k-nearest-neighbour classifier on a UniFrac embedding ([Borry et al. 2019](https://doi.org/10.21105/joss.01540)) | class probabilities | none | GPL-3 | Bioconda |

All but SourcePredict estimate what share of the sink each source contributed, plus an
unknown share. SourcePredict instead asks which single environment the sink most likely
belongs to; its numbers are probabilities, and the outputs mark them as such
(`quantity` = `probability`).

### Licences

MSTFlow itself is MIT-licensed and does not include the methods' code.

- **FEAST** may be used for non-commercial purposes only. It runs from MSTFlow's
  attribution image, which carries its licence file; with `-profile conda` it is built
  from source at the same commit (`--feast_code`).
- **SourceID-NMF** has no licence in its repository, so MSTFlow does not redistribute
  it: each run downloads its five Python files from the authors' repository at a fixed
  commit (`--sourceid_nmf_code`). On machines without internet access, clone the
  repository and point `--sourceid_nmf_code` at that directory.

## How a sink is called

For each method and sink, MSTFlow names the leading known source only when

- the unknown share is below `--attribution_max_unknown` (0.5),
- the leading source reaches `--attribution_min_top_known` (0.5), and
- it leads the runner-up by at least `--attribution_min_margin` (0.05).

Otherwise the call is `ABSTAIN`, with every reason that applied. An abstention is a
result: the sources on hand do not explain the sink well enough to name one.

`comparison.tsv` lines up the calls: `agree:<env>` when every method that made a call
named the same environment, `disagree` when they named different ones, `one_method`
when only one method made a call, and `all_abstain`.

## Cross-validation on your sources

```bash
--attribution_cv true --attribution_cv_group project
```

Each method is run again with the sources of one study held out and treated as sinks
whose environment is known, once per study. Holding out whole studies, not single
libraries, matters: libraries from one study share a laboratory, a protocol and often a
site, and a method that has seen a study's other libraries is judged too kindly.

`evaluation/evaluation.tsv` then reports, per method:

| Column | Meaning |
|---|---|
| `coverage` | share of held-out sources the method made a call for |
| `accuracy` | share of those calls that were right |
| `balanced_accuracy` | the same, averaged over environments, so one large environment cannot hide a small one that is always wrong |
| `correct_overall` | right calls over all held-out sources, abstentions counting as not right; the single number to compare methods by |
| `abstained_on_unseen_env` | held-out sources whose environment had no other source in the fold, where abstaining is the right answer |

`predictions.tsv` lists every held-out source with each method's call.

A method that scores poorly here should not be trusted on the real sinks, whatever its
paper reports.

## How well the shares hold

Cross-validation asks whether a method names the right environment for a held-out
source; whether its shares are right shows only on sinks of known composition. On the
test data's mixtures SourceTracker2 comes within 0.06 of the composition, while FEAST and
SourceID-NMF put much of each mixture in their unknown source and FastST gives the
dominant source far too much ([validation](validation.md#shares-on-sinks-of-known-composition)).
Read SourceTracker2's shares as proportions, and the other methods' calls as a check on
its call.

## Known behaviour of each method

**SourceTracker2.** `--sourcetracker2_alpha2` is the prior on the unknown source, as a
fraction of each sink's counts. MSTFlow uses 0.1, the default of the original
SourceTracker and of SourceTracker2's source code; the 2.0.1 release on Bioconda uses
0.001. With 0.001 the unknown source can take on the composition of the sink itself:
when a sink holds only a few dozen taxa and part of it is poorly represented in the
sources, the unknown share grows from that part to absorb the whole sink, including taxa
that match a source closely. In MSTFlow's test data, 0.001 put 93% of the ancient
calculus, 60% of the paleosoil and even 96% of the modern faeces into the unknown source,
and SourceTracker2 named none of the held-out sources; with 0.1, 94% of the calculus went
to the oral sources and 94% of the paleosoil to soil, and it was right on 86% of
held-out sources. Profiles with few taxa, from
small databases or shallow libraries, are where this matters; check it
with cross-validation.

**FEAST** chooses one rarefaction depth for all sinks (`--feast_coverage`; by default
the smallest sink or pooled source), since its own default would take it from the
first sink alone.

**FastST** fits only the taxa present in a sink, so a source can be scaled up freely
as long as it matches the sink on those taxa. Its calls were right on every held-out
source, but its shares of the test mixtures were far off. Treat it as experimental.

**SourceID-NMF** initialises its unknown profile from what the sources do not explain;
on sinks far from every source it tends to assign everything to the unknown source.

**SourcePredict** was built for large reference panels. Its automatic choice of
neighbours and its t-SNE embedding both need around a hundred sources; with fewer,
MSTFlow uses five neighbours and an MDS embedding. Its default count normalisation,
GMPR, fails on sinks that share few taxa with the sources, so MSTFlow subsamples
instead (`--sourcepredict_normalization`).

## Outputs

`attribution/<label>/` (`--attribution_label`, default `default`):

| Path | Contents |
|---|---|
| `inputs/` | the count table and map every method read, and in `attribution_sample_qc.tsv` which samples were left out and why |
| `<method>/` | each method's results in the shared table format |
| `proportions.tsv` | every estimate: method, sink, source, estimate, uncertainty, quantity |
| `assignments.tsv` | per method and sink: the call, the leading sources and the unknown share |
| `comparison.tsv` | per sink: every method's call and whether they agree |
| `attribution.png` | one panel per method, one bar per sink |
| `evaluation/` | with `--attribution_cv`: the folds, per-method scores, predictions and a figure |
