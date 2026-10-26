# Source attribution

Source attribution estimates where a sample's microbes came from: how much of an ancient dental calculus is oral, how much is soil or laboratory contamination. The samples to explain are *sinks*; reference samples of known environment are *sources*.

MSTFlow runs five published methods on the same count table, lines up their results and scores each on your own sources. The methods often disagree, more on real data than in their papers, so agreement between them, and accuracy on held-out sources, are part of the result.

## Samplesheet

| Column | Meaning |
|---|---|
| `source_sink` | `source` or `sink`; empty leaves the library out of attribution |
| `model_source_env` | For sources, the environment: `oral`, `gut`, `soil`, `laboratory` ... |
| `project` | The study a library comes from; cross-validation holds out one study at a time |
| `unit_id` | Optional: runs of one biological sample, merged into one unit |

Sources of one environment are pooled, so each environment needs at least one source that passes the depth checks. SourcePredict learns from each source on its own and needs at least two per environment, and many more to work well.

## Methods

```bash
--attribution_methods sourcetracker2,feast,fastst,sourceid_nmf,sourcepredict
```

| Method | Model | Reports | Licence | Comes from |
|---|---|---|---|---|
| `sourcetracker2` | Bayesian mixture, Gibbs sampling ([Knights et al. 2011](https://doi.org/10.1038/nmeth.1650)) | shares, with standard deviation over draws | BSD-3 | Bioconda |
| `feast` | Expectation-maximisation mixture ([Shenhav et al. 2019](https://doi.org/10.1038/s41592-019-0431-x)) | shares | CC BY-NC-SA 4.0: **non-commercial use only** | GitHub, fixed commit |
| `fastst` | Generalised least squares ([Choi et al. 2025](https://doi.org/10.7717/peerj.20161)); experimental | shares, with standard errors | MIT | GitHub, fixed commit |
| `sourceid_nmf` | Weighted non-negative matrix factorisation ([Huang et al. 2024](https://doi.org/10.1093/bioinformatics/btae227)) | shares | none published | the authors' repository, at run time |
| `sourcepredict` | k-nearest neighbours on a UniFrac embedding ([Borry et al. 2019](https://doi.org/10.21105/joss.01540)) | class probabilities | GPL-3 | Bioconda |

All but SourcePredict estimate the share each source contributed, plus an unknown share. SourcePredict names the most likely single environment; its numbers are probabilities (`quantity` = `probability`).

MSTFlow is MIT-licensed and includes none of the methods' code:

- **FEAST and FastST** are on GitHub only. They come in MSTFlow's attribution image, which carries their licence files, or, with `-profile conda`, are built from source at the same commits (`--feast_code`, `--fastst_code`, `--npreg_code`).
- **SourceID-NMF** has no licence, so each run downloads its five Python files at a fixed commit (`--sourceid_nmf_code`). Offline, clone the repository and point the parameter at the directory.

## How a sink is called

For each method and sink, MSTFlow names the leading known source only when

- the unknown share is below `--attribution_max_unknown` (0.5),
- the leading source reaches `--attribution_min_top_known` (0.5), and
- it leads the runner-up by at least `--attribution_min_margin` (0.05).

Otherwise the call is `ABSTAIN`, with every reason that applied. An abstention is a result: the sources on hand do not explain the sink well enough to name one.

`comparison.tsv` lines up the calls per sink: `agree:<env>`, `disagree`, `one_method` (only one method made a call) or `all_abstain`.

## Cross-validation on your sources

```bash
--attribution_cv true --attribution_cv_group project
```

Each method is run again once per study, with that study's sources held out as sinks of known environment. Whole studies are held out, not single libraries, because libraries from one study share a laboratory and a protocol, and a method that has seen the others is judged too kindly.

`evaluation/evaluation.tsv` reports, per method:

| Column | Meaning |
|---|---|
| `coverage` | share of held-out sources the method made a call for |
| `accuracy` | share of those calls that were right |
| `balanced_accuracy` | the same, averaged over environments, so a small environment that is always wrong stays visible |
| `correct_overall` | right calls over all held-out sources, abstentions counting as wrong; the number to compare methods by |
| `abstained_on_unseen_env` | held-out sources whose environment had no other source in the fold, where abstaining is right |

`predictions.tsv` lists every held-out source with each method's call. A method that scores poorly here should not be trusted on your real sinks, whatever its paper reports.

Cross-validation tests the call. Whether the shares are right shows only on sinks of known composition. On the test mixtures SourceTracker2 comes within 0.06, FEAST and SourceID-NMF put much of each mixture in their unknown source, and FastST gives the dominant source far too much ([validation](validation.md#shares-on-sinks-of-known-composition)). Read SourceTracker2's shares as proportions, and the other methods' calls as a check on its call.

## Known behaviour of each method

- **SourceTracker2.** `--sourcetracker2_alpha2` is the prior on the unknown source. MSTFlow uses 0.1, the default of the original SourceTracker; the 2.0.1 release on Bioconda uses 0.001, with which the unknown source can absorb a whole sink that clearly matches a source. This bites when sinks hold few taxa, as with small databases or shallow libraries; cross-validation shows it ([measured](validation.md#choices-the-tests-decided)).
- **FEAST** uses one rarefaction depth for all sinks (`--feast_coverage`, by default the smallest sink or pooled source), since its own default would take it from the first sink alone.
- **FastST** fits only the taxa present in a sink, so a source can be scaled up freely if it matches on those taxa. Its calls were right on every held-out source but its shares of the test mixtures were far off; treat it as experimental.
- **SourceID-NMF** starts its unknown profile from what the sources do not explain, and on sinks far from every source it tends to give everything to it.
- **SourcePredict** was built for large panels: its automatic neighbours and t-SNE embedding need about a hundred sources. With fewer, MSTFlow uses five neighbours and MDS. Its default normalisation, GMPR, fails on sinks that share few taxa with the sources, so MSTFlow subsamples instead (`--sourcepredict_normalization`).

## Outputs

In `attribution/<label>/` (`--attribution_label`, default `default`):

| Path | Contents |
|---|---|
| `inputs/` | the table and map every method read; `attribution_sample_qc.tsv` says which samples were left out and why |
| `<method>/` | each method's results in the shared format |
| `proportions.tsv` | every estimate: method, sink, source, estimate, uncertainty, quantity |
| `assignments.tsv` | per method and sink: the call, the leading sources, the unknown share |
| `comparison.tsv` | per sink: every method's call and whether they agree |
| `attribution.png` | one panel per method, one bar per sink |
| `evaluation/` | with `--attribution_cv`: folds, per-method scores, predictions and a figure |
