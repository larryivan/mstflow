/*
 * SourceTracker2 (Knights et al. 2011): a Bayesian mixture of the source environments
 * and an unknown source, fitted by Gibbs sampling. Sources of one environment are
 * pooled, and sources and sinks are rarefied to fixed depths.
 */
process SOURCETRACKER2 {
    label 'process_medium'
    tag "${meta.analysis}"
    container "${params.container_sourcetracker2}"
    publishDir "${params.outdir}/attribution/${params.attribution_label}/sourcetracker2", mode: 'copy',
        saveAs: { name -> meta.analysis == 'main' ? name : null }

    input:
    tuple val(meta), path(counts), path(mapping)
    path seeded_runner
    path results_script
    path attribution_module

    output:
    tuple val(meta), path("${meta.analysis}.sourcetracker2.attribution.tsv"), emit: results
    path "${meta.analysis}.sourcetracker2.versions.yml", emit: versions

    script:
    def seed = params.attribution_seed as int
    """
    export PYTHONHASHSEED=${seed}
    biom convert -i ${counts} -o feature_table.biom --to-hdf5 --table-type 'OTU table'

    python3 ${seeded_runner} --seed ${seed} gibbs \\
        -i feature_table.biom \\
        -m ${mapping} \\
        -o sourcetracker2.out \\
        --source_sink_column SourceSink \\
        --source_column_value source \\
        --sink_column_value sink \\
        --source_category_column Env \\
        --source_rarefaction_depth ${params.sourcetracker2_source_rarefaction_depth} \\
        --sink_rarefaction_depth ${params.sourcetracker2_sink_rarefaction_depth} \\
        --alpha1 ${params.sourcetracker2_alpha1} \\
        --alpha2 ${params.sourcetracker2_alpha2} \\
        --beta ${params.sourcetracker2_beta} \\
        --restarts ${params.sourcetracker2_restarts} \\
        --draws_per_restart ${params.sourcetracker2_draws_per_restart} \\
        --burnin ${params.sourcetracker2_burnin} \\
        --delay ${params.sourcetracker2_delay} \\
        --jobs 1

    # SourceTracker2 2.0.1 misspells some output names in some code paths.
    means=\$(ls sourcetracker2.out/mixing_propor*.txt | grep -v std | head -n 1)
    stds=\$(ls sourcetracker2.out/mixing_propor*std*.txt | head -n 1)
    python3 ${results_script} \\
        --proportions "\$means" \\
        --stds "\$stds" \\
        --map ${mapping} \\
        --analysis ${meta.analysis} \\
        --output ${meta.analysis}.sourcetracker2.attribution.tsv

    printf '"${task.process}":\\n    sourcetracker2: "%s"\\n' "\$(sourcetracker2 --version 2>&1 | head -n 1)" \\
        > ${meta.analysis}.sourcetracker2.versions.yml
    """

    stub:
    """
    python3 - <<'PY'
from mstflow_attribution import UNKNOWN, read_map, roles, write_results
envs, sinks = roles(read_map("${mapping}"))
share = 1.0 / (len(envs) + 1)
write_results("${meta.analysis}.sourcetracker2.attribution.tsv", "${meta.analysis}", "sourcetracker2", "proportion",
              {sink: {source: (share, 0.0) for source in [*envs, UNKNOWN]} for sink in sinks})
PY
    printf '"${task.process}":\\n    sourcetracker2: "stub"\\n' > ${meta.analysis}.sourcetracker2.versions.yml
    """
}
