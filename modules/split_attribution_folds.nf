include { splitList } from '../subworkflows/utils'

/*
 * Hold out the sources of one group (by default one study) at a time and turn them
 * into sinks of known origin, so every method can be scored on data it has not seen.
 */
process SPLIT_ATTRIBUTION_FOLDS {
    label 'process_low'
    tag "${params.attribution_cv_group}"
    container "${params.container_python_tools}"
    publishDir "${params.outdir}/attribution/${params.attribution_label}/evaluation", mode: 'copy',
        saveAs: { name -> name == 'folds.tsv' ? name : null }

    input:
    path counts
    path mapping
    path split_script
    path attribution_module

    output:
    path "folds/*.counts.tsv", emit: counts
    path "folds/*.map.tsv", emit: maps
    path "folds.tsv", emit: summary

    script:
    def withSourceTracker = 'sourcetracker2' in splitList(params.attribution_methods)
    def minSink = withSourceTracker ? params.sourcetracker2_sink_rarefaction_depth : 0
    """
    python3 ${split_script} \\
        --counts ${counts} \\
        --map ${mapping} \\
        --group-column '${params.attribution_cv_group}' \\
        --min-sink-count ${minSink} \\
        --out-dir folds \\
        --summary folds.tsv
    """
}
