/*
 * Score every method on the held-out sources of the cross-validation folds.
 */
process EVALUATE_ATTRIBUTION {
    label 'process_low'
    tag "${params.attribution_cv_group}"
    container "${params.container_python_tools}"
    publishDir "${params.outdir}/attribution/${params.attribution_label}/evaluation", mode: 'copy'

    input:
    path results, stageAs: 'results/*'
    path maps, stageAs: 'maps/*'
    path folds
    path evaluate_script
    path attribution_module

    output:
    path "evaluation.tsv", emit: evaluation
    path "evaluation_by_env.tsv", emit: by_env
    path "predictions.tsv", emit: predictions
    path "evaluation.png", emit: plot

    script:
    """
    export MPLCONFIGDIR=\$PWD/mpl
    python3 ${evaluate_script} \\
        --folds ${folds} \\
        --maps maps/* \\
        --results results/* \\
        --max-unknown ${params.attribution_max_unknown} \\
        --min-top-known ${params.attribution_min_top_known} \\
        --min-margin ${params.attribution_min_margin}
    """

    stub:
    """
    python3 ${evaluate_script} --no-plot \\
        --folds ${folds} \\
        --maps maps/* \\
        --results results/* \\
        --max-unknown ${params.attribution_max_unknown} \\
        --min-top-known ${params.attribution_min_top_known} \\
        --min-margin ${params.attribution_min_margin}
    touch evaluation.png
    """
}
