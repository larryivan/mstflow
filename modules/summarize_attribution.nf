/*
 * Put every method's estimates for the real sinks side by side: one table of
 * estimates, one call per method and sink, a comparison of the calls, and a figure.
 */
process SUMMARIZE_ATTRIBUTION {
    label 'process_low'
    tag "${params.attribution_label}"
    container "${params.container_python_tools}"
    publishDir "${params.outdir}/attribution/${params.attribution_label}", mode: 'copy'

    input:
    path results
    path mapping
    path summarize_script
    path attribution_module

    output:
    path "proportions.tsv", emit: proportions
    path "assignments.tsv", emit: assignments
    path "comparison.tsv", emit: comparison
    path "attribution.png", emit: plot

    script:
    """
    export MPLCONFIGDIR=\$PWD/mpl
    python3 ${summarize_script} \\
        --map ${mapping} \\
        --max-unknown ${params.attribution_max_unknown} \\
        --min-top-known ${params.attribution_min_top_known} \\
        --min-margin ${params.attribution_min_margin} \\
        --title '${params.multiqc_title}: source attribution' \\
        ${results}
    """

    stub:
    """
    python3 ${summarize_script} --no-plot \\
        --map ${mapping} \\
        --max-unknown ${params.attribution_max_unknown} \\
        --min-top-known ${params.attribution_min_top_known} \\
        --min-margin ${params.attribution_min_margin} \\
        --title '${params.multiqc_title}: source attribution' \\
        ${results}
    touch attribution.png
    """
}
