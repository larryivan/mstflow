process PLOT_SOURCETRACKER2 {
    label 'process_low'
    tag 'sourcetracker2_plot'

    publishDir "${params.outdir}/source_attribution/sourcetracker2", mode: 'copy'

    input:
    path proportions
    path mapping

    output:
    path "sourcetracker2_assignments.tsv", emit: assignments
    path "sourcetracker2_stacked_bar.png", emit: plot
    path "versions.yml", emit: versions

    script:
    """
    python3 ${projectDir}/bin/plot_sourcetracker2.py \\
        --proportions ${proportions} \\
        --mapping ${mapping} \\
        --assignments sourcetracker2_assignments.tsv \\
        --plot sourcetracker2_stacked_bar.png \\
        --title '${params.multiqc_title} SourceTracker2'

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > versions.yml
    """

    stub:
    """
    printf 'SampleID\\texpected_env\\ttop_known_env\\ttop_known_prop\\ttop_all_env\\ttop_all_prop\\ttop_known_correct\\tunknown_prop\\n' > sourcetracker2_assignments.tsv
    printf '\\x89PNG\\r\\n\\x1a\\n' > sourcetracker2_stacked_bar.png
    printf '"${task.process}":\\n    stub: "true"\\n' > versions.yml
    """
}
