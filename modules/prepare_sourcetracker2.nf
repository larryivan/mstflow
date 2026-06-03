process PREPARE_SOURCETRACKER2_INPUTS {
    label 'process_low'
    tag 'sourcetracker2_inputs'

    publishDir "${params.outdir}/source_attribution/sourcetracker2", mode: 'copy'

    input:
    path counts
    path metadata

    output:
    path "sourcetracker2_counts.tsv", emit: counts
    path "sourcetracker2_map.tsv", emit: map
    path "sourcetracker2_input_summary.tsv", emit: summary
    path "versions.yml", emit: versions

    script:
    """
    python3 ${projectDir}/bin/prepare_sourcetracker2_inputs.py \\
        --counts ${counts} \\
        --metadata ${metadata} \\
        --output-counts sourcetracker2_counts.tsv \\
        --output-map sourcetracker2_map.tsv \\
        --summary sourcetracker2_input_summary.tsv

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > versions.yml
    """
}
