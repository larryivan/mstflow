process VALIDATE_METADATA {
    label 'process_low'
    tag 'metadata'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path metadata

    output:
    path "metadata.validated.tsv", emit: metadata
    path "versions.yml", emit: versions

    script:
    """
    python3 ${projectDir}/bin/validate_metadata.py \\
        --input ${metadata} \\
        --output metadata.validated.tsv

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > versions.yml
    """

    stub:
    """
    python3 ${projectDir}/bin/validate_metadata.py \\
        --input ${metadata} \\
        --output metadata.validated.tsv
    printf '"${task.process}":\\n    stub: "true"\\n' > versions.yml
    """
}
