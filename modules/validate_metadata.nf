include { isOn } from '../subworkflows/utils'

process VALIDATE_METADATA {
    label 'process_low'
    tag 'metadata'
    container "${params.container_python_tools}"

    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path metadata
    path validator_script

    output:
    path "metadata.validated.tsv", emit: metadata
    path "validate_metadata_versions.yml", emit: versions

    script:
    // The host rules matter only when host reads are removed.
    def removeHost = isOn(params.remove_host)
    def hostArgs = [
        removeHost && isOn(params.require_explicit_host_for_nondefault) ? '--require-explicit-host-for-nondefault' : '',
        removeHost && isOn(params.require_explicit_host_action) ? '--require-explicit-host-action' : '',
    ].join(' ')
    """
    python3 "${validator_script}" \\
        --input "${metadata}" \\
        --output "metadata.validated.tsv" \\
        --default-host-group '${params.default_host_group}' \\
        ${hostArgs}

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > validate_metadata_versions.yml
    """
}
