process VALIDATE_METADATA {
    label 'process_low'
    tag 'metadata'

    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    path metadata
    path validator_script

    output:
    path "metadata.validated.tsv", emit: metadata
    path "validate_metadata_versions.yml", emit: versions

    script:
    def defaultHostGroup = params.get('default_host_group', 'human')
    def requireExplicitHost = params.get('require_explicit_host_for_nondefault', false).toString().toBoolean()
    def removeHost = params.get('remove_host', true).toString().toBoolean()
    def requireExplicitHostArg = removeHost && requireExplicitHost ? '--require-explicit-host-for-nondefault' : ''
    def configuredHostFasta = params.get('host_fasta', '') ?: ''
    def defaultHostFasta = configuredHostFasta ? file(configuredHostFasta).toString() : ''
    def defaultHostFastaArg = removeHost ? "--default-host-fasta '${defaultHostFasta}'" : ''
    """
    python3 "${validator_script}" \\
        --input "${metadata}" \\
        --output "metadata.validated.tsv" \\
        --default-host-group '${defaultHostGroup}' \\
        ${requireExplicitHostArg} ${defaultHostFastaArg}

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > validate_metadata_versions.yml
    """

    stub:
    def defaultHostGroup = params.get('default_host_group', 'human')
    def requireExplicitHost = params.get('require_explicit_host_for_nondefault', false).toString().toBoolean()
    def removeHost = params.get('remove_host', true).toString().toBoolean()
    def requireExplicitHostArg = removeHost && requireExplicitHost ? '--require-explicit-host-for-nondefault' : ''
    def configuredHostFasta = params.get('host_fasta', '') ?: ''
    def defaultHostFasta = configuredHostFasta ? file(configuredHostFasta).toString() : ''
    def defaultHostFastaArg = removeHost ? "--default-host-fasta '${defaultHostFasta}'" : ''
    """
    python3 "${validator_script}" \\
        --input "${metadata}" \\
        --output "metadata.validated.tsv" \\
        --default-host-group '${defaultHostGroup}' \\
        ${requireExplicitHostArg} ${defaultHostFastaArg}
    printf '"${task.process}":\\n    stub: "true"\\n' > validate_metadata_versions.yml
    """
}
