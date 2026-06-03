process FASTQC {
    label 'process_low'
    tag "${meta.id}"

    publishDir { "${params.outdir}/qc/${task.process.toLowerCase()}" }, mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*_fastqc.zip"), path("*_fastqc.html"), emit: qc
    path "versions.yml", emit: versions

    script:
    def readArgs = reads.collect { it.toString() }.join(' ')
    """
    fastqc \\
        --threads ${task.cpus} \\
        --outdir . \\
        ${readArgs}

    printf '"${task.process}":\\n    fastqc: "%s"\\n' "\$(fastqc --version | sed 's/FastQC v//')" > versions.yml
    """

    stub:
    def outputs = reads.collect { read ->
        def name = read.getBaseName().replaceFirst(/\\.f(ast)?q$/, '')
        """
        touch ${name}_fastqc.zip
        printf '<!doctype html><html><body>stub fastqc ${name}</body></html>\\n' > ${name}_fastqc.html
        """
    }.join('\n')
    """
    ${outputs}
    printf '"${task.process}":\\n    stub: "true"\\n' > versions.yml
    """
}
