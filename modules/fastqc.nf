process FASTQC {
    label 'process_low'
    tag "${meta.id}:${meta.qc_stage}"
    container "${params.container_fastqc}"

    publishDir { "${params.outdir}/qc/${task.process.toLowerCase()}" }, mode: 'copy'

    input:
    tuple val(meta), path(reads, stageAs: 'reads??/*')

    output:
    tuple val(meta), path("${meta.id}.${meta.qc_stage}.R*_fastqc.zip"), path("${meta.id}.${meta.qc_stage}.R*_fastqc.html"), emit: qc
    path "${meta.id}.${meta.qc_stage}.fastqc.versions.yml", emit: versions

    script:
    def id = meta.id
    def stage = meta.qc_stage
    def readList = reads instanceof List ? reads : [reads]
    def commands = readList.withIndex().collect { read, index ->
        def readNumber = index + 1
        def outputBase = "${id}.${stage}.R${readNumber}_fastqc"
        """
        mkdir -p fastqc_${readNumber}
        fastqc --threads ${task.cpus} --outdir "fastqc_${readNumber}" "${read}"
        mv "\$(find "fastqc_${readNumber}" -maxdepth 1 -type f -name '*_fastqc.zip' -print -quit)" "${outputBase}.zip"
        mv "\$(find "fastqc_${readNumber}" -maxdepth 1 -type f -name '*_fastqc.html' -print -quit)" "${outputBase}.html"
        """
    }.join('\n')
    """
    ${commands}

    printf '"${task.process}":\\n    fastqc: "%s"\\n' "\$(fastqc --version | sed 's/FastQC v//')" > ${id}.${stage}.fastqc.versions.yml
    """

    stub:
    def id = meta.id
    def stage = meta.qc_stage
    def readList = reads instanceof List ? reads : [reads]
    def outputs = readList.withIndex().collect { _read, index ->
        def name = "${id}.${stage}.R${index + 1}_fastqc"
        """
        touch ${name}.zip
        printf '<!doctype html><html><body>stub fastqc ${name}</body></html>\\n' > ${name}.html
        """
    }.join('\n')
    """
    ${outputs}
    printf '"${task.process}":\\n    stub: "true"\\n' > ${id}.${stage}.fastqc.versions.yml
    """
}
