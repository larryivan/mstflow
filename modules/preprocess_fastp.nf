process PREPROCESS_FASTP {
    label 'process_medium'
    tag "${meta.id}"

    publishDir "${params.outdir}/preprocess/fastp", mode: 'copy', saveAs: { filename ->
        (params.publish_intermediates || !filename.endsWith('.fq.gz')) ? filename : null
    }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.prehost*.fq.gz"), emit: reads
    path "*.fastp.{json,html,log}", emit: logs
    path "${meta.id}.fastp.versions.yml", emit: versions

    script:
    def id = meta.id
    def readList = (reads instanceof List ? reads : [reads]).collect { read -> read.toString() }
    def extraArgs = params.fastp_extra_args ?: ''
    if (readList.size() == 2) {
        """
        fastp \\
            --in1 "${readList[0]}" \\
            --in2 "${readList[1]}" \\
            --out1 "${id}.prehost_R1.fq.gz" \\
            --out2 "${id}.prehost_R2.fq.gz" \\
            --detect_adapter_for_pe \\
            --cut_front --cut_tail \\
            --cut_window_size 4 \\
            --cut_mean_quality 20 \\
            --length_required ${params.min_length_modern} \\
            --thread ${task.cpus} \\
            --html "${id}.fastp.html" \\
            --json "${id}.fastp.json" \\
            ${extraArgs} \\
            > "${id}.fastp.log" 2>&1

        printf '"${task.process}":\\n    fastp: "%s"\\n' "\$(fastp --version 2>&1 | sed 's/fastp //')" > ${id}.fastp.versions.yml
        """
    } else {
        """
        fastp \\
            --in1 "${readList[0]}" \\
            --out1 "${id}.prehost.fq.gz" \\
            --cut_front --cut_tail \\
            --cut_window_size 4 \\
            --cut_mean_quality 20 \\
            --length_required ${params.min_length_modern} \\
            --thread ${task.cpus} \\
            --html "${id}.fastp.html" \\
            --json "${id}.fastp.json" \\
            ${extraArgs} \\
            > "${id}.fastp.log" 2>&1

        printf '"${task.process}":\\n    fastp: "%s"\\n' "\$(fastp --version 2>&1 | sed 's/fastp //')" > ${id}.fastp.versions.yml
        """
    }

    stub:
    def id = meta.id
    def readList = (reads instanceof List ? reads : [reads]).collect { read -> read.toString() }
    def stubReads = readList.size() == 2 ? """
        printf '@${id}/1\\nACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.prehost_R1.fq.gz
        printf '@${id}/2\\nACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.prehost_R2.fq.gz
        """ : """
        printf '@${id}\\nACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.prehost.fq.gz
        """
    """
    ${stubReads}
    printf '{"stub": true}\\n' > ${id}.fastp.json
    printf '<!doctype html><html><body>stub fastp</body></html>\\n' > ${id}.fastp.html
    echo 'stub fastp' > ${id}.fastp.log
    printf '"${task.process}":\\n    stub: "true"\\n' > ${id}.fastp.versions.yml
    """
}
