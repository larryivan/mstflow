process KRAKEN2_PROFILE {
    label 'process_high'
    tag "${meta.id}"

    publishDir "${params.outdir}/kraken2", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.kraken2.report.txt"), path("*.kraken2.output.txt"), emit: reports
    path "versions.yml", emit: versions

    script:
    def id = meta.id
    def readList = reads.collect { it.toString() }
    def db = file(params.kraken2_db).toString()
    def extraArgs = params.kraken2_extra_args ?: ''
    if (readList.size() == 2) {
        """
        kraken2 \\
            --db ${db} \\
            --threads ${task.cpus} \\
            --gzip-compressed \\
            --paired \\
            --report ${id}.kraken2.report.txt \\
            --output ${id}.kraken2.output.txt \\
            ${extraArgs} \\
            ${readList[0]} ${readList[1]}

        printf '"${task.process}":\\n    kraken2: "%s"\\n' "\$(kraken2 --version 2>&1 | head -n 1 | awk '{print \$NF}')" > versions.yml
        """
    } else {
        """
        kraken2 \\
            --db ${db} \\
            --threads ${task.cpus} \\
            --gzip-compressed \\
            --report ${id}.kraken2.report.txt \\
            --output ${id}.kraken2.output.txt \\
            ${extraArgs} \\
            ${readList[0]}

        printf '"${task.process}":\\n    kraken2: "%s"\\n' "\$(kraken2 --version 2>&1 | head -n 1 | awk '{print \$NF}')" > versions.yml
        """
    }

    stub:
    def id = meta.id
    """
    printf '  100.00\\t1\\t1\\tU\\t0\\tunclassified\\n' > ${id}.kraken2.report.txt
    printf 'U\\tstub\\t0\\t32\\tA:32\\n' > ${id}.kraken2.output.txt
    printf '"${task.process}":\\n    stub: "true"\\n' > versions.yml
    """
}
