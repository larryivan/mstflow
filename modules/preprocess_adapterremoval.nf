process PREPROCESS_ADAPTERREMOVAL {
    label 'process_medium'
    tag "${meta.id}"

    publishDir "${params.outdir}/preprocess/adapterremoval", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.prehost.fq.gz"), emit: reads
    path "*.adapterremoval.log", emit: logs
    path "versions.yml", emit: versions

    script:
    def id = meta.id
    def readList = reads.collect { it.toString() }
    def extraArgs = params.adapterremoval_extra_args ?: ''
    if (readList.size() == 2) {
        """
        AdapterRemoval \\
            --file1 ${readList[0]} \\
            --file2 ${readList[1]} \\
            --basename ${id}.ar \\
            --threads ${task.cpus} \\
            --trimns \\
            --trimqualities \\
            --minquality 20 \\
            --minlength ${params.min_length_ancient} \\
            --gzip \\
            --collapse \\
            ${extraArgs} \\
            > ${id}.adapterremoval.log 2>&1

        : > ${id}.retained_parts.txt
        for f in \\
            ${id}.ar.collapsed.gz \\
            ${id}.ar.collapsed.truncated.gz \\
            ${id}.ar.singleton.truncated.gz \\
            ${id}.ar.pair1.truncated.gz \\
            ${id}.ar.pair2.truncated.gz
        do
            if [[ -s "\$f" ]]; then
                echo "\$f" >> ${id}.retained_parts.txt
            fi
        done

        if [[ ! -s ${id}.retained_parts.txt ]]; then
            echo "AdapterRemoval produced no retained reads for ${id}" >&2
            cat ${id}.adapterremoval.log >&2
            exit 1
        fi

        xargs cat < ${id}.retained_parts.txt > ${id}.prehost.fq.gz

        printf '"${task.process}":\\n    AdapterRemoval: "%s"\\n' "\$(AdapterRemoval --version 2>&1 | head -n 1 | sed 's/^AdapterRemoval ver. //')" > versions.yml
        """
    } else {
        """
        AdapterRemoval \\
            --file1 ${readList[0]} \\
            --basename ${id}.ar \\
            --threads ${task.cpus} \\
            --trimns \\
            --trimqualities \\
            --minquality 20 \\
            --minlength ${params.min_length_ancient} \\
            --gzip \\
            ${extraArgs} \\
            > ${id}.adapterremoval.log 2>&1

        : > ${id}.retained_parts.txt
        for f in \\
            ${id}.ar.truncated.gz \\
            ${id}.ar.singleton.truncated.gz
        do
            if [[ -s "\$f" ]]; then
                echo "\$f" >> ${id}.retained_parts.txt
            fi
        done

        if [[ ! -s ${id}.retained_parts.txt ]]; then
            echo "AdapterRemoval produced no retained reads for ${id}" >&2
            cat ${id}.adapterremoval.log >&2
            exit 1
        fi

        xargs cat < ${id}.retained_parts.txt > ${id}.prehost.fq.gz

        printf '"${task.process}":\\n    AdapterRemoval: "%s"\\n' "\$(AdapterRemoval --version 2>&1 | head -n 1 | sed 's/^AdapterRemoval ver. //')" > versions.yml
        """
    }

    stub:
    def id = meta.id
    """
    printf '@${id}\\nACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.prehost.fq.gz
    echo 'stub AdapterRemoval' > ${id}.adapterremoval.log
    printf '"${task.process}":\\n    stub: "true"\\n' > versions.yml
    """
}
