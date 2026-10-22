process ADAPTERREMOVAL_RAW {
    label 'process_medium'
    tag "${meta.id}"
    container "${params.container_adapterremoval}"

    publishDir "${params.outdir}/preprocess/adapterremoval", mode: 'copy',
        saveAs: { filename -> filename == 'adapterremoval_raw' ? null : filename }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("adapterremoval_raw"), emit: raw
    path "*.adapterremoval.log", emit: logs
    path "*.preprocess_adapterremoval.versions.yml", emit: versions

    script:
    def id = meta.id
    def readList = (reads instanceof List ? reads : [reads]).collect { read -> read.toString() }
    def files = readList.size() == 2
        ? "--file1 \"${readList[0]}\" --file2 \"${readList[1]}\" --collapse --collapse-deterministic"
        : "--file1 \"${readList[0]}\""
    // Recent Illumina instruments write quality 42 ('K'); AdapterRemoval rejects anything
    // above 41 unless told otherwise.
    """
    mkdir -p adapterremoval_raw
    AdapterRemoval \
        ${files} \
        --basename "adapterremoval_raw/${id}.ar" \
        --threads ${task.cpus} \
        --trimns \
        --trimqualities \
        --minquality 20 \
        --qualitymax 42 \
        --minlength ${params.min_length} \
        --gzip \
        ${params.adapterremoval_extra_args} \
        > "${id}.adapterremoval.log" 2>&1

    version=\$(AdapterRemoval --version 2>&1 | head -n 1 | sed 's/^AdapterRemoval ver. //')
    printf '"${task.process}":\n    AdapterRemoval: "%s"\n' "\$version" > "${id}.preprocess_adapterremoval.versions.yml"
    """

    stub:
    def id = meta.id
    def readList = (reads instanceof List ? reads : [reads]).collect { read -> read.toString() }
    def rawReads = readList.size() == 2 ? """
        printf '@${id}/collapsed\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > adapterremoval_raw/${id}.ar.collapsed.gz
        printf '@${id}/pair/1\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > adapterremoval_raw/${id}.ar.pair1.truncated.gz
        printf '@${id}/pair/2\\nTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > adapterremoval_raw/${id}.ar.pair2.truncated.gz
        printf '' | gzip -c > adapterremoval_raw/${id}.ar.singleton.truncated.gz
        """ : """
        printf '@${id}\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > adapterremoval_raw/${id}.ar.truncated.gz
        """
    """
    mkdir -p adapterremoval_raw
    ${rawReads}
    echo 'stub AdapterRemoval' > ${id}.adapterremoval.log
    printf '"${task.process}":\n    stub: "true"\n' > ${id}.preprocess_adapterremoval.versions.yml
    """
}
