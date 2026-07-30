process NORMALIZE_ADAPTERREMOVAL {
    label 'process_low'
    tag "${meta.id}"
    container "${params.container_python_tools}"

    publishDir "${params.outdir}/preprocess/adapterremoval", mode: 'copy',
        saveAs: { filename ->
            def publishIntermediates = (params.publish_intermediates ?: false).toString().toBoolean()
            filename.endsWith('.fq.gz') && !publishIntermediates ? null : filename
        }

    input:
    tuple val(meta), path(raw_dir)
    path normalizer_script
    path fragment_helper

    output:
    tuple val(meta), path("*.prehost*.fq.gz"), emit: reads
    path "*.adapterremoval.fragment_qc.tsv", emit: qc
    path "*.normalize_adapterremoval.versions.yml", emit: versions

    script:
    def id = meta.id
    def ancientPeMode = (params.ancient_pe_mode ?: 'collapsed_only').toString()
    def mode = meta.layout == 'PE' ? ancientPeMode : 'adapterremoval_se'
    def gzipLevel = params.fastq_gzip_level == null ? 6 : (params.fastq_gzip_level as int)
    """
    python3 "${normalizer_script}" \
        --sample-id "${id}" \
        --raw-dir "${raw_dir}" \
        --mode "${mode}" \
        --kmer-length ${params.kraken2_kmer_len ?: 35} \
        --gzip-level ${gzipLevel} \
        --fragment-helper "${fragment_helper}"

    printf '"${task.process}":\n    python: "%s"\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > "${id}.normalize_adapterremoval.versions.yml"
    """

    stub:
    def id = meta.id
    def ancientPeMode = (params.ancient_pe_mode ?: 'collapsed_only').toString()
    def stubReads = (meta.layout == 'PE' && ancientPeMode == 'all_fragments') ? """
        printf '@${id}/collapsed\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.prehost_se.fq.gz
        printf '@${id}/pair/1\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.prehost_R1.fq.gz
        printf '@${id}/pair/2\\nTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.prehost_R2.fq.gz
        """ : """
        printf '@${id}\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.prehost.fq.gz
        """
    """
    ${stubReads}
    printf 'metric\tvalue\nancient_pe_mode\t${ancientPeMode}\nfastq_gzip_level\t${params.fastq_gzip_level}\nretained_fragments\t1\n' > ${id}.adapterremoval.fragment_qc.tsv
    printf '"${task.process}":\n    stub: "true"\n' > ${id}.normalize_adapterremoval.versions.yml
    """
}
