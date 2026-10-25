process PREPARE_KRAKEN_INPUT {
    label 'process_medium'
    tag "${meta.id}"
    container "${params.container_python_tools}"

    input:
    tuple val(meta), path(reads)
    path prepare_script
    path fragment_helper

    output:
    tuple val(meta), path("*.kraken_input*.fq.gz"), path("*.kraken_input.read_length.tsv"), emit: prepared
    path "*.prepare_kraken_input.versions.yml", emit: versions

    script:
    def id = meta.id
    def readArgs = (reads instanceof List ? reads : [reads]).collect { read -> "--read \"${read}\"" }.join(' ')
    """
    python3 "${prepare_script}" \
        --sample-id "${id}" \
        ${readArgs} \
        --minimum-usable-span ${params.kraken2_min_read_len} \
        --kmer-length ${params.kraken2_kmer_len} \
        --gzip-level ${params.fastq_gzip_level} \
        --fragment-helper "${fragment_helper}"

    printf '"${task.process}":\n    python: "%s"\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > "${id}.prepare_kraken_input.versions.yml"
    """

    stub:
    def id = meta.id
    def paired = reads instanceof List && reads.size() > 1
    def stubReads = paired ? """
        printf '@${id}/1\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.kraken_input_R1.fq.gz
        printf '@${id}/2\\nTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.kraken_input_R2.fq.gz
        """ : """
        printf '@${id}\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.kraken_input.fq.gz
        """
    """
    ${stubReads}
    printf 'metric\tvalue\nretained_fragments\t1\nmedian_informative_bases\t75\nmedian_kmer_equivalent_read_length\t75\nmedian_kraken_confidence_denominator\t41\n' > ${id}.kraken_input.read_length.tsv
    printf '"${task.process}":\n    stub: "true"\n' > ${id}.prepare_kraken_input.versions.yml
    """
}
