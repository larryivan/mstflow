process PREPARE_KRAKEN_INPUT {
    label 'process_medium'
    tag "${meta.id}"
    container "${params.container_python_tools}"

    publishDir "${params.outdir}/kraken2", mode: 'copy',
        saveAs: { filename -> filename.endsWith('.fq.gz') ? null : filename }

    input:
    tuple val(meta), path(reads)
    path prepare_script
    path fragment_helper

    output:
    tuple val(meta), path("*.kraken_input*.fq.gz"), path("*.kraken_input.read_length.tsv"), path("*.kraken_input.mode.txt"), emit: prepared
    path "*.prepare_kraken_input.versions.yml", emit: versions

    script:
    def id = meta.id
    def readList = (reads instanceof List ? reads : [reads]).collect { read -> read.toString() }
    def readArgs = readList.collect { read -> "--read \"${read}\"" }.join(' ')
    def kmerLength = (params.kraken2_kmer_len ?: 35) as int
    def minimumUsableSpan = Math.max(kmerLength, (params.kraken2_min_read_len ?: kmerLength) as int)
    def gzipLevel = params.fastq_gzip_level == null ? 6 : (params.fastq_gzip_level as int)
    """
    python3 "${prepare_script}" \
        --sample-id "${id}" \
        --age-group "${meta.age_group}" \
        ${readArgs} \
        --minimum-usable-span ${minimumUsableSpan} \
        --kmer-length ${kmerLength} \
        --gzip-level ${gzipLevel} \
        --fragment-helper "${fragment_helper}"

    printf '"${task.process}":\n    python: "%s"\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > "${id}.prepare_kraken_input.versions.yml"
    """

    stub:
    def id = meta.id
    def readList = (reads instanceof List ? reads : [reads]).collect { read -> read.toString() }
    def paired = readList.size() > 1
    def stubReads = paired ? """
        printf '@${id}/1\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.kraken_input_R1.fq.gz
        printf '@${id}/2\\nTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.kraken_input_R2.fq.gz
        """ : """
        printf '@${id}\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.kraken_input.fq.gz
        """
    """
    ${stubReads}
    printf 'metric\tvalue\nretained_fragments\t1\nmedian_informative_bases\t75\nmedian_kmer_equivalent_read_length\t75\nmedian_bracken_read_length\t75\nmedian_kraken_confidence_denominator\t41\n' > ${id}.kraken_input.read_length.tsv
    printf '${paired ? 'paired' : 'single'}\n' > ${id}.kraken_input.mode.txt
    printf '"${task.process}":\n    stub: "true"\n' > ${id}.prepare_kraken_input.versions.yml
    """
}
