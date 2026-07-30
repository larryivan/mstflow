process SELECT_BRACKEN_LENGTH {
    label 'process_low'
    tag "${meta.id}"
    container "${params.container_python_tools}"

    publishDir "${params.outdir}/bracken", mode: 'copy'

    input:
    tuple val(meta), path(read_length_stats), path(bracken_db)
    path fragment_helper

    output:
    tuple val(meta), path("*.bracken.read_length_selection.tsv"), emit: selections
    path "*.select_bracken_length.versions.yml", emit: versions

    script:
    def id = meta.id
    def sampleReadLength = meta.bracken_read_len
    def hasSampleReadLength = sampleReadLength && !(sampleReadLength.toString().toUpperCase() in ['NA', 'N/A', 'AUTO'])
    def requestedReadLength = hasSampleReadLength ? sampleReadLength : (params.bracken_read_len ?: 'auto')
    """
    python3 "${fragment_helper}" select-bracken-length \
        --stats "${read_length_stats}" \
        --database "${bracken_db}" \
        --requested "${requestedReadLength}" \
        --sample-id "${id}" \
        --output "${id}.bracken.read_length_selection.tsv" \
        > "${id}.selected_read_length.txt"

    selected_read_length=\$(tr -d '\r\n' < "${id}.selected_read_length.txt")
    [[ "\$selected_read_length" =~ ^[1-9][0-9]*\$ ]] || {
        echo "Invalid selected Bracken read length: \$selected_read_length" >&2
        exit 1
    }
    printf '"${task.process}":\n    python: "%s"\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > "${id}.select_bracken_length.versions.yml"
    """

    stub:
    def id = meta.id
    def sampleReadLength = meta.bracken_read_len
    def hasSampleReadLength = sampleReadLength && !(sampleReadLength.toString().toUpperCase() in ['NA', 'N/A', 'AUTO'])
    def requestedReadLength = hasSampleReadLength ? sampleReadLength : (params.bracken_read_len ?: 'auto')
    def selectedReadLength = requestedReadLength.toString().isInteger() ? requestedReadLength : 75
    """
    printf 'sample_id\tobserved_median_bracken_read_length\tselected_bracken_read_length\tselection_mode\tdistribution_file\n${id}\t75\t${selectedReadLength}\tstub\tdatabase${selectedReadLength}mers.kmer_distrib\n' \
        > ${id}.bracken.read_length_selection.tsv
    printf '"${task.process}":\n    stub: "true"\n' > ${id}.select_bracken_length.versions.yml
    """
}
