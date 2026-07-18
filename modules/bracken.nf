process BRACKEN_PROFILE {
    label 'process_medium'
    tag "${meta.id}"

    publishDir "${params.outdir}/bracken", mode: 'copy'

    input:
    tuple val(meta), path(report), path(read_length_stats), path(bracken_db)
    path fragment_helper

    output:
    tuple val(meta), path("*.bracken.tsv"), path("*.bracken.report.txt"), path("*.bracken.read_length_selection.tsv"), emit: tables
    path "*.bracken.versions.yml", emit: versions

    script:
    def id = meta.id
    def extraArgs = params.bracken_extra_args ?: ''
    def sampleReadLength = meta.bracken_read_len
    def hasSampleReadLength = sampleReadLength && !(sampleReadLength.toString().toUpperCase() in ['NA', 'N/A', 'AUTO'])
    def requestedReadLength = hasSampleReadLength ? sampleReadLength : (params.bracken_read_len ?: 'auto')
    """
    selected_read_length=\$(python3 "${fragment_helper}" select-bracken-length \
        --stats "${read_length_stats}" \
        --database "${bracken_db}" \
        --requested "${requestedReadLength}" \
        --sample-id "${id}" \
        --output "${id}.bracken.read_length_selection.tsv")

    bracken \
        ${extraArgs} \
        -d "${bracken_db}" \
        -i "${report}" \
        -o "${id}.bracken.tsv" \
        -w "${id}.bracken.report.txt" \
        -r "\$selected_read_length" \
        -l "${params.bracken_rank}" \
        -t "${params.bracken_threshold}"

    # Bracken does not consistently implement -v/--version across releases.
    # Parse the help banner without allowing its non-zero help exit to fail a task.
    version=\$(bracken -h 2>&1 \
        | sed -nE 's/.*[Bb]racken[^0-9]*([0-9]+(\\.[0-9]+){1,2}).*/\\1/p' \
        | head -n 1 || true)
    [[ -n "\$version" ]] || version="unknown"
    printf '"${task.process}":\n    bracken: "%s"\n' "\$version" > "${id}.bracken.versions.yml"
    """

    stub:
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
        > /dev/null
    printf 'name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\nBacteroides\t816\tG\t1\t0\t1\t1.0\n' \
        > ${id}.bracken.tsv
    cp ${id}.bracken.tsv ${id}.bracken.report.txt
    printf '"${task.process}":\n    stub: "true"\n' > ${id}.bracken.versions.yml
    """
}
