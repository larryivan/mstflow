include { isOn } from '../subworkflows/utils'

process KRAKEN2_PROFILE {
    label 'process_high'
    tag "${meta.id}"
    container "${params.container_kraken2}"

    publishDir "${params.outdir}/kraken2", mode: 'copy',
        saveAs: { filename ->
            filename.endsWith('.kraken2.output.txt') && !isOn(params.kraken2_save_read_classifications) ? null : filename
        }

    input:
    tuple val(meta), path(reads), path(read_length_stats), path(kraken_db)

    output:
    tuple val(meta), path("*.kraken2.report.txt"), path("*.kraken2.read_length.tsv"), emit: reports
    path "*.kraken2.minimizer.report.txt", optional: true, emit: minimizer_reports
    path "*.kraken2.output.txt", optional: true, emit: classifications
    path "*.kraken2.versions.yml", emit: versions

    script:
    def id = meta.id
    def readList = (reads instanceof List ? reads : [reads]).collect { read -> read.toString() }
    def paired = readList.size() == 2
    def saveClassifications = isOn(params.kraken2_save_read_classifications)
    def reportMinimizers = isOn(params.kraken2_report_minimizer_data)
    def reportTarget = reportMinimizers ? "${id}.kraken2.minimizer.report.txt" : "${id}.kraken2.report.txt"
    // PREPARE_KRAKEN_INPUT names the mates; do not rely on staging order for R1 and R2.
    def orderedReads = paired
        ? [readList.find { read -> read.contains('.kraken_input_R1.') }, readList.find { read -> read.contains('.kraken_input_R2.') }]
        : readList
    """
    # Kraken2 2.1.6 ignores unknown long options while returning success, so its
    # stderr is replayed and the task fails if any requested option was ignored.
    set +e
    kraken2 \
        --db "${kraken_db}" --threads ${task.cpus} --gzip-compressed \
        ${paired ? '--paired' : ''} \
        ${params.kraken2_extra_args} \
        --confidence ${params.kraken2_confidence} --minimum-hit-groups ${params.kraken2_minimum_hit_groups} \
        ${isOn(params.kraken2_use_names) ? '--use-names' : ''} \
        ${reportMinimizers ? '--report-minimizer-data' : ''} \
        --report "${reportTarget}" \
        --output ${saveClassifications ? "\"${id}.kraken2.output.txt\"" : '/dev/null'} \
        ${orderedReads.collect { read -> "\"${read}\"" }.join(' ')} \
        2> "${id}.kraken2.stderr.log"
    kraken_exit=\$?
    set -e
    cat "${id}.kraken2.stderr.log" >&2
    if [[ "\$kraken_exit" -ne 0 ]]; then
        exit "\$kraken_exit"
    fi
    if grep -q '^Unknown option:' "${id}.kraken2.stderr.log"; then
        printf 'ERROR: Kraken2 ignored an unrecognised option; refusing to trust this profile.\\n' >&2
        exit 1
    fi

    ${reportMinimizers ? "cut -f1-3,6-8 \"${id}.kraken2.minimizer.report.txt\" > \"${id}.kraken2.report.txt\"" : ''}
    cp "${read_length_stats}" "${id}.kraken2.read_length.tsv"
    version=\$(kraken2 --version 2>&1 | awk 'NR == 1 { print \$NF }')
    printf '"${task.process}":\n    kraken2: "%s"\n' "\$version" > "${id}.kraken2.versions.yml"
    """

    stub:
    def id = meta.id
    def stubClassification = isOn(params.kraken2_save_read_classifications) ? "printf 'U\\tstub\\t0\\t40\\tA:40\\n' > ${id}.kraken2.output.txt" : ''
    def stubMinimizers = isOn(params.kraken2_report_minimizer_data) ? "printf '0.00\\t0\\t0\\t0\\t0\\tU\\t0\\tunclassified\\n100.00\\t10\\t0\\t10\\t10\\tR\\t1\\troot\\n100.00\\t10\\t10\\t10\\t10\\tG\\t1301\\t  Streptococcus\\n' > ${id}.kraken2.minimizer.report.txt" : ''
    """
    printf '  0.00\t0\t0\tU\t0\tunclassified\n  100.00\t10\t0\tR\t1\troot\n  100.00\t10\t10\tG\t1301\t    Streptococcus\n' > ${id}.kraken2.report.txt
    cp "${read_length_stats}" ${id}.kraken2.read_length.tsv
    ${stubClassification}
    ${stubMinimizers}
    printf '"${task.process}":\n    stub: "true"\n' > ${id}.kraken2.versions.yml
    """
}
