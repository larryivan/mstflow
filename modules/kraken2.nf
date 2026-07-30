process KRAKEN2_PROFILE {
    label 'process_high'
    tag "${meta.id}"
    container "${params.container_kraken2}"

    publishDir "${params.outdir}/kraken2", mode: 'copy',
        saveAs: { filename ->
            def saveClassifications = (params.kraken2_save_read_classifications ?: false).toString().toBoolean()
            filename.endsWith('.kraken2.output.txt') && !saveClassifications ? null : filename
        }

    input:
    tuple val(meta), path(reads), path(read_length_stats), path(input_mode), path(kraken_db)

    output:
    tuple val(meta), path("*.kraken2.report.txt"), path("*.kraken2.read_length.tsv"), emit: reports
    path "*.kraken2.fragment_qc.tsv", emit: qc
    path "*.kraken2.minimizer.report.txt", optional: true, emit: minimizer_reports
    path "*.kraken2.output.txt", optional: true, emit: classifications
    path "*.kraken2.versions.yml", emit: versions

    script:
    def id = meta.id
    def readList = (reads instanceof List ? reads : [reads]).collect { read -> read.toString() }
    def db = kraken_db.toString()
    def extraArgs = params.kraken2_extra_args ?: ''
    def confidence = params.kraken2_confidence == null ? 0.05 : params.kraken2_confidence
    def minimumHitGroups = params.kraken2_minimum_hit_groups == null ? 2 : params.kraken2_minimum_hit_groups
    def useNames = params.kraken2_use_names == null ? true : params.kraken2_use_names.toString().toBoolean()
    def reportMinimizers = params.kraken2_report_minimizer_data == null ? true : params.kraken2_report_minimizer_data.toString().toBoolean()
    def saveClassifications = (params.kraken2_save_read_classifications ?: false).toString().toBoolean()
    def outputArg = saveClassifications ? "--output \"${id}.kraken2.output.txt\"" : '--output /dev/null'
    def useNamesArg = useNames ? '--use-names' : ''
    def minimizerArg = reportMinimizers ? '--report-minimizer-data' : ''
    def reportTarget = reportMinimizers ? "${id}.kraken2.minimizer.report.txt" : "${id}.kraken2.report.txt"
    def standardizeReport = reportMinimizers \
        ? "cut -f1-3,6-8 \"${id}.kraken2.minimizer.report.txt\" > \"${id}.kraken2.report.txt\"" \
        : ':'
    def pairedR1 = readList.findAll { read -> read.contains('.kraken_input_R1.') }
    def pairedR2 = readList.findAll { read -> read.contains('.kraken_input_R2.') }
    if (readList.size() == 2 && (pairedR1.size() != 1 || pairedR2.size() != 1)) {
        error "Prepared paired Kraken input for ${id} must contain exactly one named R1 and one named R2 FASTQ"
    }
    // Do not make filesystem/glob ordering part of the scientific interface.
    def orderedReads = readList.size() == 2 ? [pairedR1[0], pairedR2[0]] : readList
    def readArgs = orderedReads.collect { read -> "\"${read}\"" }.join(' ')

    // Kraken2 2.1.6 ignores unknown long options while returning success.
    // Capture and replay stderr, but fail closed if its wrapper reports that
    // any requested option was ignored.
    def inspectKrakenStderr = """
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
    """.stripIndent().trim()

    """
    kraken_mode=\$(tr -d '\r\n' < "${input_mode}")
    set +e
    case "\$kraken_mode" in
        paired)
            [[ ${readList.size()} -eq 2 ]] || {
                echo 'Prepared paired Kraken input does not contain exactly two FASTQs' >&2
                exit 1
            }
            kraken2 \
                --db "${db}" --threads ${task.cpus} --gzip-compressed \
                --paired \
                ${extraArgs} \
                --confidence ${confidence} --minimum-hit-groups ${minimumHitGroups} \
                ${useNamesArg} ${minimizerArg} \
                --report "${reportTarget}" ${outputArg} \
                ${readArgs} \
                2> "${id}.kraken2.stderr.log"
            ;;
        single)
            [[ ${readList.size()} -eq 1 ]] || {
                echo 'Prepared single Kraken input does not contain exactly one FASTQ' >&2
                exit 1
            }
            kraken2 \
                --db "${db}" --threads ${task.cpus} --gzip-compressed \
                ${extraArgs} \
                --confidence ${confidence} --minimum-hit-groups ${minimumHitGroups} \
                ${useNamesArg} ${minimizerArg} \
                --report "${reportTarget}" ${outputArg} \
                ${readArgs} \
                2> "${id}.kraken2.stderr.log"
            ;;
        *)
            echo "Invalid prepared Kraken input mode: \$kraken_mode" >&2
            exit 1
            ;;
    esac
    ${inspectKrakenStderr}

    ${standardizeReport}
    cp "${read_length_stats}" "${id}.kraken2.read_length.tsv"
    cp "${read_length_stats}" "${id}.kraken2.fragment_qc.tsv"
    printf 'kraken2_confidence\t%s\nkraken2_minimum_hit_groups\t%s\nkraken2_use_names\t%s\nkraken2_report_minimizer_data\t%s\nfastq_gzip_level\t%s\n' \
        "${confidence}" "${minimumHitGroups}" "${useNames}" "${reportMinimizers}" "${params.fastq_gzip_level}" \
        >> "${id}.kraken2.fragment_qc.tsv"
    version=\$(kraken2 --version 2>&1 | awk 'NR == 1 { print \$NF }')
    printf '"${task.process}":\n    kraken2: "%s"\n' "\$version" > "${id}.kraken2.versions.yml"
    """

    stub:
    def id = meta.id
    def saveClassifications = (params.kraken2_save_read_classifications ?: false).toString().toBoolean()
    def reportMinimizers = params.kraken2_report_minimizer_data == null ? true : params.kraken2_report_minimizer_data.toString().toBoolean()
    def stubClassification = saveClassifications ? "printf 'U\\tstub\\t0\\t40\\tA:40\\n' > ${id}.kraken2.output.txt" : ':'
    def stubMinimizers = reportMinimizers ? "printf '0.00\\t0\\t0\\t0\\t0\\tU\\t0\\tunclassified\\n100.00\\t1\\t0\\t1\\t1\\tR\\t1\\troot\\n' > ${id}.kraken2.minimizer.report.txt" : ':'
    """
    printf '  0.00\t0\t0\tU\t0\tunclassified\n  100.00\t1\t0\tR\t1\troot\n' > ${id}.kraken2.report.txt
    cp "${read_length_stats}" ${id}.kraken2.read_length.tsv
    cp "${read_length_stats}" ${id}.kraken2.fragment_qc.tsv
    ${stubClassification}
    ${stubMinimizers}
    printf '"${task.process}":\n    stub: "true"\n' > ${id}.kraken2.versions.yml
    """
}
