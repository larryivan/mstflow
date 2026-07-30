process KRAKEN2_PROFILE {
    label 'process_high'
    tag "${meta.id}"

    publishDir "${params.outdir}/kraken2", mode: 'copy',
        saveAs: { filename ->
            def saveClassifications = (params.kraken2_save_read_classifications ?: false).toString().toBoolean()
            filename.endsWith('.kraken2.output.txt') && !saveClassifications ? null : filename
        }

    input:
    tuple val(meta), path(reads), path(kraken_db)
    path fragment_helper

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
    def kmerLength = (params.kraken2_kmer_len ?: 35) as int
    def requestedMinimum = (params.kraken2_min_read_len ?: kmerLength) as int
    def minimumUsableSpan = Math.max(kmerLength, requestedMinimum)
    def gzipLevel = params.fastq_gzip_level == null ? 6 : (params.fastq_gzip_level as int)
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

    // Kraken2 stderr is captured so it can be inspected, then always replayed
    // so a real failure is still visible in the task log.
    //
    // Kraken2 2.1.6 parses long options in its Perl wrapper with Getopt::Long
    // and never checks the parse result, so an unrecognised long option only
    // prints "Unknown option: <name>" and classification proceeds with exit 0.
    // Fail closed instead of silently classifying under unintended settings.
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

    def qcTail = """
    cp "${id}.kraken2.read_length.tsv" "${id}.kraken2.fragment_qc.tsv"
    printf 'kraken2_confidence\\t%s\\nkraken2_minimum_hit_groups\\t%s\\nkraken2_use_names\\t%s\\nkraken2_report_minimizer_data\\t%s\\nfastq_gzip_level\\t%s\\n' \\
        "${confidence}" "${minimumHitGroups}" "${useNames}" "${reportMinimizers}" "${gzipLevel}" >> "${id}.kraken2.fragment_qc.tsv"
    version=\$(kraken2 --version 2>&1 | awk 'NR == 1 { print \$NF }')
    printf '"${task.process}":\\n    kraken2: "%s"\\n' "\$version" > "${id}.kraken2.versions.yml"
    """.stripIndent().trim()

    if (meta.age_group == 'ancient' && readList.size() == 3) {
        def seRead = readList.find { read ->
            read.contains('.prehost_se.') || read.contains('.dehosted_se.')
        }
        def r1Read = readList.find { read ->
            read.contains('.prehost_R1.') || read.contains('.dehosted_R1.')
        }
        def r2Read = readList.find { read ->
            read.contains('.prehost_R2.') || read.contains('.dehosted_R2.')
        }
        if (!seRead || !r1Read || !r2Read) {
            error "Ancient three-file input for ${id} must contain se, R1, and R2 biological FASTQs"
        }
        """
        # ancient_pe_mode=all_fragments reaches this branch with either
        # prehost biological files (host removal skipped) or dehosted files.
        # Collapsed molecules and complete unmerged pairs are emitted as native
        # Kraken2 paired input, one classification record per molecule. A
        # collapsed fragment receives a sub-k-mer placeholder mate that
        # contributes zero minimizer positions, so its confidence denominator
        # equals native single-end classification.
        python3 "${fragment_helper}" pair-fragments \\
            --single "${seRead}" \\
            --r1 "${r1Read}" --r2 "${r2Read}" \\
            --output-r1 "${id}.kraken_R1.fq.gz" --output-r2 "${id}.kraken_R2.fq.gz" \\
            --stats "${id}.kraken2.read_length.tsv" \\
            --minimum-usable-span ${minimumUsableSpan} \\
            --kmer-length ${kmerLength} --gzip-level ${gzipLevel}

        retained=\$(awk -F '\t' '\$1 == "retained_fragments" { print \$2 }' "${id}.kraken2.read_length.tsv")
        if [[ -z "\$retained" || "\$retained" -eq 0 ]]; then
            echo "No ${id} fragments contain a usable span of at least ${minimumUsableSpan} bp" >&2
            exit 1
        fi

        set +e
        kraken2 \\
            --db "${db}" \\
            --threads ${task.cpus} \\
            --gzip-compressed \\
            --paired \\
            ${extraArgs} \\
            --confidence ${confidence} \\
            --minimum-hit-groups ${minimumHitGroups} \\
            ${useNamesArg} ${minimizerArg} \\
            --report "${reportTarget}" \\
            ${outputArg} \\
            "${id}.kraken_R1.fq.gz" "${id}.kraken_R2.fq.gz" \\
            2> "${id}.kraken2.stderr.log"

        ${inspectKrakenStderr}
        ${standardizeReport}
        ${qcTail}
        """
    } else if (readList.size() == 2) {
        """
        # Native paired classification. The helper fails closed on unequal
        # record counts, mismatched normalised mate IDs, and out-of-order mates
        # before Kraken2 runs; Kraken2 2.1.6 has no mate-name check of its own.
        python3 "${fragment_helper}" stats \\
            --r1 "${readList[0]}" --r2 "${readList[1]}" \\
            --output "${id}.kraken2.read_length.tsv" \\
            --minimum-usable-span ${minimumUsableSpan} --kmer-length ${kmerLength}
        retained=\$(awk -F '\t' '\$1 == "retained_fragments" { print \$2 }' "${id}.kraken2.read_length.tsv")
        if [[ -z "\$retained" || "\$retained" -eq 0 ]]; then
            echo "No ${id} pairs contain a usable span of at least ${minimumUsableSpan} bp" >&2
            exit 1
        fi
        kraken_r1="${readList[0]}"
        kraken_r2="${readList[1]}"
        dropped=\$(awk -F '\t' '\$1 == "dropped_below_minimum_usable_span" { print \$2 }' "${id}.kraken2.read_length.tsv")
        if [[ "\$dropped" -gt 0 ]]; then
            python3 "${fragment_helper}" prepare-paired \\
                --r1 "${readList[0]}" --r2 "${readList[1]}" \\
                --output-r1 "${id}.kraken_R1.fq.gz" --output-r2 "${id}.kraken_R2.fq.gz" \\
                --stats "${id}.kraken2.filtered.stats.tsv" \\
                --minimum-usable-span ${minimumUsableSpan} --kmer-length ${kmerLength} \\
                --gzip-level ${gzipLevel}
            kraken_r1="${id}.kraken_R1.fq.gz"
            kraken_r2="${id}.kraken_R2.fq.gz"
        fi
        set +e
        kraken2 \\
            --db "${db}" --threads ${task.cpus} --gzip-compressed \\
            --paired \\
            ${extraArgs} \\
            --confidence ${confidence} --minimum-hit-groups ${minimumHitGroups} \\
            ${useNamesArg} ${minimizerArg} \\
            --report "${reportTarget}" ${outputArg} \\
            "\$kraken_r1" "\$kraken_r2" \\
            2> "${id}.kraken2.stderr.log"
        ${inspectKrakenStderr}
        ${standardizeReport}
        ${qcTail}
        """
    } else {
        """
        python3 "${fragment_helper}" stats \\
            --single "${readList[0]}" \\
            --output "${id}.kraken2.read_length.tsv" \\
            --minimum-usable-span ${minimumUsableSpan} --kmer-length ${kmerLength}
        retained=\$(awk -F '\t' '\$1 == "retained_fragments" { print \$2 }' "${id}.kraken2.read_length.tsv")
        if [[ -z "\$retained" || "\$retained" -eq 0 ]]; then
            echo "No ${id} reads contain a usable span of at least ${minimumUsableSpan} bp" >&2
            exit 1
        fi
        kraken_input="${readList[0]}"
        dropped=\$(awk -F '\t' '\$1 == "dropped_below_minimum_usable_span" { print \$2 }' "${id}.kraken2.read_length.tsv")
        if [[ "\$dropped" -gt 0 ]]; then
            python3 "${fragment_helper}" prepare-single \\
                --single "${readList[0]}" --output "${id}.kraken_input.fq.gz" \\
                --stats "${id}.kraken2.filtered.stats.tsv" \\
                --minimum-usable-span ${minimumUsableSpan} --kmer-length ${kmerLength} \\
                --gzip-level ${gzipLevel}
            kraken_input="${id}.kraken_input.fq.gz"
        fi
        set +e
        kraken2 \\
            --db "${db}" --threads ${task.cpus} --gzip-compressed \\
            ${extraArgs} \\
            --confidence ${confidence} --minimum-hit-groups ${minimumHitGroups} \\
            ${useNamesArg} ${minimizerArg} \\
            --report "${reportTarget}" ${outputArg} \\
            "\$kraken_input" \\
            2> "${id}.kraken2.stderr.log"
        ${inspectKrakenStderr}
        ${standardizeReport}
        ${qcTail}
        """
    }

    stub:
    def id = meta.id
    def saveClassifications = (params.kraken2_save_read_classifications ?: false).toString().toBoolean()
    def reportMinimizers = params.kraken2_report_minimizer_data == null ? true : params.kraken2_report_minimizer_data.toString().toBoolean()
    def stubClassification = saveClassifications ? "printf 'U\\tstub\\t0\\t40\\tA:40\\n' > ${id}.kraken2.output.txt" : ':'
    def stubMinimizers = reportMinimizers ? "printf '0.00\\t0\\t0\\t0\\t0\\tU\\t0\\tunclassified\\n100.00\\t1\\t0\\t1\\t1\\tR\\t1\\troot\\n' > ${id}.kraken2.minimizer.report.txt" : ':'
    """
    printf '  0.00\t0\t0\tU\t0\tunclassified\n  100.00\t1\t0\tR\t1\troot\n' > ${id}.kraken2.report.txt
    printf 'metric\tvalue\nretained_fragments\t1\nmedian_informative_bases\t75\nmedian_kmer_equivalent_read_length\t75\nmedian_bracken_read_length\t75\nmedian_kraken_confidence_denominator\t41\n' \
        > ${id}.kraken2.read_length.tsv
    cp ${id}.kraken2.read_length.tsv ${id}.kraken2.fragment_qc.tsv
    ${stubClassification}
    ${stubMinimizers}
    printf '"${task.process}":\n    stub: "true"\n' > ${id}.kraken2.versions.yml
    """
}
