process HOST_REMOVE_BWA {
    label 'process_high'
    tag "${meta.id}"

    publishDir "${params.outdir}/host_removed", mode: 'copy',
        saveAs: { filename ->
            def publishIntermediates = (params.publish_intermediates ?: false).toString().toBoolean()
            filename.endsWith('.fq.gz') && !publishIntermediates ? null : filename
        }

    input:
    tuple val(meta), path(reads), path(host_bundle)
    path fragment_helper

    output:
    tuple val(meta), path("*.dehosted*.fq.gz"), emit: reads
    path "*.host_removal.log", emit: logs
    path "*.host_removal.fragment_qc.tsv", emit: qc
    path "*.host_remove_bwa.versions.yml", emit: versions

    script:
    def id = meta.id
    def readList = (reads instanceof List ? reads : [reads]).collect { read -> read.toString() }
    def hostFiles = host_bundle instanceof List ? host_bundle.collect { hostFile -> hostFile.toString() } : [host_bundle.toString()]
    if (!hostFiles) {
        error "HOST_REMOVE_BWA requires a staged host FASTA/index bundle for ${id}"
    }
    def hostRef = hostFiles[0]
    def gzipLevel = params.fastq_gzip_level == null ? 6 : (params.fastq_gzip_level as int)
    def hostMinMapq = params.host_min_mapq == null ? 0 : (params.host_min_mapq as int)
    if (hostMinMapq < 0) {
        error "host_min_mapq must be non-negative"
    }

    if (meta.age_group == 'ancient' && readList.size() == 3) {
        def seRead = readList.find { read -> read.contains('.prehost_se.') }
        def r1Read = readList.find { read -> read.contains('.prehost_R1.') }
        def r2Read = readList.find { read -> read.contains('.prehost_R2.') }
        if (!seRead || !r1Read || !r2Read) {
            error "Ancient PE input for ${id} must contain prehost_se, prehost_R1, and prehost_R2 files"
        }
        """
        : > "${id}.host_removal.log"

        bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} \
            "${hostRef}" "${seRead}" > "${id}.host_se.sai" 2>> "${id}.host_removal.log"
        bwa samse "${hostRef}" "${id}.host_se.sai" "${seRead}" \
            2>> "${id}.host_removal.log" \
            | python3 "${fragment_helper}" filter-host-sam \
                --sam - --output "${id}.dehosted_se.fq.gz" \
                --minimum-mapq ${hostMinMapq} --stats "${id}.host_se.stats.tsv" \
                --gzip-level ${gzipLevel}

        bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} \
            "${hostRef}" "${r1Read}" > "${id}.host_R1.sai" 2>> "${id}.host_removal.log"
        bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} \
            "${hostRef}" "${r2Read}" > "${id}.host_R2.sai" 2>> "${id}.host_removal.log"
        bwa sampe "${hostRef}" "${id}.host_R1.sai" "${id}.host_R2.sai" "${r1Read}" "${r2Read}" \
            2>> "${id}.host_removal.log" \
            | python3 "${fragment_helper}" filter-host-sam \
                --sam - --output-r1 "${id}.dehosted_R1.fq.gz" --output-r2 "${id}.dehosted_R2.fq.gz" \
                --minimum-mapq ${hostMinMapq} --stats "${id}.host_pair.stats.tsv" \
                --gzip-level ${gzipLevel}

        se_input=\$(awk -F '\t' '\$1 == "input_fragments" { print \$2 }' "${id}.host_se.stats.tsv")
        pe_input=\$(awk -F '\t' '\$1 == "input_fragments" { print \$2 }' "${id}.host_pair.stats.tsv")
        se_retained=\$(awk -F '\t' '\$1 == "retained_fragments" { print \$2 }' "${id}.host_se.stats.tsv")
        pe_retained=\$(awk -F '\t' '\$1 == "retained_fragments" { print \$2 }' "${id}.host_pair.stats.tsv")
        input_fragments=\$((se_input + pe_input))
        retained_fragments=\$((se_retained + pe_retained))
        printf 'metric\tvalue\nhost_min_mapq\t%s\nfastq_gzip_level\t%s\ninput_fragments\t%s\nretained_fragments\t%s\nremoved_host_fragments\t%s\nretained_single_fragments\t%s\nretained_paired_fragments\t%s\n' \
            "${hostMinMapq}" "${gzipLevel}" "\$input_fragments" "\$retained_fragments" \
            "\$((input_fragments - retained_fragments))" "\$se_retained" "\$pe_retained" \
            > "${id}.host_removal.fragment_qc.tsv"

        version=\$(bwa 2>&1 | awk '/^Version:/ { print \$2; exit }' || true)
        [[ -n "\$version" ]] || version="unknown"
        printf '"${task.process}":\n    bwa: "%s"\n' "\$version" > "${id}.host_remove_bwa.versions.yml"
        """
    } else if (meta.age_group == 'ancient' && readList.size() == 2) {
        """
        : > "${id}.host_removal.log"
        bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} \
            "${hostRef}" "${readList[0]}" > "${id}.host_R1.sai" 2>> "${id}.host_removal.log"
        bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} \
            "${hostRef}" "${readList[1]}" > "${id}.host_R2.sai" 2>> "${id}.host_removal.log"
        bwa sampe "${hostRef}" "${id}.host_R1.sai" "${id}.host_R2.sai" "${readList[0]}" "${readList[1]}" \
            2>> "${id}.host_removal.log" \
            | python3 "${fragment_helper}" filter-host-sam \
                --sam - --output-r1 "${id}.dehosted_R1.fq.gz" --output-r2 "${id}.dehosted_R2.fq.gz" \
                --minimum-mapq ${hostMinMapq} --stats "${id}.host_removal.fragment_qc.tsv" \
                --gzip-level ${gzipLevel}
        version=\$(bwa 2>&1 | awk '/^Version:/ { print \$2; exit }' || true)
        [[ -n "\$version" ]] || version="unknown"
        printf '"${task.process}":\n    bwa: "%s"\n' "\$version" > "${id}.host_remove_bwa.versions.yml"
        """
    } else if (readList.size() == 2) {
        """
        : > "${id}.host_removal.log"
        bwa mem -t ${task.cpus} "${hostRef}" "${readList[0]}" "${readList[1]}" \
            2>> "${id}.host_removal.log" \
            | python3 "${fragment_helper}" filter-host-sam \
                --sam - --output-r1 "${id}.dehosted_R1.fq.gz" --output-r2 "${id}.dehosted_R2.fq.gz" \
                --minimum-mapq ${hostMinMapq} --stats "${id}.host_removal.fragment_qc.tsv" \
                --gzip-level ${gzipLevel}
        version=\$(bwa 2>&1 | awk '/^Version:/ { print \$2; exit }' || true)
        [[ -n "\$version" ]] || version="unknown"
        printf '"${task.process}":\n    bwa: "%s"\n' "\$version" > "${id}.host_remove_bwa.versions.yml"
        """
    } else if (meta.age_group == 'ancient') {
        """
        : > "${id}.host_removal.log"
        bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} \
            "${hostRef}" "${readList[0]}" > "${id}.host.sai" 2>> "${id}.host_removal.log"
        bwa samse "${hostRef}" "${id}.host.sai" "${readList[0]}" \
            2>> "${id}.host_removal.log" \
            | python3 "${fragment_helper}" filter-host-sam \
                --sam - --output "${id}.dehosted.fq.gz" \
                --minimum-mapq ${hostMinMapq} --stats "${id}.host_removal.fragment_qc.tsv" \
                --gzip-level ${gzipLevel}
        version=\$(bwa 2>&1 | awk '/^Version:/ { print \$2; exit }' || true)
        [[ -n "\$version" ]] || version="unknown"
        printf '"${task.process}":\n    bwa: "%s"\n' "\$version" > "${id}.host_remove_bwa.versions.yml"
        """
    } else {
        """
        : > "${id}.host_removal.log"
        bwa mem -t ${task.cpus} "${hostRef}" "${readList[0]}" \
            2>> "${id}.host_removal.log" \
            | python3 "${fragment_helper}" filter-host-sam \
                --sam - --output "${id}.dehosted.fq.gz" \
                --minimum-mapq ${hostMinMapq} --stats "${id}.host_removal.fragment_qc.tsv" \
                --gzip-level ${gzipLevel}
        version=\$(bwa 2>&1 | awk '/^Version:/ { print \$2; exit }' || true)
        [[ -n "\$version" ]] || version="unknown"
        printf '"${task.process}":\n    bwa: "%s"\n' "\$version" > "${id}.host_remove_bwa.versions.yml"
        """
    }

    stub:
    def id = meta.id
    def readList = (reads instanceof List ? reads : [reads]).collect { read -> read.toString() }
    def hostMinMapq = params.host_min_mapq == null ? 0 : (params.host_min_mapq as int)
    // Keep ancient all_fragments as three biological files through clean QC.
    // The Kraken task alone creates its temporary placeholder mates.
    def stubReads = readList.size() == 3 ? """
        printf '@${id}\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.dehosted_se.fq.gz
        printf '@${id}/1\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.dehosted_R1.fq.gz
        printf '@${id}/2\\nTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.dehosted_R2.fq.gz
        """ : readList.size() == 2 ? """
        printf '@${id}/1\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.dehosted_R1.fq.gz
        printf '@${id}/2\\nTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.dehosted_R2.fq.gz
        """ : """
        printf '@${id}\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.dehosted.fq.gz
        """
    """
    ${stubReads}
    echo 'stub host removal' > ${id}.host_removal.log
    printf 'metric\tvalue\nhost_min_mapq\t${hostMinMapq}\nfastq_gzip_level\t${params.fastq_gzip_level}\ninput_fragments\t1\nretained_fragments\t1\nremoved_host_fragments\t0\n' > ${id}.host_removal.fragment_qc.tsv
    printf '"${task.process}":\n    stub: "true"\n' > ${id}.host_remove_bwa.versions.yml
    """
}
