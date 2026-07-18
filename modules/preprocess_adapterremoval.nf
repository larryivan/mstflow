process PREPROCESS_ADAPTERREMOVAL {
    label 'process_medium'
    tag "${meta.id}"

    publishDir "${params.outdir}/preprocess/adapterremoval", mode: 'copy',
        saveAs: { filename ->
            def publishIntermediates = (params.publish_intermediates ?: false).toString().toBoolean()
            filename.endsWith('.fq.gz') && !publishIntermediates ? null : filename
        }

    input:
    tuple val(meta), path(reads)
    path fragment_helper

    output:
    tuple val(meta), path("*.prehost*.fq.gz"), emit: reads
    path "*.adapterremoval.log", emit: logs
    path "*.adapterremoval.fragment_qc.tsv", emit: qc
    path "*.preprocess_adapterremoval.versions.yml", emit: versions

    script:
    def id = meta.id
    def readList = (reads instanceof List ? reads : [reads]).collect { read -> read.toString() }
    def extraArgs = params.adapterremoval_extra_args ?: ''
    def ancientPeMode = (params.ancient_pe_mode ?: 'collapsed_only').toString()
    if (!(ancientPeMode in ['collapsed_only', 'all_fragments'])) {
        error "ancient_pe_mode must be collapsed_only or all_fragments (got '${ancientPeMode}')"
    }
    if (readList.size() == 2 && ancientPeMode == 'collapsed_only') {
        """
        AdapterRemoval \
            --file1 "${readList[0]}" \
            --file2 "${readList[1]}" \
            --basename "${id}.ar" \
            --threads ${task.cpus} \
            --trimns \
            --trimqualities \
            --minquality 20 \
            --minlength ${params.min_length_ancient} \
            --gzip \
            --collapse \
            --collapse-deterministic \
            ${extraArgs} \
            > "${id}.adapterremoval.log" 2>&1

        # The default analysis unit follows steps.md: one collapsed molecule per
        # FASTQ record. Singletons and unmerged pairs are quantified but excluded.
        : > "${id}.prehost.fq.gz"
        for f in "${id}.ar.collapsed.gz" "${id}.ar.collapsed.truncated.gz"; do
            if [[ -f "\$f" ]]; then
                cat "\$f" >> "${id}.prehost.fq.gz"
            fi
        done
        if [[ ! -s "${id}.prehost.fq.gz" ]]; then
            printf '' | gzip -c > "${id}.prehost.fq.gz"
        fi

        for f in \
            "${id}.ar.singleton.truncated.gz" \
            "${id}.ar.pair1.truncated.gz" \
            "${id}.ar.pair2.truncated.gz"
        do
            if [[ ! -f "\$f" ]]; then
                printf '' | gzip -c > "\$f"
            fi
        done

        python3 "${fragment_helper}" stats \
            --single "${id}.prehost.fq.gz" \
            --kmer-length ${params.kraken2_kmer_len ?: 35} \
            --output "${id}.collapsed.stats.tsv"
        python3 "${fragment_helper}" stats \
            --single "${id}.ar.singleton.truncated.gz" \
            --kmer-length ${params.kraken2_kmer_len ?: 35} \
            --output "${id}.singleton.stats.tsv"
        python3 "${fragment_helper}" stats \
            --r1 "${id}.ar.pair1.truncated.gz" \
            --r2 "${id}.ar.pair2.truncated.gz" \
            --kmer-length ${params.kraken2_kmer_len ?: 35} \
            --output "${id}.unmerged.stats.tsv"

        collapsed=\$(awk -F '\t' '\$1 == "retained_fragments" { print \$2 }' "${id}.collapsed.stats.tsv")
        singleton=\$(awk -F '\t' '\$1 == "retained_fragments" { print \$2 }' "${id}.singleton.stats.tsv")
        unmerged=\$(awk -F '\t' '\$1 == "retained_fragments" { print \$2 }' "${id}.unmerged.stats.tsv")
        printf 'metric\tvalue\nancient_pe_mode\tcollapsed_only\nretained_collapsed_fragments\t%s\ndiscarded_singleton_fragments\t%s\ndiscarded_unmerged_fragments\t%s\nretained_fragments\t%s\n' \
            "\$collapsed" "\$singleton" "\$unmerged" "\$collapsed" \
            > "${id}.adapterremoval.fragment_qc.tsv"
        if [[ "\$collapsed" -eq 0 ]]; then
            echo "AdapterRemoval produced no collapsed molecular fragments for ${id}" >&2
            cat "${id}.adapterremoval.log" >&2
            exit 1
        fi

        version=\$(AdapterRemoval --version 2>&1 | head -n 1 | sed 's/^AdapterRemoval ver. //')
        printf '"${task.process}":\n    AdapterRemoval: "%s"\n' "\$version" > "${id}.preprocess_adapterremoval.versions.yml"
        """
    } else if (readList.size() == 2) {
        """
        AdapterRemoval \
            --file1 "${readList[0]}" \
            --file2 "${readList[1]}" \
            --basename "${id}.ar" \
            --threads ${task.cpus} \
            --trimns \
            --trimqualities \
            --minquality 20 \
            --minlength ${params.min_length_ancient} \
            --gzip \
            --collapse \
            --collapse-deterministic \
            ${extraArgs} \
            > "${id}.adapterremoval.log" 2>&1

        # Retain only complete molecular fragments, matching the modern PE
        # branch: collapsed molecules plus complete unmerged pairs. Quantify
        # orphan singletons below, but do not classify them.
        : > "${id}.prehost_se.fq.gz"
        for f in \
            "${id}.ar.collapsed.gz" \
            "${id}.ar.collapsed.truncated.gz"
        do
            if [[ -f "\$f" ]]; then
                cat "\$f" >> "${id}.prehost_se.fq.gz"
            fi
        done
        if [[ ! -s "${id}.prehost_se.fq.gz" ]]; then
            printf '' | gzip -c > "${id}.prehost_se.fq.gz"
        fi

        for mate in 1 2; do
            source_file="${id}.ar.pair\${mate}.truncated.gz"
            output_file="${id}.prehost_R\${mate}.fq.gz"
            if [[ -f "\$source_file" ]]; then
                cp "\$source_file" "\$output_file"
            else
                printf '' | gzip -c > "\$output_file"
            fi
        done

        if [[ ! -f "${id}.ar.singleton.truncated.gz" ]]; then
            printf '' | gzip -c > "${id}.ar.singleton.truncated.gz"
        fi

        python3 "${fragment_helper}" stats \
            --single "${id}.prehost_se.fq.gz" \
            --r1 "${id}.prehost_R1.fq.gz" \
            --r2 "${id}.prehost_R2.fq.gz" \
            --kmer-length ${params.kraken2_kmer_len ?: 35} \
            --output "${id}.all_fragments.stats.tsv"
        python3 "${fragment_helper}" stats \
            --single "${id}.ar.singleton.truncated.gz" \
            --kmer-length ${params.kraken2_kmer_len ?: 35} \
            --output "${id}.singleton.stats.tsv"

        {
            head -n 1 "${id}.all_fragments.stats.tsv"
            printf 'ancient_pe_mode\tall_fragments\n'
            tail -n +2 "${id}.all_fragments.stats.tsv"
            singleton=\$(awk -F '\t' '\$1 == "retained_fragments" { print \$2 }' "${id}.singleton.stats.tsv")
            printf 'discarded_singleton_fragments\t%s\n' "\$singleton"
        } > "${id}.adapterremoval.fragment_qc.tsv"

        retained=\$(awk -F '\t' '\$1 == "retained_fragments" { print \$2 }' "${id}.adapterremoval.fragment_qc.tsv")
        if [[ -z "\$retained" || "\$retained" -eq 0 ]]; then
            echo "AdapterRemoval produced no retained molecular fragments for ${id}" >&2
            cat "${id}.adapterremoval.log" >&2
            exit 1
        fi

        version=\$(AdapterRemoval --version 2>&1 | head -n 1 | sed 's/^AdapterRemoval ver. //')
        printf '"${task.process}":\n    AdapterRemoval: "%s"\n' "\$version" > "${id}.preprocess_adapterremoval.versions.yml"
        """
    } else {
        """
        AdapterRemoval \
            --file1 "${readList[0]}" \
            --basename "${id}.ar" \
            --threads ${task.cpus} \
            --trimns \
            --trimqualities \
            --minquality 20 \
            --minlength ${params.min_length_ancient} \
            --gzip \
            ${extraArgs} \
            > "${id}.adapterremoval.log" 2>&1

        : > "${id}.prehost.fq.gz"
        for f in "${id}.ar.truncated.gz" "${id}.ar.singleton.truncated.gz"; do
            if [[ -f "\$f" ]]; then
                cat "\$f" >> "${id}.prehost.fq.gz"
            fi
        done
        if [[ ! -s "${id}.prehost.fq.gz" ]]; then
            printf '' | gzip -c > "${id}.prehost.fq.gz"
        fi

        python3 "${fragment_helper}" stats \
            --single "${id}.prehost.fq.gz" \
            --kmer-length ${params.kraken2_kmer_len ?: 35} \
            --output "${id}.adapterremoval.fragment_qc.tsv"

        retained=\$(awk -F '\t' '\$1 == "retained_fragments" { print \$2 }' "${id}.adapterremoval.fragment_qc.tsv")
        if [[ -z "\$retained" || "\$retained" -eq 0 ]]; then
            echo "AdapterRemoval produced no retained molecular fragments for ${id}" >&2
            cat "${id}.adapterremoval.log" >&2
            exit 1
        fi

        version=\$(AdapterRemoval --version 2>&1 | head -n 1 | sed 's/^AdapterRemoval ver. //')
        printf '"${task.process}":\n    AdapterRemoval: "%s"\n' "\$version" > "${id}.preprocess_adapterremoval.versions.yml"
        """
    }

    stub:
    def id = meta.id
    def readList = (reads instanceof List ? reads : [reads]).collect { read -> read.toString() }
    def ancientPeMode = (params.ancient_pe_mode ?: 'collapsed_only').toString()
    def stubReads = (readList.size() == 2 && ancientPeMode == 'all_fragments') ? """
        printf '@${id}/collapsed\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.prehost_se.fq.gz
        printf '@${id}/pair/1\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.prehost_R1.fq.gz
        printf '@${id}/pair/2\\nTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.prehost_R2.fq.gz
        """ : """
        printf '@${id}\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.prehost.fq.gz
        """
    """
    ${stubReads}
    echo 'stub AdapterRemoval' > ${id}.adapterremoval.log
    printf 'metric\tvalue\nancient_pe_mode\t${ancientPeMode}\nretained_fragments\t1\n' > ${id}.adapterremoval.fragment_qc.tsv
    printf '"${task.process}":\n    stub: "true"\n' > ${id}.preprocess_adapterremoval.versions.yml
    """
}
