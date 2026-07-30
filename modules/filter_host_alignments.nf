process FILTER_HOST_ALIGNMENTS {
    label 'process_medium'
    tag "${meta.id}"
    container "${params.container_python_tools}"

    publishDir "${params.outdir}/host_removed", mode: 'copy',
        saveAs: { filename ->
            def publishIntermediates = (params.publish_intermediates ?: false).toString().toBoolean()
            filename.endsWith('.fq.gz') && !publishIntermediates ? null : filename
        }

    input:
    tuple val(meta), path(alignment_dir)
    path filter_script
    path fragment_helper

    output:
    tuple val(meta), path("*.dehosted*.fq.gz"), emit: reads
    path "*.host_removal.fragment_qc.tsv", emit: qc
    path "*.filter_host_alignments.versions.yml", emit: versions

    script:
    def id = meta.id
    def hostMinMapq = params.host_min_mapq == null ? 0 : (params.host_min_mapq as int)
    def gzipLevel = params.fastq_gzip_level == null ? 6 : (params.fastq_gzip_level as int)
    """
    python3 "${filter_script}" \
        --sample-id "${id}" \
        --alignment-dir "${alignment_dir}" \
        --minimum-mapq ${hostMinMapq} \
        --gzip-level ${gzipLevel} \
        --fragment-helper "${fragment_helper}"

    printf '"${task.process}":\n    python: "%s"\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > "${id}.filter_host_alignments.versions.yml"
    """

    stub:
    def id = meta.id
    def hostMinMapq = params.host_min_mapq == null ? 0 : (params.host_min_mapq as int)
    def gzipLevel = params.fastq_gzip_level == null ? 6 : (params.fastq_gzip_level as int)
    """
    mode=\$(tr -d '\r\n' < "${alignment_dir}/mode.txt")
    case "\$mode" in
        ancient_mixed)
            printf '@${id}\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.dehosted_se.fq.gz
            printf '@${id}/1\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.dehosted_R1.fq.gz
            printf '@${id}/2\\nTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.dehosted_R2.fq.gz
            ;;
        ancient_pe|modern_pe)
            printf '@${id}/1\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.dehosted_R1.fq.gz
            printf '@${id}/2\\nTGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.dehosted_R2.fq.gz
            ;;
        *)
            printf '@${id}\\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.dehosted.fq.gz
            ;;
    esac
    printf 'metric\tvalue\nhost_min_mapq\t${hostMinMapq}\nfastq_gzip_level\t${gzipLevel}\ninput_fragments\t1\nretained_fragments\t1\nremoved_host_fragments\t0\n' > ${id}.host_removal.fragment_qc.tsv
    printf '"${task.process}":\n    stub: "true"\n' > ${id}.filter_host_alignments.versions.yml
    """
}
