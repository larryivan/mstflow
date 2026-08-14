process SUMMARIZE_SAMPLES {
    label 'process_low'
    tag 'summary'
    container "${params.container_python_tools}"

    publishDir "${params.outdir}/summary", mode: 'copy', pattern: 'sample_summary.tsv'

    input:
    path metadata
    path sample_qc
    path qc_files, stageAs: 'qc/*'
    path summarize_script

    output:
    path "sample_summary.tsv", emit: summary
    path "sample_summary_mqc.tsv", emit: multiqc
    path "summarize_samples_versions.yml", emit: versions

    script:
    """
    python3 "${summarize_script}" \\
        --metadata "${metadata}" \\
        --sample-qc "${sample_qc}" \\
        --taxon-filter ${params.taxon_filter} \\
        --rank ${params.bracken_rank} \\
        --min-classified ${params.qc_min_classified_reads} \\
        --min-taxa ${params.qc_min_taxa} \\
        qc/*

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > summarize_samples_versions.yml
    """

    stub:
    """
    printf 'sample_id\\tflags\\n' > sample_summary.tsv
    cp sample_summary.tsv sample_summary_mqc.tsv
    printf '"${task.process}":\\n    stub: "true"\\n' > summarize_samples_versions.yml
    """
}
