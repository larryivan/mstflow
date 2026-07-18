process MAKE_BIOM {
    label 'process_low'
    tag 'feature_table'

    publishDir "${params.outdir}/source_attribution/sourcetracker2/${params.get('sourcetracker2_analysis_label', 'default')}", mode: 'copy'

    input:
    path counts

    output:
    path "feature_table.biom", emit: biom
    path "make_biom_versions.yml", emit: versions

    script:
    """
    biom convert \\
        -i "${counts}" \\
        -o "feature_table.biom" \\
        --to-hdf5 \\
        --table-type 'OTU table'

    printf '"${task.process}":\\n    biom: "%s"\\n' "\$(biom --version | awk '{print \$NF}')" > make_biom_versions.yml
    """

    stub:
    """
    printf '{"id": "stub", "format": "Biological Observation Matrix 1.0"}\\n' > feature_table.biom
    printf '"${task.process}":\\n    stub: "true"\\n' > make_biom_versions.yml
    """
}
