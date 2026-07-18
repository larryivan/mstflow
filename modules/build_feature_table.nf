process BUILD_FEATURE_TABLE {
    label 'process_low'
    tag "${params.bracken_rank}:${params.taxon_filter}"

    publishDir "${params.outdir}/tables", mode: 'copy'

    input:
    path metadata
    path bracken_files
    path taxdump_nodes
    path build_feature_table_script

    output:
    path "bracken_${params.bracken_rank}_counts.tsv", emit: raw_counts
    path "bracken_${params.bracken_rank}_relative_abundance.tsv", emit: raw_relative
    path "bracken_${params.bracken_rank}_features.tsv", emit: raw_features
    path "bracken_${params.bracken_rank}_${params.taxon_filter}_counts.tsv", emit: counts
    path "bracken_${params.bracken_rank}_${params.taxon_filter}_relative_abundance.tsv", emit: relative
    path "bracken_${params.bracken_rank}_${params.taxon_filter}_features.tsv", emit: features
    path "bracken_${params.bracken_rank}_${params.taxon_filter}_filter_qc.tsv", emit: qc
    path "bracken_${params.bracken_rank}_${params.taxon_filter}_sample_qc.tsv", emit: sample_qc
    path "build_feature_table_versions.yml", emit: versions

    script:
    def taxdumpArg = params.taxon_filter == 'prokaryote' ? "--taxdump \"${taxdump_nodes}\"" : ''
    def brackenArgs = bracken_files.collect { table -> "\"${table}\"" }.join(' ')
    """
    python3 "${build_feature_table_script}" \\
        --metadata "${metadata}" \\
        --prefix bracken_${params.bracken_rank} \\
        --taxon-filter ${params.taxon_filter} \\
        --drop-taxids '${params.drop_taxids}' \\
        ${taxdumpArg} \\
        ${brackenArgs}

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > build_feature_table_versions.yml
    """
}


/*
 * Strict scientific-QC entry point.  It preserves BUILD_FEATURE_TABLE outputs
 * while requiring one corresponding Kraken report per metadata sample.
 */
process BUILD_FEATURE_TABLE_WITH_KRAKEN_QC {
    label 'process_low'
    tag "${params.bracken_rank}:${params.taxon_filter}:kraken_qc"

    publishDir "${params.outdir}/tables", mode: 'copy'

    input:
    path metadata
    path bracken_files
    path kraken_reports
    path taxdump_nodes
    path build_feature_table_script

    output:
    path "bracken_${params.bracken_rank}_counts.tsv", emit: raw_counts
    path "bracken_${params.bracken_rank}_relative_abundance.tsv", emit: raw_relative
    path "bracken_${params.bracken_rank}_features.tsv", emit: raw_features
    path "bracken_${params.bracken_rank}_${params.taxon_filter}_counts.tsv", emit: counts
    path "bracken_${params.bracken_rank}_${params.taxon_filter}_relative_abundance.tsv", emit: relative
    path "bracken_${params.bracken_rank}_${params.taxon_filter}_features.tsv", emit: features
    path "bracken_${params.bracken_rank}_${params.taxon_filter}_filter_qc.tsv", emit: qc
    path "bracken_${params.bracken_rank}_${params.taxon_filter}_sample_qc.tsv", emit: sample_qc
    path "build_feature_table_versions.yml", emit: versions

    script:
    def taxdumpArg = params.taxon_filter == 'prokaryote' ? "--taxdump \"${taxdump_nodes}\"" : ''
    def brackenArgs = bracken_files.collect { table -> "\"${table}\"" }.join(' ')
    def krakenArgs = kraken_reports.collect { report -> "\"${report}\"" }.join(' ')
    """
    python3 "${build_feature_table_script}" \\
        --metadata "${metadata}" \\
        --prefix bracken_${params.bracken_rank} \\
        --taxon-filter ${params.taxon_filter} \\
        --drop-taxids '${params.drop_taxids}' \\
        ${taxdumpArg} \\
        ${brackenArgs} \\
        --kraken-reports ${krakenArgs}

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > build_feature_table_versions.yml
    """
}
