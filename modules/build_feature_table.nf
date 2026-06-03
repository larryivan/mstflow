process BUILD_FEATURE_TABLE {
    label 'process_low'
    tag "${params.bracken_rank}:${params.taxon_filter}"

    publishDir "${params.outdir}/tables", mode: 'copy'

    input:
    path metadata
    path bracken_files

    output:
    path "bracken_${params.bracken_rank}_counts.tsv", emit: raw_counts
    path "bracken_${params.bracken_rank}_relative_abundance.tsv", emit: raw_relative
    path "bracken_${params.bracken_rank}_features.tsv", emit: raw_features
    path "bracken_${params.bracken_rank}_${params.taxon_filter}_counts.tsv", emit: counts
    path "bracken_${params.bracken_rank}_${params.taxon_filter}_relative_abundance.tsv", emit: relative
    path "bracken_${params.bracken_rank}_${params.taxon_filter}_features.tsv", emit: features
    path "bracken_${params.bracken_rank}_${params.taxon_filter}_filter_qc.tsv", emit: qc
    path "versions.yml", emit: versions

    script:
    def taxdumpArg = params.taxdump ? "--taxdump ${file(params.taxdump)}" : ''
    """
    python3 ${projectDir}/bin/build_feature_table.py \\
        --metadata ${metadata} \\
        --prefix bracken_${params.bracken_rank} \\
        --taxon-filter ${params.taxon_filter} \\
        --drop-taxids '${params.drop_taxids}' \\
        ${taxdumpArg} \\
        ${bracken_files.join(' ')}

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > versions.yml
    """
}
