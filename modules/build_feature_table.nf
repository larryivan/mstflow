/*
 * Merge the classifier's per-sample clade counts at one rank into one count table,
 * before and after the taxon filter, and record for every sample how many reads each
 * step kept. Taxa with too few distinct minimizers (Kraken2) or unique k-mers
 * (KrakenUniq) are removed sample by sample first.
 */
process BUILD_FEATURE_TABLE {
    label 'process_low'
    tag "${params.rank}:${params.taxon_filter}"
    container "${params.container_python_tools}"
    publishDir "${params.outdir}/tables", mode: 'copy'

    input:
    path metadata
    path reports, stageAs: 'reports/*'                // Kraken2 or KrakenUniq reports
    path minimizer_reports, stageAs: 'minimizers/*'   // Kraken2 only
    path taxdump_nodes                                // nodes.dmp for the prokaryote filter, else empty
    path build_feature_table_script

    output:
    // The prefix names the classifier and the rank: kraken2_G or krakenuniq_G.
    path "${params.classifier}_${params.rank}_counts.tsv", emit: raw_counts
    path "${params.classifier}_${params.rank}_relative_abundance.tsv", emit: raw_relative
    path "${params.classifier}_${params.rank}_features.tsv", emit: raw_features
    path "${params.classifier}_${params.rank}_kmer_filter.tsv", emit: kmer_filter
    path "${params.classifier}_${params.rank}_${params.taxon_filter}_counts.tsv", emit: counts
    path "${params.classifier}_${params.rank}_${params.taxon_filter}_relative_abundance.tsv", emit: relative
    path "${params.classifier}_${params.rank}_${params.taxon_filter}_features.tsv", emit: features
    path "${params.classifier}_${params.rank}_${params.taxon_filter}_filter_qc.tsv", emit: qc
    path "${params.classifier}_${params.rank}_${params.taxon_filter}_sample_qc.tsv", emit: sample_qc
    path "build_feature_table_versions.yml", emit: versions

    script:
    def krakenUniq = params.classifier == 'krakenuniq'
    def minDistinct = (krakenUniq ? params.krakenuniq_min_kmers : params.kraken2_min_distinct_minimizers) as int
    def taxdumpArg = params.taxon_filter == 'prokaryote' ? "--taxdump \"${taxdump_nodes}\"" : ''
    // Options taking several files go last, so none of them swallows the reports.
    def minimizerArgs = !krakenUniq && minDistinct > 0 ? '--minimizer-reports minimizers/*' : ''
    """
    python3 "${build_feature_table_script}" \\
        --metadata "${metadata}" \\
        --prefix ${params.classifier}_${params.rank} \\
        --rank ${params.rank} \\
        --taxon-filter ${params.taxon_filter} \\
        --drop-taxids '${params.drop_taxids}' \\
        --min-distinct-kmers ${minDistinct} \\
        ${taxdumpArg} \\
        reports/* \\
        ${minimizerArgs}

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > build_feature_table_versions.yml
    """
}
