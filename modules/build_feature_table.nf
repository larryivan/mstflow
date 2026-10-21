/*
 * Merge the per-sample Bracken estimates (or, without Bracken, the classifier's clade
 * counts) into one count table, before and after the taxon filter, and record for
 * every sample how many reads each step kept. Taxa backed by too few distinct
 * minimizers (Kraken2) or unique k-mers (KrakenUniq) are removed sample by sample
 * first, if asked.
 */
process BUILD_FEATURE_TABLE {
    label 'process_low'
    tag "${params.bracken_rank}:${params.taxon_filter}"
    container "${params.container_python_tools}"
    publishDir "${params.outdir}/tables", mode: 'copy'

    input:
    path metadata
    path bracken_files, stageAs: 'bracken/*'          // empty without Bracken
    path reports, stageAs: 'reports/*'                // Kraken2 or KrakenUniq reports
    path minimizer_reports, stageAs: 'minimizers/*'   // Kraken2 only
    path taxdump_nodes
    path build_feature_table_script

    output:
    // The prefix names what made the counts: bracken_, kraken2_ or krakenuniq_<rank>.
    path "*_${params.bracken_rank}_counts.tsv", emit: raw_counts
    path "*_${params.bracken_rank}_relative_abundance.tsv", emit: raw_relative
    path "*_${params.bracken_rank}_features.tsv", emit: raw_features
    path "*_${params.bracken_rank}_kmer_filter.tsv", emit: kmer_filter
    path "*_${params.bracken_rank}_${params.taxon_filter}_counts.tsv", emit: counts
    path "*_${params.bracken_rank}_${params.taxon_filter}_relative_abundance.tsv", emit: relative
    path "*_${params.bracken_rank}_${params.taxon_filter}_features.tsv", emit: features
    path "*_${params.bracken_rank}_${params.taxon_filter}_filter_qc.tsv", emit: qc
    path "*_${params.bracken_rank}_${params.taxon_filter}_sample_qc.tsv", emit: sample_qc
    path "build_feature_table_versions.yml", emit: versions

    script:
    def krakenUniq = params.classifier == 'krakenuniq'
    def bracken = params.run_bracken.toString().toBoolean()
    def prefix = "${bracken ? 'bracken' : params.classifier}_${params.bracken_rank}"
    def minDistinct = (krakenUniq ? params.krakenuniq_min_kmers : params.kraken2_min_distinct_minimizers) as int
    def minReads = (krakenUniq ? params.krakenuniq_min_reads : 0) as int
    def taxdumpArg = params.taxon_filter == 'prokaryote' ? "--taxdump \"${taxdump_nodes}\"" : ''
    // Options taking several files go last, so none of them swallows another's files.
    def brackenArgs = bracken ? 'bracken/*' : ''
    def minimizerArgs = !krakenUniq && minDistinct > 0 ? '--minimizer-reports minimizers/*' : ''
    """
    python3 "${build_feature_table_script}" \\
        --metadata "${metadata}" \\
        --prefix ${prefix} \\
        --rank ${params.bracken_rank} \\
        --taxon-filter ${params.taxon_filter} \\
        --drop-taxids '${params.drop_taxids}' \\
        --min-distinct-kmers ${minDistinct} \\
        --min-taxon-reads ${minReads} \\
        ${taxdumpArg} \\
        ${brackenArgs} \\
        --kraken-reports reports/* \\
        ${minimizerArgs}

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > build_feature_table_versions.yml
    """
}

