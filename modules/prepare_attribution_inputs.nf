/*
 * Choose the analysis units for source attribution: drop excluded roles and
 * environments and units too shallow to use, merge runs of one biological sample,
 * and write the count table and map that every method reads.
 */
process PREPARE_ATTRIBUTION_INPUTS {
    label 'process_low'
    tag "${params.attribution_label}"
    container "${params.container_python_tools}"
    publishDir "${params.outdir}/attribution/${params.attribution_label}/inputs", mode: 'copy'

    input:
    path counts
    path metadata
    path feature_sample_qc
    path prepare_script

    output:
    path "attribution_counts.tsv", emit: counts
    path "attribution_map.tsv", emit: map
    path "attribution_input_summary.tsv", emit: summary
    path "attribution_excluded_samples.tsv", emit: exclusions
    path "attribution_sample_qc.tsv", emit: sample_qc
    path "prepare_attribution_inputs_versions.yml", emit: versions

    script:
    // SourceTracker2 cannot rarefy a unit shallower than its depths, so with it
    // those units are left out for every method, keeping the methods comparable.
    def withSourceTracker = 'sourcetracker2' in params.attribution_methods.toString().tokenize(',')*.trim()
    def minSource = withSourceTracker ? params.sourcetracker2_source_rarefaction_depth : 0
    def minSink = withSourceTracker ? params.sourcetracker2_sink_rarefaction_depth : 0
    def patternArgs = [
        'include-panel-role': params.attribution_include_panel_roles,
        'exclude-panel-role': params.attribution_exclude_panel_roles,
        'exclude-source-env': params.attribution_exclude_source_envs,
        'exclude-source-class': params.attribution_exclude_source_classes,
    ].findAll { _flag, value -> value?.toString()?.trim() }
        .collect { flag, value -> "--${flag} '${value}'" }
        .join(' ')
    """
    python3 ${prepare_script} \\
        --counts ${counts} \\
        --metadata ${metadata} \\
        --feature-sample-qc ${feature_sample_qc} \\
        --output-counts attribution_counts.tsv \\
        --output-map attribution_map.tsv \\
        --summary attribution_input_summary.tsv \\
        --exclusions attribution_excluded_samples.tsv \\
        --sample-qc attribution_sample_qc.tsv \\
        --aggregate-by '${params.attribution_aggregate_by}' \\
        --analysis-label '${params.attribution_label}' \\
        --taxon-scope '${params.attribution_taxon_scope}' \\
        --min-source-count ${minSource} \\
        --min-sink-count ${minSink} \\
        --min-kraken-classified-fraction ${params.attribution_min_classified_fraction} \\
        --min-taxon-filter-retained-reads ${params.attribution_min_retained_reads} \\
        --min-taxon-filter-retained-fraction ${params.attribution_min_retained_fraction} \\
        ${patternArgs}

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > prepare_attribution_inputs_versions.yml
    """
}
