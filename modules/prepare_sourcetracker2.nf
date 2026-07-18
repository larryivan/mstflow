process PREPARE_SOURCETRACKER2_INPUTS {
    label 'process_low'
    tag "sourcetracker2_inputs:${params.get('sourcetracker2_analysis_label', 'default')}"

    publishDir "${params.outdir}/source_attribution/sourcetracker2/${params.get('sourcetracker2_analysis_label', 'default')}", mode: 'copy'

    input:
    path counts
    path metadata
    path feature_sample_qc
    path prepare_sourcetracker2_script

    output:
    path "sourcetracker2_counts.tsv", emit: counts
    path "sourcetracker2_map.tsv", emit: map
    path "sourcetracker2_input_summary.tsv", emit: summary
    path "sourcetracker2_excluded_samples.tsv", emit: exclusions
    path "sourcetracker2_sample_qc.tsv", emit: sample_qc
    path "prepare_sourcetracker2_versions.yml", emit: versions

    script:
    def aggregateBy = params.get('sourcetracker2_aggregate_by', 'sample_id')
    def analysisLabel = params.get('sourcetracker2_analysis_label', 'default')
    def taxonScope = params.get('sourcetracker2_taxon_scope', 'filtered')
    def includeRoles = params.get('sourcetracker2_include_panel_roles', '').toString().trim()
    def excludeRoles = params.get('sourcetracker2_exclude_panel_roles', 'control,holdout').toString().trim()
    def excludeEnvs = params.get('sourcetracker2_exclude_source_envs', '').toString().trim()
    def excludeClasses = params.get('sourcetracker2_exclude_source_classes', '').toString().trim()
    def includeRolesArg = includeRoles ? "--include-panel-role '${includeRoles}'" : ''
    def excludeRolesArg = excludeRoles ? "--exclude-panel-role '${excludeRoles}'" : ''
    def excludeEnvsArg = excludeEnvs ? "--exclude-source-env '${excludeEnvs}'" : ''
    def excludeClassesArg = excludeClasses ? "--exclude-source-class '${excludeClasses}'" : ''
    def looArg = params.get('sourcetracker2_loo', false).toString().toBoolean() ? '--loo' : ''
    def minClassified = params.get('sourcetracker2_min_kraken_classified_fraction', 0.001)
    def minRetainedReads = params.get('sourcetracker2_min_taxon_filter_retained_reads', 1000)
    def minRetainedFraction = params.get('sourcetracker2_min_taxon_filter_retained_fraction', 0.001)
    """
    python3 ${prepare_sourcetracker2_script} \\
        --counts ${counts} \\
        --metadata ${metadata} \\
        --feature-sample-qc ${feature_sample_qc} \\
        --output-counts sourcetracker2_counts.tsv \\
        --output-map sourcetracker2_map.tsv \\
        --summary sourcetracker2_input_summary.tsv \\
        --exclusions sourcetracker2_excluded_samples.tsv \\
        --sample-qc sourcetracker2_sample_qc.tsv \\
        --aggregate-by '${aggregateBy}' \\
        --analysis-label '${analysisLabel}' \\
        --random-seed ${params.sourcetracker2_seed} \\
        --taxon-scope '${taxonScope}' \\
        --source-rarefaction-depth ${params.sourcetracker2_source_rarefaction_depth} \\
        --sink-rarefaction-depth ${params.sourcetracker2_sink_rarefaction_depth} \\
        --min-kraken-classified-fraction ${minClassified} \\
        --min-taxon-filter-retained-reads ${minRetainedReads} \\
        --min-taxon-filter-retained-fraction ${minRetainedFraction} \\
        ${includeRolesArg} \\
        ${excludeRolesArg} \\
        ${excludeEnvsArg} \\
        ${excludeClassesArg} \\
        ${looArg}

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > prepare_sourcetracker2_versions.yml
    """
}
