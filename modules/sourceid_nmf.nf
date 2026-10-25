/*
 * SourceID-NMF (Huang et al. 2024): a weighted non-negative matrix factorisation of
 * each sink into the source profiles and one unknown profile. Its repository has no
 * licence, so its code is not distributed with MSTFlow: `code` is staged from the
 * authors' repository at the commit set by --sourceid_nmf_code.
 */
process SOURCEID_NMF {
    label 'process_medium'
    tag "${meta.analysis}"
    container "${params.container_attribution}"
    publishDir "${params.outdir}/attribution/${params.attribution_label}/sourceid_nmf", mode: 'copy',
        saveAs: { name -> meta.analysis == 'main' ? name : null }

    input:
    tuple val(meta), path(counts), path(mapping)
    path code, stageAs: 'sourceid_nmf/*'
    path run_script
    path attribution_module

    output:
    tuple val(meta), path("${meta.analysis}.sourceid_nmf.attribution.tsv"), emit: results
    path "${meta.analysis}.sourceid_nmf.versions.yml", emit: versions

    script:
    """
    python3 ${run_script} \\
        --counts ${counts} \\
        --map ${mapping} \\
        --analysis ${meta.analysis} \\
        --code-dir sourceid_nmf \\
        --iterations ${params.sourceid_nmf_iterations} \\
        --rho ${params.sourceid_nmf_rho} \\
        --weight ${params.sourceid_nmf_weight} \\
        --tolerance ${params.sourceid_nmf_tolerance} \\
        --threads ${task.cpus} \\
        --output ${meta.analysis}.sourceid_nmf.attribution.tsv

    printf '"${task.process}":\\n    SourceID-NMF: "%s"\\n' "${params.sourceid_nmf_code}" > ${meta.analysis}.sourceid_nmf.versions.yml
    """

    stub:
    """
    python3 - <<'PY'
from mstflow_attribution import UNKNOWN, read_map, roles, write_results
envs, sinks = roles(read_map("${mapping}"))
share = 1.0 / (len(envs) + 1)
write_results("${meta.analysis}.sourceid_nmf.attribution.tsv", "${meta.analysis}", "sourceid_nmf", "proportion",
              {sink: {source: (share, None) for source in [*envs, UNKNOWN]} for sink in sinks})
PY
    printf '"${task.process}":\\n    SourceID-NMF: "stub"\\n' > ${meta.analysis}.sourceid_nmf.versions.yml
    """
}
