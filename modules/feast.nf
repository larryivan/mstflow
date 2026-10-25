/*
 * FEAST (Shenhav et al. 2019): an expectation-maximisation mixture of the source
 * environments and an unknown source. FEAST is licensed for non-commercial use only.
 */
process FEAST {
    label 'process_low'
    tag "${meta.analysis}"
    container "${params.container_attribution}"
    publishDir "${params.outdir}/attribution/${params.attribution_label}/feast", mode: 'copy',
        saveAs: { name -> meta.analysis == 'main' ? name : null }

    input:
    tuple val(meta), path(counts), path(mapping)
    path run_script
    path attribution_helpers   // mstflow_attribution.R, and .py for the stub

    output:
    tuple val(meta), path("${meta.analysis}.feast.attribution.tsv"), emit: results
    path "${meta.analysis}.feast.versions.yml", emit: versions

    script:
    """
    Rscript ${run_script} \\
        --counts ${counts} \\
        --map ${mapping} \\
        --analysis ${meta.analysis} \\
        --iterations ${params.feast_em_iterations} \\
        --coverage ${params.feast_coverage} \\
        --seed ${params.attribution_seed} \\
        --output ${meta.analysis}.feast.attribution.tsv

    printf '"${task.process}":\\n    FEAST: "%s"\\n    R: "%s"\\n' \\
        "\$(cat /usr/local/share/licenses/FEAST/COMMIT)" "\$(Rscript -e 'cat(R.version\$major, R.version\$minor, sep=".")')" \\
        > ${meta.analysis}.feast.versions.yml
    """

    stub:
    """
    python3 - <<'PY'
from mstflow_attribution import UNKNOWN, read_map, roles, write_results
envs, sinks = roles(read_map("${mapping}"))
share = 1.0 / (len(envs) + 1)
write_results("${meta.analysis}.feast.attribution.tsv", "${meta.analysis}", "feast", "proportion",
              {sink: {source: (share, None) for source in [*envs, UNKNOWN]} for sink in sinks})
PY
    printf '"${task.process}":\\n    FEAST: "stub"\\n' > ${meta.analysis}.feast.versions.yml
    """
}
