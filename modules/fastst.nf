/*
 * FastST (Choi et al. 2025): a generalised-least-squares mixture of the source
 * environments and an unknown source, with standard errors. Experimental: it fits
 * only the taxa present in a sink, which misleads it on sinks with few taxa.
 */
process FASTST {
    label 'process_low'
    tag "${meta.analysis}"
    container "${params.container_attribution}"
    publishDir "${params.outdir}/attribution/${params.attribution_label}/fastst", mode: 'copy',
        saveAs: { name -> meta.analysis == 'main' ? name : null }

    input:
    tuple val(meta), path(counts), path(mapping)
    path run_script
    path attribution_helpers   // mstflow_attribution.R, and .py for the stub

    output:
    tuple val(meta), path("${meta.analysis}.fastst.attribution.tsv"), emit: results
    path "${meta.analysis}.fastst.versions.yml", emit: versions

    script:
    """
    Rscript ${run_script} \\
        --counts ${counts} \\
        --map ${mapping} \\
        --analysis ${meta.analysis} \\
        --output ${meta.analysis}.fastst.attribution.tsv

    printf '"${task.process}":\\n    FastST: "%s"\\n' "\$(cat /usr/local/share/licenses/FastST/COMMIT)" \\
        > ${meta.analysis}.fastst.versions.yml
    """

    stub:
    """
    python3 - <<'PY'
from mstflow_attribution import UNKNOWN, read_map, roles, write_results
envs, sinks = roles(read_map("${mapping}"))
share = 1.0 / (len(envs) + 1)
write_results("${meta.analysis}.fastst.attribution.tsv", "${meta.analysis}", "fastst", "proportion",
              {sink: {source: (share, None) for source in [*envs, UNKNOWN]} for sink in sinks})
PY
    printf '"${task.process}":\\n    FastST: "stub"\\n' > ${meta.analysis}.fastst.versions.yml
    """
}
