/*
 * SourcePredict (Borry et al. 2019) needs the NCBI taxonomy as an ete3 database to
 * compute UniFrac distances. It is built once from the same taxdump as the
 * classifier database, so the two never disagree about a taxon.
 */
process PREPARE_SOURCEPREDICT_TAXONOMY {
    label 'process_low'
    tag "ete3"
    container "${params.container_sourcepredict}"

    input:
    path dumps, stageAs: 'taxdump/*'

    output:
    path "taxa.sqlite*", emit: database

    script:
    """
    [[ -e taxdump/merged.dmp ]] || : > taxdump/merged.dmp
    tar -czhf taxdump.tar.gz -C taxdump nodes.dmp names.dmp merged.dmp
    export HOME=\$PWD
    python3 - <<'PY'
import warnings
warnings.filterwarnings("ignore")
from ete3 import NCBITaxa
NCBITaxa(dbfile="taxa.sqlite", taxdump_file="taxdump.tar.gz", update=False)
PY
    """

    stub:
    """
    touch taxa.sqlite taxa.sqlite.traverse.pkl
    """
}

/*
 * SourcePredict: a k-nearest-neighbour classifier on a UniFrac embedding of the
 * sources, trained on every source library separately. Its estimates are class
 * probabilities, not shares of the sink.
 */
process SOURCEPREDICT {
    label 'process_medium'
    tag "${meta.analysis}"
    container "${params.container_sourcepredict}"
    publishDir "${params.outdir}/attribution/${params.attribution_label}/sourcepredict", mode: 'copy',
        saveAs: { name -> meta.analysis == 'main' ? name : null }

    input:
    tuple val(meta), path(counts), path(mapping)
    path taxonomy, stageAs: '.etetoolkit/*'
    path run_script
    path attribution_module

    output:
    tuple val(meta), path("${meta.analysis}.sourcepredict.attribution.tsv"), emit: results
    path "${meta.analysis}.sourcepredict.versions.yml", emit: versions

    script:
    def neighbours = params.sourcepredict_neighbors as int
    """
    export HOME=\$PWD NUMBA_CACHE_DIR=\$PWD/numba_cache MPLCONFIGDIR=\$PWD/mpl PYTHONWARNINGS=ignore
    python3 ${run_script} \\
        --counts ${counts} \\
        --map ${mapping} \\
        --analysis ${meta.analysis} \\
        --rank ${params.bracken_rank} \\
        --neighbors ${neighbours} \\
        --unknown-alpha ${params.sourcepredict_unknown_alpha} \\
        --normalization ${params.sourcepredict_normalization} \\
        --distance ${params.sourcepredict_distance} \\
        --embedding ${params.sourcepredict_embedding} \\
        --seed ${params.attribution_seed} \\
        --threads ${task.cpus} \\
        ${meta.analysis == 'main' ? '' : '--drop-lonely-envs'} \\
        --output ${meta.analysis}.sourcepredict.attribution.tsv

    printf '"${task.process}":\\n    sourcepredict: "%s"\\n' \\
        "\$(python3 -c 'import sourcepredict; print(sourcepredict.__version__)' 2>/dev/null)" \\
        > ${meta.analysis}.sourcepredict.versions.yml
    """

    stub:
    """
    python3 - <<'PY'
from mstflow_attribution import UNKNOWN, read_map, roles, write_results
envs, sinks = roles(read_map("${mapping}"))
share = 1.0 / (len(envs) + 1)
write_results("${meta.analysis}.sourcepredict.attribution.tsv", "${meta.analysis}", "sourcepredict", "probability",
              {sink: {source: (share, None) for source in [*envs, UNKNOWN]} for sink in sinks})
PY
    printf '"${task.process}":\\n    sourcepredict: "stub"\\n' > ${meta.analysis}.sourcepredict.versions.yml
    """
}
