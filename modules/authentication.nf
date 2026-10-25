// The most abundant taxa of each ancient library, and the reference sequences of each.
process SELECT_AUTHENTICATION_TARGETS {
    label 'process_low'
    tag 'authentication'
    container "${params.container_python_tools}"
    publishDir "${params.outdir}/authentication", mode: 'copy', pattern: 'targets.tsv'

    input:
    path counts
    path metadata
    path genomes
    path nodes
    path select_script

    output:
    path "targets.tsv", emit: targets
    path "references.tsv", emit: references

    script:
    """
    python3 "${select_script}" \\
        --counts "${counts}" \\
        --metadata "${metadata}" \\
        --fasta "${genomes}" \\
        --nodes "${nodes}" \\
        --top ${params.authentication_top_taxa} \\
        --min-reads ${params.authentication_min_reads}
    """

    stub:
    """
    printf 'sample_id\\ttaxid\\tname\\ttable_reads\\thas_reference\\n' > targets.tsv
    printf 'seqid\\ttaxid\\n' > references.tsv
    """
}


// PyDamage (Borry et al. 2021) on the reads of each target taxon, with how much of
// its reference genomes they cover.
process AUTHENTICATE {
    label 'process_medium'
    tag "${meta.id}"
    container "${params.container_pydamage}"

    input:
    tuple val(meta), path(bam)
    path targets
    path references
    path authenticate_script

    output:
    path "${meta.id}.authentication.tsv", emit: results
    path "${meta.id}.pydamage.versions.yml", emit: versions

    script:
    // PyDamage compiles with Numba, which caches beside its code unless told otherwise;
    // inside a read-only image that fails, so both caches go to the task directory.
    """
    export NUMBA_CACHE_DIR=\$PWD/.cache MPLCONFIGDIR=\$PWD/.cache
    python3 "${authenticate_script}" \\
        --bam "${bam}" \\
        --sample "${meta.id}" \\
        --targets "${targets}" \\
        --references "${references}" \\
        --min-reads ${params.authentication_min_reads} \\
        --threads ${task.cpus} \\
        --output ${meta.id}.authentication.tsv

    printf '"${task.process}":\\n    pydamage: "%s"\\n' "\$(pydamage --version | awk '{ print \$NF }')" > ${meta.id}.pydamage.versions.yml
    """

    stub:
    """
    printf 'sample_id\\ttaxid\\tverdict\\n${meta.id}\\t1301\\tdamaged\\n' > ${meta.id}.authentication.tsv
    printf '"${task.process}":\\n    pydamage: "stub"\\n' > ${meta.id}.pydamage.versions.yml
    """
}
