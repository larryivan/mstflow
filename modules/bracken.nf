process BRACKEN_PROFILE {
    label 'process_medium'
    tag "${meta.id}"

    publishDir "${params.outdir}/bracken", mode: 'copy'

    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path("*.bracken.tsv"), path("*.bracken.report.txt"), emit: tables
    path "versions.yml", emit: versions

    script:
    def id = meta.id
    def extraArgs = params.bracken_extra_args ?: ''
    """
    bracken \\
        -d ${file(params.bracken_db)} \\
        -i ${report} \\
        -o ${id}.bracken.tsv \\
        -w ${id}.bracken.report.txt \\
        -r ${params.bracken_read_len} \\
        -l ${params.bracken_rank} \\
        -t ${params.bracken_threshold} \\
        ${extraArgs}

    printf '"${task.process}":\\n    bracken: "%s"\\n' "\$(bracken -v 2>&1 | head -n 1 | sed 's/^.*Bracken version //')" > versions.yml
    """

    stub:
    def id = meta.id
    """
    cat > ${id}.bracken.tsv <<-EOF
    name	taxonomy_id	taxonomy_lvl	kraken_assigned_reads	added_reads	new_est_reads	fraction_total_reads
    Bacteroides	816	G	80	20	100	0.6667
    Pseudomonas	286	G	40	10	50	0.3333
    EOF
    cp ${id}.bracken.tsv ${id}.bracken.report.txt
    printf '"${task.process}":\\n    stub: "true"\\n' > versions.yml
    """
}
