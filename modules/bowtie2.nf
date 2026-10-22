/*
 * Bowtie2 (Langmead and Salzberg 2012) for authentication: the reads of an ancient
 * library are aligned end to end to the reference genomes of the taxa to authenticate.
 */
process BOWTIE2_BUILD {
    label 'process_medium'
    tag "${fasta.name}"
    container "${params.container_bowtie2}"

    input:
    path fasta

    output:
    path "${fasta.name}.*.bt2*", emit: index

    script:
    """
    bowtie2-build --threads ${task.cpus} ${fasta} ${fasta.name}
    """

    stub:
    """
    touch ${fasta.name}.{1,2,3,4,rev.1,rev.2}.bt2
    """
}


// Every read file is aligned as single reads: damage shows as C-to-T at the 5' end of
// any read, whichever mate it is. Only reads placed uniquely (MAPQ 30 or more) are
// kept, so a read shared by related genomes counts for none of them.
process BOWTIE2_ALIGN {
    label 'process_medium'
    tag "${meta.id}"
    container "${params.container_bowtie2}"
    publishDir "${params.outdir}/authentication/bowtie2", mode: 'copy'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: bam
    path "${meta.id}.bowtie2.log", emit: log
    path "${meta.id}.bowtie2.versions.yml", emit: versions

    script:
    def prefix = index[0].name.replaceFirst(/\.(rev\.)?\d\.bt2l?$/, '')
    """
    bowtie2 --very-sensitive --threads ${task.cpus} -x ${prefix} -U ${[reads].flatten().join(',')} 2> ${meta.id}.bowtie2.log \\
        | samtools view -b -q 30 -F 4 \\
        | samtools sort -o ${meta.id}.bam

    printf '"${task.process}":\\n    bowtie2: "%s"\\n    samtools: "%s"\\n' \\
        "\$(bowtie2 --version | awk 'NR == 1 { print \$NF }')" "\$(samtools --version | awk 'NR == 1 { print \$2 }')" \\
        > ${meta.id}.bowtie2.versions.yml
    """

    stub:
    """
    touch ${meta.id}.bam ${meta.id}.bowtie2.log
    printf '"${task.process}":\\n    bowtie2: "stub"\\n' > ${meta.id}.bowtie2.versions.yml
    """
}
