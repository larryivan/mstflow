process HOST_REMOVE_BWA {
    label 'process_high'
    tag "${meta.id}"

    publishDir "${params.outdir}/host_removed", mode: 'copy'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.dehosted*.fq.gz"), emit: reads
    path "*.host_removal.log", emit: logs
    path "versions.yml", emit: versions

    script:
    def id = meta.id
    def readList = reads.collect { it.toString() }
    def sampleHostRef = meta.host_ref ?: 'default'
    def configuredHostRef = (sampleHostRef in ['default', 'NA', '']) ? params.host_fasta : sampleHostRef
    if (!configuredHostRef) {
        error "HOST_REMOVE_BWA requires params.host_fasta or a per-sample host_ref path for ${id}"
    }
    def hostRef = file(configuredHostRef).toString()
    if (meta.age_group == 'ancient') {
        """
        bwa aln \\
            -l 1024 \\
            -n 0.01 \\
            -o 2 \\
            -t ${task.cpus} \\
            ${hostRef} \\
            ${readList[0]} \\
            > ${id}.host.sai \\
            2> ${id}.host_removal.log

        bwa samse \\
            ${hostRef} \\
            ${id}.host.sai \\
            ${readList[0]} \\
            2>> ${id}.host_removal.log \\
            | samtools view -@ ${task.cpus} -b -f 4 -o ${id}.unmapped.bam -

        samtools fastq -@ ${task.cpus} ${id}.unmapped.bam \\
            | gzip -c > ${id}.dehosted.fq.gz

        printf '"${task.process}":\\n    bwa: "%s"\\n    samtools: "%s"\\n' "\$(bwa 2>&1 | grep '^Version' | awk '{print \$2}')" "\$(samtools --version | head -n 1 | awk '{print \$2}')" > versions.yml
        """
    } else if (readList.size() == 2) {
        """
        bwa mem \\
            -t ${task.cpus} \\
            ${hostRef} \\
            ${readList[0]} \\
            ${readList[1]} \\
            2> ${id}.host_removal.log \\
            | samtools view -@ ${task.cpus} -b -f 12 -F 2304 -o ${id}.both_unmapped.bam -

        samtools fastq \\
            -@ ${task.cpus} \\
            -1 ${id}.dehosted_R1.fq.gz \\
            -2 ${id}.dehosted_R2.fq.gz \\
            -0 /dev/null \\
            -s /dev/null \\
            -n \\
            ${id}.both_unmapped.bam

        printf '"${task.process}":\\n    bwa: "%s"\\n    samtools: "%s"\\n' "\$(bwa 2>&1 | grep '^Version' | awk '{print \$2}')" "\$(samtools --version | head -n 1 | awk '{print \$2}')" > versions.yml
        """
    } else {
        """
        bwa mem \\
            -t ${task.cpus} \\
            ${hostRef} \\
            ${readList[0]} \\
            2> ${id}.host_removal.log \\
            | samtools view -@ ${task.cpus} -b -f 4 -o ${id}.unmapped.bam -

        samtools fastq -@ ${task.cpus} ${id}.unmapped.bam \\
            | gzip -c > ${id}.dehosted.fq.gz

        printf '"${task.process}":\\n    bwa: "%s"\\n    samtools: "%s"\\n' "\$(bwa 2>&1 | grep '^Version' | awk '{print \$2}')" "\$(samtools --version | head -n 1 | awk '{print \$2}')" > versions.yml
        """
    }

    stub:
    def id = meta.id
    def readList = reads.collect { it.toString() }
    def stubReads = (meta.age_group == 'modern' && readList.size() == 2) ? """
        printf '@${id}/1\\nACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.dehosted_R1.fq.gz
        printf '@${id}/2\\nACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.dehosted_R2.fq.gz
        """ : """
        printf '@${id}\\nACGTACGTACGTACGTACGTACGTACGTACGT\\n+\\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\\n' | gzip -c > ${id}.dehosted.fq.gz
        """
    """
    ${stubReads}
    echo 'stub host removal' > ${id}.host_removal.log
    printf '"${task.process}":\\n    stub: "true"\\n' > versions.yml
    """
}
