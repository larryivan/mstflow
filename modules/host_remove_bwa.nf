process HOST_ALIGN_BWA {
    label 'process_high'
    tag "${meta.id}"
    container "${params.container_bwa}"

    publishDir "${params.outdir}/host_removed", mode: 'copy',
        saveAs: { filename -> filename == 'host_alignments' ? null : filename }

    input:
    tuple val(meta), path(reads), path(host_bundle)

    output:
    tuple val(meta), path("host_alignments"), emit: alignments
    path "*.host_removal.log", emit: logs
    path "*.host_remove_bwa.versions.yml", emit: versions

    script:
    def id = meta.id
    def readList = (reads instanceof List ? reads : [reads]).collect { read -> read.toString() }
    def hostRef = host_bundle[0]
    def ancient = meta.age_group == 'ancient'
    // Short, damaged ancient reads go through bwa aln with seeding disabled; modern
    // reads through bwa mem. mode.txt tells FILTER_HOST_ALIGNMENTS which files exist.
    def aln = "bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} \"${hostRef}\""
    def mode
    def commands
    if (ancient && readList.size() == 3) {
        def se = readList.find { read -> read.contains('.prehost_se.') }
        def r1 = readList.find { read -> read.contains('.prehost_R1.') }
        def r2 = readList.find { read -> read.contains('.prehost_R2.') }
        mode = 'ancient_mixed'
        commands = """
        ${aln} "${se}" > se.sai
        bwa samse "${hostRef}" se.sai "${se}" | gzip -1 -c > host_alignments/se.sam.gz
        ${aln} "${r1}" > R1.sai
        ${aln} "${r2}" > R2.sai
        bwa sampe "${hostRef}" R1.sai R2.sai "${r1}" "${r2}" | gzip -1 -c > host_alignments/pair.sam.gz
        """
    } else if (ancient) {
        mode = 'ancient_se'
        commands = """
        ${aln} "${readList[0]}" > single.sai
        bwa samse "${hostRef}" single.sai "${readList[0]}" | gzip -1 -c > host_alignments/single.sam.gz
        """
    } else if (readList.size() == 2) {
        mode = 'modern_pe'
        commands = """
        bwa mem -t ${task.cpus} "${hostRef}" "${readList[0]}" "${readList[1]}" | gzip -1 -c > host_alignments/pair.sam.gz
        """
    } else {
        mode = 'modern_se'
        commands = """
        bwa mem -t ${task.cpus} "${hostRef}" "${readList[0]}" | gzip -1 -c > host_alignments/single.sam.gz
        """
    }
    """
    mkdir -p host_alignments
    printf '${mode}\\n' > host_alignments/mode.txt
    {
    ${commands}
    } 2> "${id}.host_removal.log"

    version=\$(bwa 2>&1 | awk '/^Version:/ { print \$2; exit }' || true)
    printf '"${task.process}":\\n    bwa: "%s"\\n' "\$version" > "${id}.host_remove_bwa.versions.yml"
    """

    stub:
    def id = meta.id
    def readList = reads instanceof List ? reads : [reads]
    def ancient = meta.age_group == 'ancient'
    def mode = ancient ? (readList.size() == 3 ? 'ancient_mixed' : 'ancient_se') : (readList.size() == 2 ? 'modern_pe' : 'modern_se')
    def samFiles = [ancient_mixed: 'se pair', ancient_se: 'single', modern_pe: 'pair', modern_se: 'single'][mode]
    """
    mkdir -p host_alignments
    printf '${mode}\\n' > host_alignments/mode.txt
    for name in ${samFiles}; do touch host_alignments/\$name.sam.gz; done
    echo 'stub host alignment' > ${id}.host_removal.log
    printf '"${task.process}":\\n    stub: "true"\\n' > ${id}.host_remove_bwa.versions.yml
    """
}
