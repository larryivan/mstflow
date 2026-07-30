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
    def hostFiles = host_bundle instanceof List ? host_bundle.collect { hostFile -> hostFile.toString() } : [host_bundle.toString()]
    if (!hostFiles) {
        error "HOST_ALIGN_BWA requires a staged host FASTA/index bundle for ${id}"
    }
    def hostRef = hostFiles[0]

    if (meta.age_group == 'ancient' && readList.size() == 3) {
        def seRead = readList.find { read -> read.contains('.prehost_se.') }
        def r1Read = readList.find { read -> read.contains('.prehost_R1.') }
        def r2Read = readList.find { read -> read.contains('.prehost_R2.') }
        if (!seRead || !r1Read || !r2Read) {
            error "Ancient PE input for ${id} must contain prehost_se, prehost_R1, and prehost_R2 files"
        }
        """
        mkdir -p host_alignments
        printf 'ancient_mixed\n' > host_alignments/mode.txt
        : > "${id}.host_removal.log"

        bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} \
            "${hostRef}" "${seRead}" > "${id}.host_se.sai" 2>> "${id}.host_removal.log"
        bwa samse "${hostRef}" "${id}.host_se.sai" "${seRead}" \
            2>> "${id}.host_removal.log" | gzip -1 -c > host_alignments/se.sam.gz

        bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} \
            "${hostRef}" "${r1Read}" > "${id}.host_R1.sai" 2>> "${id}.host_removal.log"
        bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} \
            "${hostRef}" "${r2Read}" > "${id}.host_R2.sai" 2>> "${id}.host_removal.log"
        bwa sampe "${hostRef}" "${id}.host_R1.sai" "${id}.host_R2.sai" "${r1Read}" "${r2Read}" \
            2>> "${id}.host_removal.log" | gzip -1 -c > host_alignments/pair.sam.gz

        version=\$(bwa 2>&1 | awk '/^Version:/ { print \$2; exit }' || true)
        [[ -n "\$version" ]] || version="unknown"
        printf '"${task.process}":\n    bwa: "%s"\n' "\$version" > "${id}.host_remove_bwa.versions.yml"
        """
    } else if (meta.age_group == 'ancient' && readList.size() == 2) {
        """
        mkdir -p host_alignments
        printf 'ancient_pe\n' > host_alignments/mode.txt
        : > "${id}.host_removal.log"
        bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} \
            "${hostRef}" "${readList[0]}" > "${id}.host_R1.sai" 2>> "${id}.host_removal.log"
        bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} \
            "${hostRef}" "${readList[1]}" > "${id}.host_R2.sai" 2>> "${id}.host_removal.log"
        bwa sampe "${hostRef}" "${id}.host_R1.sai" "${id}.host_R2.sai" "${readList[0]}" "${readList[1]}" \
            2>> "${id}.host_removal.log" | gzip -1 -c > host_alignments/pair.sam.gz
        version=\$(bwa 2>&1 | awk '/^Version:/ { print \$2; exit }' || true)
        [[ -n "\$version" ]] || version="unknown"
        printf '"${task.process}":\n    bwa: "%s"\n' "\$version" > "${id}.host_remove_bwa.versions.yml"
        """
    } else if (readList.size() == 2) {
        """
        mkdir -p host_alignments
        printf 'modern_pe\n' > host_alignments/mode.txt
        : > "${id}.host_removal.log"
        bwa mem -t ${task.cpus} "${hostRef}" "${readList[0]}" "${readList[1]}" \
            2>> "${id}.host_removal.log" | gzip -1 -c > host_alignments/pair.sam.gz
        version=\$(bwa 2>&1 | awk '/^Version:/ { print \$2; exit }' || true)
        [[ -n "\$version" ]] || version="unknown"
        printf '"${task.process}":\n    bwa: "%s"\n' "\$version" > "${id}.host_remove_bwa.versions.yml"
        """
    } else if (meta.age_group == 'ancient') {
        """
        mkdir -p host_alignments
        printf 'ancient_se\n' > host_alignments/mode.txt
        : > "${id}.host_removal.log"
        bwa aln -l 1024 -n 0.01 -o 2 -t ${task.cpus} \
            "${hostRef}" "${readList[0]}" > "${id}.host.sai" 2>> "${id}.host_removal.log"
        bwa samse "${hostRef}" "${id}.host.sai" "${readList[0]}" \
            2>> "${id}.host_removal.log" | gzip -1 -c > host_alignments/single.sam.gz
        version=\$(bwa 2>&1 | awk '/^Version:/ { print \$2; exit }' || true)
        [[ -n "\$version" ]] || version="unknown"
        printf '"${task.process}":\n    bwa: "%s"\n' "\$version" > "${id}.host_remove_bwa.versions.yml"
        """
    } else {
        """
        mkdir -p host_alignments
        printf 'modern_se\n' > host_alignments/mode.txt
        : > "${id}.host_removal.log"
        bwa mem -t ${task.cpus} "${hostRef}" "${readList[0]}" \
            2>> "${id}.host_removal.log" | gzip -1 -c > host_alignments/single.sam.gz
        version=\$(bwa 2>&1 | awk '/^Version:/ { print \$2; exit }' || true)
        [[ -n "\$version" ]] || version="unknown"
        printf '"${task.process}":\n    bwa: "%s"\n' "\$version" > "${id}.host_remove_bwa.versions.yml"
        """
    }

    stub:
    def id = meta.id
    def readList = (reads instanceof List ? reads : [reads]).collect { read -> read.toString() }
    def mode = meta.age_group == 'ancient' && readList.size() == 3 ? 'ancient_mixed' :
        (meta.age_group == 'ancient' && readList.size() == 2 ? 'ancient_pe' :
        (readList.size() == 2 ? 'modern_pe' :
        (meta.age_group == 'ancient' ? 'ancient_se' : 'modern_se')))
    def samFiles = mode == 'ancient_mixed' ? 'touch host_alignments/se.sam.gz host_alignments/pair.sam.gz' :
        (mode.endsWith('_pe') ? 'touch host_alignments/pair.sam.gz' : 'touch host_alignments/single.sam.gz')
    """
    mkdir -p host_alignments
    printf '${mode}\n' > host_alignments/mode.txt
    ${samFiles}
    echo 'stub host alignment' > ${id}.host_removal.log
    printf '"${task.process}":\n    stub: "true"\n' > ${id}.host_remove_bwa.versions.yml
    """
}
