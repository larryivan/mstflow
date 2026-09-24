/*
 * KrakenUniq (Breitwieser et al. 2018): exact k-mer classification that also counts,
 * for every taxon, the unique k-mers its reads cover. Its databases may be larger
 * than memory; --krakenuniq_preload_size loads them in chunks of that size.
 */
process KRAKENUNIQ_PROFILE {
    label 'process_high'
    tag "${meta.id}"
    container "${params.container_krakenuniq}"
    publishDir "${params.outdir}/krakenuniq", mode: 'copy'

    input:
    tuple val(meta), path(reads), path(read_length_stats), path(input_mode), path(database)

    output:
    tuple val(meta), path("${meta.id}.krakenuniq.report.txt"), emit: reports
    path "${meta.id}.krakenuniq.fragment_qc.tsv", emit: qc
    path "${meta.id}.krakenuniq.versions.yml", emit: versions

    script:
    def readList = (reads instanceof List ? reads : [reads]).collect { read -> read.toString() }
    def pairedR1 = readList.findAll { read -> read.contains('.kraken_input_R1.') }
    def pairedR2 = readList.findAll { read -> read.contains('.kraken_input_R2.') }
    def orderedReads = readList.size() == 2 ? [pairedR1[0], pairedR2[0]] : readList
    def preload = params.krakenuniq_preload_size ? "--preload-size ${params.krakenuniq_preload_size}" : '--preload'
    """
    case "\$(tr -d '\\r\\n' < ${input_mode})" in
        paired) paired='--paired' ;;
        single) paired='' ;;
        *) echo "Invalid prepared input mode in ${input_mode}" >&2; exit 1 ;;
    esac
    krakenuniq \\
        --db ${database} \\
        --threads ${task.cpus} \\
        ${preload} \\
        \$paired \\
        --report-file ${meta.id}.krakenuniq.report.txt \\
        --output off \\
        ${orderedReads.collect { read -> "\"${read}\"" }.join(' ')}

    # KrakenUniq's wrapper exits 0 even when its classifier is killed, typically for
    # memory: the database index alone (database.idx) is loaded whole. A report
    # without an unclassified or root row means the classification never finished.
    if ! awk -F'\t' '\$7 == "0" || \$7 == "1" { found = 1 } END { exit !found }' ${meta.id}.krakenuniq.report.txt; then
        echo "KrakenUniq wrote no results for ${meta.id}. It most likely ran out of memory: the database index is loaded whole, and the database too unless --krakenuniq_preload_size is set. Raise --kraken2_memory." >&2
        exit 1
    fi

    cp ${read_length_stats} ${meta.id}.krakenuniq.fragment_qc.tsv
    printf '"${task.process}":\\n    krakenuniq: "%s"\\n' "\$(krakenuniq --version 2>&1 | awk 'NR == 1 { print \$NF }')" \\
        > ${meta.id}.krakenuniq.versions.yml
    """

    stub:
    """
    printf '%%\\treads\\ttaxReads\\tkmers\\tdup\\tcov\\ttaxID\\trank\\ttaxName\\n0\\t0\\t0\\t0\\t0\\tNA\\t0\\tno rank\\tunclassified\\n100\\t10\\t0\\t10\\t1\\tNA\\t1\\tno rank\\troot\\n100\\t10\\t10\\t10\\t1\\tNA\\t1301\\tgenus\\t  Streptococcus\\n' \\
        > ${meta.id}.krakenuniq.report.txt
    cp ${read_length_stats} ${meta.id}.krakenuniq.fragment_qc.tsv
    printf '"${task.process}":\\n    krakenuniq: "stub"\\n' > ${meta.id}.krakenuniq.versions.yml
    """
}
