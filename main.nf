#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VALIDATE_METADATA } from './modules/validate_metadata'
include { FASTQC as FASTQC_RAW } from './modules/fastqc'
include { FASTQC as FASTQC_CLEAN } from './modules/fastqc'
include { PREPROCESS_FASTP } from './modules/preprocess_fastp'
include { PREPROCESS_ADAPTERREMOVAL } from './modules/preprocess_adapterremoval'
include { HOST_REMOVE_BWA } from './modules/host_remove_bwa'
include { KRAKEN2_PROFILE } from './modules/kraken2'
include { BRACKEN_PROFILE } from './modules/bracken'
include { BUILD_FEATURE_TABLE } from './modules/build_feature_table'
include { PREPARE_SOURCETRACKER2_INPUTS } from './modules/prepare_sourcetracker2'
include { MAKE_BIOM } from './modules/make_biom'
include { SOURCETRACKER2_GIBBS } from './modules/sourcetracker2'
include { PLOT_SOURCETRACKER2 } from './modules/plot_sourcetracker2'
include { MULTIQC_REPORT } from './modules/multiqc'

def asBool(value) {
    if (value == null) {
        return false
    }
    if (value instanceof Boolean) {
        return value
    }
    return value.toString().toLowerCase() in ['true', '1', 'yes', 'y']
}

def readListFromRow(row) {
    def fastq1 = row.fastq_1.toString()
    def fastq2 = row.fastq_2.toString()
    if (fastq2 == 'NA' || fastq2 == '') {
        return [file(fastq1)]
    }
    return [file(fastq1), file(fastq2)]
}

def metaFromRow(row) {
    return [
        id             : row.sample_id.toString(),
        layout         : row.layout.toString(),
        age_group      : row.age_group.toString(),
        preprocess_mode: row.preprocess_mode.toString(),
        host_ref       : row.host_ref.toString()
    ]
}

workflow {
    if (asBool(params.run_bracken) && !asBool(params.run_kraken2)) {
        error 'run_bracken requires run_kraken2'
    }
    if (asBool(params.run_sourcetracker2) && !asBool(params.run_bracken)) {
        error 'run_sourcetracker2 requires run_bracken'
    }

    input_ch = Channel.fromPath(params.input, checkIfExists: true)

    VALIDATE_METADATA(input_ch)

    samples_ch = VALIDATE_METADATA.out.metadata
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(metaFromRow(row), readListFromRow(row)) }

    sample_views = samples_ch.multiMap { meta, reads ->
        raw: tuple(meta, reads)
        preprocess: tuple(meta, reads)
    }

    FASTQC_RAW(sample_views.raw)

    preprocess_ch = sample_views.preprocess.branch { meta, reads ->
        fastp: meta.preprocess_mode.startsWith('fastp')
        adapterremoval: meta.preprocess_mode.startsWith('adapterremoval')
    }

    PREPROCESS_FASTP(preprocess_ch.fastp)
    PREPROCESS_ADAPTERREMOVAL(preprocess_ch.adapterremoval)

    prehost_ch = PREPROCESS_FASTP.out.reads.mix(PREPROCESS_ADAPTERREMOVAL.out.reads)

    if (asBool(params.remove_host)) {
        HOST_REMOVE_BWA(prehost_ch)
        clean_ch = HOST_REMOVE_BWA.out.reads
        host_qc_ch = HOST_REMOVE_BWA.out.logs
    } else {
        clean_ch = prehost_ch
        host_qc_ch = Channel.empty()
    }

    FASTQC_CLEAN(clean_ch)

    if (asBool(params.run_kraken2)) {
        KRAKEN2_PROFILE(clean_ch)
        kraken_qc_ch = KRAKEN2_PROFILE.out.reports.map { meta, report, output -> [report, output] }.flatten()
    } else {
        kraken_qc_ch = Channel.empty()
    }

    if (asBool(params.run_bracken)) {
        BRACKEN_PROFILE(KRAKEN2_PROFILE.out.reports.map { meta, report, output -> tuple(meta, report) })
        bracken_files_ch = BRACKEN_PROFILE.out.tables.map { meta, bracken, report -> bracken }.collect()
        BUILD_FEATURE_TABLE(VALIDATE_METADATA.out.metadata, bracken_files_ch)
        bracken_qc_ch = BRACKEN_PROFILE.out.tables.map { meta, bracken, report -> report }

        if (asBool(params.run_sourcetracker2)) {
            PREPARE_SOURCETRACKER2_INPUTS(BUILD_FEATURE_TABLE.out.counts, VALIDATE_METADATA.out.metadata)
            MAKE_BIOM(PREPARE_SOURCETRACKER2_INPUTS.out.counts)
            SOURCETRACKER2_GIBBS(MAKE_BIOM.out.biom, PREPARE_SOURCETRACKER2_INPUTS.out.map)
            PLOT_SOURCETRACKER2(SOURCETRACKER2_GIBBS.out.proportions, PREPARE_SOURCETRACKER2_INPUTS.out.map)
        }
    } else {
        bracken_qc_ch = Channel.empty()
    }

    raw_fastqc_files = FASTQC_RAW.out.qc.map { meta, zip, html -> [zip, html] }.flatten()
    clean_fastqc_files = FASTQC_CLEAN.out.qc.map { meta, zip, html -> [zip, html] }.flatten()
    fastp_qc_files = PREPROCESS_FASTP.out.logs.flatten()
    adapterremoval_qc_files = PREPROCESS_ADAPTERREMOVAL.out.logs.flatten()
    host_qc_files = host_qc_ch.flatten()

    multiqc_files = raw_fastqc_files
        .mix(clean_fastqc_files)
        .mix(fastp_qc_files)
        .mix(adapterremoval_qc_files)
        .mix(host_qc_files)
        .mix(kraken_qc_ch)
        .mix(bracken_qc_ch)
        .collect()

    MULTIQC_REPORT(multiqc_files)
}
