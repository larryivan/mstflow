#!/usr/bin/env nextflow

/*
 * MSTFlow: metagenomic profiling of ancient and modern DNA that counts molecules,
 * not reads, followed by optional source attribution with several methods.
 *
 *   samplesheet -> PREPROCESS -> PROFILE -> ATTRIBUTION
 *                  trimming,     Kraken2 and Bracken,    SourceTracker2, FEAST, FastST,
 *                  host removal  or KrakenUniq           SourceID-NMF, SourcePredict
 *                                          \-> AUTHENTICATION
 *                                              Bowtie2 and PyDamage
 */

include { VALIDATE_METADATA } from './modules/validate_metadata'
include { SUMMARIZE_SAMPLES } from './modules/summarize_samples'
include { MULTIQC_REPORT } from './modules/multiqc'
include { PREPROCESS } from './subworkflows/preprocess'
include { PROFILE } from './subworkflows/profile'
include { ATTRIBUTION } from './subworkflows/attribution'
include { AUTHENTICATION } from './subworkflows/authentication'
include {
    binScript; isArchive; isOn; metaFromRow; readListFromRow; rejectReservedArgs;
    requireBrackenDistributions; requireDirectory; splitList; taxdumpFile
} from './subworkflows/utils'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


workflow {
    // Types, ranges, choices and file existence come from nextflow_schema.json.
    validateParameters()
    log.info paramsSummaryLog(workflow)
    checkParameters()

    def methods = splitList(params.attribution_methods)

    def samplesheet = file(params.input, checkIfExists: true)
    if (samplesheet instanceof List) {
        error "The input parameter must identify exactly one samplesheet; matched ${samplesheet.size()} files"
    }
    VALIDATE_METADATA(samplesheet, binScript('validate_metadata.py'))
    metadata = VALIDATE_METADATA.out.metadata

    // The validator runs in a container and treats paths as text; existence is
    // checked here, against the samplesheet's own directory first.
    samples = metadata
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(metaFromRow(row), readListFromRow(row, samplesheet.parent)) }

    PREPROCESS(samples)
    qc_files = PREPROCESS.out.qc

    if (isOn(params.run_kraken2)) {
        def taxdumpNodes = params.taxon_filter == 'prokaryote'
            ? taxdumpFile(params.taxdump, 'nodes.dmp', 'taxon_filter=prokaryote')
            : []
        PROFILE(PREPROCESS.out.reads, metadata, channel.value(resolveDatabase()), channel.value(taxdumpNodes))
        qc_files = qc_files.mix(PROFILE.out.qc)

        // One row per sample: molecules kept by each step, and depth flags.
        SUMMARIZE_SAMPLES(
            metadata,
            PROFILE.out.sample_qc,
            PREPROCESS.out.step_counts.toSortedList { a, b -> a.name <=> b.name },
            binScript('summarize_samples.py')
        )
        qc_files = qc_files.mix(SUMMARIZE_SAMPLES.out.multiqc)

        if (params.authentication_genomes) {
            AUTHENTICATION(PREPROCESS.out.reads, PROFILE.out.filtered_counts, metadata)
        }

        if (methods) {
            ATTRIBUTION(
                params.attribution_taxon_scope == 'all' ? PROFILE.out.all_counts : PROFILE.out.filtered_counts,
                metadata,
                PROFILE.out.sample_qc
            )
        }
    }

    // Sorted so MultiQC always stages the same file into the same qcNN directory,
    // which keeps its inputs, and so its cache key, independent of task order.
    MULTIQC_REPORT(qc_files.toSortedList { a, b -> a.name <=> b.name })

    // Printed once the run ends: where the results are and which samples need care.
    // The handler runs outside the workflow's scope, so it gets plain values.
    def run = workflow
    def report = log
    def summaryFile = file("${params.outdir}/summary/sample_summary.tsv")
    def outputs = [
        Summary: 'summary/sample_summary.tsv',
        Tables: 'tables/',
        Attribution: "attribution/${params.attribution_label}/",
        Authentication: 'authentication/authentication.tsv',
        Report: 'qc/multiqc/multiqc_report.html',
    ].collect { name, path -> ["  ${name.padRight(15)} ${params.outdir}/${path}".toString(), file("${params.outdir}/${path}")] }
    run.onComplete {
        def lines = ['']
        if (run.success) {
            lines << "MSTFlow finished in ${run.duration}.".toString()
            if (summaryFile.exists()) {
                def rows = summaryFile.splitCsv(header: true, sep: '\t')
                def flagged = rows.findAll { row -> row.flags != 'ok' }
                lines << "  ${rows.size()} samples, ${flagged.size()} flagged as too shallow to interpret".toString()
                flagged.each { row -> lines << "    ${row.sample_id}: ${row.flags}".toString() }
            }
            outputs.findAll { _line, output -> output.exists() }.each { line, _output -> lines << line }
        } else {
            lines << 'MSTFlow stopped with an error. The message above says why; .nextflow.log has the details.'
        }
        report.info lines.join('\n')
    }
}


// Constraints between parameters, which the schema cannot express.
def checkParameters() {
    ['min_length_ancient', 'min_length_modern', 'kraken2_min_read_len'].each { name ->
        if ((params[name] as int) < (params.kraken2_kmer_len as int)) {
            error "Parameter '${name}' (${params[name]}) must be at least kraken2_kmer_len (${params.kraken2_kmer_len})"
        }
    }

    // '--check-names' is deliberately rejected: it is not a Kraken2 2.1.6 option.
    // Its Perl wrapper parses long options with Getopt::Long and ignores the parse
    // result, so an unknown long option only warns and classification continues
    // with exit 0. KRAKEN2_PROFILE fails closed on that warning, and mate identity
    // is verified by the fragment helper before Kraken2 runs.
    rejectReservedArgs('kraken2_extra_args', params.kraken2_extra_args, [
        '--db', '--threads', '--paired', '--report', '--output', '--confidence',
        '--minimum-hit-groups', '--use-names', '--report-minimizer-data',
        '--use-mpa-style', '--bzip2-compressed', '--check-names'
    ])
    rejectReservedArgs('fastp_extra_args', params.fastp_extra_args, [
        '-i', '--in1', '-I', '--in2', '-o', '--out1', '-O', '--out2',
        '--detect_adapter_for_pe', '--cut_front', '--cut_tail', '--cut_window_size',
        '--cut_mean_quality', '-l', '--length_required', '-w', '--thread',
        '-h', '--html', '-j', '--json'
    ])
    rejectReservedArgs('adapterremoval_extra_args', params.adapterremoval_extra_args, [
        '--file1', '--file2', '--basename', '--threads', '--trimns',
        '--trimqualities', '--minquality', '--minlength', '--gzip', '--collapse',
        '--collapse-deterministic', '--collapse-conservatively', '--seed',
        '--output1', '--output2', '--outputcollapsed', '--outputcollapsedtruncated',
        '--outputsingleton', '--outputsingletontruncated', '--outputdiscarded',
        '--singleton', '--discarded', '--settings', '--combined-output', '--interleaved',
        '--interleaved-input', '--interleaved-output', '--demultiplex-only', '--bzip2'
    ])
    rejectReservedArgs('bracken_extra_args', params.bracken_extra_args, ['-d', '-i', '-o', '-w', '-r', '-l', '-t'])

    if (isOn(params.run_bracken) && !isOn(params.run_kraken2)) {
        error 'run_bracken requires classification (run_kraken2)'
    }
    if ((params.kraken2_min_distinct_minimizers as int) > 0 && !isOn(params.kraken2_report_minimizer_data)) {
        error 'kraken2_min_distinct_minimizers needs kraken2_report_minimizer_data = true'
    }
    def methods = splitList(params.attribution_methods)
    if ((methods || params.authentication_genomes) && !isOn(params.run_kraken2)) {
        error 'Source attribution and authentication need the count table, so classification must run (run_kraken2)'
    }
    // FEAST and FastST come from GitHub and are only packaged in the attribution image.
    def imageOnly = methods.intersect(['feast', 'fastst'])
    if (imageOnly && !workflow.containerEngine && !workflow.stubRun) {
        error "${imageOnly.join(' and ')} run only from MSTFlow's attribution container; use -profile docker, singularity or apptainer"
    }
    if (isOn(params.attribution_cv) && !methods) {
        error 'attribution_cv needs at least one method in attribution_methods'
    }
}


// A database is a directory, or a .tar.gz archive (local or a URL) that the workflow
// unpacks once. Bracken reads the classifier's database.
def resolveDatabase() {
    if (params.classifier == 'krakenuniq') {
        if (isArchive(params.krakenuniq_db)) {
            return file(params.krakenuniq_db, checkIfExists: true)
        }
        def directory = requireDirectory('krakenuniq_db', params.krakenuniq_db, ['database.kdb', 'database.idx', 'taxDB'],
            'It should point at a KrakenUniq database, such as one built with krakenuniq-build.')
        if (isOn(params.run_bracken)) {
            requireBrackenDistributions(directory, 'krakenuniq')
        }
        return directory
    }
    if (isArchive(params.kraken2_db)) {
        return file(params.kraken2_db, checkIfExists: true)
    }
    def directory = requireDirectory('kraken2_db', params.kraken2_db, ['hash.k2d', 'opts.k2d', 'taxo.k2d'],
        'It should point at a Kraken2 database; prebuilt ones are listed at https://benlangmead.github.io/aws-indexes/k2')
    if (isOn(params.run_bracken)) {
        requireBrackenDistributions(directory, 'kraken2')
    }
    return directory
}
