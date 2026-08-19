#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { VALIDATE_METADATA } from './modules/validate_metadata'
include { FASTQC as FASTQC_RAW } from './modules/fastqc'
include { FASTQC as FASTQC_CLEAN } from './modules/fastqc'
include { PREPROCESS_FASTP } from './modules/preprocess_fastp'
include { ADAPTERREMOVAL_RAW } from './modules/preprocess_adapterremoval'
include { NORMALIZE_ADAPTERREMOVAL } from './modules/normalize_adapterremoval'
include { HOST_ALIGN_BWA } from './modules/host_remove_bwa'
include { FILTER_HOST_ALIGNMENTS } from './modules/filter_host_alignments'
include { PREPARE_KRAKEN_INPUT } from './modules/prepare_kraken_input'
include { KRAKEN2_PROFILE } from './modules/kraken2'
include { SELECT_BRACKEN_LENGTH } from './modules/select_bracken_length'
include { BRACKEN_PROFILE } from './modules/bracken'
include { BUILD_FEATURE_TABLE_WITH_KRAKEN_QC } from './modules/build_feature_table'
include { PREPARE_SOURCETRACKER2_INPUTS } from './modules/prepare_sourcetracker2'
include { MAKE_BIOM } from './modules/make_biom'
include { SOURCETRACKER2_GIBBS } from './modules/sourcetracker2'
include { PLOT_SOURCETRACKER2_WITH_UNCERTAINTY } from './modules/plot_sourcetracker2'
include { MULTIQC_REPORT } from './modules/multiqc'
include { SUMMARIZE_SAMPLES } from './modules/summarize_samples'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'


def rejectReservedArgs(name, value, reservedFlags) {
    def text = value == null ? '' : value.toString()
    reservedFlags.each { flag ->
        def isShortOption = flag ==~ /^-[^-]$/
        def optionEnd = isShortOption ? '' : '(?==|[\\s\\x27\\x22;|&]|$)'
        def pattern = java.util.regex.Pattern.compile(
            '(^|[\\s\\x27\\x22;|&])' + java.util.regex.Pattern.quote(flag) +
                optionEnd
        )
        if (pattern.matcher(text).find()) {
            throw new IllegalArgumentException(
                "Parameter '${name}' must not override workflow-controlled option '${flag}'"
            )
        }
    }
}


def resolveReadPath(sampleId, column, value, sheetDir) {
    // Relative paths are read relative to the samplesheet, so a sheet and its
    // reads can move together. Paths relative to the launch directory are still
    // accepted for older sheets; URLs and absolute paths are used as given.
    def text = value.toString().trim()
    if (text ==~ /^[A-Za-z][A-Za-z0-9+.-]*:\/\/.+/ || java.nio.file.Paths.get(text).isAbsolute()) {
        return file(text, checkIfExists: true)
    }
    def besideSheet = sheetDir.resolve(text)
    if (besideSheet.exists()) {
        return besideSheet
    }
    def fromLaunchDir = file(text)
    if (fromLaunchDir.exists()) {
        return fromLaunchDir
    }
    throw new IllegalArgumentException(
        "Sample '${sampleId}': ${column} '${text}' was not found. Relative paths are " +
        "looked up next to the samplesheet (${besideSheet}) and then in the launch " +
        "directory (${fromLaunchDir})."
    )
}


def readListFromRow(row, sheetDir) {
    def read1 = resolveReadPath(row.sample_id, 'fastq_1', row.fastq_1, sheetDir)
    def fastq2 = row.fastq_2.toString()
    if (fastq2 == 'NA' || fastq2 == '') {
        return [read1]
    }
    return [read1, resolveReadPath(row.sample_id, 'fastq_2', fastq2, sheetDir)]
}


def metaFromRow(row) {
    def meta = row.collectEntries { key, value ->
        [(key.toString()): value == null ? '' : value.toString()]
    }
    meta.id = meta.sample_id
    meta.layout = meta.layout.toUpperCase()
    meta.age_group = meta.age_group.toLowerCase()
    return meta
}


def requireDirectory(name, value, requiredFiles = []) {
    if (!value) {
        throw new IllegalArgumentException("Parameter '${name}' is required")
    }
    def directory = file(value, checkIfExists: true)
    if (!java.nio.file.Files.isDirectory(directory)) {
        throw new IllegalArgumentException("Parameter '${name}' must be a directory: ${directory}")
    }
    requiredFiles.each { filename ->
        def required = directory.resolve(filename)
        if (!java.nio.file.Files.isRegularFile(required)) {
            throw new IllegalArgumentException(
                "${name} is missing required file: ${required}. It should point at a Kraken2 " +
                "database directory; prebuilt ones are listed at https://benlangmead.github.io/aws-indexes/k2"
            )
        }
    }
    return directory
}


def requireBrackenDistributions(directory) {
    def matches = java.nio.file.Files.newDirectoryStream(directory, 'database*mers.kmer_distrib')
    def foundDistribution = matches.iterator().hasNext()
    matches.close()
    if (!foundDistribution) {
        throw new IllegalArgumentException(
            "bracken_db contains no database<N>mers.kmer_distrib files: ${directory}. " +
            "Build one per read length with `bracken-build -d ${directory} -l 100`, " +
            "or use a prebuilt index, which ships them."
        )
    }
}


def resolveTaxdumpNodes(value) {
    if (!value) {
        throw new IllegalArgumentException("Parameter 'taxdump' is required for taxon_filter=prokaryote")
    }
    def candidate = file(value, checkIfExists: true)
    def nodes = java.nio.file.Files.isDirectory(candidate) ? candidate.resolve('nodes.dmp') : candidate
    if (!java.nio.file.Files.isRegularFile(nodes)) {
        throw new IllegalArgumentException("taxdump nodes.dmp does not exist: ${nodes}")
    }
    return nodes
}


def hostReferenceFor(meta, defaultReference) {
    def requested = (meta.host_ref ?: '').toString().trim()
    if (requested in ['', 'NA', 'N/A', 'default']) {
        if (!defaultReference) {
            throw new IllegalArgumentException(
                "Sample '${meta.id}' uses host_ref=default but parameter 'host_fasta' is not set"
            )
        }
        return defaultReference.toString()
    }
    return requested
}


def hostBundleFor(reference) {
    def fasta = file(reference, checkIfExists: true)
    if (!java.nio.file.Files.isRegularFile(fasta)) {
        throw new IllegalArgumentException("Host FASTA is not a regular file: ${fasta}")
    }
    def bundle = [fasta]
    ['.amb', '.ann', '.bwt', '.pac', '.sa'].each { suffix ->
        def index = file("${fasta}${suffix}", checkIfExists: true)
        if (!java.nio.file.Files.isRegularFile(index)) {
            throw new IllegalArgumentException("Host FASTA is missing BWA index: ${index}")
        }
        bundle << index
    }
    def fai = file("${fasta}.fai")
    if (java.nio.file.Files.isRegularFile(fai)) {
        bundle << fai
    }
    return bundle
}


workflow {
    // Types, ranges, choices and file existence come from nextflow_schema.json.
    validateParameters()
    log.info paramsSummaryLog(workflow)

    // Constraints between parameters, which the schema cannot express.
    ['min_length_ancient', 'min_length_modern', 'kraken2_min_read_len', 'fragment_pair_gap'].each { name ->
        if ((params[name] as int) < (params.kraken2_kmer_len as int)) {
            error "Parameter '${name}' (${params[name]}) must be at least kraken2_kmer_len (${params.kraken2_kmer_len})"
        }
    }

    def removeHost = params.remove_host
    def runKraken2 = params.run_kraken2
    def runBracken = params.run_bracken
    def runSourceTracker2 = params.run_sourcetracker2
    if (!(params.sourcetracker2_analysis_label ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) {
        error "Parameter 'sourcetracker2_analysis_label' must be filesystem-safe"
    }

    // '--check-names' is deliberately rejected: it is not a Kraken2 2.1.6 option.
    // Its Perl wrapper parses long options with Getopt::Long and ignores the
    // parse result, so an unknown long option only warns and classification
    // continues with exit 0. KRAKEN2_PROFILE fails closed on that warning, and
    // mate identity is verified by the fragment helper before Kraken2 runs.
    rejectReservedArgs(
        'kraken2_extra_args',
        params.kraken2_extra_args,
        [
            '--db', '--threads', '--paired', '--report', '--output', '--confidence',
            '--minimum-hit-groups', '--use-names', '--report-minimizer-data',
            '--use-mpa-style', '--bzip2-compressed', '--check-names'
        ]
    )
    rejectReservedArgs(
        'fastp_extra_args',
        params.fastp_extra_args,
        [
            '-i', '--in1', '-I', '--in2', '-o', '--out1', '-O', '--out2',
            '--detect_adapter_for_pe', '--cut_front', '--cut_tail', '--cut_window_size',
            '--cut_mean_quality', '-l', '--length_required', '-w', '--thread',
            '-h', '--html', '-j', '--json'
        ]
    )
    rejectReservedArgs(
        'adapterremoval_extra_args',
        params.adapterremoval_extra_args,
        [
            '--file1', '--file2', '--basename', '--threads', '--trimns',
            '--trimqualities', '--minquality', '--minlength', '--gzip', '--collapse',
            '--collapse-deterministic', '--collapse-conservatively', '--seed',
            '--output1', '--output2', '--outputcollapsed', '--outputcollapsedtruncated',
            '--outputsingleton', '--outputsingletontruncated', '--outputdiscarded',
            '--singleton', '--discarded', '--settings', '--combined-output', '--interleaved',
            '--interleaved-input', '--interleaved-output', '--demultiplex-only', '--bzip2'
        ]
    )
    rejectReservedArgs(
        'bracken_extra_args',
        params.bracken_extra_args,
        ['-d', '-i', '-o', '-w', '-r', '-l', '-t']
    )
    rejectReservedArgs(
        'sourcetracker2_extra_args',
        params.sourcetracker2_extra_args,
        [
            '-i', '-m', '-o', '--loo', '--jobs', '--source_rarefaction_depth',
            '--sink_rarefaction_depth', '--source_sink_column', '--source_category_column',
            '--source_column_value', '--sink_column_value', '--seed'
        ]
    )

    if (runBracken && !runKraken2) {
        error 'run_bracken requires run_kraken2'
    }
    if (runSourceTracker2 && !runBracken) {
        error 'run_sourcetracker2 requires run_bracken'
    }
    if (
        runSourceTracker2 &&
        params.sourcetracker2_loo &&
        params.sourcetracker2_source_rarefaction_depth != 0
    ) {
        error 'SourceTracker2 2.0.1 LOO requires sourcetracker2_source_rarefaction_depth=0; its CLI subtracts unrarefied held-out counts from already-rarefied source counts'
    }
    if (!params.input) {
        error 'The input parameter is required and must identify exactly one metadata TSV file'
    }

    def krakenDbPath = null
    def brackenDbPath = null
    if (runKraken2) {
        krakenDbPath = requireDirectory(
            'kraken2_db', params.kraken2_db, ['hash.k2d', 'opts.k2d', 'taxo.k2d']
        )
    }
    if (runBracken) {
        brackenDbPath = requireDirectory(
            'bracken_db', params.bracken_db ?: params.kraken2_db, ['hash.k2d', 'opts.k2d', 'taxo.k2d']
        )
        requireBrackenDistributions(brackenDbPath)
        if (krakenDbPath.toRealPath() != brackenDbPath.toRealPath()) {
            error 'kraken2_db and bracken_db must resolve to the same database snapshot'
        }
    }

    def noTaxdumpSentinel = file("${projectDir}/assets/no_taxdump.nodes.dmp", checkIfExists: true)
    def taxdumpNodesPath = (runBracken && params.taxon_filter == 'prokaryote') \
        ? resolveTaxdumpNodes(params.taxdump) \
        : noTaxdumpSentinel

    def validatorScript = file("${projectDir}/bin/validate_metadata.py", checkIfExists: true)
    def fragmentHelper = file("${projectDir}/bin/mstflow_fragment_tools.py", checkIfExists: true)
    def adapterRemovalNormalizer = file(
        "${projectDir}/bin/normalize_adapterremoval.py", checkIfExists: true
    )
    def hostAlignmentFilter = file(
        "${projectDir}/bin/filter_host_alignments.py", checkIfExists: true
    )
    def krakenInputPreparer = file(
        "${projectDir}/bin/prepare_kraken_input.py", checkIfExists: true
    )
    def buildFeatureTableScript = file("${projectDir}/bin/build_feature_table.py", checkIfExists: true)
    def prepareSourceTrackerScript = file(
        "${projectDir}/bin/prepare_sourcetracker2_inputs.py", checkIfExists: true
    )
    def sourceTrackerRunner = file(
        "${projectDir}/bin/run_sourcetracker2_seeded.py", checkIfExists: true
    )
    def plotSourceTrackerScript = file("${projectDir}/bin/plot_sourcetracker2.py", checkIfExists: true)
    def summarizeScript = file("${projectDir}/bin/summarize_samples.py", checkIfExists: true)

    validator_script_ch = channel.value(validatorScript)
    fragment_helper_ch = channel.value(fragmentHelper)
    adapterremoval_normalizer_ch = channel.value(adapterRemovalNormalizer)
    host_alignment_filter_ch = channel.value(hostAlignmentFilter)
    kraken_input_preparer_ch = channel.value(krakenInputPreparer)
    build_feature_table_script_ch = channel.value(buildFeatureTableScript)
    prepare_sourcetracker2_script_ch = channel.value(prepareSourceTrackerScript)
    sourcetracker2_runner_ch = channel.value(sourceTrackerRunner)
    plot_sourcetracker2_script_ch = channel.value(plotSourceTrackerScript)
    taxdump_nodes_ch = channel.value(taxdumpNodesPath)
    summarize_script_ch = channel.value(summarizeScript)

    def inputMatches = file(params.input, checkIfExists: true)
    inputMatches = inputMatches instanceof List ? inputMatches : [inputMatches]
    if (inputMatches.size() != 1) {
        error "The input parameter must identify exactly one samplesheet; matched ${inputMatches.size()} files"
    }
    def samplesheetDir = inputMatches.first().parent

    input_ch = channel.fromPath(params.input, checkIfExists: true)
        .collect()
        .map { matches ->
            if (matches.size() != 1) {
                error "The input parameter must identify exactly one metadata TSV file; matched ${matches.size()} files"
            }
            matches.first()
        }

    VALIDATE_METADATA(input_ch, validator_script_ch)

    // The validator runs inside a portable Python container and therefore
    // treats filesystem paths embedded in TSV text as opaque. The Nextflow
    // controller performs the authoritative existence checks here, and
    // hostBundleFor() independently checks each selected FASTA/index bundle.
    samples_ch = VALIDATE_METADATA.out.metadata
        .splitCsv(header: true, sep: '\t')
        .map { row -> tuple(metaFromRow(row), readListFromRow(row, samplesheetDir)) }

    host_bundle_ch = removeHost \
        ? samples_ch
            .filter { meta, _reads -> meta.host_action == 'remove' }
            .map { meta, _reads ->
                tuple(meta.id, hostBundleFor(hostReferenceFor(meta, params.host_fasta)))
            }
        : channel.empty()

    sample_views = samples_ch.multiMap { meta, reads ->
        raw: tuple(meta, reads)
        preprocess: tuple(meta, reads)
    }

    raw_qc_ch = sample_views.raw.map { meta, reads ->
        tuple(meta + [qc_stage: 'raw'], reads)
    }
    FASTQC_RAW(raw_qc_ch)

    preprocess_ch = sample_views.preprocess.branch { meta, _reads ->
        fastp: meta.preprocess_mode.startsWith('fastp')
        adapterremoval: meta.preprocess_mode.startsWith('adapterremoval')
    }

    PREPROCESS_FASTP(preprocess_ch.fastp)
    ADAPTERREMOVAL_RAW(preprocess_ch.adapterremoval)
    NORMALIZE_ADAPTERREMOVAL(
        ADAPTERREMOVAL_RAW.out.raw,
        adapterremoval_normalizer_ch,
        fragment_helper_ch
    )

    prehost_ch = PREPROCESS_FASTP.out.reads.mix(NORMALIZE_ADAPTERREMOVAL.out.reads)

    if (removeHost) {
        // Fail closed rather than silently dropping samples. `host_action` is
        // the one branch key whose input column is optional, so a validated
        // manifest written before host policies existed would otherwise lose
        // every row here without any error.
        host_views = prehost_ch.branch { meta, _reads ->
            remove: meta.host_action == 'remove'
            skip: meta.host_action == 'skip'
            undeclared: true
        }
        undeclared_host_action_ch = host_views.undeclared.map { meta, _reads ->
            error "Sample '${meta.id}' has host_action '${meta.host_action}'; validated metadata must declare remove or skip"
        }
        host_input_ch = host_views.remove
            .map { meta, reads -> tuple(meta.id, meta, reads) }
            .join(host_bundle_ch)
            .map { _id, meta, reads, hostBundle -> tuple(meta, reads, hostBundle) }
        HOST_ALIGN_BWA(host_input_ch)
        FILTER_HOST_ALIGNMENTS(
            HOST_ALIGN_BWA.out.alignments,
            host_alignment_filter_ch,
            fragment_helper_ch
        )
        clean_ch = FILTER_HOST_ALIGNMENTS.out.reads
            .mix(host_views.skip)
            .mix(undeclared_host_action_ch)
        host_qc_ch = HOST_ALIGN_BWA.out.logs.mix(FILTER_HOST_ALIGNMENTS.out.qc)
    } else {
        clean_ch = prehost_ch
        host_qc_ch = channel.empty()
    }

    clean_qc_ch = clean_ch.map { meta, reads ->
        tuple(meta + [qc_stage: 'clean'], reads)
    }
    FASTQC_CLEAN(clean_qc_ch)

    if (runKraken2) {
        PREPARE_KRAKEN_INPUT(clean_ch, kraken_input_preparer_ch, fragment_helper_ch)
        kraken_input_ch = PREPARE_KRAKEN_INPUT.out.prepared.map {
            meta, reads, readLength, inputMode ->
                tuple(meta, reads, readLength, inputMode, krakenDbPath)
        }
        KRAKEN2_PROFILE(kraken_input_ch)
        kraken_qc_ch = KRAKEN2_PROFILE.out.reports
            .map { _meta, report, readLength -> [report, readLength] }
            .flatten()
            .mix(KRAKEN2_PROFILE.out.qc)
            .mix(KRAKEN2_PROFILE.out.minimizer_reports)
    } else {
        kraken_qc_ch = channel.empty()
    }

    if (runBracken) {
        bracken_selection_input_ch = KRAKEN2_PROFILE.out.reports.map {
            meta, _report, readLength ->
                tuple(meta, readLength, brackenDbPath)
        }
        SELECT_BRACKEN_LENGTH(bracken_selection_input_ch, fragment_helper_ch)
        bracken_report_keyed_ch = KRAKEN2_PROFILE.out.reports.map {
            meta, report, _readLength -> tuple(meta.id, meta, report)
        }
        bracken_selection_keyed_ch = SELECT_BRACKEN_LENGTH.out.selections.map {
            meta, selection -> tuple(meta.id, selection)
        }
        bracken_input_ch = bracken_report_keyed_ch
            .join(
                bracken_selection_keyed_ch,
                failOnDuplicate: true,
                failOnMismatch: true
            )
            .map { _id, meta, report, selection ->
                tuple(meta, report, selection, brackenDbPath)
            }
        BRACKEN_PROFILE(bracken_input_ch)

        // `collect()` preserves task completion order, so the same inputs can
        // reach a terminal process in a different order on every run. Sorting
        // by file name makes the staged list, the rendered command line, and
        // therefore the task cache key depend only on the inputs. The tables
        // themselves are order-invariant, so results do not change.
        bracken_files_ch = BRACKEN_PROFILE.out.tables
            .map { _meta, bracken, _report -> bracken }
            .toSortedList { left, right -> left.name <=> right.name }
        kraken_reports_ch = KRAKEN2_PROFILE.out.reports
            .map { _meta, report, _readLength -> report }
            .toSortedList { left, right -> left.name <=> right.name }

        BUILD_FEATURE_TABLE_WITH_KRAKEN_QC(
            VALIDATE_METADATA.out.metadata,
            bracken_files_ch,
            kraken_reports_ch,
            taxdump_nodes_ch,
            build_feature_table_script_ch
        )

        bracken_qc_ch = BRACKEN_PROFILE.out.tables
            .map { _meta, _bracken, report -> report }
            .mix(SELECT_BRACKEN_LENGTH.out.selections.map { _meta, selection -> selection })
            .mix(BRACKEN_PROFILE.out.status)
        table_qc_ch = BUILD_FEATURE_TABLE_WITH_KRAKEN_QC.out.qc
            .mix(BUILD_FEATURE_TABLE_WITH_KRAKEN_QC.out.sample_qc)

        if (runSourceTracker2) {
            st_counts_ch = params.sourcetracker2_taxon_scope == 'all' \
                ? BUILD_FEATURE_TABLE_WITH_KRAKEN_QC.out.raw_counts \
                : BUILD_FEATURE_TABLE_WITH_KRAKEN_QC.out.counts
            PREPARE_SOURCETRACKER2_INPUTS(
                st_counts_ch,
                VALIDATE_METADATA.out.metadata,
                BUILD_FEATURE_TABLE_WITH_KRAKEN_QC.out.sample_qc,
                prepare_sourcetracker2_script_ch
            )
            MAKE_BIOM(PREPARE_SOURCETRACKER2_INPUTS.out.counts)
            SOURCETRACKER2_GIBBS(
                MAKE_BIOM.out.biom,
                PREPARE_SOURCETRACKER2_INPUTS.out.map,
                sourcetracker2_runner_ch
            )
            PLOT_SOURCETRACKER2_WITH_UNCERTAINTY(
                SOURCETRACKER2_GIBBS.out.proportions,
                SOURCETRACKER2_GIBBS.out.uncertainty,
                PREPARE_SOURCETRACKER2_INPUTS.out.map,
                plot_sourcetracker2_script_ch
            )
            st_qc_ch = PREPARE_SOURCETRACKER2_INPUTS.out.summary
                .mix(PREPARE_SOURCETRACKER2_INPUTS.out.exclusions)
                .mix(PREPARE_SOURCETRACKER2_INPUTS.out.sample_qc)
        } else {
            st_qc_ch = channel.empty()
        }
    } else {
        bracken_qc_ch = channel.empty()
        table_qc_ch = channel.empty()
        st_qc_ch = channel.empty()
    }

    raw_fastqc_files = FASTQC_RAW.out.qc.map { _meta, zip, html -> [zip, html] }.flatten()
    clean_fastqc_files = FASTQC_CLEAN.out.qc.map { _meta, zip, html -> [zip, html] }.flatten()
    fastp_qc_files = PREPROCESS_FASTP.out.logs.flatten()
    adapterremoval_qc_files = ADAPTERREMOVAL_RAW.out.logs
        .flatten()
        .mix(NORMALIZE_ADAPTERREMOVAL.out.qc)

    // One row per sample: molecules surviving each step, and depth flags.
    if (runBracken) {
        summary_inputs_ch = fastp_qc_files
            .mix(adapterremoval_qc_files)
            .mix(host_qc_ch)
            .filter { qcFile ->
                qcFile.name ==~ /.+\.(fastp\.json|adapterremoval\.log|adapterremoval\.fragment_qc\.tsv|host_removal\.fragment_qc\.tsv)/
            }
            .toSortedList { left, right -> left.name <=> right.name }
        SUMMARIZE_SAMPLES(
            VALIDATE_METADATA.out.metadata,
            BUILD_FEATURE_TABLE_WITH_KRAKEN_QC.out.sample_qc,
            summary_inputs_ch,
            summarize_script_ch
        )
        summary_mqc_ch = SUMMARIZE_SAMPLES.out.multiqc
    } else {
        summary_mqc_ch = channel.empty()
    }

    // Sorted so MultiQC always stages the same file into the same qcNN
    // directory, which keeps the report input set and the task cache key
    // independent of task completion order.
    multiqc_files = raw_fastqc_files
        .mix(clean_fastqc_files)
        .mix(fastp_qc_files)
        .mix(adapterremoval_qc_files)
        .mix(host_qc_ch)
        .mix(kraken_qc_ch)
        .mix(bracken_qc_ch)
        .mix(table_qc_ch)
        .mix(st_qc_ch)
        .mix(summary_mqc_ch)
        .toSortedList { left, right -> left.name <=> right.name }

    MULTIQC_REPORT(multiqc_files)

    // Printed once the run ends: where the results are and which samples need care.
    // The handler runs outside the workflow's scope, so it gets plain values.
    def run = workflow
    def report = log
    def outdir = params.outdir
    def summaryFile = file("${params.outdir}/summary/sample_summary.tsv")
    run.onComplete {
        def lines = ['']
        if (run.success) {
            lines << "MSTFlow finished in ${run.duration}.".toString()
            if (summaryFile.exists()) {
                def rows = summaryFile.splitCsv(header: true, sep: '\t')
                def flagged = rows.findAll { row -> row.flags != 'ok' }
                lines << "  ${rows.size()} samples, ${flagged.size()} flagged as too shallow to interpret".toString()
                flagged.each { row -> lines << "    ${row.sample_id}: ${row.flags}".toString() }
                lines << "  Summary  ${summaryFile}".toString()
            }
            lines << "  Tables   ${outdir}/tables/".toString()
            lines << "  Report   ${outdir}/qc/multiqc/multiqc_report.html".toString()
        } else {
            lines << 'MSTFlow stopped with an error. The message above says why; .nextflow.log has the details.'
        }
        report.info lines.join('\n')
    }
}
