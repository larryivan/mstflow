include { UNTAR_DATABASE } from '../modules/untar_database'
include { PREPARE_KRAKEN_INPUT } from '../modules/prepare_kraken_input'
include { KRAKEN2_PROFILE } from '../modules/kraken2'
include { KRAKENUNIQ_PROFILE } from '../modules/krakenuniq'
include { SELECT_BRACKEN_LENGTH } from '../modules/select_bracken_length'
include { BRACKEN_PROFILE } from '../modules/bracken'
include { BUILD_FEATURE_TABLE } from '../modules/build_feature_table'
include { binScript; isArchive; isOn } from './utils'

/*
 * Classify every molecule with Kraken2 or KrakenUniq (--classifier), re-estimate
 * abundances at one rank with Bracken, and merge the samples into one count table.
 * Without Bracken (--run_bracken false), the classifier's own clade counts at that
 * rank are used. Taxa with too few distinct minimizers (Kraken2) or unique k-mers
 * (KrakenUniq) can be removed sample by sample when the table is built.
 */
workflow PROFILE {
    take:
    reads           // tuple(meta, [reads]) after host removal
    metadata        // validated samplesheet
    database        // classifier database directory, or a .tar.gz archive
    taxdump_nodes   // nodes.dmp for the prokaryote filter, or []

    main:
    def krakenUniq = params.classifier == 'krakenuniq'
    def fragmentHelper = channel.value(binScript('mstflow_fragment_tools.py'))

    if (isArchive(krakenUniq ? params.krakenuniq_db : params.kraken2_db)) {
        UNTAR_DATABASE(database, "${params.classifier}${isOn(params.run_bracken) ? '+bracken' : ''}")
        db = UNTAR_DATABASE.out.db.first()
    } else {
        db = database
    }

    PREPARE_KRAKEN_INPUT(reads, binScript('prepare_kraken_input.py'), fragmentHelper)
    classifier_input = PREPARE_KRAKEN_INPUT.out.prepared.combine(db)

    // Both classifiers give tuple(meta, report, read-length statistics).
    if (krakenUniq) {
        KRAKENUNIQ_PROFILE(classifier_input)
        reports = KRAKENUNIQ_PROFILE.out.reports
        minimizer_reports = channel.value([])
    } else {
        KRAKEN2_PROFILE(classifier_input)
        reports = KRAKEN2_PROFILE.out.reports
        minimizer_reports = KRAKEN2_PROFILE.out.minimizer_reports.toSortedList { a, b -> a.name <=> b.name }
    }

    if (isOn(params.run_bracken)) {
        SELECT_BRACKEN_LENGTH(reports.map { meta, _report, readLength -> tuple(meta, readLength) }.combine(db), fragmentHelper)
        bracken_input = reports
            .map { meta, report, _readLength -> tuple(meta.id, meta, report) }
            .join(
                SELECT_BRACKEN_LENGTH.out.selections.map { meta, selection -> tuple(meta.id, selection) },
                failOnDuplicate: true,
                failOnMismatch: true
            )
            .map { _id, meta, report, selection -> tuple(meta, report, selection) }
            .combine(db)
        BRACKEN_PROFILE(bracken_input)
        // Per-sample files are sorted by name before the merge, so the staged list and
        // so the cache key depend only on the inputs, not on which sample finished first.
        bracken_tables = BRACKEN_PROFILE.out.tables
            .map { _meta, table, _report -> table }
            .toSortedList { a, b -> a.name <=> b.name }
        bracken_reports = BRACKEN_PROFILE.out.tables.map { _meta, _table, report -> report }
    } else {
        bracken_tables = channel.value([])
        bracken_reports = channel.empty()
    }

    BUILD_FEATURE_TABLE(
        metadata,
        bracken_tables,
        reports.map { _meta, report, _readLength -> report }.toSortedList { a, b -> a.name <=> b.name },
        minimizer_reports,
        taxdump_nodes,
        binScript('build_feature_table.py')
    )

    emit:
    all_counts = BUILD_FEATURE_TABLE.out.raw_counts       // one row per taxon, before the taxon filter
    filtered_counts = BUILD_FEATURE_TABLE.out.counts      // after the taxon filter
    sample_qc = BUILD_FEATURE_TABLE.out.sample_qc         // per-sample classified and retained reads
    qc = reports.map { _meta, report, _readLength -> report }.mix(bracken_reports)   // Kraken-style reports, for MultiQC
}
