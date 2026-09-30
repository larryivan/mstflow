include { UNTAR_DATABASE } from '../modules/untar_database'
include { PREPARE_KRAKEN_INPUT } from '../modules/prepare_kraken_input'
include { KRAKEN2_PROFILE } from '../modules/kraken2'
include { SELECT_BRACKEN_LENGTH } from '../modules/select_bracken_length'
include { BRACKEN_PROFILE } from '../modules/bracken'
include { KRAKENUNIQ_PROFILE } from '../modules/krakenuniq'
include { BUILD_FEATURE_TABLE } from '../modules/build_feature_table'
include { binScript; isArchive } from './utils'

/*
 * Classify every molecule and merge the samples into one count table at one rank.
 *
 *   kraken2     Kraken2, then Bracken re-estimates abundances at the rank (default)
 *   krakenuniq  KrakenUniq; its clade counts at the rank are used as they are, and
 *               taxa covering too few unique k-mers are removed
 */
workflow PROFILE {
    take:
    reads           // tuple(meta, [reads]) after host removal
    metadata        // validated samplesheet
    database        // classifier database directory, or a .tar.gz archive
    taxdump_nodes   // nodes.dmp for the prokaryote filter, or a placeholder

    main:
    def krakenUniq = params.classifier == 'krakenuniq'
    def fragmentHelper = channel.value(binScript('mstflow_fragment_tools.py'))

    if (isArchive(krakenUniq ? params.krakenuniq_db : params.kraken2_db)) {
        UNTAR_DATABASE(database, krakenUniq ? 'krakenuniq' : (params.run_bracken ? 'kraken2+bracken' : 'kraken2'))
        db = UNTAR_DATABASE.out.db.first()
    } else {
        db = database
    }

    PREPARE_KRAKEN_INPUT(reads, binScript('prepare_kraken_input.py'), fragmentHelper)
    classifier_input = PREPARE_KRAKEN_INPUT.out.prepared.combine(db)

    // Per-sample files are sorted by name before the merge, so the staged list and
    // so the cache key depend only on the inputs, not on which sample finished first.
    if (krakenUniq) {
        KRAKENUNIQ_PROFILE(classifier_input)
        reports = KRAKENUNIQ_PROFILE.out.reports
            .map { _meta, report -> report }
            .toSortedList { a, b -> a.name <=> b.name }
        bracken_tables = channel.value([])
        minimizer_reports = channel.value([])
        classifier_qc = KRAKENUNIQ_PROFILE.out.reports
            .map { _meta, report -> report }
            .mix(KRAKENUNIQ_PROFILE.out.qc)
    } else {
        KRAKEN2_PROFILE(classifier_input)
        reports = KRAKEN2_PROFILE.out.reports
            .map { _meta, report, _readLength -> report }
            .toSortedList { a, b -> a.name <=> b.name }
        minimizer_reports = KRAKEN2_PROFILE.out.minimizer_reports.toSortedList { a, b -> a.name <=> b.name }
        classifier_qc = KRAKEN2_PROFILE.out.reports
            .map { _meta, report, readLength -> [report, readLength] }
            .flatten()
            .mix(KRAKEN2_PROFILE.out.qc)
            .mix(KRAKEN2_PROFILE.out.minimizer_reports)

        if (params.run_bracken) {
            SELECT_BRACKEN_LENGTH(
                KRAKEN2_PROFILE.out.reports.map { meta, _report, readLength -> tuple(meta, readLength) }.combine(db),
                fragmentHelper
            )
            bracken_input = KRAKEN2_PROFILE.out.reports
                .map { meta, report, _readLength -> tuple(meta.id, meta, report) }
                .join(
                    SELECT_BRACKEN_LENGTH.out.selections.map { meta, selection -> tuple(meta.id, selection) },
                    failOnDuplicate: true,
                    failOnMismatch: true
                )
                .map { _id, meta, report, selection -> tuple(meta, report, selection) }
                .combine(db)
            BRACKEN_PROFILE(bracken_input)
            bracken_tables = BRACKEN_PROFILE.out.tables
                .map { _meta, table, _report -> table }
                .toSortedList { a, b -> a.name <=> b.name }
            classifier_qc = classifier_qc
                .mix(BRACKEN_PROFILE.out.tables.map { _meta, _table, report -> report })
                .mix(SELECT_BRACKEN_LENGTH.out.selections.map { _meta, selection -> selection })
                .mix(BRACKEN_PROFILE.out.status)
        } else {
            bracken_tables = channel.empty()
        }
    }

    // Kraken2 reports alone make no count table; Bracken or KrakenUniq do.
    if (krakenUniq || params.run_bracken) {
        BUILD_FEATURE_TABLE(
            metadata, bracken_tables, reports, minimizer_reports, taxdump_nodes, binScript('build_feature_table.py')
        )
        table_qc = BUILD_FEATURE_TABLE.out.qc
            .mix(BUILD_FEATURE_TABLE.out.sample_qc)
            .mix(BUILD_FEATURE_TABLE.out.kmer_filter)
        all_counts = BUILD_FEATURE_TABLE.out.raw_counts
        filtered_counts = BUILD_FEATURE_TABLE.out.counts
        sample_qc = BUILD_FEATURE_TABLE.out.sample_qc
    } else {
        table_qc = channel.empty()
        all_counts = channel.empty()
        filtered_counts = channel.empty()
        sample_qc = channel.empty()
    }

    emit:
    all_counts                            // one row per taxon, before the taxon filter
    filtered_counts                       // after the taxon filter
    sample_qc                             // per-sample classified and retained reads
    qc = classifier_qc.mix(table_qc)      // for MultiQC
}
