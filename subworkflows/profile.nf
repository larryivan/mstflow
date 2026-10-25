include { UNTAR_DATABASE } from '../modules/untar_database'
include { PREPARE_KRAKEN_INPUT } from '../modules/prepare_kraken_input'
include { KRAKEN2_PROFILE } from '../modules/kraken2'
include { KRAKENUNIQ_PROFILE } from '../modules/krakenuniq'
include { BUILD_FEATURE_TABLE } from '../modules/build_feature_table'
include { binScript; isArchive } from './utils'

/*
 * Classify every fragment with Kraken2 or KrakenUniq (--classifier) and merge the
 * samples into one count table: for each taxon at one rank (--rank), the reads
 * classified to it or below it. Taxa with too few distinct minimizers (Kraken2) or
 * unique k-mers (KrakenUniq) can be removed sample by sample.
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
        UNTAR_DATABASE(database, params.classifier)
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

    BUILD_FEATURE_TABLE(
        metadata,
        reports.map { _meta, report, _readLength -> report }.toSortedList { a, b -> a.name <=> b.name },
        minimizer_reports,
        taxdump_nodes,
        binScript('build_feature_table.py')
    )

    emit:
    all_counts = BUILD_FEATURE_TABLE.out.raw_counts       // one row per taxon, before the taxon filter
    filtered_counts = BUILD_FEATURE_TABLE.out.counts      // after the taxon filter
    sample_qc = BUILD_FEATURE_TABLE.out.sample_qc         // per-sample classified and retained reads
    qc = reports.map { _meta, report, _readLength -> report }   // Kraken-style reports, for MultiQC
}
