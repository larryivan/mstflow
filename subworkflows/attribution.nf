include { PREPARE_ATTRIBUTION_INPUTS } from '../modules/prepare_attribution_inputs'
include { SPLIT_ATTRIBUTION_FOLDS } from '../modules/split_attribution_folds'
include { SOURCETRACKER2 } from '../modules/sourcetracker2'
include { FEAST } from '../modules/feast'
include { FASTST } from '../modules/fastst'
include { R_METHODS_INSTALL } from '../modules/r_methods_install'
include { SOURCEID_NMF } from '../modules/sourceid_nmf'
include { PREPARE_SOURCEPREDICT_TAXONOMY; SOURCEPREDICT } from '../modules/sourcepredict'
include { SUMMARIZE_ATTRIBUTION } from '../modules/summarize_attribution'
include { EVALUATE_ATTRIBUTION } from '../modules/evaluate_attribution'
include { binScript; isOn; splitList; taxdumpFile; taxdumpPath } from './utils'

/*
 * Estimate where each sink's microbes came from with every requested method.
 *
 * All methods read the same inputs and write the same result table, so they can be
 * compared directly. With --attribution_cv, each is also run on cross-validation
 * folds in which whole groups of sources (by default, studies) are held out as
 * sinks of known origin, and scored on them.
 */
workflow ATTRIBUTION {
    take:
    counts        // taxon count table (all taxa or after the taxon filter)
    metadata      // validated samplesheet
    sample_qc     // per-sample classification QC from the profiling step

    main:
    def methods = splitList(params.attribution_methods)
    def attributionModule = binScript('mstflow_attribution.py')
    def attributionHelpersR = [binScript('mstflow_attribution.R'), attributionModule]

    PREPARE_ATTRIBUTION_INPUTS(counts, metadata, sample_qc, binScript('prepare_attribution_inputs.py'))
    main_analysis = PREPARE_ATTRIBUTION_INPUTS.out.counts
        .combine(PREPARE_ATTRIBUTION_INPUTS.out.map)
        .map { table, mapping -> tuple([analysis: 'main'], table, mapping) }

    if (isOn(params.attribution_cv)) {
        SPLIT_ATTRIBUTION_FOLDS(
            PREPARE_ATTRIBUTION_INPUTS.out.counts,
            PREPARE_ATTRIBUTION_INPUTS.out.map,
            binScript('split_attribution_folds.py'),
            attributionModule
        )
        fold_counts = SPLIT_ATTRIBUTION_FOLDS.out.counts.flatten()
            .map { table -> tuple(table.name - ~/\.counts\.tsv$/, table) }
        fold_maps = SPLIT_ATTRIBUTION_FOLDS.out.maps.flatten()
            .map { mapping -> tuple(mapping.name - ~/\.map\.tsv$/, mapping) }
        folds = fold_counts
            .join(fold_maps, failOnMismatch: true, failOnDuplicate: true)
            .map { fold, table, mapping -> tuple([analysis: fold], table, mapping) }
        analyses = main_analysis.mix(folds)
    } else {
        analyses = main_analysis
    }

    // The attribution image holds FEAST and FastST. Without it, their R packages are built
    // from source, once, and handed to the tasks as a library.
    if (isOn(params.r_methods_from_source) && ('feast' in methods || 'fastst' in methods)) {
        R_METHODS_INSTALL(file(params.feast_code, checkIfExists: true), file(params.npreg_code, checkIfExists: true))
        r_library = R_METHODS_INSTALL.out.library
        fastst_code = channel.value(file(params.fastst_code, checkIfExists: true))
    } else {
        r_library = channel.value([])
        fastst_code = channel.value([])
    }

    results = channel.empty()
    if ('sourcetracker2' in methods) {
        SOURCETRACKER2(analyses, binScript('run_sourcetracker2_seeded.py'), binScript('sourcetracker2_results.py'),
                       attributionModule)
        results = results.mix(SOURCETRACKER2.out.results)
    }
    if ('feast' in methods) {
        FEAST(analyses, binScript('run_feast.R'), attributionHelpersR, r_library)
        results = results.mix(FEAST.out.results)
    }
    if ('fastst' in methods) {
        FASTST(analyses, binScript('run_fastst.R'), attributionHelpersR, fastst_code, r_library)
        results = results.mix(FASTST.out.results)
    }
    if ('sourceid_nmf' in methods) {
        SOURCEID_NMF(analyses, sourceidNmfCode(params.sourceid_nmf_code), binScript('run_sourceid_nmf.py'),
                     attributionModule)
        results = results.mix(SOURCEID_NMF.out.results)
    }
    if ('sourcepredict' in methods) {
        def purpose = 'SourcePredict, which needs the NCBI taxonomy'
        PREPARE_SOURCEPREDICT_TAXONOMY(
            ['nodes.dmp', 'names.dmp'].collect { name -> taxdumpFile(params.taxdump, name, purpose) } +
                [taxdumpPath(params.taxdump, 'merged.dmp')].findAll { merged -> merged.exists() }
        )
        SOURCEPREDICT(analyses, PREPARE_SOURCEPREDICT_TAXONOMY.out.database.collect(), binScript('run_sourcepredict.py'),
                      attributionModule)
        results = results.mix(SOURCEPREDICT.out.results)
    }

    // Sorted by name so the staged inputs, and so the cache key, do not depend on
    // which task finished first.
    by_analysis = results.branch { meta, _result ->
        sinks: meta.analysis == 'main'
        folds: true
    }
    SUMMARIZE_ATTRIBUTION(
        by_analysis.sinks.map { _meta, result -> result }.toSortedList { a, b -> a.name <=> b.name },
        PREPARE_ATTRIBUTION_INPUTS.out.map,
        binScript('summarize_attribution.py'),
        attributionModule
    )
    if (isOn(params.attribution_cv)) {
        EVALUATE_ATTRIBUTION(
            by_analysis.folds.map { _meta, result -> result }.toSortedList { a, b -> a.name <=> b.name },
            SPLIT_ATTRIBUTION_FOLDS.out.maps.flatten().toSortedList { a, b -> a.name <=> b.name },
            SPLIT_ATTRIBUTION_FOLDS.out.summary,
            binScript('evaluate_attribution.py'),
            attributionModule
        )
    }
}


def sourceidNmfCode(location) {
    // The code is staged file by file from a local directory or, by default, from
    // the authors' repository at a pinned commit, so a run never picks up changes.
    def names = ['SourceID-NMF.py', 'admm.py', 'data_clustering.py', 'data_estimation.py', 'data_initialization.py']
    def base = location.toString().replaceAll('/+$', '')
    return names.collect { name -> file("${base}/${name}", checkIfExists: true) }
}
