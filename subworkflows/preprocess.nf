include { FASTQC as FASTQC_RAW } from '../modules/fastqc'
include { FASTQC as FASTQC_CLEAN } from '../modules/fastqc'
include { PREPROCESS_FASTP } from '../modules/preprocess_fastp'
include { ADAPTERREMOVAL_RAW } from '../modules/preprocess_adapterremoval'
include { NORMALIZE_ADAPTERREMOVAL } from '../modules/normalize_adapterremoval'
include { HOST_ALIGN_BWA } from '../modules/host_remove_bwa'
include { FILTER_HOST_ALIGNMENTS } from '../modules/filter_host_alignments'
include { binScript; hostBundleFor; hostReferenceFor; isOn } from './utils'

/*
 * Trim every library, collapse overlapping ancient read pairs into the molecules they
 * came from, and remove host reads, with FastQC before and after. Modern libraries go
 * through fastp and ancient ones through AdapterRemoval; both leave one FASTQ per
 * fragment class, counted in molecules.
 */
workflow PREPROCESS {
    take:
    samples      // tuple(meta, [reads])

    main:
    def fragmentHelper = channel.value(binScript('mstflow_fragment_tools.py'))

    FASTQC_RAW(samples.map { meta, reads -> tuple(meta + [qc_stage: 'raw'], reads) })

    trimming = samples.branch { meta, _reads ->
        fastp: meta.preprocess_mode.startsWith('fastp')
        adapterremoval: meta.preprocess_mode.startsWith('adapterremoval')
    }
    PREPROCESS_FASTP(trimming.fastp)
    ADAPTERREMOVAL_RAW(trimming.adapterremoval)
    NORMALIZE_ADAPTERREMOVAL(
        ADAPTERREMOVAL_RAW.out.raw,
        binScript('normalize_adapterremoval.py'),
        fragmentHelper
    )
    trimmed = PREPROCESS_FASTP.out.reads.mix(NORMALIZE_ADAPTERREMOVAL.out.reads)

    if (isOn(params.remove_host)) {
        // Fail closed rather than silently dropping samples. `host_action` is the
        // one branch key whose input column is optional, so a validated manifest
        // written before host policies existed would otherwise lose every row here.
        host_views = trimmed.branch { meta, _reads ->
            remove: meta.host_action == 'remove'
            skip: meta.host_action == 'skip'
            undeclared: true
        }
        undeclared = host_views.undeclared.map { meta, _reads ->
            error "Sample '${meta.id}' has host_action '${meta.host_action}'; validated metadata must declare remove or skip"
        }
        host_bundles = samples
            .filter { meta, _reads -> meta.host_action == 'remove' }
            .map { meta, _reads -> tuple(meta.id, hostBundleFor(hostReferenceFor(meta, params.host_fasta))) }
        host_input = host_views.remove
            .map { meta, reads -> tuple(meta.id, meta, reads) }
            .join(host_bundles)
            .map { _id, meta, reads, bundle -> tuple(meta, reads, bundle) }
        HOST_ALIGN_BWA(host_input)
        FILTER_HOST_ALIGNMENTS(
            HOST_ALIGN_BWA.out.alignments,
            binScript('filter_host_alignments.py'),
            fragmentHelper
        )
        clean = FILTER_HOST_ALIGNMENTS.out.reads.mix(host_views.skip).mix(undeclared)
        host_qc = HOST_ALIGN_BWA.out.logs.mix(FILTER_HOST_ALIGNMENTS.out.qc)
    } else {
        clean = trimmed
        host_qc = channel.empty()
    }

    FASTQC_CLEAN(clean.map { meta, reads -> tuple(meta + [qc_stage: 'clean'], reads) })

    trimming_qc = PREPROCESS_FASTP.out.logs.flatten()
        .mix(ADAPTERREMOVAL_RAW.out.logs.flatten())
        .mix(NORMALIZE_ADAPTERREMOVAL.out.qc)
    fastqc = FASTQC_RAW.out.qc.mix(FASTQC_CLEAN.out.qc)
        .map { _meta, zip, html -> [zip, html] }
        .flatten()

    emit:
    reads       = clean                                  // tuple(meta, [reads]) after host removal
    qc          = fastqc.mix(trimming_qc).mix(host_qc)   // every QC file, for MultiQC
    step_counts = trimming_qc.mix(host_qc)               // molecules kept by each step
}
