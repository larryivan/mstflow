include { BOWTIE2_BUILD; BOWTIE2_ALIGN } from '../modules/bowtie2'
include { SELECT_AUTHENTICATION_TARGETS; AUTHENTICATE } from '../modules/authentication'
include { binScript; taxdumpFile } from './utils'

/*
 * Check that the most abundant taxa of each ancient library are themselves ancient.
 * Its reads are aligned to the taxa's reference genomes (--authentication_genomes),
 * and each taxon is tested for the C-to-T damage at read ends that ancient DNA carries
 * and modern contamination does not.
 */
workflow AUTHENTICATION {
    take:
    reads       // tuple(meta, [reads]) after host removal
    counts      // filtered count table
    metadata    // validated samplesheet

    main:
    def genomes = file(params.authentication_genomes, checkIfExists: true)
    def index = bowtie2Index(genomes)
    if (index) {
        index_files = channel.value(index)
    } else {
        BOWTIE2_BUILD(genomes)
        index_files = BOWTIE2_BUILD.out.index
    }

    SELECT_AUTHENTICATION_TARGETS(
        counts,
        metadata,
        genomes,
        taxdumpFile(params.taxdump, 'nodes.dmp', 'authentication'),
        binScript('select_authentication_targets.py')
    )
    BOWTIE2_ALIGN(reads.filter { meta, _reads -> meta.age_group == 'ancient' }, index_files)
    AUTHENTICATE(
        BOWTIE2_ALIGN.out.bam,
        SELECT_AUTHENTICATION_TARGETS.out.targets,
        SELECT_AUTHENTICATION_TARGETS.out.references,
        binScript('authenticate_taxa.py')
    )
    AUTHENTICATE.out.results.collectFile(
        name: 'authentication.tsv', keepHeader: true, sort: { result -> result.name },
        storeDir: "${params.outdir}/authentication"
    )
}


// The Bowtie2 index beside the FASTA (<fasta>.1.bt2 and the rest, or .bt2l for large
// genomes), or null when there is none and it has to be built.
def bowtie2Index(fasta) {
    return ['bt2', 'bt2l'].findResult { extension ->
        def files = ['1', '2', '3', '4', 'rev.1', 'rev.2'].collect { part ->
            fasta.resolveSibling("${fasta.name}.${part}.${extension}")
        }
        files.every { index -> index.exists() } ? files : null
    }
}
