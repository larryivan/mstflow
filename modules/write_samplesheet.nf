/*
 * Write a checked samplesheet for the steps that read one. It runs in the controller, so
 * a table is rewritten, and the steps reading it rerun, only when its content changes.
 */
process WRITE_SAMPLESHEET {
    tag "${name}"
    publishDir "${params.outdir}/pipeline_info", mode: 'copy', saveAs: { filename -> filename == 'samplesheet.tsv' ? filename : null }

    input:
    tuple val(name), val(text)

    output:
    tuple val(name), path(name)

    exec:
    task.workDir.resolve(name).text = text
}
