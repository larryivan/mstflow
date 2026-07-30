process MULTIQC_REPORT {
    label 'process_low'
    tag 'multiqc'
    container "${params.container_multiqc}"

    publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

    input:
    path files, stageAs: 'qc??/*'

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data", emit: data
    path "multiqc.versions.yml", emit: versions

    script:
    """
    set +e
    multiqc . \\
        --filename multiqc_report.html \\
        --outdir . \\
        --title '${params.multiqc_title}' \\
        > multiqc.stdout.log 2> multiqc.stderr.log
    status=\$?
    set -e

    if [[ "\$status" -ne 0 ]]; then
        cat multiqc.stdout.log >&2
        cat multiqc.stderr.log >&2
        exit "\$status"
    fi

    mkdir -p multiqc_data
    if [[ ! -s multiqc_report.html ]]; then
        printf '<!doctype html>\\n<html><body><h1>mstflow MultiQC</h1><p>No recognized QC files were found.</p></body></html>\\n' > multiqc_report.html
        echo -e "metric\\tvalue\\nrecognized_results\\t0" > multiqc_data/mstflow_placeholder.tsv
    fi

    printf '"${task.process}":\\n    multiqc: "%s"\\n' "\$(multiqc --version | awk '{print \$NF}')" > multiqc.versions.yml
    """

    stub:
    """
    mkdir -p multiqc_data
    printf '<!doctype html><html><body>stub multiqc</body></html>\\n' > multiqc_report.html
    echo -e "metric\\tvalue\\nstub\\ttrue" > multiqc_data/stub.tsv
    printf '"${task.process}":\\n    stub: "true"\\n' > multiqc.versions.yml
    """
}
