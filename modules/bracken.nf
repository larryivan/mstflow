process BRACKEN_PROFILE {
    label 'process_medium'
    tag "${meta.id}"
    container "${params.container_bracken}"

    publishDir "${params.outdir}/bracken", mode: 'copy'

    input:
    tuple val(meta), path(report), path(read_length_selection), path(bracken_db)

    output:
    tuple val(meta), path("*.bracken.tsv"), path("*.bracken.report.txt"), emit: tables
    path "*.bracken.status.tsv", emit: status
    path "*.bracken.versions.yml", emit: versions

    script:
    def id = meta.id
    """
    selected_read_length=\$(awk -F '\t' 'NR == 2 { print \$3 }' "${read_length_selection}")

    set +e
    bracken \
        ${params.bracken_extra_args} \
        -d "${bracken_db}" \
        -i "${report}" \
        -o "${id}.bracken.tsv" \
        -w "${id}.bracken.report.txt" \
        -r "\$selected_read_length" \
        -l "${params.bracken_rank}" \
        -t "${params.bracken_threshold}" \
        > "${id}.bracken.stdout.log" \
        2> "${id}.bracken.stderr.log"
    bracken_exit=\$?
    set -e

    cat "${id}.bracken.stdout.log"

    if [[ "\$bracken_exit" -eq 0 ]]; then
        cat "${id}.bracken.stderr.log" >&2
        expected_header=\$'name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads'
        actual_header=\$(head -n 1 "${id}.bracken.tsv" 2>/dev/null || true)
        if [[ "\$actual_header" != "\$expected_header" ]]; then
            printf 'ERROR: Bracken returned success but %s has an invalid header.\n' \
                "${id}.bracken.tsv" >&2
            exit 1
        fi
        if ! awk 'NR > 1 && NF { found=1 } END { exit(found ? 0 : 1) }' \
            "${id}.bracken.tsv"; then
            printf 'ERROR: Bracken returned success but %s has no estimate rows.\n' \
                "${id}.bracken.tsv" >&2
            exit 1
        fi
        printf 'sample_id\tstatus\treason\tbracken_rank\tthreshold\tselected_read_length\tbracken_exit_status\n%s\testimated\tabundance_estimated\t%s\t%s\t%s\t0\n' \
            "${id}" "${params.bracken_rank}" "${params.bracken_threshold}" "\$selected_read_length" \
            > "${id}.bracken.status.tsv"
    elif [[ "\$bracken_exit" -eq 1 ]] && awk '
        /^>> Checking report file: .+\$/ { next }
        /^Error: no reads found[.] Please check your Kraken report\$/ { matched += 1; next }
        /^[[:space:]]*\$/ { next }
        { unexpected=1 }
        END { exit(matched == 1 && !unexpected ? 0 : 1) }
    ' "${id}.bracken.stderr.log"; then
        # Preserve valid low-biomass/control samples as explicit zero-estimate
        # columns. Never lower the declared scientific threshold.
        printf 'name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n' \
            > "${id}.bracken.tsv"
        printf '0.00\t0\t0\tR\t1\troot\n' > "${id}.bracken.report.txt"
        printf 'sample_id\tstatus\treason\tbracken_rank\tthreshold\tselected_read_length\tbracken_exit_status\n%s\tzero_estimates\tno_estimates_at_requested_rank_threshold\t%s\t%s\t%s\t%s\n' \
            "${id}" "${params.bracken_rank}" "${params.bracken_threshold}" "\$selected_read_length" "\$bracken_exit" \
            > "${id}.bracken.status.tsv"
        printf 'NOTICE: Bracken produced no estimates for sample %s at rank %s with threshold %s; retaining a zero-estimate sample column.\n' \
            "${id}" "${params.bracken_rank}" "${params.bracken_threshold}" >&2
    else
        cat "${id}.bracken.stderr.log" >&2
        exit "\$bracken_exit"
    fi

    version=\$(bracken -v 2>&1 | sed -nE 's/.*[Bb]racken[^0-9]*([0-9]+(\\.[0-9]+){1,2}).*/\\1/p' | head -n 1 || true)
    printf '"${task.process}":\n    bracken: "%s"\n' "\$version" > "${id}.bracken.versions.yml"
    """

    stub:
    def id = meta.id
    """
    selected_read_length=\$(awk -F '\t' 'NR == 2 { print \$3 }' "${read_length_selection}")
    printf 'name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\nBacteroides\t816\tG\t1\t0\t1\t1.0\n' \
        > ${id}.bracken.tsv
    cp ${id}.bracken.tsv ${id}.bracken.report.txt
    printf 'sample_id\tstatus\treason\tbracken_rank\tthreshold\tselected_read_length\tbracken_exit_status\n${id}\tstub\tstub_run\t${params.bracken_rank}\t${params.bracken_threshold}\t%s\t0\n' \
        "\$selected_read_length" \
        > ${id}.bracken.status.tsv
    printf '"${task.process}":\n    stub: "true"\n' > ${id}.bracken.versions.yml
    """
}
