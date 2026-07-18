process PLOT_SOURCETRACKER2 {
    label 'process_low'
    tag 'sourcetracker2_plot'

    publishDir "${params.outdir}/source_attribution/sourcetracker2/${params.get('sourcetracker2_analysis_label', 'default')}", mode: 'copy'

    input:
    path proportions
    path mapping
    path plot_sourcetracker2_script

    output:
    path "sourcetracker2_assignments.tsv", emit: assignments
    path "sourcetracker2_stacked_bar.png", emit: plot
    path "plot_sourcetracker2_versions.yml", emit: versions

    script:
    def maxUnknown = params.get('sourcetracker2_max_unknown_prop', 0.5)
    def minKnown = params.get('sourcetracker2_min_top_known_prop', 0.5)
    def minMargin = params.get('sourcetracker2_min_top_known_margin', 0.05)
    def looArg = params.get('sourcetracker2_loo', false).toString().toBoolean() ? '--loo' : ''
    """
    python3 ${plot_sourcetracker2_script} \\
        --proportions ${proportions} \\
        --mapping ${mapping} \\
        --assignments sourcetracker2_assignments.tsv \\
        --plot sourcetracker2_stacked_bar.png \\
        --max-unknown-prop ${maxUnknown} \\
        --min-top-known-prop ${minKnown} \\
        --min-top-known-margin ${minMargin} \\
        ${looArg} \\
        --title '${params.multiqc_title} SourceTracker2'

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > plot_sourcetracker2_versions.yml
    """

    stub:
    """
    printf 'SampleID\\texpected_env\\ttop_known_env\\ttop_known_prop\\tassignment\\tabstained\\tabstention_reason\\n' > sourcetracker2_assignments.tsv
    python3 - <<'PY'
import base64
from pathlib import Path
Path('sourcetracker2_stacked_bar.png').write_bytes(base64.b64decode(
    'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII='
))
PY
    printf '"${task.process}":\\n    stub: "true"\\n' > plot_sourcetracker2_versions.yml
    """
}


process PLOT_SOURCETRACKER2_WITH_UNCERTAINTY {
    label 'process_low'
    tag 'sourcetracker2_plot_with_uncertainty'

    publishDir "${params.outdir}/source_attribution/sourcetracker2/${params.get('sourcetracker2_analysis_label', 'default')}", mode: 'copy'

    input:
    path proportions
    path uncertainty
    path mapping
    path plot_sourcetracker2_script

    output:
    path "sourcetracker2_assignments.tsv", emit: assignments
    path "sourcetracker2_stacked_bar.png", emit: plot
    path "plot_sourcetracker2_versions.yml", emit: versions

    script:
    def maxUnknown = params.get('sourcetracker2_max_unknown_prop', 0.5)
    def minKnown = params.get('sourcetracker2_min_top_known_prop', 0.5)
    def minMargin = params.get('sourcetracker2_min_top_known_margin', 0.05)
    def looArg = params.get('sourcetracker2_loo', false).toString().toBoolean() ? '--loo' : ''
    """
    python3 ${plot_sourcetracker2_script} \\
        --proportions ${proportions} \\
        --uncertainty ${uncertainty} \\
        --mapping ${mapping} \\
        --assignments sourcetracker2_assignments.tsv \\
        --plot sourcetracker2_stacked_bar.png \\
        --max-unknown-prop ${maxUnknown} \\
        --min-top-known-prop ${minKnown} \\
        --min-top-known-margin ${minMargin} \\
        ${looArg} \\
        --title '${params.multiqc_title} SourceTracker2'

    printf '"${task.process}":\\n    python: "%s"\\n' "\$(python3 --version 2>&1 | sed 's/Python //')" > plot_sourcetracker2_versions.yml
    """

    stub:
    """
    printf 'SampleID\\texpected_env\\ttop_known_env\\ttop_known_prop\\tassignment\\tabstained\\tabstention_reason\\n' > sourcetracker2_assignments.tsv
    python3 - <<'PY'
import base64
from pathlib import Path
Path('sourcetracker2_stacked_bar.png').write_bytes(base64.b64decode(
    'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNk+A8AAQUBAScY42YAAAAASUVORK5CYII='
))
PY
    printf '"${task.process}":\\n    stub: "true"\\n' > plot_sourcetracker2_versions.yml
    """
}
