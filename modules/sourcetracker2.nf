process SOURCETRACKER2_GIBBS {
    label 'process_medium'
    tag 'sourcetracker2'
    conda "${projectDir}/envs/sourcetracker2.yml"

    publishDir "${params.outdir}/source_attribution/sourcetracker2", mode: 'copy'

    input:
    path biom
    path mapping

    output:
    path "sourcetracker2.out", emit: outdir
    path "mixing_proportions.txt", emit: proportions
    path "versions.yml", emit: versions

    script:
    def extraArgs = params.sourcetracker2_extra_args ?: ''
    """
    rm -rf sourcetracker2.out

    sourcetracker2 gibbs \\
        -i ${biom} \\
        -m ${mapping} \\
        -o sourcetracker2.out \\
        --source_sink_column SourceSink \\
        --source_column_value source \\
        --sink_column_value sink \\
        --source_category_column Env \\
        --source_rarefaction_depth ${params.sourcetracker2_source_rarefaction_depth} \\
        --sink_rarefaction_depth ${params.sourcetracker2_sink_rarefaction_depth} \\
        --jobs ${task.cpus} \\
        ${extraArgs}

    if [[ -s sourcetracker2.out/mixing_proportions.txt ]]; then
        cp sourcetracker2.out/mixing_proportions.txt mixing_proportions.txt
    elif [[ -s sourcetracker2.out/mixing_proporitions.txt ]]; then
        cp sourcetracker2.out/mixing_proporitions.txt mixing_proportions.txt
    else
        echo "SourceTracker2 did not produce mixing proportions" >&2
        exit 1
    fi

    printf '"${task.process}":\\n    sourcetracker2: "%s"\\n' "\$(sourcetracker2 --version 2>&1 | head -n 1 || true)" > versions.yml
    """

    stub:
    """
    mkdir -p sourcetracker2.out
    python3 - <<'PY'
import csv
sources = []
sinks = []
with open('${mapping}', newline='') as handle:
    for row in csv.DictReader(handle, delimiter='\\t'):
        if row['SourceSink'] == 'source' and row['Env'] not in sources:
            sources.append(row['Env'])
        elif row['SourceSink'] == 'sink':
            sinks.append(row['#SampleID'])
sources.append('Unknown')
with open('mixing_proportions.txt', 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\\t', lineterminator='\\n')
    writer.writerow(['SampleID', *sources])
    for sink in sinks:
        values = [0.0] * len(sources)
        values[-1] = 1.0
        writer.writerow([sink, *values])
PY
    cp mixing_proportions.txt sourcetracker2.out/mixing_proportions.txt
    printf '"${task.process}":\\n    stub: "true"\\n' > versions.yml
    """
}
