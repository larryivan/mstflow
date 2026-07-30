process SOURCETRACKER2_GIBBS {
    label 'process_medium'
    tag "sourcetracker2:${params.get('sourcetracker2_analysis_label', 'default')}"
    container "${params.container_sourcetracker2}"
    publishDir "${params.outdir}/source_attribution/sourcetracker2/${params.get('sourcetracker2_analysis_label', 'default')}", mode: 'copy'

    input:
    path biom
    path mapping
    path seeded_runner

    output:
    path "sourcetracker2.out", emit: outdir
    path "mixing_proportions.txt", emit: proportions
    path "mixing_proportions_stds.txt", emit: uncertainty
    path "sourcetracker2_versions.yml", emit: versions

    script:
    def extraArgs = params.sourcetracker2_extra_args ?: ''
    def looEnabled = params.get('sourcetracker2_loo', false).toString().toBoolean()
    def looArg = looEnabled ? '--loo' : ''
    def seed = params.sourcetracker2_seed as int
    def sourceDepth = looEnabled ? 0 : params.sourcetracker2_source_rarefaction_depth
    def sinkDepth = looEnabled ? 0 : params.sourcetracker2_sink_rarefaction_depth
    """
    rm -rf sourcetracker2.out
    export PYTHONHASHSEED=${seed}

    python3 ${seeded_runner} --seed ${seed} gibbs \\
        -i ${biom} \\
        -m ${mapping} \\
        -o sourcetracker2.out \\
        --source_sink_column SourceSink \\
        --source_column_value source \\
        --sink_column_value sink \\
        --source_category_column Env \\
        --source_rarefaction_depth ${sourceDepth} \\
        --sink_rarefaction_depth ${sinkDepth} \\
        --jobs 1 \\
        ${looArg} \\
        ${extraArgs}

    mean_file=''
    for candidate in \\
        sourcetracker2.out/mixing_proportions.txt \\
        sourcetracker2.out/mixing_proporitions.txt
    do
        if [[ -s "\$candidate" ]]; then
            mean_file="\$candidate"
            break
        fi
    done
    if [[ -z "\$mean_file" ]]; then
        echo "SourceTracker2 did not produce mixing proportions" >&2
        exit 1
    fi
    cp "\$mean_file" mixing_proportions.txt

    std_file=''
    for candidate in \\
        sourcetracker2.out/mixing_proportions_stds.txt \\
        sourcetracker2.out/mixing_proportion_stds.txt \\
        sourcetracker2.out/mixing_proporitions_stds.txt \\
        sourcetracker2.out/mixing_proporition_stds.txt
    do
        if [[ -s "\$candidate" ]]; then
            std_file="\$candidate"
            break
        fi
    done
    if [[ -z "\$std_file" ]]; then
        echo "SourceTracker2 did not produce mixing-proportion standard deviations" >&2
        exit 1
    fi
    cp "\$std_file" mixing_proportions_stds.txt

    printf '"${task.process}":\\n    sourcetracker2: "%s"\\n    seed: "${seed}"\\n    jobs: "1"\\n    source_rarefaction_depth: "${sourceDepth}"\\n    sink_rarefaction_depth: "${sinkDepth}"\\n' "\$(sourcetracker2 --version 2>&1 | head -n 1 || true)" > sourcetracker2_versions.yml
    """

    stub:
    def looEnabled = params.get('sourcetracker2_loo', false).toString().toBoolean()
    def seed = params.sourcetracker2_seed as int
    """
    mkdir -p sourcetracker2.out
    python3 - <<'PY'
import csv

loo = '${looEnabled}' == 'true'
sources = []
source_rows = []
sinks = []
with open('${mapping}', newline='') as handle:
    for row in csv.DictReader(handle, delimiter='\\t'):
        if row['SourceSink'] == 'source':
            if row['Env'] not in sources:
                sources.append(row['Env'])
            source_rows.append(row)
        elif row['SourceSink'] == 'sink':
            sinks.append(row)
targets = source_rows if loo else sinks
columns = [*sources, 'Unknown']
with open('mixing_proportions.txt', 'w', newline='') as means, open(
    'mixing_proportions_stds.txt', 'w', newline=''
) as stds:
    mean_writer = csv.writer(means, delimiter='\\t', lineterminator='\\n')
    std_writer = csv.writer(stds, delimiter='\\t', lineterminator='\\n')
    mean_writer.writerow(['SampleID', *columns])
    std_writer.writerow(['SampleID', *columns])
    for target in targets:
        values = [0.0] * len(columns)
        if loo and target['Env'] in sources:
            values[sources.index(target['Env'])] = 1.0
        else:
            values[-1] = 1.0
        mean_writer.writerow([target['#SampleID'], *values])
        std_writer.writerow([target['#SampleID'], *([0.0] * len(columns))])
PY
    cp mixing_proportions.txt sourcetracker2.out/mixing_proportions.txt
    cp mixing_proportions_stds.txt sourcetracker2.out/mixing_proportions_stds.txt
    printf '"${task.process}":\\n    stub: "true"\\n    seed: "${seed}"\\n    jobs: "1"\\n' > sourcetracker2_versions.yml
    """
}
