/*
 * What the user gives: the samples to trace (--input), their own reference samples
 * (--sources) and the sources of an earlier run (--panel). All of it is read and checked
 * here, in the controller, before any task starts, and every problem is reported at once.
 *
 * Which file a library comes from decides its role, so the sheets carry no role column:
 * every --input row is a sample to trace, every --sources row a source.
 */

include { isOn } from './utils'


def isNa(value) {
    return value == null || value.toString().trim().toUpperCase() in ['', 'NA', 'N/A', 'NONE', 'NULL']
}


// A sheet cell, or '' when it holds no value.
def cell(row, column) {
    return isNa(row[column]) ? '' : row[column]
}


def isSafeId(value) {
    return value ==~ /^[A-Za-z0-9][A-Za-z0-9._-]*$/
}


// Environments are compared without case or surrounding spaces: Oral and oral are one.
def environmentName(value) {
    return value.toString().trim().toLowerCase()
}


// The columns the workflow writes for every library, in this order. Columns the user
// adds follow them and reach the attribution map.
def libraryColumns() {
    return ['sample_id', 'fastq_1', 'fastq_2', 'layout', 'age_group', 'preprocess_mode',
            'host_action', 'host_ref', 'host_group']
}


def labelColumns() {
    return ['source_sink', 'model_source_env', 'unit_id', 'project']
}


// Modern libraries are trimmed by fastp; ancient ones by AdapterRemoval, which also
// collapses overlapping read pairs into the fragments they came from.
def preprocessMode(ageGroup, layout) {
    return [
        'modern PE': 'fastp_pe', 'modern SE': 'fastp_se',
        'ancient PE': 'adapterremoval_pe_collapse', 'ancient SE': 'adapterremoval_se',
    ]["${ageGroup} ${layout}".toString()] ?: 'NA'
}


/*
 * Reads --input, --sources and --panel. Returns the libraries, each as
 * [fields: <map of column to text>, reads: [<FASTQ paths>]], and the panel, as
 * [counts: <file>, sources: <file>, rows: [<map per source>]] or null. Fails with every
 * problem found, one per line.
 */
def readInputs() {
    def errors = []
    def libraries = []
    if (params.input) {
        libraries += readSheet(params.input, 'sink', errors)
    }
    if (params.sources) {
        libraries += readSheet(params.sources, 'source', errors)
    }
    def panel = params.panel ? readPanel(params.panel, errors) : null

    def repeated = libraries.groupBy { library -> library.fields.sample_id }.findAll { _id, copies -> copies.size() > 1 }
    repeated.each { id, _copies -> errors << "sample_id '${id}' is used more than once across --input and --sources".toString() }
    libraries.findAll { library -> !(library.fields.sample_id in repeated) }.groupBy { library -> library.fields.unit_id }.each { unit, members ->
        def labels = members.collect { library -> [library.fields.source_sink, library.fields.model_source_env] }.unique()
        if (labels.size() > 1) {
            errors << "unit_id '${unit}' joins libraries with different roles or environments: ${members*.fields*.sample_id.join(', ')}".toString()
        }
    }
    // A run with --sources writes its panel to <outdir>/panel, which must not be the one it reads.
    if (panel && params.sources && panel.sources.parent.normalize() == file("${params.outdir}/panel").normalize()) {
        errors << "--panel ${params.panel} would be overwritten by this run's panel; copy it elsewhere or choose another --outdir".toString()
    }
    if (panel) {
        def taken = (libraries*.fields*.sample_id + libraries*.fields*.unit_id) as Set
        def clashes = panel.rows*.source_id.findAll { id -> id in taken }
        if (clashes) {
            errors << "the panel already has source(s) named ${clashes.join(', ')}; rename them in your sheets".toString()
        }
    }

    if (!params.input && !params.sources) {
        errors << 'nothing to analyse: give --input (samples to trace), --sources (reference samples), or both'
    }
    if (errors) {
        error((['MSTFlow cannot start:'] + errors.collect { problem -> "  ${problem}".toString() }).join('\n'))
    }
    return [libraries: libraries, panel: panel]
}


// The rows of a TSV or CSV file (told apart by its header line), as maps of trimmed text,
// each with its line number. Blank lines are skipped.
def tableRows(sheet) {
    def header = sheet.readLines().find { line -> line.trim() } ?: ''
    def separator = header.contains('\t') ? '\t' : ','
    def rows = sheet.splitCsv(header: true, sep: separator, quote: '"')
    return rows.withIndex().collect { row, index ->
        def fields = row.collectEntries { key, value ->
            [(key.toString().trim()): value == null ? '' : value.toString().replaceAll(/[\t\r\n]+/, ' ').trim()]
        }
        [line: index + 2, fields: fields]
    }.findAll { row -> row.fields.values().any { value -> value } }
}


def readSheet(location, role, errors) {
    def option = role == 'sink' ? '--input' : '--sources'
    def sheet = file(location)
    if (sheet instanceof List || !sheet.exists()) {
        errors << "${option} ${location}: no such file".toString()
        return []
    }
    if (sheet.name ==~ /(?i).+\.(xlsx|xls|ods)$/) {
        errors << "${option} ${sheet.name} is a spreadsheet; save it as tab- or comma-separated text".toString()
        return []
    }
    def rows = tableRows(sheet)
    def columns = rows ? rows[0].fields.keySet() : (sheet.readLines().find { line -> line.trim() } ?: '').split(/[\t,]/)*.trim()
    def problems = []
    def missing = ['sample_id', 'fastq_1', 'age_group'] + (role == 'source' ? ['environment'] : [])
    missing = missing.findAll { column -> !(column in columns) }
    if (missing) {
        problems << "missing column(s) ${missing.join(', ')}".toString()
    }
    if (role == 'sink' && columns.any { column -> column in ['source_sink', 'model_source_env', 'environment'] }) {
        problems << ('lists samples to trace, but has source labels (source_sink, model_source_env or environment); ' +
            'give the reference samples in their own sheet with --sources')
    }
    if (!rows && !missing) {
        problems << 'has no rows'
    }
    if (problems) {
        problems.each { problem -> errors << "${option} ${sheet.name}: ${problem}".toString() }
        return []
    }
    return rows.collect { row -> readLibrary(row.fields, role, sheet, "${option} ${sheet.name} line ${row.line}", errors) }
        .findAll { library -> library }
}


// One samplesheet row as a library, with the defaults filled in, or null after adding
// its problems to errors.
def readLibrary(row, role, sheet, where, errors) {
    def problems = []

    def sampleId = cell(row, 'sample_id')
    if (!sampleId) {
        problems << 'sample_id is empty'
    } else if (!isSafeId(sampleId)) {
        problems << "sample_id '${sampleId}' may use only letters, digits, '.', '_' and '-', and must start with a letter or digit"
    }
    def ageGroup = cell(row, 'age_group').toLowerCase()
    if (!(ageGroup in ['ancient', 'modern'])) {
        problems << "age_group must be ancient or modern, not '${row.age_group}'"
    }

    def fastq1 = cell(row, 'fastq_1')
    def fastq2 = cell(row, 'fastq_2')
    def layout = cell(row, 'layout').toUpperCase() ?: (fastq2 ? 'PE' : 'SE')
    if (!(layout in ['SE', 'PE'])) {
        problems << "layout must be SE or PE, or empty to infer it from fastq_2, not '${row.layout}'"
    }
    def reads = []
    if (!fastq1) {
        problems << 'fastq_1 is empty'
    } else {
        reads << readFile(fastq1, sheet, problems)
    }
    if (layout == 'PE' && !fastq2) {
        problems << 'layout is PE but fastq_2 is empty'
    } else if (layout == 'SE' && fastq2) {
        problems << 'layout is SE but fastq_2 is set'
    } else if (fastq2) {
        reads << readFile(fastq2, sheet, problems)
    }

    def hostAction = cell(row, 'host_action').toLowerCase()
    if (!hostAction && isOn(params.remove_host) && isOn(params.require_explicit_host_action)) {
        problems << 'host_action must be set to remove or skip'
    }
    hostAction = hostAction ?: 'remove'
    if (!(hostAction in ['remove', 'skip'])) {
        problems << "host_action must be remove or skip, not '${row.host_action}'"
    }
    def hostRef = cell(row, 'host_ref')
    def hostGroup = cell(row, 'host_group').toLowerCase()
    def defaultGroup = params.default_host_group.toString().trim().toLowerCase()
    if (hostAction == 'skip') {
        if (hostRef) {
            problems << 'host_action is skip, so host_ref must be empty'
        }
        hostRef = 'NA'
        hostGroup = hostGroup ?: 'NA'
    } else {
        hostRef = hostRef ?: 'default'
        hostGroup = hostGroup ?: defaultGroup
        if (isOn(params.remove_host) && isOn(params.require_explicit_host_for_nondefault)
                && hostGroup != defaultGroup && hostRef == 'default') {
            problems << "host_group '${hostGroup}' needs its own host_ref; --host_fasta is the ${defaultGroup} genome"
        }
    }

    def environment = role == 'source' ? environmentName(cell(row, 'environment')) : 'NA'
    if (role == 'source' && !environment) {
        problems << 'environment is empty'
    } else if (environment in ['unknown', 'sampleid']) {
        problems << "environment '${environment}' is reserved; choose another name"
    }
    def unitId = cell(row, 'unit_id')
    if (unitId && !isSafeId(unitId)) {
        problems << "unit_id '${unitId}' may use only letters, digits, '.', '_' and '-'"
    }
    unitId = unitId ?: sampleId

    if (problems) {
        problems.each { problem -> errors << "${where}: ${problem}".toString() }
        return null
    }
    def known = libraryColumns() + labelColumns() + ['environment']
    def fields = [
        sample_id: sampleId, fastq_1: fastq1, fastq_2: fastq2 ?: 'NA', layout: layout, age_group: ageGroup,
        preprocess_mode: preprocessMode(ageGroup, layout), host_action: hostAction, host_ref: hostRef,
        host_group: hostGroup, source_sink: role, model_source_env: environment, unit_id: unitId,
        project: cell(row, 'project') ?: 'NA',
    ] + row.findAll { column, _value -> !(column in known) }
    return [fields: fields, reads: reads]
}


// A read file named in a sheet. Relative paths are read relative to the sheet, so a sheet
// and its reads can move together, then relative to the launch directory; URLs and
// absolute paths are used as given.
def readFile(text, sheet, problems) {
    if (text ==~ /^[A-Za-z][A-Za-z0-9+.-]*:\/\/.+/ || java.nio.file.Paths.get(text).isAbsolute()) {
        def given = file(text)
        if (!given.exists()) {
            problems << "${text} was not found"
        }
        return given
    }
    def besideSheet = sheet.parent.resolve(text)
    if (besideSheet.exists()) {
        return besideSheet
    }
    def fromLaunchDir = file(text)
    if (fromLaunchDir.exists()) {
        return fromLaunchDir
    }
    problems << "${text} was not found next to the sheet or in the launch directory"
    return besideSheet
}


// A panel is a directory, or a URL, holding counts.tsv and sources.tsv as MSTFlow writes
// them to results/panel.
def readPanel(location, errors) {
    def base = location.toString().replaceAll('/+$', '')
    def counts = file("${base}/counts.tsv")
    def sources = file("${base}/sources.tsv")
    def absent = [counts, sources].findAll { panelFile -> !panelFile.exists() }
    if (absent) {
        errors << "--panel ${location} has no ${absent*.name.join(' or ')}; a panel is the panel/ directory of an earlier run".toString()
        return null
    }
    def rows = tableRows(sources)
    def missing = ['source_id', 'environment'].findAll { column -> !rows || !(column in rows[0].fields.keySet()) }
    if (missing) {
        errors << "--panel ${location}: sources.tsv is missing column(s) ${missing.join(', ')}".toString()
        return null
    }
    return [
        counts: counts,
        sources: sources,
        rows: rows.collect { row ->
            [source_id: row.fields.source_id, environment: environmentName(row.fields.environment),
             project: isNa(row.fields.project) ? 'NA' : row.fields.project]
        },
    ]
}


// Only the fields that decide how reads are processed go into a task's meta, and so into
// its cache key: relabelling a source reruns source attribution, not trimming or Kraken2.
def metaFromLibrary(library) {
    def fields = library.fields
    return [id: fields.sample_id, layout: fields.layout, age_group: fields.age_group,
            preprocess_mode: fields.preprocess_mode, host_action: fields.host_action, host_ref: fields.host_ref]
}


def tableText(columns, rows) {
    def lines = [columns.join('\t')] + rows.collect { row -> columns.collect { column -> row[column] ?: '' }.join('\t') }
    return lines.join('\n') + '\n'
}


// The libraries as the processing steps read them, without any label.
def librariesText(libraries) {
    return tableText(libraryColumns(), libraries*.fields)
}


// The libraries with their roles, labels and any column the user added.
def samplesheetText(libraries) {
    def columns = libraryColumns() + labelColumns()
    libraries.each { library -> library.fields.keySet().each { column -> if (!(column in columns)) { columns << column } } }
    return tableText(columns, libraries*.fields)
}


// What MSTFlow understood, printed before the run starts.
def describeInputs(inputs, attribute) {
    def libraries = inputs.libraries*.fields
    def samples = libraries.findAll { library -> library.source_sink == 'sink' }
    def sources = libraries.findAll { library -> library.source_sink == 'source' }
    def lines = ['']
    if (samples) {
        def ancient = samples.count { library -> library.age_group == 'ancient' }
        lines << "Samples to trace  ${samples.size()} (${ancient} ancient, ${samples.size() - ancient} modern)  ${params.input}".toString()
    }
    def references = sources.collect { library -> [id: library[params.attribution_aggregate_by], environment: library.model_source_env, project: library.project] }.unique()
    if (sources) {
        lines << "Sources           ${references.size()}  ${params.sources}".toString()
    }
    if (inputs.panel) {
        lines << "Panel sources     ${inputs.panel.rows.size()}  ${params.panel}".toString()
        references += inputs.panel.rows.collect { row -> [id: row.source_id, environment: row.environment, project: row.project] }
    }
    if (references) {
        def width = Math.max(11, references*.environment*.length().max())
        lines << "  ${'environment'.padRight(width)}  sources  studies".toString()
        references.groupBy { reference -> reference.environment }.sort().each { environment, members ->
            def studies = members*.project.findAll { project -> !isNa(project) }.unique().size()
            lines << "  ${environment.padRight(width)}  ${members.size().toString().padLeft(7)}  ${(studies ?: '-').toString().padLeft(7)}".toString()
        }
    }
    def methods = params.attribution_methods.toString().split(',')*.trim().findAll { method -> method }
    def scoring = isOn(params.attribution_cv) ? ", scored by holding out one ${params.attribution_cv_group} at a time" : ''
    if (!attribute) {
        lines << 'Source attribution is off: give --sources or --panel to trace where the samples came from.'
    } else if (!samples) {
        lines << "No samples to trace: the sources are profiled${scoring}.".toString()
    } else {
        lines << "Source attribution with ${methods.join(', ') ?: 'no method'}${scoring}.".toString()
    }
    if (params.sources) {
        lines << "The sources used will be saved to ${params.outdir}/panel, for reuse with --panel.".toString()
    }
    return (lines + ['']).join('\n')
}
