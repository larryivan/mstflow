/*
 * Helpers for turning parameters and samplesheet rows into files and channels.
 * They run in the Nextflow controller, before any task starts, so a wrong path or
 * a missing index fails the run at once with a message that says how to fix it.
 */

def rejectReservedArgs(name, value, reservedFlags) {
    def text = value == null ? '' : value.toString()
    reservedFlags.each { flag ->
        def isShortOption = flag ==~ /^-[^-]$/
        def optionEnd = isShortOption ? '' : '(?==|[\\s\\x27\\x22;|&]|$)'
        def pattern = java.util.regex.Pattern.compile(
            '(^|[\\s\\x27\\x22;|&])' + java.util.regex.Pattern.quote(flag) +
                optionEnd
        )
        if (pattern.matcher(text).find()) {
            throw new IllegalArgumentException(
                "Parameter '${name}' must not override workflow-controlled option '${flag}'"
            )
        }
    }
}


// Values given on the command line arrive as text, so `--remove_host false` is the
// string "false", which Groovy treats as true. Every on/off parameter is read here.
def isOn(value) {
    return value instanceof Boolean ? value : value?.toString()?.trim()?.toLowerCase() == 'true'
}


def splitList(value) {
    return (value ?: '').toString().split(',')*.trim().findAll { item -> item }
}


def isArchive(value) {
    return value.toString() ==~ /.+\.(tar\.gz|tgz)/
}


def isLocal(path) {
    return path.toUri().scheme in [null, 'file']
}


def requireDirectory(name, value, requiredFiles, hint) {
    if (!value) {
        throw new IllegalArgumentException("Parameter '${name}' is required")
    }
    def directory = file(value, checkIfExists: true)
    if (!java.nio.file.Files.isDirectory(directory)) {
        throw new IllegalArgumentException("Parameter '${name}' must be a directory: ${directory}")
    }
    requiredFiles.each { filename ->
        def required = directory.resolve(filename)
        if (!java.nio.file.Files.isRegularFile(required)) {
            throw new IllegalArgumentException("${name} is missing required file: ${required}. ${hint}")
        }
    }
    return directory
}


// `taxdump` is a directory holding the NCBI .dmp files, or one of those files, in
// which case the others are looked up beside it. Either may be a URL.
def taxdumpPath(value, name) {
    def taxdump = file(value, checkIfExists: true)
    def isDirectory = isLocal(taxdump) && java.nio.file.Files.isDirectory(taxdump)
    return isDirectory ? taxdump.resolve(name) : (taxdump.name == name ? taxdump : taxdump.resolveSibling(name))
}


def taxdumpFile(value, name, purpose) {
    if (!value) {
        throw new IllegalArgumentException("Parameter 'taxdump' is required for ${purpose}")
    }
    def candidate = taxdumpPath(value, name)
    if (!candidate.exists()) {
        throw new IllegalArgumentException("taxdump has no ${name}: ${candidate.toUriString()}")
    }
    return candidate
}


def hostReferenceFor(meta, defaultReference) {
    def requested = (meta.host_ref ?: '').toString().trim()
    if (requested in ['', 'NA', 'N/A', 'default']) {
        if (!defaultReference) {
            throw new IllegalArgumentException(
                "Sample '${meta.id}' uses host_ref=default but parameter 'host_fasta' is not set"
            )
        }
        return defaultReference.toString()
    }
    return requested
}


def hostBundleFor(reference) {
    // Index files are looked up next to the FASTA, so a remote host reference works
    // as long as its .amb .ann .bwt .pac .sa files sit beside it.
    def fasta = file(reference, checkIfExists: true)
    if (isLocal(fasta) && !java.nio.file.Files.isRegularFile(fasta)) {
        throw new IllegalArgumentException("Host FASTA is not a regular file: ${fasta}")
    }
    def bundle = [fasta]
    ['.amb', '.ann', '.bwt', '.pac', '.sa'].each { suffix ->
        def index = fasta.resolveSibling(fasta.name + suffix)
        if (!index.exists()) {
            throw new IllegalArgumentException(
                "Host FASTA is missing BWA index: ${index.toUriString()}. Build it with `bwa index ${fasta.name}`."
            )
        }
        bundle << index
    }
    return bundle
}


def binScript(name) {
    return file("${projectDir}/bin/${name}", checkIfExists: true)
}
