/*
 * Unpack a classifier database supplied as a .tar.gz archive, local or remote, such
 * as the prebuilt Kraken2 indexes at https://benlangmead.github.io/aws-indexes/k2.
 * Files may sit at the archive root or inside one top-level directory. `kind` is the
 * classifier, kraken2 or krakenuniq, followed by +bracken when Bracken's read-length
 * models must be there too.
 */
process UNTAR_DATABASE {
    label 'process_low'
    tag "${archive.name}"
    container "${params.container_python_tools}"

    input:
    path archive
    val kind

    output:
    path "database", emit: db

    script:
    def krakenUniq = kind.startsWith('krakenuniq')
    def marker = krakenUniq ? 'database.kdb' : 'hash.k2d'
    def required = krakenUniq ? 'database.kdb database.idx taxDB' : 'hash.k2d opts.k2d taxo.k2d'
    """
    mkdir -p unpacked
    tar -xzf "${archive}" -C unpacked

    root=unpacked
    if [[ ! -f unpacked/${marker} ]]; then
        mapfile -t found < <(find unpacked -mindepth 2 -maxdepth 2 -name ${marker})
        if [[ \${#found[@]} -ne 1 ]]; then
            echo "${archive.name} does not contain a ${krakenUniq ? 'KrakenUniq' : 'Kraken2'} database (${required})" >&2
            exit 1
        fi
        root=\$(dirname "\${found[0]}")
    fi
    for name in ${required}; do
        [[ -s "\$root/\$name" ]] || { echo "${archive.name} is missing \$name" >&2; exit 1; }
    done
    if [[ "${kind}" == *+bracken ]] && ! compgen -G "\$root/database*mers.kmer_distrib" > /dev/null; then
        echo "${archive.name} has no Bracken database<N>mers.kmer_distrib files; build them with bracken-build, or set --run_bracken false" >&2
        exit 1
    fi
    mv "\$root" database
    """

    stub:
    """
    mkdir -p database
    touch database/hash.k2d database/opts.k2d database/taxo.k2d database/database100mers.kmer_distrib
    touch database/database.kdb database/database.idx database/taxDB
    """
}
