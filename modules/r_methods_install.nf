/*
 * FEAST, and the npreg package that FastST uses, are not on conda-forge. Without the
 * attribution image they are compiled once per run, from the sources at --feast_code and
 * --npreg_code, into a library that the FEAST and FastST tasks add to R_LIBS.
 */
process R_METHODS_INSTALL {
    label 'process_medium'
    container "${params.container_attribution}"

    input:
    path feast_source
    path npreg_source

    output:
    path 'rlib', emit: library

    script:
    """
    mkdir -p rlib feast
    tar -xzf ${feast_source} -C feast --strip-components=1
    # The FEAST repository ships compiled objects; remove them so R builds its own.
    rm -f feast/src/*.o feast/src/*.so
    export MAKEFLAGS="-j${task.cpus}"
    R CMD INSTALL --no-docs --library=rlib feast
    R CMD INSTALL --no-docs --library=rlib ${npreg_source}

    # What was built, from the archive's top directory (FEAST-<commit>) and the tarball name.
    printf '%s\\n' "\$(tar -tzf ${feast_source} | sed -n '1s#/.*##p' | sed 's/^FEAST-//')" > rlib/FEAST.commit
    """

    stub:
    """
    mkdir -p rlib
    printf 'stub\\n' > rlib/FEAST.commit
    """
}
