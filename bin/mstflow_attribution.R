# R counterpart of bin/mstflow_attribution.py: reads the attribution inputs and writes
# the standard result table, with the same checks, for the methods written in R.

read_attribution_inputs <- function(counts_path, map_path) {
    counts <- as.matrix(read.delim(counts_path, check.names = FALSE, comment.char = "",
                                   row.names = 1, colClasses = c("character")))
    storage.mode(counts) <- "numeric"
    if (any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts))) {
        stop(counts_path, " must hold non-negative whole counts")
    }
    units <- read.delim(map_path, check.names = FALSE, comment.char = "", colClasses = "character")
    if (!all(c("#SampleID", "SourceSink", "Env") %in% names(units))) {
        stop(map_path, " must have #SampleID, SourceSink and Env columns")
    }
    if (!setequal(units[["#SampleID"]], colnames(counts))) {
        stop("the attribution count table and map list different units")
    }
    role <- tolower(units$SourceSink)
    sources <- units[["#SampleID"]][role == "source"]
    sinks <- units[["#SampleID"]][role == "sink"]
    if (length(sources) == 0 || length(sinks) == 0) stop("attribution needs at least one source and one sink")
    env <- setNames(units$Env, units[["#SampleID"]])
    list(counts = counts, sources = sources, sinks = sinks, env = env,
         envs = sort(unique(env[sources])))
}

pool_sources <- function(inputs) {
    pooled <- sapply(inputs$envs, function(e) {
        rowSums(inputs$counts[, inputs$sources[inputs$env[inputs$sources] == e], drop = FALSE])
    })
    matrix(pooled, nrow = nrow(inputs$counts), dimnames = list(rownames(inputs$counts), inputs$envs))
}

# `estimates` (and `sds`, if given) are sinks x sources matrices with an Unknown column.
write_attribution_results <- function(path, analysis, method, quantity, estimates, sds = NULL,
                                      tolerance = 1e-3) {
    if (!"Unknown" %in% colnames(estimates)) stop(method, " reported no Unknown share")
    sources <- c(sort(setdiff(colnames(estimates), "Unknown")), "Unknown")
    rows <- list()
    for (sink in rownames(estimates)) {
        values <- pmax(estimates[sink, sources], 0)
        if (any(!is.finite(values))) stop(method, ": non-finite estimate for ", sink)
        total <- sum(values)
        if (abs(total - 1) > tolerance) {
            stop(sprintf("%s: estimates for %s sum to %.6g, not 1", method, sink, total))
        }
        sd <- if (is.null(sds)) rep(NA_real_, length(sources)) else sds[sink, sources] / total
        rows[[sink]] <- data.frame(analysis = analysis, method = method, sink = sink, source = sources,
                                   estimate = signif(values / total, 6),
                                   sd = ifelse(is.finite(sd), signif(sd, 6), NA),
                                   quantity = quantity)
    }
    write.table(do.call(rbind, rows), path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
