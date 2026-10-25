#!/usr/bin/env Rscript
# FEAST (Shenhav et al. 2019, Nature Methods): fits each sink as a mixture of the
# source environments and an unknown source by expectation-maximisation. Sources of
# one environment are pooled first, as SourceTracker2 does, so both methods estimate
# the same quantities. Writes the standard attribution table (bin/mstflow_attribution.py).

suppressPackageStartupMessages({
    library(optparse)
    library(FEAST)
})

opts <- parse_args(OptionParser(option_list = list(
    make_option("--counts"),
    make_option("--map"),
    make_option("--analysis"),
    make_option("--output"),
    make_option("--iterations", type = "integer", default = 1000L),
    make_option("--coverage", type = "integer", default = 0L,
                help = "Rarefaction depth; 0 uses the smallest sink or pooled source"),
    make_option("--seed", type = "integer", default = 1L)
)))

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "mstflow_attribution.R"))

inputs <- read_attribution_inputs(opts$counts, opts$map)
pooled <- pool_sources(inputs)
sinks <- inputs$counts[, inputs$sinks, drop = FALSE]

# FEAST takes samples in rows. Pooled sources get neutral row names so they can never
# collide with a sink's name; results are mapped back by position.
C <- t(cbind(pooled, sinks))
rownames(C) <- c(sprintf("source%03d", seq_len(ncol(pooled))), sprintf("sink%05d", seq_along(inputs$sinks)))
metadata <- data.frame(
    Env = c(colnames(pooled), rep("sink", length(inputs$sinks))),
    SourceSink = c(rep("Source", ncol(pooled)), rep("Sink", length(inputs$sinks))),
    id = c(rep(NA, ncol(pooled)), seq_along(inputs$sinks)),
    row.names = rownames(C)
)

# FEAST fixes the depth from the first sink when none is given and reuses it for the
# others, so the depth is chosen here, once, for all sinks.
coverage <- if (opts$coverage > 0) opts$coverage else min(rowSums(C))
if (min(rowSums(C)) < coverage) {
    stop(sprintf("FEAST coverage %d exceeds the smallest sink or pooled source (%d counts)",
                 coverage, min(rowSums(C))))
}

set.seed(opts$seed)
invisible(capture.output(FEAST(C = C, metadata = metadata, different_sources_flag = 0,
                               EM_iterations = opts$iterations, COVERAGE = coverage,
                               dir_path = getwd(), outfile = "feast")))
result <- as.matrix(read.table("feast_source_contributions_matrix.txt", header = TRUE, sep = "\t",
                               check.names = FALSE))
if (!identical(dim(result), c(length(inputs$sinks), ncol(pooled) + 1L))) {
    stop("FEAST returned a result of unexpected shape")
}
colnames(result) <- c(colnames(pooled), "Unknown")
rownames(result) <- inputs$sinks

write_attribution_results(opts$output, opts$analysis, "feast", "proportion", result)
