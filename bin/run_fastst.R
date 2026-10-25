#!/usr/bin/env Rscript
# FastST (Choi et al. 2025, PeerJ): estimates each sink as a mixture of the source
# environments and an unknown source by generalised least squares, with standard
# errors. Sources of one environment are pooled first, as for the other methods.
# Writes the standard attribution table (bin/mstflow_attribution.py).

suppressPackageStartupMessages({
    library(optparse)
    library(npreg)
})

opts <- parse_args(OptionParser(option_list = list(
    make_option("--counts"),
    make_option("--map"),
    make_option("--analysis"),
    make_option("--output")
)))

source(file.path(Sys.getenv("FASTST_HOME", "/opt/fastst"), "funST.R"))
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))),
                 "mstflow_attribution.R"))

inputs <- read_attribution_inputs(opts$counts, opts$map)
pooled <- pool_sources(inputs)
columns <- c("Unknown", colnames(pooled))   # glsest() reports the unknown source first
estimates <- matrix(NA_real_, length(inputs$sinks), length(columns), dimnames = list(inputs$sinks, columns))
sds <- estimates

for (sink in inputs$sinks) {
    fit <- suppressWarnings(glsest(inputs$counts[, sink], t(pooled)))
    raw <- fit[[1]]
    if (!all(is.finite(raw)) || sum(raw) <= 0) stop("FastST found no solution for sink ", sink)
    # As FastST's own command-line tool does: shares are rescaled to sum to one, and
    # the standard errors by the same factor.
    estimates[sink, ] <- raw / sum(raw)
    sds[sink, ] <- fit[[2]] / sum(raw)
}

write_attribution_results(opts$output, opts$analysis, "fastst", "proportion", estimates, sds)
