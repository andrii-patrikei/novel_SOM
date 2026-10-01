#!/usr/bin/env Rscript
# Precomputes the cache of som_next_steps.Rmd outside the knit, so that the knit itself takes a minute:
#   Rscript next_steps/run_benchmark.R [seeds = 5] [passes = 10] [cores = all but two]
# Run it from the project root. The knit uses the same key (code + settings) and finds the results.
args  <- commandArgs(trailingOnly = TRUE)
seeds <- seq_len(if (length(args) >= 1) as.integer(args[1]) else 5L)
rlen  <- if (length(args) >= 2) as.integer(args[2]) else 10L
cores <- if (length(args) >= 3) as.integer(args[3]) else max(1L, parallel::detectCores() - 2L)
if (!file.exists("next_steps/ns_core.cpp")) stop("run this from the project root (the folder that holds next_steps/)")
suppressMessages(library(Rcpp))
Sys.setenv(PKG_LIBS = "-lquadmath")
sourceCpp("next_steps/ns_core.cpp", cacheDir = "next_steps/.cpp_cache")
for (f in c("ns_data.R", "ns_maps.R", "ns_quality.R", "ns_bench.R")) source(file.path("next_steps/R", f))
datasets <- load_datasets()
cat(sprintf("%d seeds, %d passes, %d cores; data: %s\n", length(seeds), rlen, cores,
            paste(sapply(datasets, function(d) sprintf("%s (%s)", d$key, d$source)), collapse = ", ")))
t0 <- Sys.time()
res <- compute_all(datasets, VARIANTS, c(5, 7, 10, 14), seeds, rlen, cores)
cat(sprintf("%d result rows, %d seed-1 maps, %d detection rows, in %s\n", nrow(res$bench), length(res$fits), nrow(res$detection),
            format(round(Sys.time() - t0, 1))))
