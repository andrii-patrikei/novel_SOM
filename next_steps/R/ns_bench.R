# ns_bench.R: the jobs, the parallel runs and the cache of som_next_steps.Rmd

NS_CODE_FILES <- c("next_steps/ns_core.cpp", "next_steps/R/ns_data.R", "next_steps/R/ns_maps.R", "next_steps/R/ns_quality.R")
NS_SLOW <- c("softdtw", "dtw", "dtw_anneal")

# a key from the code and the settings: change either and the results are computed again
bench_key <- function(..., files = NS_CODE_FILES) {
  code <- lapply(files, function(f) if (file.exists(f)) xfun::read_utf8(f) else NA)
  substr(digest::digest(list(code, list(...)), algo = "sha1"), 1, 12)
}

cached <- function(prefix, key, fn, dir = "next_steps/cache") {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(dir, sprintf("%s_%s.rds", prefix, key))
  if (file.exists(path)) return(readRDS(path))
  value <- fn()
  saveRDS(value, path)
  value
}

# ----------------------------------------------------------------------------------------------- the jobs
# one row per map: data set, variant, map size (units per side), seed
make_jobs <- function(datasets, variants, sizes, seeds, skip = NULL) {
  jobs <- expand.grid(dataset = names(datasets), variant = names(variants), size = sizes, seed = seeds, stringsAsFactors = FALSE)
  if (!is.null(skip)) jobs <- jobs[!skip(jobs), ]
  jobs <- jobs[order(-jobs$size, -(jobs$variant %in% NS_SLOW)), ]       # the slow ones first, so the cores stay busy
  rownames(jobs) <- NULL
  jobs
}
# soft-DTW aligns every unit twice per step and costs minutes per map on the largest size: it stops at 10 x 10
skip_softdtw_large <- function(jobs) jobs$variant == "softdtw" & jobs$size > 10

# train one map and judge it; returns one row per judge, and for seed 1 the map itself (for the pictures)
run_job <- function(job, datasets, variants, rlen, with_trust = TRUE, keep_fit = job$seed == 1) {
  ds <- datasets[[job$dataset]]; spec <- variants[[job$variant]]
  fit <- train_map(ds, spec, job$size, job$seed, rlen = rlen)
  m <- measures(ds, fit, spec, with_trust = with_trust)
  rows <- cbind(data.frame(dataset = job$dataset, variant = job$variant, size = job$size, units = fit$grid$K, seed = job$seed,
                           seconds = fit$seconds, comparisons = fit$comparisons,
                           not_nearest = mean(fit$win != fit$nearest),      # the soft and conscience winners
                           gap_median = median(fit$gap, na.rm = TRUE), gap_tiny = mean(fit$gap < 1e-6, na.rm = TRUE),
                           stringsAsFactors = FALSE), m, row.names = NULL)
  kept <- NULL
  if (keep_fit) {
    bmu <- map_items(cross_dist(ds$X, fit$M, ds, judge_args(ds, "own", spec)))$bmu
    kept <- list(dataset = job$dataset, variant = job$variant, size = job$size, M = fit$M, bmu = bmu, counts = fit$counts)
  }
  list(rows = rows, fit = kept)
}

run_benchmark <- function(jobs, datasets, variants, rlen, cores, with_trust = TRUE) {
  res <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
    tryCatch(run_job(jobs[i, ], datasets, variants, rlen, with_trust),
             error = function(e) { message("job ", i, " failed: ", conditionMessage(e)); NULL })
  }, mc.cores = cores, mc.preschedule = FALSE)
  res <- Filter(Negate(is.null), res)
  fits <- Filter(Negate(is.null), lapply(res, `[[`, "fit"))
  names(fits) <- sapply(fits, function(f) paste(f$dataset, f$variant, f$size, sep = "|"))
  list(rows = do.call(rbind, lapply(res, `[[`, "rows")), fits = fits)
}

# everything the report needs, computed once and cached: the benchmark (with the seed-1 maps) and the
# detection runs. The key holds the code of the core, the data, the maps and the measures, and the settings.
compute_all <- function(datasets, variants, sizes, seeds, rlen, cores, det_size = 7, dir = "next_steps/cache") {
  data_id <- lapply(datasets, function(d) list(d$X, d$cls, d$L, d$nch, d$band_share))
  jobs <- make_jobs(datasets, variants, sizes, seeds, skip = skip_softdtw_large)
  bench <- cached("bench", bench_key(names(variants), sizes, seeds, rlen, data_id),
                  function() run_benchmark(jobs, datasets, variants, rlen, cores), dir)
  det <- cached("detect", bench_key(names(variants), det_size, seeds, rlen, data_id),
                function() run_detection(datasets, variants, det_size, seeds, rlen, cores), dir)
  list(bench = bench$rows, fits = bench$fits, detection = det, jobs = jobs)
}

# ------------------------------------------------------------------------------------- detection as novelty
# train on the background class only, score every item by its distance to the winner, and ask how well
# that score tells the other classes from the background (AUROC); one row per map
run_detection <- function(datasets, variants, sizes, seeds, rlen, cores, background = c(pd = "noise")) {
  jobs <- expand.grid(dataset = names(background), variant = names(variants), size = sizes, seed = seeds, stringsAsFactors = FALSE)
  rows <- parallel::mclapply(seq_len(nrow(jobs)), function(i) {
    job <- jobs[i, ]; ds <- datasets[[job$dataset]]; spec <- variants[[job$variant]]
    keep <- ds$cls == background[[job$dataset]]
    bg <- ds; bg$X <- ds$X[keep, , drop = FALSE]; bg$cls <- droplevels(ds$cls[keep])
    fit <- train_map(bg, spec, job$size, job$seed, rlen = rlen)
    ja <- judge_args(bg, "own", spec)                                 # the scales the map was trained with
    score <- map_items(cross_dist(ds$X, fit$M, ds, ja))$dist
    pos <- !keep
    per_class <- sapply(setdiff(levels(ds$cls), background[[job$dataset]]), function(cl)
      auroc(score[keep | ds$cls == cl], (ds$cls == cl)[keep | ds$cls == cl]))
    cbind(data.frame(dataset = job$dataset, variant = job$variant, size = job$size, seed = job$seed,
                     AUROC = auroc(score, pos), stringsAsFactors = FALSE), as.list(setNames(per_class, paste0("AUROC_", names(per_class)))))
  }, mc.cores = cores, mc.preschedule = FALSE)
  do.call(rbind, rows)
}
