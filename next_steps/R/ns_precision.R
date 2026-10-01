# ns_precision.R: the precision ladder of som_next_steps.Rmd
#
# The same map, the same seed, trained in fp32, fp64, 80-bit extended and 128-bit quad precision by the
# C++ core of ns_core.cpp, and in MPFR arbitrary precision (default 256 bits) by ns_mpfr.cpp, which
# repeats the core's operations one by one: distances to every prototype, the first unit with the
# smallest distance, a Gaussian pull with units below 1e-10 skipped, and m += (alpha h) (x - m). The R
# loop train_mpfr() does the same with Rmpfr, as an independent check of the C++ MPFR code (it is about
# a hundred times slower). Then: at which step the winners first differ from the most precise run, how
# many differ in all, how far the prototypes drift, and whether any measure moves.

# the data as fp32: every value rounded to the nearest single-precision number, then used in fp64
round_fp32 <- function(x) {
  v <- readBin(writeBin(as.numeric(x), raw(), size = 4), "double", n = length(x), size = 4)
  if (is.matrix(x)) matrix(v, nrow(x), ncol(x)) else v
}

# the online Euclidean SOM in Rmpfr; X, M0: matrices; gd: map distances; returns the prototypes as
# doubles, the winners, and snapshots at the steps in snap_at
train_mpfr <- function(X, M0, gd, pick, alpha, radius, bits = 256, snap_at = integer(0)) {
  if (!requireNamespace("Rmpfr", quietly = TRUE)) stop("Rmpfr is needed for the R check of the MPFR reference")
  K <- nrow(M0); D <- ncol(M0); S <- length(pick)
  Xm <- Rmpfr::mpfr(X, bits); M <- Rmpfr::mpfr(M0, bits)
  G2 <- Rmpfr::mpfr(gd, bits); G2 <- G2 * G2
  two <- Rmpfr::mpfr(2, bits); thr <- Rmpfr::mpfr(1e-10, bits); ones <- Rmpfr::mpfr(rep(1, D), bits)
  win <- integer(S); snaps <- list()
  last_r <- NA; H <- NULL
  for (s in seq_len(S)) {
    a <- Rmpfr::mpfr(alpha[s], bits); r <- radius[s]
    if (!identical(r, last_r)) {                                  # the neighbourhood matrix, once per radius value
      rm <- Rmpfr::mpfr(r, bits)
      H <- exp(-G2 / (two * rm * rm))
      last_r <- r
    }
    diff <- Xm[rep(pick[s], K), ] - M                             # x - m, one row per unit
    d <- sqrt((diff * diff) %*% ones)                             # the distances, K x 1
    bmu <- which(as.vector(d == min(d)))[1]
    h <- H[bmu, ]
    h[h <= thr] <- 0
    step <- a * h                                                 # alpha h, one per unit
    M <- M + diff * step                                          # (alpha h) (x - m): step[k] recycles down row k
    win[s] <- bmu
    if (s %in% snap_at) snaps[[length(snaps) + 1]] <- matrix(as.numeric(M), K, D)
  }
  list(M = matrix(as.numeric(M), K, D), win = win, snapshots = snaps)
}

# run one setting through the ladder; returns the summary table and the drift curves
precision_ladder <- function(ds, size, seed, rlen, modes = c("euclid", "dtw"), mpfr_bits = 256, n_snaps = 40,
                             with_mpfr = TRUE, with_rmpfr = FALSE, fp32_data = TRUE) {
  grid <- make_grid(size, size)
  N <- nrow(ds$X); K <- grid$K; S <- rlen * N
  set.seed(seed); start <- sample(N, K); pick <- as.integer(floor(N * runif(S)) + 1)
  M0 <- ds$X[start, , drop = FALSE]
  X32 <- round_fp32(ds$X); M0_32 <- X32[start, , drop = FALSE]       # the fp32 data: items and start alike
  sch <- schedule(S, grid); sc <- scales_of(ds)
  snap_at <- unique(round(S * (seq_len(n_snaps) / n_snaps)^1.5)); snap_at <- snap_at[snap_at >= 1]
  rows <- list(); curves <- list(); gaps <- list()
  cpp <- function(X, prec, mode, kahan = FALSE, M_start = M0)
    som_train_cpp(X, M_start, grid$dist, pick, sch$alpha, sch$radius, rep(0, S), NS_MODES[[mode]], ds$L, ds$nch, sc$band, 2, 1, sc$maxshift,
                  0, 0, 0, prec, kahan, as.integer(snap_at))
  timed <- function(f) { t0 <- proc.time()[["elapsed"]]; r <- f(); r$seconds <- proc.time()[["elapsed"]] - t0; r }
  mp_name <- sprintf("MPFR %d-bit", mpfr_bits)
  for (mode in modes) {
    r <- list(fp32 = timed(function() cpp(ds$X, 0L, mode)), fp64 = timed(function() cpp(ds$X, 1L, mode)),
              `fp64 + Kahan` = timed(function() cpp(ds$X, 1L, mode, TRUE)), fp80 = timed(function() cpp(ds$X, 2L, mode)),
              fp128 = timed(function() cpp(ds$X, 3L, mode)))
    if (fp32_data) r[["fp32 data, fp64 arithmetic"]] <- timed(function() cpp(X32, 1L, mode, M_start = M0_32))
    if (with_mpfr)
      r[[mp_name]] <- timed(function() som_train_mpfr_cpp(ds$X, M0, grid$dist, pick, sch$alpha, sch$radius, NS_MODES[[mode]], ds$L, ds$nch, sc$band,
                                                          as.integer(mpfr_bits), as.integer(snap_at)))
    if (with_rmpfr && mode == "euclid")
      r[[paste(mp_name, "(Rmpfr check)")]] <- timed(function() train_mpfr(ds$X, M0, grid$dist, pick, sch$alpha, sch$radius, mpfr_bits, snap_at))
    ref_name <- if (with_mpfr) mp_name else "fp128"
    ref <- r[[ref_name]]
    for (nm in names(r)) {
      f <- r[[nm]]
      flips <- which(f$win != ref$win)
      Dm <- cross_dist(ds$X, f$M, ds, judge_args(ds, mode)); mp <- map_items(Dm)
      rows[[length(rows) + 1]] <- data.frame(mode = mode, precision = nm, reference = ref_name,
        first_flip = if (length(flips)) flips[1] else NA, flips = length(flips), flip_share = length(flips) / S,
        drift = max(abs(f$M - ref$M)), QE = mean(mp$dist), TE = mean(grid$dist[cbind(mp$bmu, mp$bmu2)] > 1.01),
        seconds = f$seconds, stringsAsFactors = FALSE)
      drift <- sapply(seq_along(snap_at), function(i) max(abs(f$snapshots[[i]] - ref$snapshots[[i]])))
      curves[[length(curves) + 1]] <- data.frame(mode = mode, precision = nm, step = snap_at, drift = drift, stringsAsFactors = FALSE)
    }
    gaps[[mode]] <- r$fp64$gap
  }
  list(table = do.call(rbind, rows), curves = do.call(rbind, curves), S = S, K = K, gaps = gaps,
       gap = gaps[[modes[1]]], gap_dtw = if ("dtw" %in% modes) gaps[["dtw"]] else NULL)
}

# the ladder of the report, computed once and cached (the same key as the report's); 'settings' names the
# size and the passes of each of the two runs
compute_ladder <- function(datasets, mpfr_bits = 256,
                           settings = list(cbf = list(size = 7, rlen = 20), pd = list(size = 7, rlen = 4), pd_raw = list(size = 7, rlen = 4)),
                           dir = "next_steps/cache") {
  small_cbf <- make_cbf(30, 40, 0.5, seed = 1)
  small_cbf$title <- "small CBF (90 series of 40 points, as the tutorial)"
  pd_raw <- datasets$pd                                             # the same pulses as raw pressure in mPa: four decades, no asinh
  pd_raw$X <- sinh(datasets$pd$X); pd_raw$key <- "pd_raw"; pd_raw$title <- "PD stand-in as raw pressure (mPa, no asinh)"
  key <- bench_key(mpfr_bits, settings, lapply(datasets$pd, function(x) x),
                   files = c(NS_CODE_FILES[c(1, 3, 4)], "next_steps/ns_mpfr.cpp", "next_steps/R/ns_precision.R"))
  ladder <- cached("precision", key, dir = dir, fn = function() list(
    cbf = precision_ladder(small_cbf, settings$cbf$size, 1, rlen = settings$cbf$rlen, modes = c("euclid", "dtw"), mpfr_bits = mpfr_bits, with_rmpfr = TRUE),
    pd  = precision_ladder(datasets$pd, settings$pd$size, 1, rlen = settings$pd$rlen, modes = c("euclid", "dtw"), mpfr_bits = mpfr_bits),
    pd_raw = precision_ladder(pd_raw, settings$pd_raw$size, 1, rlen = settings$pd_raw$rlen, modes = c("euclid", "dtw"), mpfr_bits = mpfr_bits)))
  ladder$small_cbf <- small_cbf
  ladder
}
