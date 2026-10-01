# ns_maps.R: lattices, schedules, the variants and one training function
#
# The conventions are the project's: a hexagonal lattice, a Gaussian neighbourhood, learning rate 0.5 to
# 0.01 and radius from kohonen's default start (the two-thirds quantile of the distances between units)
# to 0, used as 0.5 below 1; the start is a draw of items and the order of items a draw per step, both
# from the seed, so that every variant with the same seed sees the same start and the same order.

NS_MODES <- c(euclid = 0L, manhattan = 1L, chebyshev = 2L, minkowski = 3L, cosine = 4L, correlation = 5L,
              dtw = 6L, softdtw = 7L, sbd = 8L)
NS_PREC  <- c(fp32 = 0L, fp64 = 1L, fp80 = 2L, fp128 = 3L)

# --------------------------------------------------------------------------------------------- lattices
make_grid <- function(xdim, ydim = xdim, topo = c("hex", "rect"), torus = FALSE) {
  topo <- match.arg(topo)
  g <- as.matrix(expand.grid(x = 1:xdim, y = 1:ydim))
  if (topo == "hex") {
    if (torus) stop("the torus is built on the rectangular lattice")
    g[, "x"] <- g[, "x"] + 0.5 * (g[, "y"] %% 2)
    g[, "y"] <- g[, "y"] * sqrt(3) / 2
    d <- as.matrix(dist(g))
  } else if (!torus) {
    d <- as.matrix(dist(g))
  } else {                                                        # wrap both ways: the shortest of the images
    dx <- abs(outer(g[, "x"], g[, "x"], "-")); dx <- pmin(dx, xdim - dx)
    dy <- abs(outer(g[, "y"], g[, "y"], "-")); dy <- pmin(dy, ydim - dy)
    d <- sqrt(dx^2 + dy^2)
  }
  dimnames(d) <- NULL
  list(pos = g, dist = d, links = (d > 0 & d < 1.01) * 1L, xdim = xdim, ydim = ydim, K = nrow(g), topo = topo, torus = torus)
}

radius_start <- function(grid) unname(quantile(grid$dist[upper.tri(grid$dist)], 2/3))   # kohonen's default

decay <- function(start, end, step, n_steps) start - (start - end) * (step - 1) / n_steps

schedule <- function(S, grid, alpha = c(0.5, 0.01), radius0 = radius_start(grid)) {
  r <- decay(radius0, 0, 1:S, S); r[r < 1] <- 0.5
  list(alpha = decay(alpha[1], alpha[2], 1:S, S), radius = r)
}

# --------------------------------------------------------------------------------------------- variants
# Each variant is a list of the arguments of som_train_cpp that differ from the online SOM. 'family'
# names the plain control that the variant must beat seed by seed.
variant <- function(id, label, mode = "euclid", family = "euclid", topo = "hex", torus = FALSE, p = 2, gamma = 0.1,
                    temp = 0, conscience = 0, huber = 0, momentum = 0, prec = "fp64", kahan = FALSE, what = "", group = "distance")
  list(id = id, label = label, mode = mode, family = family, topo = topo, torus = torus, p = p, gamma = gamma, temp = temp,
       conscience = conscience, huber = huber, momentum = momentum, prec = prec, kahan = kahan, what = what, group = group)

VARIANTS <- list(
  variant("euclid",      "Online SOM (Euclidean)", what = "The project's online SOM: the plain control of the Euclidean family.", group = "control"),
  variant("manhattan",   "Manhattan (L1)",          mode = "manhattan",   what = "The sum of absolute differences chooses the winner; the update is unchanged (a straight line is an L1 geodesic too). Less swayed by one loud sample."),
  variant("chebyshev",   "Chebyshev (L-inf)",       mode = "chebyshev",   what = "Only the largest difference counts: the winner is the prototype whose worst sample is least wrong."),
  variant("minkowski3",  "Minkowski p = 3",         mode = "minkowski", p = 3, what = "Between Euclidean and Chebyshev: large differences weigh more than in L2."),
  variant("cosine",      "Cosine",                  mode = "cosine",      what = "The angle between item and prototype: amplitude drops out, shape stays. For PD pulses this removes the loudness."),
  variant("correlation", "Correlation (1 - r)",     mode = "correlation", what = "Cosine after removing the mean: level and amplitude drop out."),
  variant("sbd",         "Shape-based (SBD)",       mode = "sbd",         what = "k-Shape's distance: 1 - the best normalised cross-correlation over shifts of up to the DTW window; the item is shifted by the best lag before the update. The rigid cousin of DTW."),
  variant("dtw",         "DTW-SOM",                 mode = "dtw", family = "dtw", what = "The project's DTW-SOM: DTW with a 10% band chooses the winner and warps the item onto each prototype before the update. The control of the elastic family.", group = "control"),
  variant("softdtw",     "Soft-DTW SOM",            mode = "softdtw", family = "dtw", gamma = 0.1, what = "Soft-DTW (Cuturi & Blondel 2017): a free energy over all alignments at temperature gamma instead of the single best path; the update averages the item with the expected alignment (the soft barycentre step)."),
  variant("huber",       "Euclidean, clipped pull", huber = 1, group = "ingredient", what = "Every component of the pull x - m is clipped at one standard deviation of the data: a robust, L1-like step that stops one loud item from dragging a prototype."),
  variant("l1_huber",    "Manhattan, clipped pull", mode = "manhattan", huber = 1, group = "ingredient", what = "The L1 winner with the clipped pull: the pairing that an L1 objective suggests."),
  variant("conscience",  "Conscience (frequency-sensitive)", conscience = 0.25, group = "ingredient", what = "DeSieno's conscience: a unit that wins too often is handicapped, so no unit stays dead and no unit takes all the loud items."),
  variant("anneal",      "Deterministic annealing", temp = 1, group = "ingredient", what = "A soft winner: the item pulls every unit with a Boltzmann weight at a temperature that falls to zero during training (the free-energy SOM of Graepel, Burger & Obermayer). Chaos was a clock without a reason; a temperature is the reason."),
  variant("momentum",    "Momentum (heavy ball)",   momentum = 0.5, group = "ingredient", what = "Each prototype keeps a velocity: half of the last step is added to the next. The physicist's inertia; the optimiser's heavy ball."),
  variant("dtw_anneal",  "DTW-SOM, annealed",       mode = "dtw", family = "dtw", temp = 1, group = "ingredient", what = "The soft winner on the DTW-SOM."),
  variant("rect",        "Rectangular sheet",       topo = "rect", family = "euclid", group = "lattice", what = "The online SOM on a square lattice (4 side neighbours at distance 1): the control for the torus."),
  variant("torus",       "Torus",                   topo = "rect", torus = TRUE, family = "rect", group = "lattice", what = "The square lattice wrapped both ways: no rim, no corner units, every unit has four neighbours. For cyclic data (the phase of the mains cycle in PD) and against rim effects.")
)
names(VARIANTS) <- sapply(VARIANTS, `[[`, "id")

# ------------------------------------------------------------------------------------ data-dependent scales
# the DTW band in samples, the SBD shift, soft-DTW's gamma, the Huber clip and the annealing temperature
# all come from the data set, so that one setting of a variant means the same thing on every set
scales_of <- function(ds) {
  v <- mean(apply(ds$X, 2, var))                                  # a typical squared difference of one sample
  list(band = max(1L, as.integer(round(ds$band_share * ds$L))), maxshift = max(1L, as.integer(round(ds$band_share * ds$L))),
       gamma_unit = v, huber_unit = sqrt(v))
}

# squared scale of the data under one distance: the mean squared distance from the items to their mean,
# the temperature unit of the annealed variants (as the chaotic-units page used s for its scale)
scale2_under <- function(ds, mode, sc, p = 2, gamma = 1) {
  Dm <- cross_dist_cpp(ds$X, rbind(colMeans(ds$X)), NS_MODES[[mode]], ds$L, ds$nch, sc$band, p, gamma, sc$maxshift)
  mean(Dm^2)
}

# ------------------------------------------------------------------------------------------ one training
train_map <- function(ds, spec, size, seed, rlen = 10, snap_at = integer(0)) {
  grid <- make_grid(size, size, spec$topo, spec$torus)
  N <- nrow(ds$X); K <- grid$K; S <- rlen * N
  set.seed(seed)
  M0 <- ds$X[sample(N, K), , drop = FALSE]
  pick <- as.integer(floor(N * runif(S)) + 1)
  sch <- schedule(S, grid)
  sc <- scales_of(ds)
  gamma <- spec$gamma * sc$gamma_unit
  temp <- rep(0, S)
  if (spec$temp > 0) {                                            # geometric cooling over the first 90%, then the hard winner
    s2 <- scale2_under(ds, spec$mode, sc, spec$p, gamma)
    cool <- floor(0.9 * S)
    temp[1:cool] <- spec$temp * s2 * exp(log(1e-3) * (0:(cool - 1)) / cool)
  }
  t0 <- proc.time()[["elapsed"]]
  r <- som_train_cpp(ds$X, M0, grid$dist, pick, sch$alpha, sch$radius, temp, NS_MODES[[spec$mode]], ds$L, ds$nch,
                     sc$band, spec$p, gamma, sc$maxshift, spec$conscience, spec$huber * sc$huber_unit, spec$momentum,
                     NS_PREC[[spec$prec]], spec$kahan, as.integer(snap_at))
  r$seconds <- proc.time()[["elapsed"]] - t0
  r$grid <- grid; r$M0 <- M0; r$pick <- pick; r$sched <- sch; r$temp <- temp; r$scales <- sc; r$gamma <- gamma
  r
}
