# ns_data.R: the data of som_next_steps.Rmd
#
# Three synthetic stand-ins, made to look like the project's hard cases, and a loader that swaps in the
# project's real sets (som_bench/data/datasets.rds) when that file is present.
#
# A data set is a list: X (one item per row; a series of L points in nch channels, channel after channel),
# cls (factor), L, nch, key, title, about, source ("synthetic" or "project"), band_share (the DTW window
# as a share of L), and x_label.

# ---------------------------------------------------------------------------------------- CBF (Saito 1994)
make_cbf <- function(n_per_class = 100, L = 64, noise = 1, seed = 1) {
  set.seed(seed)
  one <- function(shape) {
    t <- 1:L; a <- runif(1, 1/8, 1/4) * L; b <- a + runif(1, 1/4, 3/4) * L; on <- t >= a & t <= b
    amp <- 6 + rnorm(1)
    y <- switch(shape, cylinder = amp * on, bell = amp * on * (t - a) / (b - a), funnel = amp * on * (b - t) / (b - a))
    y + rnorm(L, sd = noise)
  }
  cls <- factor(rep(c("cylinder", "bell", "funnel"), each = n_per_class), levels = c("cylinder", "bell", "funnel"))
  X <- t(sapply(as.character(cls), one)); rownames(X) <- NULL
  list(X = X, cls = cls, L = L, nch = 1L, key = "cbf", title = "CBF", source = "synthetic", band_share = 0.1,
       x_label = "time", about = sprintf(
         "Cylinder-Bell-Funnel (Saito, 1994): %d series of %d points, noise %g. One event per series that starts and ends at random times: made for DTW.",
         nrow(X), L, noise))
}

# ------------------------------------------------------------------------- partial-discharge-like pulses
# 125 samples at 500 kS/s (2 us per sample, 250 us in all), as the project's acoustic PD data, triggered
# near the onset. Four kinds of waveform, each a damped oscillation with its own ringing frequency, decay
# and loudness, or band-limited noise. Loudness spans four decades, so the pressure p (in mPa) is used as
# asinh(p), as the project does. The four amplitude ranges overlap on purpose: amplitude alone cannot
# tell the kinds apart; the ringing frequency, the decay, a reflection and the rise time can.
make_pd <- function(n_per_class = 150, L = 125, seed = 1, jitter_us = 4) {
  set.seed(seed)
  t <- (0:(L - 1)) * 2                                            # microseconds
  pulse <- function(A, f_khz, tau, t0, rise = 0, reflect = NULL) {
    s <- pmax(t - t0, 0); on <- t >= t0
    env <- exp(-s / tau) * (if (rise > 0) (1 - exp(-s / rise)) else 1)
    y <- A * env * sin(2 * pi * f_khz * s / 1000) * on
    if (!is.null(reflect)) {
      s2 <- pmax(t - t0 - reflect$delay, 0); on2 <- t >= t0 + reflect$delay
      y <- y + reflect$share * A * exp(-s2 / tau) * sin(2 * pi * f_khz * s2 / 1000) * on2
    }
    y
  }
  one <- function(kind) {
    t0 <- 40 + rnorm(1, 0, jitter_us)
    p <- switch(kind,
      corona  = pulse(10^runif(1, -0.3, 0.9), runif(1, 55, 95), runif(1, 8, 16), t0),
      air_gap = pulse(10^runif(1, 1.6, 3.4), runif(1, 22, 42), runif(1, 30, 60), t0,
                      reflect = list(delay = runif(1, 50, 90), share = runif(1, 0.25, 0.55))),
      void    = pulse(10^runif(1, 0.4, 2.1), runif(1, 38, 62), runif(1, 20, 40), t0, rise = runif(1, 4, 10)),
      noise   = { sd <- 10^runif(1, -0.8, 0.6); e <- rnorm(L); for (i in 2:L) e[i] <- 0.7 * e[i - 1] + e[i]; sd * e / sd(e) })
    asinh(p + rnorm(L, sd = 0.2))                                 # 0.2 mPa of background noise, then the asinh scale
  }
  levels <- c("corona", "air_gap", "void", "noise")
  cls <- factor(rep(levels, each = n_per_class), levels = levels)
  X <- t(sapply(as.character(cls), one)); rownames(X) <- NULL
  list(X = X, cls = cls, L = L, nch = 1L, key = "pd", title = "Partial discharges (synthetic)", source = "synthetic",
       band_share = 0.1, x_label = "time (2 us per sample)", about = sprintf(
         "Synthetic stand-in for the project's airborne-ultrasound partial discharges: %d waveforms of %d samples at 500 kS/s, four kinds (corona: fast ringing, quick decay, quiet; air gap: slow ringing, long decay, loud, with a reflection; void: a slower rise; noise: band-limited noise), loudness over four decades, pressure as asinh(p / 1 mPa), onsets jittered by %g us around the trigger.",
         nrow(X), L, jitter_us))
}

# --------------------------------------------------------------------------- gesture-like 3-channel series
# Four hand gestures seen by a 3-axis sensor: a template per class in normalised time, then for every
# instance a random speed profile (so the same shape comes at different times, as in the UCR gesture
# sets), a random amplitude, a small random rotation of the three axes (the sensor is never held the same
# way twice), a small offset per axis (gravity) and noise.
make_gestures <- function(n_per_class = 100, L = 64, seed = 1, rotation_deg = 35, noise = 0.25, speed_var = 0.6) {
  set.seed(seed)
  bump <- function(u, c, w) exp(-((u - c) / w)^2 / 2)
  templates <- list(                                              # flick and twist share their first axis: the hard pair
    circle = function(u) cbind(sin(2 * pi * u), cos(2 * pi * u), 0.3 * sin(4 * pi * u)),
    shake  = function(u) { b <- bump(u, 0.5, 0.22); cbind(sin(8 * pi * u) * b, 0.4 * cos(8 * pi * u) * b, 0.2 * sin(8 * pi * u) * b) },
    flick  = function(u) cbind(bump(u, 0.35, 0.06) - 0.8 * bump(u, 0.5, 0.08), 0.6 * bump(u, 0.42, 0.07), -0.3 * bump(u, 0.4, 0.1)),
    twist  = function(u) cbind(bump(u, 0.35, 0.06) - 0.8 * bump(u, 0.5, 0.08), 0.5 * tanh(6 * (u - 0.5)), 0.4 * bump(u, 0.55, 0.12)))
  rotation <- function(deg) {                                     # a random rotation by 'deg' degrees about a random axis
    a <- rnorm(3); a <- a / sqrt(sum(a^2)); th <- deg * pi / 180
    Kx <- matrix(c(0, a[3], -a[2], -a[3], 0, a[1], a[2], -a[1], 0), 3)
    diag(3) + sin(th) * Kx + (1 - cos(th)) * Kx %*% Kx
  }
  one <- function(kind) {
    u0 <- seq(0, 1, length.out = L)
    speed <- 1 + speed_var * (sin(2 * pi * u0 + runif(1, 0, 2 * pi)) * runif(1) + 0.5 * sin(4 * pi * u0 + runif(1, 0, 2 * pi)) * runif(1))
    u <- cumsum(pmax(speed, 0.2)); u <- (u - u[1]) / (u[L] - u[1])   # a monotone warp of time
    Y <- templates[[kind]](u) * exp(rnorm(1, 0, 0.35))
    Y <- Y %*% rotation(rnorm(1, 0, rotation_deg))
    Y <- Y + matrix(rnorm(3, 0, 0.25), L, 3, byrow = TRUE) + matrix(rnorm(3 * L, 0, noise), L, 3)
    as.vector(Y)                                                  # channel after channel
  }
  levels <- names(templates)
  cls <- factor(rep(levels, each = n_per_class), levels = levels)
  X <- t(sapply(as.character(cls), one)); rownames(X) <- NULL
  list(X = X, cls = cls, L = L, nch = 3L, key = "gestures", title = "Gestures (synthetic, 3 axes)", source = "synthetic",
       band_share = 0.1, x_label = "time", about = sprintf(
         "Synthetic stand-in for movement data: %d gestures of %d samples on three axes, four kinds (circle, shake, flick, twist; flick and twist share one axis), each instance with its own speed profile, amplitude, a rotation of the axes of about %g degrees, an offset per axis (gravity) and noise of %g. Stored as one row per gesture, channel after channel; DTW, soft-DTW and SBD align all three channels together.",
         nrow(X), L, rotation_deg, noise))
}

# -------------------------------------------------------------------------------------------- the loader
# The project's real sets replace the synthetic ones when som_bench/data/datasets.rds is found (CBF, the
# acoustic partial discharges, and GunPoint as the single-channel movement set); the synthetic gestures
# stay, because the project has no multi-channel set.
load_datasets <- function(project_rds = "som_bench/data/datasets.rds", n_per_class = c(cbf = 100, pd = 150, gestures = 100)) {
  ds <- list(cbf = make_cbf(n_per_class[["cbf"]]), pd = make_pd(n_per_class[["pd"]]), gestures = make_gestures(n_per_class[["gestures"]]))
  if (file.exists(project_rds)) {
    real <- readRDS(project_rds)
    take <- function(key, new_key, title) {
      d <- real[[key]]
      if (is.null(d)) return(NULL)
      list(X = d$X, cls = d$cls, L = ncol(d$X), nch = 1L, key = new_key, title = title, source = "project",
           band_share = if (!is.null(d$window) && !is.na(d$window)) d$window / ncol(d$X) else 0.1,
           x_label = if (is.null(d$x_label)) "time" else d$x_label, about = paste("The project's own data set:", d$title))
    }
    for (k in list(c("cbf", "cbf", "CBF (project)"), c("pd_acoustic", "pd", "Partial discharges (project)"),
                   c("gunpoint", "gunpoint", "GunPoint (project)"))) {
      r <- take(k[1], k[2], k[3])
      if (!is.null(r)) ds[[k[2]]] <- r
    }
  }
  ds
}
