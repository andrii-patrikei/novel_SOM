# Checks of the compiled core; run from the project root: Rscript next_steps/tests/test_core.R
if (!file.exists("next_steps/ns_core.cpp")) stop("run this from the project root")
suppressMessages(library(Rcpp))
Sys.setenv(PKG_LIBS = "-lquadmath")
t0 <- Sys.time()
sourceCpp("next_steps/ns_core.cpp", cacheDir = "next_steps/.cpp_cache")
cat("compiled in", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")

# --- 1. DTW against a plain R implementation (squared local cost, band) ---------------------------------
dtw_r <- function(x, y, band) {
  L <- length(x); D <- matrix(Inf, L + 1, L + 1); D[1, 1] <- 0
  for (i in 1:L) for (j in max(1, i - band):min(L, i + band))
    D[i + 1, j + 1] <- (x[i] - y[j])^2 + min(D[i, j], D[i, j + 1], D[i + 1, j])
  sqrt(D[L + 1, L + 1])
}
set.seed(1)
ok <- replicate(200, {
  L <- sample(10:40, 1); x <- rnorm(L); y <- rnorm(L); b <- sample(0:L, 1)
  a <- dtw_r(x, y, b); c <- cross_dist_cpp(rbind(x), rbind(y), 6L, L, 1L, b, 2, 1, -1L)[1, 1]
  abs(a - c) < 1e-12
})
cat("DTW == R reference in", sum(ok), "of 200\n")

# --- 2. soft-DTW -> DTW^2 as gamma -> 0; and E-warp sums to a convex combination -------------------------
x <- rnorm(30); y <- rnorm(30)
d_hard <- cross_dist_cpp(rbind(x), rbind(y), 6L, 30L, 1L, 5L, 2, 1, -1L)[1, 1]^2
for (g in c(1, 0.1, 0.01, 0.001)) {
  d_soft <- cross_dist_cpp(rbind(x), rbind(y), 7L, 30L, 1L, 5L, 2, g, -1L)[1, 1]
  cat(sprintf("gamma %6.3f: softDTW %.5f  DTW^2 %.5f\n", g, d_soft, d_hard))
}
xw_hard <- item_seen_cpp(x, y, 6L, 30L, 1L, 5L, 2, 1, -1L)
xw_soft <- item_seen_cpp(x, y, 7L, 30L, 1L, 5L, 2, 0.001, -1L)
cat("soft warp ~ hard warp at gamma 0.001: max diff", max(abs(xw_hard - xw_soft)), "\n")
cat("warped item within range of x:", all(xw_soft >= min(x) - 1e-9 & xw_soft <= max(x) + 1e-9), "\n")

# --- 3. SBD: a shifted copy has distance ~0 and is shifted back -----------------------------------------
s <- c(rep(0, 10), dnorm(1:20, 10, 3), rep(0, 10)); s2 <- c(rep(0, 4), s)[1:40]
d_sbd <- cross_dist_cpp(rbind(s2), rbind(s), 8L, 40L, 1L, 0L, 2, 1, 10L)[1, 1]
back <- item_seen_cpp(s2, s, 8L, 40L, 1L, 0L, 2, 1, 10L)
cat(sprintf("SBD of a 4-step shifted copy: %.4f; shifted back max diff %.4f\n", d_sbd, max(abs(back - s))))

# --- 4. multichannel dependent DTW equals single-channel DTW on stacked channels when L matches ----------
x3 <- c(x, x * 2, x * 0.5); y3 <- c(y, y * 2, y * 0.5)
d3 <- cross_dist_cpp(rbind(x3), rbind(y3), 6L, 30L, 3L, 5L, 2, 1, -1L)[1, 1]
cat("3-channel DTW (should be sqrt(1+4+0.25) x single):", d3 / cross_dist_cpp(rbind(x), rbind(y), 6L, 30L, 1L, 5L, 2, 1, -1L)[1, 1], "vs", sqrt(5.25), "\n")

# --- 5. pointwise distances against R ------------------------------------------------------------------
X <- matrix(rnorm(20 * 8), 20); M <- matrix(rnorm(5 * 8), 5)
chk <- function(mode, ref) max(abs(cross_dist_cpp(X, M, mode, 8L, 1L, 0L, 3, 1, 0L) - ref))
cat("euclid", chk(0L, as.matrix(dist(rbind(X, M)))[1:20, 21:25]),
    " manhattan", chk(1L, as.matrix(dist(rbind(X, M), "manhattan"))[1:20, 21:25]),
    " chebyshev", chk(2L, as.matrix(dist(rbind(X, M), "maximum"))[1:20, 21:25]),
    " minkowski3", chk(3L, as.matrix(dist(rbind(X, M), "minkowski", p = 3))[1:20, 21:25]),
    " cosine", chk(4L, 1 - (X %*% t(M)) / outer(sqrt(rowSums(X^2)), sqrt(rowSums(M^2)))),
    " correlation", chk(5L, 1 - cor(t(X), t(M))), "\n")

# --- 6. training: four precisions on a small CBF-like problem -------------------------------------------
make_cbf <- function(shape, L, noise) {
  t <- 1:L; a <- runif(1, 1/8, 1/4) * L; b <- a + runif(1, 1/4, 3/4) * L; on <- t >= a & t <= b; amp <- 6 + rnorm(1)
  y <- switch(shape, cylinder = amp * on, bell = amp * on * (t - a) / (b - a), funnel = amp * on * (b - t) / (b - a))
  y + rnorm(L, sd = noise)
}
set.seed(1); Xc <- t(sapply(rep(c("cylinder", "bell", "funnel"), each = 30), make_cbf, L = 40, noise = 0.5))
grid <- as.matrix(expand.grid(x = 1:5, y = 1:5)); grid[, 1] <- grid[, 1] + 0.5 * (grid[, 2] %% 2); grid[, 2] <- grid[, 2] * sqrt(3) / 2
gd <- as.matrix(dist(grid)); K <- 25; S <- 180
set.seed(1); M0 <- Xc[sample(90, K), ]; pick <- as.integer(floor(90 * runif(S)) + 1)
decay <- function(a, b, s, n) a - (a - b) * (s - 1) / n
alpha <- decay(0.5, 0.01, 1:S, S); radius <- decay(3, 0, 1:S, S); radius[radius < 1] <- 0.5
run <- function(prec, mode = 0L, band = 8L, ...) som_train_cpp(Xc, M0, gd, pick, alpha, radius, rep(0, S), mode, 40L, 1L, band, 2, 1, -1L,
                                                                0, 0, 0, prec, FALSE, integer(0))
for (mode in c(0L, 6L)) {
  res <- lapply(0:3, run, mode = mode)
  cat(sprintf("mode %d: max |M_prec - M_quad| float %.2e double %.2e longdouble %.2e; winner flips vs quad: %d %d %d\n", mode,
              max(abs(res[[1]]$M - res[[4]]$M)), max(abs(res[[2]]$M - res[[4]]$M)), max(abs(res[[3]]$M - res[[4]]$M)),
              sum(res[[1]]$win != res[[4]]$win), sum(res[[2]]$win != res[[4]]$win), sum(res[[3]]$win != res[[4]]$win)))
}
# the plain SOM in R, the same operations: must equal the double run to the last bit (or nearly)
M <- M0
for (s in 1:S) {
  x <- Xc[pick[s], ]; d <- sqrt(colSums((t(M) - x)^2)); bmu <- which.min(d)
  h <- exp(-gd[bmu, ]^2 / (2 * radius[s]^2)); h[h <= 1e-10] <- 0
  M <- M + h * alpha[s] * (matrix(x, K, 40, byrow = TRUE) - M)
}
cat("R loop vs C++ double: max diff", max(abs(M - run(1L)$M)), "\n")

# --- 7. options run and change something ----------------------------------------------------------------
base <- run(1L)$M
opt <- function(...) som_train_cpp(Xc, M0, gd, pick, alpha, radius, ..., 0L, 40L, 1L, 8L, 2, 1, -1L, 0, 0, 0, 1L, FALSE, integer(0))
temp <- rep(0, S); temp2 <- seq(5, 0.01, length.out = S)
cat("annealing changes M:", max(abs(som_train_cpp(Xc, M0, gd, pick, alpha, radius, temp2, 0L, 40L, 1L, 8L, 2, 1, -1L, 0, 0, 0, 1L, FALSE, integer(0))$M - base)) > 0,
    " conscience:", max(abs(som_train_cpp(Xc, M0, gd, pick, alpha, radius, temp, 0L, 40L, 1L, 8L, 2, 1, -1L, 0.5, 0, 0, 1L, FALSE, integer(0))$M - base)) > 0,
    " huber:", max(abs(som_train_cpp(Xc, M0, gd, pick, alpha, radius, temp, 0L, 40L, 1L, 8L, 2, 1, -1L, 0, 1, 0, 1L, FALSE, integer(0))$M - base)) > 0,
    " momentum:", max(abs(som_train_cpp(Xc, M0, gd, pick, alpha, radius, temp, 0L, 40L, 1L, 8L, 2, 1, -1L, 0, 0, 0.5, 1L, FALSE, integer(0))$M - base)) > 0,
    " kahan == plain:", max(abs(som_train_cpp(Xc, M0, gd, pick, alpha, radius, temp, 0L, 40L, 1L, 8L, 2, 1, -1L, 0, 0, 0, 1L, TRUE, integer(0))$M - base)), "\n")
t1 <- Sys.time(); r7 <- som_train_cpp(Xc, M0, gd, pick, alpha, radius, temp, 7L, 40L, 1L, 8L, 2, 0.5, -1L, 0, 0, 0, 1L, FALSE, integer(0))
cat("softDTW training 180 steps 25 units:", round(as.numeric(Sys.time() - t1, units = "secs"), 2), "s; finite:", all(is.finite(r7$M)), "\n")
t1 <- Sys.time(); r8 <- som_train_cpp(Xc, M0, gd, pick, alpha, radius, temp, 8L, 40L, 1L, 8L, 2, 0.5, 8L, 0, 0, 0, 1L, FALSE, integer(0))
cat("SBD training:", round(as.numeric(Sys.time() - t1, units = "secs"), 2), "s; finite:", all(is.finite(r8$M)), "\n")

# --- 8. topographic product and map ways -----------------------------------------------------------------
A <- as.matrix(dist(cbind(1:10, 0))); Ml <- cbind(1:10, 0) * 2; set.seed(1); Ms <- Ml[sample(10), ]
cat(sprintf("TP ordered chain %.4f, scrambled %.4f\n", topographic_product_cpp(as.matrix(dist(Ml)), A), topographic_product_cpp(as.matrix(dist(Ms)), A)))
links <- (A > 0 & A < 1.01) * 1L
P <- map_ways_cpp(as.matrix(dist(Ml)), links); cat("way 1->3 ordered chain (should be 4):", P[1, 3], "\n")
