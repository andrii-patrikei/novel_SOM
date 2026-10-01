# ns_quality.R: the measures of som_next_steps.Rmd, the project's and a few new ones
#
# Every map is judged with three distances: the Euclidean distance, DTW with the data set's band, and the
# distance it was trained with ("own"). From the table of winners come QE, EV, TE, KL, TP, the share of
# units used, purity, ARI and NMI, as in som_metrics.Rmd and the benchmark; new here are the
# trustworthiness and continuity of the map (Venna & Kaski 2001) and, for detection, the AUROC of the
# distance to the winner as a novelty score.

# ------------------------------------------------------------------------------------- distances to judge
judge_args <- function(ds, judge, spec = NULL) {
  sc <- scales_of(ds)
  if (judge == "own") {
    list(mode = NS_MODES[[spec$mode]], p = spec$p, gamma = spec$gamma * sc$gamma_unit, band = sc$band, maxshift = sc$maxshift)
  } else list(mode = NS_MODES[[judge]], p = 2, gamma = 1, band = sc$band, maxshift = sc$maxshift)
}
cross_dist <- function(A, B, ds, ja)
  cross_dist_cpp(A, B, ja$mode, ds$L, ds$nch, ja$band, ja$p, ja$gamma, ja$maxshift)

# winner, its distance, runner-up (ties: the first unit), from an items x units table of distances
map_items <- function(Dm) {
  bmu <- max.col(-Dm, ties.method = "first")
  D2 <- Dm; D2[cbind(seq_len(nrow(Dm)), bmu)] <- Inf
  bmu2 <- max.col(-D2, ties.method = "first")
  data.frame(bmu = bmu, dist = Dm[cbind(seq_len(nrow(Dm)), bmu)], bmu2 = bmu2)
}

# -------------------------------------------------------------------------------------- class agreement
purity <- function(bmu, cls) sum(tapply(cls, bmu, function(v) max(table(v)))) / length(cls)
ari <- function(a, b) {
  pairs <- function(n) n * (n - 1) / 2
  tab <- table(a, b); both <- sum(pairs(tab)); in_a <- sum(pairs(rowSums(tab))); in_b <- sum(pairs(colSums(tab)))
  chance <- in_a * in_b / pairs(length(a))
  den <- (in_a + in_b) / 2 - chance
  if (den == 0) return(0)
  (both - chance) / den
}
nmi <- function(a, b) {
  tab <- table(a, b) / length(a)
  H <- function(p) { p <- p[p > 0]; -sum(p * log(p)) }
  ha <- H(rowSums(tab)); hb <- H(colSums(tab))
  if (ha == 0 || hb == 0) return(0)
  pij <- tab[tab > 0]; pi_ <- (rowSums(tab) %o% colSums(tab))[tab > 0]
  2 * sum(pij * log(pij / pi_)) / (ha + hb)
}
# the prototypes of the used units cut into as many groups as classes (Ward), every item in its winner's group
groups_of <- function(dV, bmu, n_groups) {
  used <- sort(unique(bmu))
  if (length(used) <= n_groups) { g <- rep(NA, nrow(dV)); g[used] <- seq_along(used); return(g[bmu]) }
  d <- dV[used, used]; d[d < 0] <- 0; d <- pmax(d, t(d))
  tree <- hclust(as.dist(d), method = "ward.D2")
  g <- rep(NA, nrow(dV)); g[used] <- cutree(tree, n_groups)
  g[bmu]
}

# ---------------------------------------------------------------------------- trustworthiness, continuity
# Venna & Kaski (2001), on the units: for every item, its k nearest items in the data against its k
# nearest items on the map (map distance of their winners, ties by data distance). Trustworthiness falls
# when far-off items land next to each other on the map; continuity falls when data neighbours are torn
# apart on the map. Both in [0, 1], 1 is best. Computed on a sample of items when there are many.
trust_cont <- function(Dx, bmu, grid_dist, k = 12, max_items = 400) {
  n <- nrow(Dx)
  idx <- if (n > max_items) sort(sample.int(n, max_items)) else seq_len(n)
  m <- length(idx); if (m <= k + 1) return(c(trust = NA, cont = NA))
  Dmap <- grid_dist[bmu[idx], bmu[idx]]
  Dd <- Dx[idx, idx]
  trust <- cont <- 0
  for (i in seq_len(m)) {
    rd <- order(Dd[i, -i]); rm <- order(Dmap[i, -i], Dd[i, -i])    # ranks without the item itself
    rank_data <- integer(m - 1); rank_data[rd] <- seq_len(m - 1)
    rank_map  <- integer(m - 1); rank_map[rm]  <- seq_len(m - 1)
    nd <- rd[1:k]; nm <- rm[1:k]
    trust <- trust + sum(pmax(rank_data[setdiff(nm, nd)] - k, 0))
    cont  <- cont  + sum(pmax(rank_map[setdiff(nd, nm)] - k, 0))
  }
  norm <- if (2 * k < m) 2 / (m * k * (2 * m - 3 * k - 1)) else 2 / (m * (m - k) * (m - k - 1))   # Venna & Kaski's two cases
  c(trust = 1 - norm * trust, cont = 1 - norm * cont)
}

# ----------------------------------------------------------------------------------------------- AUROC
# the probability that a random positive scores above a random negative (Mann-Whitney)
auroc <- function(score, positive) {
  r <- rank(score); n1 <- sum(positive); n0 <- sum(!positive)
  if (n1 == 0 || n0 == 0) return(NA)
  (sum(r[positive]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# ------------------------------------------------------------------------------------ all measures at once
# returns a one-row data frame per judge
measures <- function(ds, fit, spec, judges = c("euclid", "dtw", "own"), with_trust = TRUE) {
  M <- fit$M; grid <- fit$grid; K <- grid$K; n_cls <- nlevels(ds$cls)
  if (ds$nch > 1 && "euclid" %in% judges) judges <- judges            # Euclidean on the stacked channels is fine
  rows <- lapply(judges, function(judge) {
    ja <- judge_args(ds, judge, spec)
    Dm <- cross_dist(ds$X, M, ds, ja)
    mp <- map_items(Dm)
    dV <- cross_dist(M, M, ds, ja); diag(dV) <- 0
    soft_own <- judge == "own" && spec$mode == "softdtw"            # soft-DTW values can be negative: no EV, KL, TP
    d_mean <- cross_dist(ds$X, rbind(colMeans(ds$X)), ds, ja)[, 1]
    EV <- if (soft_own) NA else 1 - mean(mp$dist^2) / mean(d_mean^2)
    TE <- mean(grid$dist[cbind(mp$bmu, mp$bmu2)] > 1.01)
    KL <- if (soft_own) NA else { P <- map_ways_cpp(pmax(dV, 0), grid$links); mean(mp$dist + P[cbind(mp$bmu, mp$bmu2)]) }
    # (two prototypes that coincide give a zero ratio, which the C++ counts as no contribution)
    TP <- if (soft_own) NA else topographic_product_cpp(pmax(dV, 0), grid$dist)
    g  <- groups_of(dV, mp$bmu, n_cls)
    # trustworthiness and continuity against the Euclidean neighbourhoods of the items (one yardstick for
    # every map; an N x N DTW table per map would cost more than the map)
    tc <- if (with_trust && judge == "euclid") { set.seed(1); trust_cont(cross_dist(ds$X, ds$X, ds, ja), mp$bmu, grid$dist, max_items = 300) }
          else c(trust = NA, cont = NA)
    data.frame(judge = judge, QE = mean(mp$dist), EV = EV, TE = TE, KL = KL, TP = TP,
               used = length(unique(mp$bmu)) / K, dead = 1 - length(unique(mp$bmu)) / K,
               purity = purity(mp$bmu, ds$cls), ARI = ari(g, ds$cls), NMI = nmi(g, ds$cls),
               trust = unname(tc["trust"]), cont = unname(tc["cont"]), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}
