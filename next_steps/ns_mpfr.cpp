// ns_mpfr.cpp: the online SOM in arbitrary precision (the MPFR library), the reference of the precision
// ladder of som_next_steps.Rmd. Euclidean or DTW distance, the Gaussian neighbourhood, the same steps in
// the same order as ns_core.cpp: distances to every unit, the first unit with the smallest distance,
// units pulled below 1e-10 skipped, m += (alpha h) (x - m), with the item warped onto the prototype
// along the DTW path when the distance is DTW. Every operation rounds to nearest at 'bits' bits.
// Compile with PKG_LIBS = "-lmpfr -lgmp".
//
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <mpfr.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;

struct MpVec {                                        // a vector of MPFR numbers, freed when it goes out of scope
  std::vector<__mpfr_struct> v;
  MpVec(size_t n, mpfr_prec_t bits) : v(n) { for (size_t i = 0; i < n; i++) mpfr_init2(&v[i], bits); }
  ~MpVec() { for (auto& x : v) mpfr_clear(&x); }
  MpVec(const MpVec&) = delete; MpVec& operator=(const MpVec&) = delete;
  mpfr_ptr operator[](size_t i) { return &v[i]; }
};
static const mpfr_rnd_t RN = MPFR_RNDN;

// the Euclidean distance between the item at x and the prototype at m (D values)
static void euclid_mp(MpVec& V, size_t x, size_t m, int D, mpfr_ptr out, mpfr_ptr t) {
  mpfr_set_zero(out, 1);
  for (int i = 0; i < D; i++) { mpfr_sub(t, V[x + i], V[m + i], RN); mpfr_sqr(t, t, RN); mpfr_add(out, out, t, RN); }
  mpfr_sqrt(out, out, RN);
}

// the DTW table of the item at x against the prototype at m, (L + 1) x (L + 1) in R; returns sqrt(cost) in out
static void dtw_mp(MpVec& V, size_t x, size_t m, int L, int nch, int band, MpVec& R, mpfr_ptr out, mpfr_ptr t, mpfr_ptr c) {
  const int W = L + 1;
  for (int i = 0; i < W * W; i++) mpfr_set_inf(R[i], 1);
  mpfr_set_zero(R[0], 1);
  for (int i = 1; i <= L; i++) {
    int lo = std::max(1, i - band), hi = std::min(L, i + band);
    for (int j = lo; j <= hi; j++) {
      mpfr_set_zero(c, 1);
      for (int ch = 0; ch < nch; ch++) { mpfr_sub(t, V[x + ch * L + i - 1], V[m + ch * L + j - 1], RN); mpfr_sqr(t, t, RN); mpfr_add(c, c, t, RN); }
      mpfr_ptr best = R[(i - 1) * W + (j - 1)];
      if (mpfr_cmp(R[(i - 1) * W + j], best) < 0) best = R[(i - 1) * W + j];
      if (mpfr_cmp(R[i * W + (j - 1)], best) < 0) best = R[i * W + (j - 1)];
      mpfr_add(R[i * W + j], c, best, RN);
    }
  }
  mpfr_sqrt(out, R[L * W + L], RN);
}

// walk the path back and average the item's points matched with every prototype point, into xw
static void dtw_warp_mp(MpVec& V, size_t x, int L, int nch, MpVec& R, MpVec& xw, std::vector<int>& cnt) {
  const int W = L + 1;
  std::fill(cnt.begin(), cnt.end(), 0);
  for (int k = 0; k < L * nch; k++) mpfr_set_zero(xw[k], 1);
  int i = L, j = L;
  while (i > 0 && j > 0) {
    for (int ch = 0; ch < nch; ch++) mpfr_add(xw[ch * L + j - 1], xw[ch * L + j - 1], V[x + ch * L + i - 1], RN);
    cnt[j - 1]++;
    if (i == 1 && j == 1) break;
    mpfr_ptr d = R[(i - 1) * W + (j - 1)], u = R[(i - 1) * W + j], l = R[i * W + (j - 1)];
    if (mpfr_cmp(d, u) <= 0 && mpfr_cmp(d, l) <= 0) { i--; j--; }
    else if (mpfr_cmp(u, l) <= 0) i--;
    else j--;
  }
  for (int j2 = 0; j2 < L; j2++) if (cnt[j2] > 0) for (int ch = 0; ch < nch; ch++) mpfr_div_ui(xw[ch * L + j2], xw[ch * L + j2], cnt[j2], RN);
}

// [[Rcpp::export]]
List som_train_mpfr_cpp(NumericMatrix X, NumericMatrix M0, NumericMatrix grid_dist, IntegerVector pick,
                        NumericVector alpha, NumericVector radius, int mode, int L, int nch, int band, int bits,
                        IntegerVector snap_at) {
  const int N = X.nrow(), D = X.ncol(), K = M0.nrow(), S = pick.size();
  if (D != L * nch || M0.ncol() != D) stop("X and M0 must have L * nch columns");
  if (mode != 0 && mode != 6) stop("mode must be 0 (Euclidean) or 6 (DTW)");
  if (band < 0) band = L;
  const mpfr_prec_t P = bits;
  MpVec V((size_t) (N + K) * D, P);                   // the items, then the prototypes
  for (int n = 0; n < N; n++) for (int d = 0; d < D; d++) mpfr_set_d(V[(size_t) n * D + d], X(n, d), RN);
  const size_t MOFF = (size_t) N * D;
  for (int k = 0; k < K; k++) for (int d = 0; d < D; d++) mpfr_set_d(V[MOFF + (size_t) k * D + d], M0(k, d), RN);
  MpVec G2((size_t) K * K, P), H((size_t) K * K, P), dist(K, P), xw(D, P), R((size_t) (L + 1) * (L + 1), P), tmp(8, P);
  mpfr_ptr t = tmp[0], c = tmp[1], a = tmp[2], r = tmp[3], step = tmp[4], delta = tmp[5], two_r2 = tmp[6];
  for (int k = 0; k < K; k++) for (int j = 0; j < K; j++) { mpfr_set_d(G2[(size_t) k * K + j], grid_dist(k, j), RN); mpfr_sqr(G2[(size_t) k * K + j], G2[(size_t) k * K + j], RN); }
  std::vector<int> win(S), cnt(L), snaps(snap_at.begin(), snap_at.end());
  List snapshots;
  double last_r = -1;
  for (int s = 0; s < S; s++) {
    const size_t x = (size_t) (pick[s] - 1) * D;
    mpfr_set_d(a, alpha[s], RN);
    if (radius[s] != last_r) {                        // the neighbourhood matrix, once per radius value
      mpfr_set_d(r, radius[s], RN); mpfr_sqr(two_r2, r, RN); mpfr_mul_ui(two_r2, two_r2, 2, RN);
      for (size_t i = 0; i < (size_t) K * K; i++) { mpfr_div(t, G2[i], two_r2, RN); mpfr_neg(t, t, RN); mpfr_exp(H[i], t, RN); }
      last_r = radius[s];
    }
    int near = 0;
    for (int k = 0; k < K; k++) {
      if (mode == 0) euclid_mp(V, x, MOFF + (size_t) k * D, D, dist[k], t);
      else dtw_mp(V, x, MOFF + (size_t) k * D, L, nch, band, R, dist[k], t, c);
      if (mpfr_cmp(dist[k], dist[near]) < 0) near = k;
    }
    win[s] = near + 1;
    for (int j = 0; j < K; j++) {
      mpfr_ptr h = H[(size_t) near * K + j];
      if (mpfr_cmp_d(h, 1e-10) <= 0) continue;
      const size_t m = MOFF + (size_t) j * D;
      mpfr_mul(step, a, h, RN);
      if (mode == 6) { dtw_mp(V, x, m, L, nch, band, R, t, t, c); dtw_warp_mp(V, x, L, nch, R, xw, cnt); }
      for (int d = 0; d < D; d++) {
        mpfr_sub(delta, mode == 6 ? xw[d] : V[x + d], V[m + d], RN);
        mpfr_mul(t, step, delta, RN);
        mpfr_add(V[m + d], V[m + d], t, RN);
      }
    }
    if (!snaps.empty() && std::find(snaps.begin(), snaps.end(), s + 1) != snaps.end()) {
      NumericMatrix Ms(K, D);
      for (int k = 0; k < K; k++) for (int d = 0; d < D; d++) Ms(k, d) = mpfr_get_d(V[MOFF + (size_t) k * D + d], RN);
      snapshots.push_back(Ms);
    }
  }
  NumericMatrix Mout(K, D);
  for (int k = 0; k < K; k++) for (int d = 0; d < D; d++) Mout(k, d) = mpfr_get_d(V[MOFF + (size_t) k * D + d], RN);
  return List::create(_["M"] = Mout, _["win"] = IntegerVector(win.begin(), win.end()), _["snapshots"] = snapshots, _["bits"] = bits);
}
