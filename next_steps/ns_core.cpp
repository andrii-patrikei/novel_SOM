// ns_core.cpp: the compiled core of som_next_steps.Rmd
//
// One online SOM, written once as a template over the scalar type, so that the same loop can run in
// fp32 (float), fp64 (double), 80-bit extended (long double) and 128-bit quad (__float128); the
// 256-bit reference runs in R with Rmpfr (next_steps/R/ns_precision.R) and repeats these operations.
//
// Every row of X is one item: a series of L points in nch channels, stored channel by channel
// (value at time i of channel c is x[c * L + i]); a feature vector is a series with nch = 1.
//
// Distances (the 'mode' argument), each with the step that moves a prototype towards an item:
//   0 EUCLID      sqrt(sum (x - m)^2)                 m += a h (x - m)
//   1 MANHATTAN   sum |x - m|                         same step (a straight line is a geodesic of L1 too)
//   2 CHEBYSHEV   max |x - m|                         same step
//   3 MINKOWSKI   (sum |x - m|^p)^(1/p)               same step
//   4 COSINE      1 - <x, m> / (|x| |m|)              same step
//   5 CORRELATION 1 - Pearson(x, m)                   same step
//   6 DTW         sqrt(DTW cost, Sakoe-Chiba band)    x is warped onto m's time axis along the path first
//   7 SOFTDTW     soft-DTW value (Cuturi & Blondel)   x is warped by the expected alignment E (soft barycentre step)
//   8 SBD         1 - max normalised cross-correlation  x is shifted by the best lag first (k-Shape)
// Multichannel series use the dependent form of DTW, soft-DTW and SBD: one alignment for all channels.
//
// Options of the training loop (all off = the online SOM of the project, Gaussian neighbourhood):
//   temp[s] > 0   soft winner (deterministic annealing, Graepel, Burger & Obermayer): P(k | x) is a
//                 softmax of -(sum_j h_kj d_j^2) / temp over the units, and unit j moves with weight
//                 sum_k P(k) h_kj instead of h_{win, j}
//   conscience    frequency-sensitive winner (Ahalt et al. 1990 / DeSieno 1988): the winner minimises
//                 d_k (1 + beta (K p_k - 1)), with p_k the running share of wins of unit k
//   huber > 0     every component of the pull x - m is clipped to [-huber, huber] (a robust, L1-like step)
//   momentum > 0  heavy-ball update: v = mu v + a h (x - m); m += v
//   kahan         Neumaier-compensated sums in the pointwise distances (fp64 and below)
//
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <quadmath.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
using namespace Rcpp;

enum { EUCLID = 0, MANHATTAN = 1, CHEBYSHEV = 2, MINKOWSKI = 3, COSINE = 4, CORRELATION = 5,
       DTW = 6, SOFTDTW = 7, SBD = 8 };
enum { P_FLOAT = 0, P_DOUBLE = 1, P_LONGDOUBLE = 2, P_QUAD = 3 };

// ---------------------------------------------------------------------------------------- scalar helpers
template <typename T> inline T t_sqrt(T x) { return std::sqrt(x); }
template <typename T> inline T t_exp(T x)  { return std::exp(x); }
template <typename T> inline T t_log(T x)  { return std::log(x); }
template <typename T> inline T t_abs(T x)  { return std::fabs(x); }
template <typename T> inline T t_pow(T x, T y) { return std::pow(x, y); }
template <> inline __float128 t_sqrt(__float128 x) { return sqrtq(x); }
template <> inline __float128 t_exp(__float128 x)  { return expq(x); }
template <> inline __float128 t_log(__float128 x)  { return logq(x); }
template <> inline __float128 t_abs(__float128 x)  { return fabsq(x); }
template <> inline __float128 t_pow(__float128 x, __float128 y) { return powq(x, y); }

template <typename T> inline T big() { return T(1e300); }   // 'infinity' inside the DP tables

struct Opts {
  int mode, L, nch, band, maxshift;
  double p, gamma;
  bool kahan;
};

// ------------------------------------------------------------------------------------- pointwise distances
template <typename T>
struct Sum {                                      // plain or Neumaier-compensated running sum
  T s = 0, c = 0; bool kahan;
  explicit Sum(bool k) : kahan(k) {}
  inline void add(T v) {
    if (!kahan) { s += v; return; }
    T t = s + v;
    if (t_abs(s) >= t_abs(v)) c += (s - t) + v; else c += (v - t) + s;
    s = t;
  }
  inline T value() const { return kahan ? s + c : s; }
};

template <typename T>
T pointwise(const T* x, const T* m, int D, const Opts& o) {
  switch (o.mode) {
    case EUCLID: { Sum<T> s(o.kahan); for (int i = 0; i < D; i++) { T d = x[i] - m[i]; s.add(d * d); } return t_sqrt(s.value()); }
    case MANHATTAN: { Sum<T> s(o.kahan); for (int i = 0; i < D; i++) s.add(t_abs(x[i] - m[i])); return s.value(); }
    case CHEBYSHEV: { T mx = 0; for (int i = 0; i < D; i++) { T d = t_abs(x[i] - m[i]); if (d > mx) mx = d; } return mx; }
    case MINKOWSKI: { Sum<T> s(o.kahan); T p = T(o.p); for (int i = 0; i < D; i++) s.add(t_pow(t_abs(x[i] - m[i]), p));
                      return t_pow(s.value(), T(1) / p); }
    case COSINE: { Sum<T> xy(o.kahan), xx(o.kahan), mm(o.kahan);
                   for (int i = 0; i < D; i++) { xy.add(x[i] * m[i]); xx.add(x[i] * x[i]); mm.add(m[i] * m[i]); }
                   T den = t_sqrt(xx.value()) * t_sqrt(mm.value());
                   return den > 0 ? T(1) - xy.value() / den : T(1); }
    case CORRELATION: { T mx = 0, mm_ = 0; for (int i = 0; i < D; i++) { mx += x[i]; mm_ += m[i]; } mx /= D; mm_ /= D;
                        Sum<T> xy(o.kahan), xx(o.kahan), mm(o.kahan);
                        for (int i = 0; i < D; i++) { T a = x[i] - mx, b = m[i] - mm_; xy.add(a * b); xx.add(a * a); mm.add(b * b); }
                        T den = t_sqrt(xx.value()) * t_sqrt(mm.value());
                        return den > 0 ? T(1) - xy.value() / den : T(1); }
  }
  return 0;
}

// ---------------------------------------------------------------------------------------------- DTW
// local cost between time i of x and time j of m: the squared distance over the channels
template <typename T>
inline T local_cost(const T* x, const T* m, int i, int j, int L, int nch) {
  T c = 0;
  for (int ch = 0; ch < nch; ch++) { T d = x[ch * L + i] - m[ch * L + j]; c += d * d; }
  return c;
}

// the DP table of DTW: (L + 1) x (L + 1), row 0 and column 0 are the border; returns sqrt(cost)
template <typename T>
T dtw_table(const T* x, const T* m, const Opts& o, std::vector<T>& R) {
  const int L = o.L, W = L + 1, b = o.band;
  R.assign((size_t) W * W, big<T>());
  R[0] = 0;
  for (int i = 1; i <= L; i++) {
    int lo = std::max(1, i - b), hi = std::min(L, i + b);
    for (int j = lo; j <= hi; j++) {
      T best = R[(size_t) (i - 1) * W + (j - 1)];
      T up = R[(size_t) (i - 1) * W + j], left = R[(size_t) i * W + (j - 1)];
      if (up < best) best = up;
      if (left < best) best = left;
      R[(size_t) i * W + j] = local_cost(x, m, i - 1, j - 1, L, o.nch) + best;
    }
  }
  return t_sqrt(R[(size_t) L * W + L]);
}

// walk the path back and average, for every prototype point j, the item points matched with it
template <typename T>
void dtw_warp(const T* x, const Opts& o, const std::vector<T>& R, T* xw) {
  const int L = o.L, W = L + 1;
  std::vector<int> cnt(L, 0);
  for (int k = 0; k < L * o.nch; k++) xw[k] = 0;
  int i = L, j = L;
  while (i > 0 && j > 0) {
    for (int ch = 0; ch < o.nch; ch++) xw[ch * L + (j - 1)] += x[ch * L + (i - 1)];
    cnt[j - 1]++;
    if (i == 1 && j == 1) break;
    T d = R[(size_t) (i - 1) * W + (j - 1)], u = R[(size_t) (i - 1) * W + j], l = R[(size_t) i * W + (j - 1)];
    if (d <= u && d <= l) { i--; j--; }            // ties: the diagonal first, then the item's time, then the prototype's
    else if (u <= l) i--;
    else j--;
  }
  for (int j2 = 0; j2 < L; j2++) if (cnt[j2] > 0) for (int ch = 0; ch < o.nch; ch++) xw[ch * L + j2] /= T(cnt[j2]);
}

// ------------------------------------------------------------------------------------------ soft-DTW
template <typename T>
inline T softmin3(T a, T b, T c, T g) {
  T m = a; if (b < m) m = b; if (c < m) m = c;
  if (m >= big<T>() / 2) return big<T>();
  T s = t_exp(-(a - m) / g) + t_exp(-(b - m) / g) + t_exp(-(c - m) / g);
  return m - g * t_log(s);
}

// forward pass; R is (L + 2) x (L + 2) so that the backward pass has its padding row and column
template <typename T>
T softdtw_table(const T* x, const T* m, const Opts& o, std::vector<T>& R) {
  const int L = o.L, W = L + 2, b = o.band;
  const T g = T(o.gamma);
  R.assign((size_t) W * W, big<T>());
  R[0] = 0;
  for (int i = 1; i <= L; i++) {
    int lo = std::max(1, i - b), hi = std::min(L, i + b);
    for (int j = lo; j <= hi; j++)
      R[(size_t) i * W + j] = local_cost(x, m, i - 1, j - 1, L, o.nch)
        + softmin3(R[(size_t) (i - 1) * W + (j - 1)], R[(size_t) (i - 1) * W + j], R[(size_t) i * W + (j - 1)], g);
  }
  return R[(size_t) L * W + L];
}

// backward pass: E[i, j] is the expected number of times (i, j) lies on the alignment under the Gibbs
// distribution over paths at temperature gamma; xw[j] = sum_i E[i, j] x[i] / sum_i E[i, j]
template <typename T>
void softdtw_warp(const T* x, const T* m, const Opts& o, std::vector<T>& R, T* xw) {
  const int L = o.L, W = L + 2, b = o.band;
  const T g = T(o.gamma), BIG = big<T>();
  std::vector<T> E((size_t) W * W, T(0)), delta((size_t) W * W, T(0));
  for (int i = 1; i <= L; i++) for (int j = 1; j <= L; j++)
    if (std::abs(i - j) <= b) delta[(size_t) i * W + j] = local_cost(x, m, i - 1, j - 1, L, o.nch);
  // padding as in Cuturi & Blondel (2017), Algorithm 2
  for (int k = 0; k <= L + 1; k++) { R[(size_t) k * W + (L + 1)] = -BIG; R[(size_t) (L + 1) * W + k] = -BIG; }
  R[(size_t) (L + 1) * W + (L + 1)] = R[(size_t) L * W + L];
  E[(size_t) (L + 1) * W + (L + 1)] = 1;
  for (int j = L; j >= 1; j--) {
    for (int i = L; i >= 1; i--) {
      T r = R[(size_t) i * W + j];
      if (r >= BIG / 2 || std::abs(i - j) > b) continue;          // outside the band: no path passes here
      T a = 0, bb = 0, c = 0;
      T r1 = R[(size_t) (i + 1) * W + j], r2 = R[(size_t) i * W + (j + 1)], r3 = R[(size_t) (i + 1) * W + (j + 1)];
      if (r1 < BIG / 2 && r1 > -BIG / 2) a  = t_exp((r1 - r - delta[(size_t) (i + 1) * W + j]) / g);
      if (r2 < BIG / 2 && r2 > -BIG / 2) bb = t_exp((r2 - r - delta[(size_t) i * W + (j + 1)]) / g);
      if (r3 < BIG / 2 && r3 > -BIG / 2) c  = t_exp((r3 - r - delta[(size_t) (i + 1) * W + (j + 1)]) / g);
      E[(size_t) i * W + j] = a * E[(size_t) (i + 1) * W + j] + bb * E[(size_t) i * W + (j + 1)] + c * E[(size_t) (i + 1) * W + (j + 1)];
    }
  }
  for (int k = 0; k < L * o.nch; k++) xw[k] = 0;
  for (int j = 1; j <= L; j++) {
    T tot = 0;
    for (int i = 1; i <= L; i++) {
      T e = E[(size_t) i * W + j];
      if (e <= 0) continue;
      tot += e;
      for (int ch = 0; ch < o.nch; ch++) xw[ch * L + (j - 1)] += e * x[ch * L + (i - 1)];
    }
    if (tot > 0) for (int ch = 0; ch < o.nch; ch++) xw[ch * L + (j - 1)] /= tot;
    else for (int ch = 0; ch < o.nch; ch++) xw[ch * L + (j - 1)] = m[ch * L + (j - 1)];   // cannot happen inside the band
  }
}

// ----------------------------------------------------------------------------------------------- SBD
// cross-correlation over lags -maxshift..maxshift, summed over the channels; returns 1 - max NCC and the lag
template <typename T>
T sbd_dist(const T* x, const T* m, const Opts& o, int* best_lag) {
  const int L = o.L;
  T xx = 0, mm = 0;
  for (int k = 0; k < L * o.nch; k++) { xx += x[k] * x[k]; mm += m[k] * m[k]; }
  T den = t_sqrt(xx) * t_sqrt(mm);
  if (den <= 0) { *best_lag = 0; return T(1); }
  T best = -2; int lag = 0;
  for (int w = -o.maxshift; w <= o.maxshift; w++) {
    T cc = 0;
    for (int ch = 0; ch < o.nch; ch++)
      for (int i = 0; i < L; i++) { int ix = i - w; if (ix >= 0 && ix < L) cc += x[ch * L + ix] * m[ch * L + i]; }
    cc /= den;
    if (cc > best) { best = cc; lag = w; }
  }
  *best_lag = lag;
  return T(1) - best;
}

template <typename T>
void sbd_shift(const T* x, const Opts& o, int lag, T* xw) {   // xw[i] = x[i - lag], zero outside
  const int L = o.L;
  for (int ch = 0; ch < o.nch; ch++)
    for (int i = 0; i < L; i++) { int ix = i - lag; xw[ch * L + i] = (ix >= 0 && ix < L) ? x[ch * L + ix] : T(0); }
}

// ---------------------------------------------------------------------------------- one distance, any mode
template <typename T>
T distance_of(const T* x, const T* m, const Opts& o, std::vector<T>& scratch, int* lag) {
  switch (o.mode) {
    case DTW:     return dtw_table(x, m, o, scratch);
    case SOFTDTW: return softdtw_table(x, m, o, scratch);
    case SBD:     return sbd_dist(x, m, o, lag);
    default:      return pointwise(x, m, o.L * o.nch, o);
  }
}

// the item as the update sees it from prototype m: warped, shifted, or itself
template <typename T>
void item_for(const T* x, const T* m, const Opts& o, std::vector<T>& scratch, int lag, T* xw) {
  const int D = o.L * o.nch;
  switch (o.mode) {
    case DTW:     dtw_table(x, m, o, scratch); dtw_warp(x, o, scratch, xw); break;
    case SOFTDTW: softdtw_table(x, m, o, scratch); softdtw_warp(x, m, o, scratch, xw); break;
    case SBD:     sbd_shift(x, o, lag, xw); break;
    default:      for (int k = 0; k < D; k++) xw[k] = x[k];
  }
}

static Opts make_opts(int mode, int L, int nch, int band, double p, double gamma, int maxshift, bool kahan) {
  Opts o; o.mode = mode; o.L = L; o.nch = nch; o.band = band < 0 ? L : band; o.p = p; o.gamma = gamma;
  o.maxshift = maxshift < 0 ? L - 1 : maxshift; o.kahan = kahan;
  return o;
}

// -------------------------------------------------------------------------------------- the training loop
template <typename T>
List train_T(const NumericMatrix& X, const NumericMatrix& M0, const NumericMatrix& grid_dist,
             const IntegerVector& pick, const NumericVector& alpha, const NumericVector& radius,
             const NumericVector& temp, const Opts& o, double conscience, double huber, double momentum,
             const IntegerVector& snap_at) {
  const int N = X.nrow(), D = X.ncol(), K = M0.nrow(), S = pick.size();
  std::vector<T> Xt((size_t) N * D), M((size_t) K * D), V, G((size_t) K * K);
  for (int n = 0; n < N; n++) for (int d = 0; d < D; d++) Xt[(size_t) n * D + d] = T(X(n, d));
  for (int k = 0; k < K; k++) for (int d = 0; d < D; d++) M[(size_t) k * D + d] = T(M0(k, d));
  for (int k = 0; k < K; k++) for (int j = 0; j < K; j++) G[(size_t) k * K + j] = T(grid_dist(k, j));
  if (momentum > 0) V.assign((size_t) K * D, T(0));
  std::vector<T> d(K), h((size_t) K * K), w(K), P(K), e(K), xw(D), scratch, pfreq(K, T(1) / T(K));
  std::vector<int> lag(K, 0), win(S), nearest(S), counts(K, 0);
  std::vector<int> snaps(snap_at.begin(), snap_at.end());
  std::vector<double> gap(S);                                    // (d_2nd - d_1st) / d_1st: how close the contest was
  List snapshots;
  long comparisons = 0;
  T last_r = -1;
  const T mu = T(momentum), hub = T(huber), beta = T(conscience);
  for (int s = 0; s < S; s++) {
    const T* x = &Xt[(size_t) (pick[s] - 1) * D];
    const T a = T(alpha[s]), r = T(radius[s]), tmp = T(temp[s]);
    if (r != last_r) {                                           // the neighbourhood matrix, once per radius value
      for (int k = 0; k < K; k++) for (int j = 0; j < K; j++) {
        T g = G[(size_t) k * K + j];
        h[(size_t) k * K + j] = t_exp(-g * g / (T(2) * r * r));
      }
      last_r = r;
    }
    // 1. distances to every prototype
    int near = 0;
    for (int k = 0; k < K; k++) {
      d[k] = distance_of(x, &M[(size_t) k * D], o, scratch, &lag[k]);
      if (d[k] < d[near]) near = k;
    }
    comparisons += K;
    nearest[s] = near + 1;
    {
      T second = big<T>();
      for (int k = 0; k < K; k++) if (k != near && d[k] < second) second = d[k];
      T dn = d[near];
      gap[s] = (K > 1 && second < big<T>() / 2) ? (double) ((second - dn) / (dn > 0 ? dn : T(1))) : NA_REAL;
    }
    // 2. the winner, or the soft assignment
    int winner = near;
    if (beta > 0) {
      T best = d[0] * (T(1) + beta * (T(K) * pfreq[0] - T(1)));
      winner = 0;
      for (int k = 1; k < K; k++) {
        T v = d[k] * (T(1) + beta * (T(K) * pfreq[k] - T(1)));
        if (v < best) { best = v; winner = k; }
      }
      for (int k = 0; k < K; k++) pfreq[k] += T(0.0001) * ((k == winner ? T(1) : T(0)) - pfreq[k]);
    }
    if (tmp > 0) {                                               // deterministic annealing: a softmax over the units
      T emin = big<T>();
      for (int k = 0; k < K; k++) {
        T ek = 0;
        for (int j = 0; j < K; j++) ek += h[(size_t) k * K + j] * d[j] * d[j];
        e[k] = ek;
        if (ek < emin) emin = ek;
      }
      T Z = 0;
      for (int k = 0; k < K; k++) { P[k] = t_exp(-(e[k] - emin) / tmp); Z += P[k]; }
      int argmax = 0;
      for (int k = 0; k < K; k++) { P[k] /= Z; if (P[k] > P[argmax]) argmax = k; }
      winner = argmax;                                           // reported winner: the most probable unit
      for (int j = 0; j < K; j++) { T wj = 0; for (int k = 0; k < K; k++) wj += P[k] * h[(size_t) k * K + j]; w[j] = wj; }
    } else {
      for (int j = 0; j < K; j++) w[j] = h[(size_t) winner * K + j];
    }
    win[s] = winner + 1;
    counts[winner]++;
    // 3. the update of every unit the item pulls
    for (int j = 0; j < K; j++) {
      if (w[j] <= T(1e-10)) continue;
      T* m = &M[(size_t) j * D];
      item_for(x, m, o, scratch, lag[j], xw.data());
      T step = a * w[j];
      for (int k = 0; k < D; k++) {
        T delta = xw[k] - m[k];
        if (hub > 0) { if (delta > hub) delta = hub; else if (delta < -hub) delta = -hub; }
        if (mu > 0) { T* v = &V[(size_t) j * D + k]; *v = mu * *v + step * delta; m[k] += *v; }
        else m[k] += step * delta;
      }
    }
    if (!snaps.empty() && std::find(snaps.begin(), snaps.end(), s + 1) != snaps.end()) {
      NumericMatrix Ms(K, D);
      for (int k = 0; k < K; k++) for (int dd = 0; dd < D; dd++) Ms(k, dd) = (double) M[(size_t) k * D + dd];
      snapshots.push_back(Ms);
    }
  }
  NumericMatrix Mout(K, D);
  for (int k = 0; k < K; k++) for (int dd = 0; dd < D; dd++) Mout(k, dd) = (double) M[(size_t) k * D + dd];
  return List::create(_["M"] = Mout, _["win"] = IntegerVector(win.begin(), win.end()),
                      _["nearest"] = IntegerVector(nearest.begin(), nearest.end()),
                      _["counts"] = IntegerVector(counts.begin(), counts.end()),
                      _["gap"] = NumericVector(gap.begin(), gap.end()),
                      _["comparisons"] = (double) comparisons, _["snapshots"] = snapshots);
}

// [[Rcpp::export]]
List som_train_cpp(NumericMatrix X, NumericMatrix M0, NumericMatrix grid_dist, IntegerVector pick,
                   NumericVector alpha, NumericVector radius, NumericVector temp,
                   int mode, int L, int nch, int band, double p, double gamma, int maxshift,
                   double conscience, double huber, double momentum, int prec, bool kahan,
                   IntegerVector snap_at) {
  if (X.ncol() != L * nch) stop("X has %d columns but L * nch = %d", X.ncol(), L * nch);
  if (pick.size() != alpha.size() || pick.size() != radius.size() || pick.size() != temp.size())
    stop("pick, alpha, radius and temp must have one value per step");
  Opts o = make_opts(mode, L, nch, band, p, gamma, maxshift, kahan);
  switch (prec) {
    case P_FLOAT:      return train_T<float>(X, M0, grid_dist, pick, alpha, radius, temp, o, conscience, huber, momentum, snap_at);
    case P_LONGDOUBLE: return train_T<long double>(X, M0, grid_dist, pick, alpha, radius, temp, o, conscience, huber, momentum, snap_at);
    case P_QUAD:       return train_T<__float128>(X, M0, grid_dist, pick, alpha, radius, temp, o, conscience, huber, momentum, snap_at);
    default:           return train_T<double>(X, M0, grid_dist, pick, alpha, radius, temp, o, conscience, huber, momentum, snap_at);
  }
}

// ----------------------------------------------------------------------------------- distances for judging
// every item of X against every row of M, in double: the table from which the winners and the measures come
// [[Rcpp::export]]
NumericMatrix cross_dist_cpp(NumericMatrix X, NumericMatrix M, int mode, int L, int nch, int band,
                             double p, double gamma, int maxshift) {
  const int N = X.nrow(), K = M.nrow(), D = X.ncol();
  if (D != L * nch || M.ncol() != D) stop("X and M must have L * nch columns");
  Opts o = make_opts(mode, L, nch, band, p, gamma, maxshift, false);
  std::vector<double> x(D), m((size_t) K * D), scratch;
  for (int k = 0; k < K; k++) for (int d = 0; d < D; d++) m[(size_t) k * D + d] = M(k, d);
  NumericMatrix out(N, K);
  int lag;
  for (int n = 0; n < N; n++) {
    for (int d = 0; d < D; d++) x[d] = X(n, d);
    for (int k = 0; k < K; k++) out(n, k) = distance_of(x.data(), &m[(size_t) k * D], o, scratch, &lag);
  }
  return out;
}

// the item warped, shifted or copied as the update of prototype m would see it (for pictures and checks)
// [[Rcpp::export]]
NumericVector item_seen_cpp(NumericVector x, NumericVector m, int mode, int L, int nch, int band,
                            double p, double gamma, int maxshift) {
  Opts o = make_opts(mode, L, nch, band, p, gamma, maxshift, false);
  std::vector<double> xx(x.begin(), x.end()), mm(m.begin(), m.end()), scratch, xw(x.size());
  int lag = 0;
  distance_of(xx.data(), mm.data(), o, scratch, &lag);
  item_for(xx.data(), mm.data(), o, scratch, lag, xw.data());
  return NumericVector(xw.begin(), xw.end());
}

// ----------------------------------------------------------------------------------- topographic product
// Bauer & Pawelzik (1992): P = 1 / (K (K - 1)) sum_j sum_k log P3(j, k). dV: distances between prototypes
// in the data space; dA: distances on the map. Ties on the map are broken by the data distance.
// [[Rcpp::export]]
double topographic_product_cpp(NumericMatrix dV, NumericMatrix dA) {
  const int K = dV.nrow();
  if (K < 3) return 0;
  double total = 0;
  std::vector<int> nV(K - 1), nA(K - 1), idx(K - 1);
  for (int j = 0; j < K; j++) {
    int c = 0;
    for (int k = 0; k < K; k++) if (k != j) idx[c++] = k;
    std::vector<int> oV = idx, oA = idx;
    std::sort(oV.begin(), oV.end(), [&](int a, int b) {
      if (dV(j, a) != dV(j, b)) return dV(j, a) < dV(j, b);
      if (dA(j, a) != dA(j, b)) return dA(j, a) < dA(j, b);
      return a < b; });
    std::sort(oA.begin(), oA.end(), [&](int a, int b) {
      if (dA(j, a) != dA(j, b)) return dA(j, a) < dA(j, b);
      if (dV(j, a) != dV(j, b)) return dV(j, a) < dV(j, b);
      return a < b; });
    double logprod = 0;
    for (int k = 1; k <= K - 1; k++) {
      double q1 = dV(j, oA[k - 1]) / dV(j, oV[k - 1]);
      double q2 = dA(j, oA[k - 1]) / dA(j, oV[k - 1]);
      if (!(q1 > 0) || !(q2 > 0) || !std::isfinite(q1) || !std::isfinite(q2)) { q1 = q2 = 1; }
      logprod += std::log(q1) + std::log(q2);
      total += logprod / (2.0 * k);
    }
  }
  return total / ((double) K * (K - 1));
}

// shortest ways along the map between all units, each step as long as the distance between the two
// prototypes (for the Kaski-Lagus error); 'links' is 1 where two units are direct neighbours
// [[Rcpp::export]]
NumericMatrix map_ways_cpp(NumericMatrix dV, IntegerMatrix links) {
  const int K = dV.nrow();
  const double INF = std::numeric_limits<double>::infinity();
  NumericMatrix P(K, K);
  for (int a = 0; a < K; a++) for (int b = 0; b < K; b++) P(a, b) = (a == b) ? 0 : (links(a, b) ? dV(a, b) : INF);
  for (int k = 0; k < K; k++) for (int a = 0; a < K; a++) {
    double pak = P(a, k);
    if (pak == INF) continue;
    for (int b = 0; b < K; b++) { double v = pak + P(k, b); if (v < P(a, b)) P(a, b) = v; }
  }
  return P;
}
