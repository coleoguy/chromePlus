// Sparse-uniformization Felsenstein pruning for mkn-style chromosome models.
// Drop-in inner kernel for chromePlus's likelihood at large k.

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

using namespace Rcpp;

// Compute exp(Q*t) %*% v via uniformization on a sparse Q stored as
// off-diagonal triples (qi, qj, qx) with diagonals equal to -row sums.
//
// Mathematical identity:
//   exp(Q t) v = sum_{n=0..inf} e^{-alpha t} (alpha t)^n / n! * P^n v,
// where alpha = max_i (-Q[i,i]) and P = I + Q/alpha is row-stochastic.
//
// Numerically: Pn_v starts at v, multiplied by P each iteration; the
// Poisson-weighted sum converges fast for small (alpha t).
static inline void expmv_uniform(
    int k,
    const std::vector<int>& qi,        // row indices, off-diagonal
    const std::vector<int>& qj,        // col indices, off-diagonal
    const std::vector<double>& qx,     // values, off-diagonal
    const std::vector<double>& qdiag,  // -Q[i,i], i.e. row sums of off-diag
    double alpha,
    double t,
    const double* v,
    double* result,
    std::vector<double>& Pn_v_buf,
    std::vector<double>& tmp_buf,
    double tol = 1e-14,
    int max_terms = 5000) {
  if (alpha <= 0.0 || t <= 0.0) {
    std::copy(v, v + k, result);
    return;
  }
  double at = alpha * t;
  double weight = std::exp(-at);

  // result = exp(-at) * v, Pn_v = v
  for (int i = 0; i < k; ++i) {
    result[i] = weight * v[i];
    Pn_v_buf[i] = v[i];
  }

  const double inv_alpha = 1.0 / alpha;
  const std::size_t nnz = qi.size();

  for (int n = 1; n <= max_terms; ++n) {
    // tmp = P * Pn_v = Pn_v + Q*Pn_v / alpha
    // Q*Pn_v[i] = -qdiag[i] * Pn_v[i] + sum over off-diag (i,j): qx * Pn_v[j]
    for (int i = 0; i < k; ++i) {
      tmp_buf[i] = Pn_v_buf[i] - qdiag[i] * Pn_v_buf[i] * inv_alpha;
    }
    for (std::size_t e = 0; e < nnz; ++e) {
      tmp_buf[qi[e]] += qx[e] * Pn_v_buf[qj[e]] * inv_alpha;
    }
    Pn_v_buf.swap(tmp_buf);

    weight = weight * at / n;
    for (int i = 0; i < k; ++i) {
      result[i] += weight * Pn_v_buf[i];
    }
    if (weight < tol && static_cast<double>(n) > at) break;
  }
}

// Felsenstein pruning. Tree given as postorder edge list; tip states as a
// dense matrix (n_tips x k); root prior selected via root_mode.
//
// root_mode:
//   0 = flat   (uniform prior; matches root="flat" in diversitree)
//   1 = obs    (proportional to D at root; matches root=ROOT.OBS, the default)
//   2 = given  (use root_prior vector)
// [[Rcpp::export]]
double mkn_loglik_sparse_cpp(
    int k,
    int n_tips,
    int n_internal,
    IntegerVector edge_parent,    // 0-indexed, postorder
    IntegerVector edge_child,     // 0-indexed
    NumericVector edge_length,
    IntegerVector qi,             // Q sparsity: row indices (off-diagonal)
    IntegerVector qj,             // Q sparsity: col indices
    NumericVector qx,             // Q values (off-diagonal)
    NumericMatrix tip_states,     // n_tips x k
    int root_mode,                // 0=flat, 1=obs, 2=given
    NumericVector root_prior      // used only when root_mode == 2
) {
  const int n_total = n_tips + n_internal;
  const int n_edges = edge_parent.size();
  const std::size_t nnz = qi.size();

  std::vector<int> qi_v(qi.begin(), qi.end());
  std::vector<int> qj_v(qj.begin(), qj.end());
  std::vector<double> qx_v(qx.begin(), qx.end());

  // Compute diagonal magnitudes (qdiag[i] = sum of off-diag entries in row i)
  std::vector<double> qdiag(k, 0.0);
  for (std::size_t e = 0; e < nnz; ++e) qdiag[qi_v[e]] += qx_v[e];
  double alpha = 0.0;
  for (int i = 0; i < k; ++i) if (qdiag[i] > alpha) alpha = qdiag[i];

  // D[node] is k-vector, stored row-major in flat array
  std::vector<double> D(static_cast<std::size_t>(n_total) * k, 0.0);
  for (int t = 0; t < n_tips; ++t) {
    for (int j = 0; j < k; ++j) {
      D[static_cast<std::size_t>(t) * k + j] = tip_states(t, j);
    }
  }
  std::vector<double> log_scale(n_total, 0.0);
  std::vector<bool> initialized(n_total, false);
  for (int t = 0; t < n_tips; ++t) initialized[t] = true;

  std::vector<double> contrib(k);
  std::vector<double> Pn_v_buf(k), tmp_buf(k);

  const double NEG_INF = -std::numeric_limits<double>::infinity();

  for (int e = 0; e < n_edges; ++e) {
    int parent = edge_parent[e];
    int child  = edge_child[e];
    const double* D_child = &D[static_cast<std::size_t>(child) * k];

    expmv_uniform(k, qi_v, qj_v, qx_v, qdiag, alpha,
                  edge_length[e], D_child, contrib.data(),
                  Pn_v_buf, tmp_buf);

    // Normalize contrib (should sum to ~1 since exp(Qt) is row-stochastic
    // and D_child is normalized; rounding may drift it slightly)
    double s_c = 0.0;
    for (int i = 0; i < k; ++i) s_c += contrib[i];
    if (!(s_c > 0.0)) return NEG_INF;
    double inv_sc = 1.0 / s_c;
    for (int i = 0; i < k; ++i) contrib[i] *= inv_sc;

    double child_scale = log_scale[child] + std::log(s_c);

    double* D_parent = &D[static_cast<std::size_t>(parent) * k];
    if (!initialized[parent]) {
      std::copy(contrib.begin(), contrib.end(), D_parent);
      log_scale[parent] = child_scale;
      initialized[parent] = true;
    } else {
      for (int i = 0; i < k; ++i) D_parent[i] *= contrib[i];
      double s = 0.0;
      for (int i = 0; i < k; ++i) s += D_parent[i];
      if (!(s > 0.0)) return NEG_INF;
      double inv_s = 1.0 / s;
      for (int i = 0; i < k; ++i) D_parent[i] *= inv_s;
      log_scale[parent] = log_scale[parent] + std::log(s) + child_scale;
    }
  }

  const int root = n_tips;  // first internal node, by ape convention (0-idx)
  const double* D_root = &D[static_cast<std::size_t>(root) * k];

  double loglik_extra = 0.0;
  if (root_mode == 0) {           // flat
    double s = 0.0;
    for (int i = 0; i < k; ++i) s += D_root[i];
    if (!(s > 0.0)) return NEG_INF;
    loglik_extra = std::log(s) - std::log(static_cast<double>(k));
  } else if (root_mode == 1) {    // obs — D_root is normalized to sum 1
    double s2 = 0.0;
    for (int i = 0; i < k; ++i) s2 += D_root[i] * D_root[i];
    if (!(s2 > 0.0)) return NEG_INF;
    loglik_extra = std::log(s2);
  } else {                        // given
    double s = 0.0;
    for (int i = 0; i < k; ++i) s += root_prior[i] * D_root[i];
    if (!(s > 0.0)) return NEG_INF;
    loglik_extra = std::log(s);
  }
  return log_scale[root] + loglik_extra;
}
