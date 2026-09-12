#include "bif/spectrum.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace bif {

namespace {
std::complex<double> eigen2(double a, double b, double c, double d) {
  const double tr = a + d;
  const double det = a * d - b * c;
  return {0.5 * tr, 0.5 * std::sqrt(std::max(0.0, 4.0 * det - tr * tr))};
}
}

SpectrumResult DenseMapEigensolver::solve(const StateLayout& layout, ConstState q,
                                         const Parameters &p, const FlowMapSpec &spec,
                                         MapDerivative &derivative, int nev) {
  SpectrumResult result;
  const std::size_t n = q.size();
  if (n == 0 || nev <= 0 || layout.localSize != n) {
    result.status = OperationStatus::error(StatusCode::invalid_state, "empty spectrum request");
    return result;
  }
  if (n > 256) {
    result.status = OperationStatus::error(StatusCode::not_supported,
                                           "dense eigensolver is only for small examples; use Anasazi for PDE states");
    return result;
  }
  std::vector<std::vector<double>> a(n, std::vector<double>(n));
  for (std::size_t j = 0; j < n; ++j) {
    State e(n, 0.0), column;
    e[j] = 1.0;
    auto status = derivative.apply(q, p, spec, e, column);
    ++result.mapApplications;
    if (!status) { result.status = status; return result; }
    if (column.size() != n) {
      result.status = OperationStatus::error(StatusCode::solver_failure, "operator returned wrong dimension");
      return result;
    }
    for (std::size_t i = 0; i < n; ++i) a[i][j] = column[i];
  }

  // Unshifted real QR is intentionally limited to the tiny reference examples.
  // Keeping this implementation here prevents a dense global PDE matrix from
  // accidentally becoming part of the production path.
  for (int iteration = 0; iteration < 400; ++iteration) {
    std::vector<std::vector<double>> qmat(n, std::vector<double>(n));
    std::vector<std::vector<double>> r(n, std::vector<double>(n));
    for (std::size_t j = 0; j < n; ++j) {
      std::vector<double> v(n);
      for (std::size_t i = 0; i < n; ++i) v[i] = a[i][j];
      for (std::size_t k = 0; k < j; ++k) {
        for (std::size_t i = 0; i < n; ++i) r[k][j] += qmat[i][k] * v[i];
        for (std::size_t i = 0; i < n; ++i) v[i] -= r[k][j] * qmat[i][k];
      }
      double norm = 0.0; for (double x : v) norm += x * x;
      r[j][j] = std::sqrt(norm);
      if (r[j][j] < 1e-14) { qmat[j][j] = 1.0; continue; }
      for (std::size_t i = 0; i < n; ++i) qmat[i][j] = v[i] / r[j][j];
    }
    std::vector<std::vector<double>> next(n, std::vector<double>(n));
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j)
        for (std::size_t k = 0; k < n; ++k) next[i][j] += r[i][k] * qmat[k][j];
    a.swap(next);
  }

  const int count = std::min<int>(nev, static_cast<int>(n));
  for (int i = 0; i < count; ++i) {
    std::complex<double> lambda(a[i][i], 0.0);
    if (i + 1 < static_cast<int>(n) && std::abs(a[i + 1][i]) > 1e-7) {
      lambda = eigen2(a[i][i], a[i][i + 1], a[i + 1][i], a[i + 1][i + 1]);
      ++i;
    }
    Eigenpair pair; pair.pairId = static_cast<int>(result.eigenpairs.size()); pair.value = lambda; pair.vector.assign(n, 0.0);
    pair.vector[std::min<std::size_t>(i, n - 1)] = 1.0;
    result.eigenpairs.push_back(std::move(pair));
  }
  std::sort(result.eigenpairs.begin(), result.eigenpairs.end(),
            [](const Eigenpair &x, const Eigenpair &y) { return std::abs(x.value) > std::abs(y.value); });
  result.status = OperationStatus::ok();
  return result;
}

PairMatch matchPair(const SpectrumResult &current, std::complex<double> previous,
                    ConstState previousVector, const StateLayout &layout) {
  PairMatch match;
  double best = -1.0;
  for (const auto &candidate : current.eigenpairs) {
    if (candidate.vector.size() != previousVector.size()) continue;
    const double denom = layout.norm(candidate.vector) * layout.norm(previousVector);
    const double overlap = denom > 0.0 ? std::abs(layout.dot(candidate.vector, previousVector)) / denom : 0.0;
    const double distance = std::abs(candidate.value - previous);
    const double score = overlap / (1.0 + distance);
    if (score > best) { best = score; match.value = candidate.value; match.overlap = overlap; match.pairId = candidate.pairId; match.matched = true; }
  }
  return match;
}

} // namespace bif
