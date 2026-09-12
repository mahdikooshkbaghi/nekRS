#include "bif/workflow.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <functional>

namespace bif {
namespace {
std::vector<double> solveDense(std::vector<std::vector<double>> a, std::vector<double> b, bool &ok) {
  const std::size_t n = b.size(); ok = true;
  for (std::size_t k = 0; k < n; ++k) {
    std::size_t pivot = k; for (std::size_t i = k + 1; i < n; ++i) if (std::abs(a[i][k]) > std::abs(a[pivot][k])) pivot = i;
    if (std::abs(a[pivot][k]) < 1e-14) { ok = false; return {}; }
    std::swap(a[pivot], a[k]); std::swap(b[pivot], b[k]);
    for (std::size_t i = k + 1; i < n; ++i) { const double f = a[i][k] / a[k][k]; for (std::size_t j = k; j < n; ++j) a[i][j] -= f * a[k][j]; b[i] -= f * b[k]; }
  }
  std::vector<double> x(n); for (std::size_t i = n; i-- > 0;) { x[i] = b[i]; for (std::size_t j = i + 1; j < n; ++j) x[i] -= a[i][j] * x[j]; x[i] /= a[i][i]; }
  return x;
}

double norm(const State &x) { double s = 0; for (double v : x) s += v * v; return std::sqrt(s); }

OperationStatus gmres(const StateLayout& layout,
                      const std::function<OperationStatus(ConstState, State&)>& apply,
                      ConstState rhs, State& solution, int restart, int cycles)
{
  solution.assign(rhs.size(), 0.0);
  for (int cycle = 0; cycle < cycles; ++cycle) {
    State Ax, residual = rhs;
    auto status = apply(solution, Ax); if (!status) return status;
    for (std::size_t i = 0; i < residual.size(); ++i) residual[i] -= Ax[i];
    const double beta = layout.norm(residual);
    if (beta < 1.0e-12) return OperationStatus::ok();
    std::vector<State> basis(static_cast<std::size_t>(restart + 1), State(rhs.size(), 0.0));
    for (std::size_t i = 0; i < rhs.size(); ++i) basis[0][i] = residual[i] / beta;
    std::vector<std::vector<double>> h(restart + 1, std::vector<double>(restart, 0.0));
    State best = solution;
    for (int j = 0; j < restart; ++j) {
      State w; status = apply(basis[j], w); if (!status) return status;
      for (int i = 0; i <= j; ++i) {
        h[i][j] = layout.dot(basis[i], w);
        for (std::size_t k = 0; k < w.size(); ++k) w[k] -= h[i][j] * basis[i][k];
      }
      h[j + 1][j] = layout.norm(w);
      if (h[j + 1][j] > 1.0e-14) for (std::size_t k = 0; k < w.size(); ++k) basis[j + 1][k] = w[k] / h[j + 1][j];
      const int m = j + 1;
      std::vector<std::vector<double>> normal(m, std::vector<double>(m, 0.0));
      std::vector<double> normalRhs(m, 0.0);
      for (int row = 0; row < m; ++row) {
        for (int col = 0; col < m; ++col) for (int k = 0; k <= j + 1; ++k) normal[row][col] += h[k][row] * h[k][col];
        normalRhs[row] = h[0][row] * beta;
      }
      bool solved = false; auto coefficients = solveDense(std::move(normal), std::move(normalRhs), solved);
      if (!solved) break;
      best = solution; for (int col = 0; col < m; ++col) for (std::size_t k = 0; k < best.size(); ++k) best[k] += coefficients[col] * basis[col][k];
      State trialImage; status = apply(best, trialImage); if (!status) return status;
      double trialNorm = 0.0; for (std::size_t k = 0; k < rhs.size(); ++k) { const double value = rhs[k] - trialImage[k]; trialNorm += value * value; }
      if (std::sqrt(trialNorm) / std::max(1.0, layout.norm(rhs)) < 1.0e-7) { solution = best; return OperationStatus::ok(); }
    }
    solution = best;
  }
  return OperationStatus::error(StatusCode::not_converged, "matrix-free GMRES reached its iteration limit");
}
}

CorrectionResult correctEquilibrium(FlowMapProvider &provider, ConstState initial,
                                    const Parameters &parameters, const FlowMapSpec &spec,
                                    const NewtonOptions &options) {
  CorrectionResult result; result.state = initial;
  auto valid = provider.validate(initial); if (!valid) { result.status = valid; return result; }
  auto snapshot = provider.capture();
  if (!snapshot) {
    result.status = OperationStatus::error(StatusCode::solver_failure, "unable to capture provider state before correction");
    return result;
  }
  FiniteDifferenceMapDerivative mapDerivative(provider, {options.finiteDifferenceStep, options.finiteDifferenceStep});
  FlowMapResidual residual(provider, mapDerivative);
  const auto layout = provider.layout();
  for (int iteration = 0; iteration < options.maxIterations; ++iteration) {
    State f; auto status = residual.residual(result.state, parameters, spec, f); ++result.mapEvaluations;
    if (!status) { result.status = status; break; }
    result.residual = layout.norm(f) / std::max(1.0, layout.norm(result.state)); result.iterations = iteration;
    if (result.residual <= options.tolerance) { result.status = OperationStatus::ok(); break; }
    const std::size_t n = result.state.size();
    State rhs = f; for (double& value : rhs) value = -value;
    State step;
    if (n <= 256) {
      std::vector<std::vector<double>> jac(n, std::vector<double>(n));
      for (std::size_t j = 0; j < n; ++j) { State e(n, 0.0), column; e[j] = 1.0; status = residual.apply(result.state, parameters, spec, e, column); result.mapEvaluations += 2; if (!status) break; for (std::size_t i = 0; i < n; ++i) jac[i][j] = column[i]; }
      if (!status) { result.status = status; break; }
      bool solved = false; step = solveDense(std::move(jac), std::move(rhs), solved);
      if (!solved) { result.status = OperationStatus::error(StatusCode::solver_failure, "singular flow-map Jacobian"); break; }
    } else {
      auto apply = [&](ConstState vector, State& output) {
        result.mapEvaluations += 2;
        return residual.apply(result.state, parameters, spec, vector, output);
      };
      status = gmres(layout, apply, rhs, step, 40, 10);
      if (!status) { result.status = status; break; }
    }
    for (std::size_t i = 0; i < n; ++i) result.state[i] += step[i];
    if (iteration + 1 == options.maxIterations)
      result.status = OperationStatus::error(StatusCode::not_converged, "Newton correction limit reached");
  }
  const auto restoreStatus = provider.restore(*snapshot);
  if (!restoreStatus) result.status = restoreStatus;
  return result;
}

LocalizationResult localizeHopf(FlowMapProvider &provider, MapEigensolver &eigensolver,
                                MapDerivative &derivative, State loState,
                                Parameters loParameters, State hiState,
                                Parameters hiParameters, const FlowMapSpec &spec,
                                int maxIterations, double bracketTolerance) {
  LocalizationResult result;
  auto growth = [&](ConstState state, const Parameters &parameters, double &g) -> OperationStatus {
    const auto spectrum = eigensolver.solve(provider.layout(), state, parameters, spec, derivative, 8);
    ++result.evaluations;
    if (!spectrum.status) return spectrum.status;
    const Eigenpair *selected = nullptr;
    for (const auto &pair : spectrum.eigenpairs)
      if (std::abs(pair.value.imag()) > 1e-8 && (!selected || std::abs(pair.value) > std::abs(selected->value))) selected = &pair;
    if (!selected) return OperationStatus::error(StatusCode::numerical_failure, "no non-real multiplier in leading spectrum");
    g = std::log(std::abs(selected->value)) / spec.horizon;
    return OperationStatus::ok();
  };
  double glo = 0, ghi = 0; auto s = growth(loState, loParameters, glo); if (!s) { result.status = s; return result; }
  s = growth(hiState, hiParameters, ghi); if (!s) { result.status = s; return result; }
  if (glo * ghi >= 0.0) { result.status = OperationStatus::error(StatusCode::invalid_parameter, "growth-rate endpoints do not bracket a crossing"); return result; }
  result.event = EventStatus::bracketed; result.bracketLo = loParameters.reynolds; result.bracketHi = hiParameters.reynolds;
  for (int iteration = 0; iteration < maxIterations && std::abs(result.bracketHi - result.bracketLo) > bracketTolerance; ++iteration) {
    const double fraction = -glo / (ghi - glo); const double parameter = loParameters.reynolds + fraction * (hiParameters.reynolds - loParameters.reynolds);
    Parameters middleParameters = loParameters; middleParameters.reynolds = parameter;
    State middle(loState.size()); for (std::size_t i = 0; i < middle.size(); ++i) middle[i] = loState[i] + fraction * (hiState[i] - loState[i]);
    auto corrected = correctEquilibrium(provider, middle, middleParameters, spec); if (!corrected.status) { result.status = corrected.status; return result; }
    double gm = 0; s = growth(corrected.state, middleParameters, gm); if (!s) { result.status = s; return result; }
    if (glo * gm <= 0) { hiState = corrected.state; hiParameters = middleParameters; ghi = gm; }
    else { loState = corrected.state; loParameters = middleParameters; glo = gm; }
    result.bracketLo = loParameters.reynolds; result.bracketHi = hiParameters.reynolds; result.parameter = middleParameters.reynolds; result.growthRate = gm;
  }
  result.parameter = 0.5 * (result.bracketLo + result.bracketHi); result.event = EventStatus::hopfCandidate; result.status = OperationStatus::ok();
  return result;
}

} // namespace bif
