#include "bif/state.hpp"
#include "bif/status.hpp"
#include "bif/task.hpp"
#include <algorithm>
#include <cmath>
#include <utility>

namespace bif {

const char *statusName(StatusCode code) {
  switch (code) {
  case StatusCode::success: return "success";
  case StatusCode::invalid_state: return "invalid_state";
  case StatusCode::invalid_parameter: return "invalid_parameter";
  case StatusCode::solver_failure: return "solver_failure";
  case StatusCode::not_converged: return "not_converged";
  case StatusCode::not_supported: return "not_supported";
  case StatusCode::io_error: return "io_error";
  case StatusCode::collective_failure: return "collective_failure";
  case StatusCode::numerical_failure: return "numerical_failure";
  }
  return "unknown";
}

const char *eventStatusName(EventStatus status) {
  switch (status) {
  case EventStatus::none: return "none";
  case EventStatus::bracketed: return "bracketed";
  case EventStatus::hopfCandidate: return "HOPF_CANDIDATE";
  case EventStatus::hopf: return "HOPF";
  case EventStatus::ambiguous: return "ambiguous";
  case EventStatus::rejected: return "rejected";
  }
  return "unknown";
}

const char *taskKindName(TaskKind kind) {
  switch (kind) {
  case TaskKind::seed: return "seed";
  case TaskKind::continueBranch: return "continue";
  case TaskKind::resume: return "resume";
  case TaskKind::stability: return "stability";
  case TaskKind::dns: return "dns";
  }
  return "unknown";
}

Scalar StateLayout::dot(ConstState x, ConstState y) const {
  if (x.size() != localSize || y.size() != localSize) return NAN;
  Scalar result = 0.0;
  for (std::size_t i = 0; i < localSize; ++i)
    result += metricWeights.empty() ? x[i] * y[i] : metricWeights[i] * x[i] * y[i];
  return communicator.globalSum ? communicator.globalSum(result) : result;
}

Scalar StateLayout::norm(ConstState x) const {
  const auto value = dot(x, x);
  return value >= 0.0 ? std::sqrt(value) : NAN;
}

double Parameters::get(const std::string &name, double fallback) const {
  if (name == "Re") return reynolds;
  for (std::size_t i = 0; i < names.size() && i < values.size(); ++i)
    if (names[i] == name) return values[i];
  return fallback;
}

void Parameters::set(const std::string &name, double value) {
  if (name == "Re") { reynolds = value; return; }
  for (std::size_t i = 0; i < names.size() && i < values.size(); ++i) {
    if (names[i] == name) { values[i] = value; return; }
  }
  names.push_back(name); values.push_back(value);
}

OperationStatus FlowMapProvider::validate(ConstState q) const {
  const auto l = layout();
  if (!l.valid() || q.size() != l.localSize)
    return OperationStatus::error(StatusCode::invalid_state, "state does not match provider layout");
  for (std::size_t i = 0; i < q.size(); ++i)
    if (!std::isfinite(q[i])) return OperationStatus::error(StatusCode::invalid_state, "state contains NaN or Inf");
  return OperationStatus::ok();
}

} // namespace bif
