#include "bif/derivatives.hpp"
#include <algorithm>
#include <cmath>

namespace bif {
namespace {
OperationStatus evaluateRestored(FlowMapProvider &provider, ConstState q,
                                  const Parameters &p, const FlowMapSpec &spec,
                                  State &out) {
  auto snapshot = provider.capture();
  if (!snapshot)
    return OperationStatus::error(StatusCode::solver_failure, "unable to capture provider state before derivative evaluation");
  auto status = provider.evaluate(q, p, spec, out);
  const auto restoreStatus = provider.restore(*snapshot);
  // A derivative evaluation is transactional even when the map itself fails;
  // report a failed restore because the provider is no longer trustworthy.
  if (!restoreStatus) return restoreStatus;
  return status;
}
}

FiniteDifferenceMapDerivative::FiniteDifferenceMapDerivative(FlowMapProvider &provider)
    : provider_(provider), options_{} {}

FiniteDifferenceMapDerivative::FiniteDifferenceMapDerivative(FlowMapProvider &provider,
                                                             Options options)
    : provider_(provider), options_(options) {}

OperationStatus FiniteDifferenceMapDerivative::apply(ConstState q, const Parameters &p,
                                                      const FlowMapSpec &spec,
                                                      ConstState direction, State &out) {
  auto valid = provider_.validate(q);
  if (!valid) return valid;
  if (direction.size() != q.size())
    return OperationStatus::error(StatusCode::invalid_state, "direction has wrong size");
  const auto layout = provider_.layout();
  const double qnorm = std::max(1.0, layout.norm(q));
  const double vnorm = layout.norm(direction);
  if (!std::isfinite(vnorm))
    return OperationStatus::error(StatusCode::invalid_state, "direction norm is not finite");
  out.assign(q.size(), 0.0);
  if (vnorm == 0.0) return OperationStatus::ok();

  const double h = options_.relativeStep * qnorm / vnorm;
  State qp(q), qm(q), fp, fm;
  for (std::size_t i = 0; i < q.size(); ++i) {
    qp[i] += h * direction[i];
    qm[i] -= h * direction[i];
  }
  auto plus = evaluateRestored(provider_, qp, p, spec, fp);
  if (!plus) return plus;
  auto minus = evaluateRestored(provider_, qm, p, spec, fm);
  if (!minus) return minus;
  if (fp.size() != q.size() || fm.size() != q.size())
    return OperationStatus::error(StatusCode::solver_failure, "flow map returned wrong state size");
  for (std::size_t i = 0; i < q.size(); ++i) out[i] = (fp[i] - fm[i]) / (2.0 * h);
  return OperationStatus::ok();
}

OperationStatus FiniteDifferenceMapDerivative::parameterAction(ConstState q, const Parameters &p,
                                                               const FlowMapSpec &spec,
                                                               const std::string &parameter,
                                                               State &out) {
  const double center = p.get(parameter, NAN);
  if (!std::isfinite(center))
    return OperationStatus::error(StatusCode::invalid_parameter, "unknown parameter: " + parameter);
  const double h = options_.parameterStep * std::max(1.0, std::abs(center));
  if (parameter == "Re" && center - h <= 0.0)
    return OperationStatus::error(StatusCode::invalid_parameter, "centered Reynolds perturbation is non-positive");
  Parameters pp = p, pm = p;
  pp.set(parameter, center + h); pm.set(parameter, center - h);
  State fp, fm;
  auto plus = evaluateRestored(provider_, q, pp, spec, fp);
  if (!plus) return plus;
  auto minus = evaluateRestored(provider_, q, pm, spec, fm);
  if (!minus) return minus;
  if (fp.size() != fm.size())
    return OperationStatus::error(StatusCode::solver_failure, "parameter action returned inconsistent sizes");
  out.resize(fp.size());
  for (std::size_t i = 0; i < out.size(); ++i) out[i] = (fp[i] - fm[i]) / (2.0 * h);
  return OperationStatus::ok();
}

OperationStatus FlowMapResidual::residual(ConstState q, const Parameters &p,
                                          const FlowMapSpec &spec, State &out) {
  auto valid = provider_.validate(q);
  if (!valid) return valid;
  State image;
  auto status = provider_.evaluate(q, p, spec, image);
  if (!status) return status;
  if (image.size() != q.size())
    return OperationStatus::error(StatusCode::solver_failure, "flow map returned wrong state size");
  out.resize(q.size());
  for (std::size_t i = 0; i < q.size(); ++i) out[i] = q[i] - image[i];
  return status;
}

OperationStatus FlowMapResidual::apply(ConstState q, const Parameters &p,
                                       const FlowMapSpec &spec, ConstState direction,
                                       State &out) {
  auto status = mapDerivative_.apply(q, p, spec, direction, out);
  if (!status) return status;
  if (out.size() != direction.size())
    return OperationStatus::error(StatusCode::solver_failure, "derivative returned wrong size");
  for (std::size_t i = 0; i < out.size(); ++i) out[i] = direction[i] - out[i];
  return status;
}

OperationStatus FlowMapResidual::parameterAction(ConstState q, const Parameters &p,
                                                 const FlowMapSpec &spec,
                                                 const std::string &parameter,
                                                 State &out) {
  auto status = mapDerivative_.parameterAction(q, p, spec, parameter, out);
  if (!status) return status;
  for (auto &value : out) value = -value;
  return status;
}

} // namespace bif
