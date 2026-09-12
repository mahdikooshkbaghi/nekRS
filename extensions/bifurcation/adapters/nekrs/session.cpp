#include "session.hpp"
#include <utility>

namespace bif::nekrs {

Session::Session(StateLayout layout, SessionCallbacks callbacks)
    : layout_(std::move(layout)), callbacks_(std::move(callbacks)) {}

OperationStatus Session::install(ConstState q, const Parameters &p) {
  if (!callbacks_.installState || !callbacks_.setParameters)
    return OperationStatus::error(StatusCode::not_supported, "nekRS state callbacks are incomplete");
  // Install the trial phase state first; applying the continuation parameter
  // afterwards prevents a serialized baseline Reynolds value from overriding
  // the parameter supplied to this evaluation.
  auto status = callbacks_.installState(q); if (!status) return status;
  return callbacks_.setParameters(p);
}

OperationStatus Session::evaluate(ConstState q, const Parameters &p,
                                  const FlowMapSpec &spec, State &out) {
  if (!spec.valid()) return OperationStatus::error(StatusCode::invalid_parameter, "invalid flow-map specification");
  if (q.size() != layout_.localSize)
    return OperationStatus::error(StatusCode::invalid_state, "state does not match nekRS layout");
  if (!callbacks_.capture || !callbacks_.restore || !callbacks_.advanceMap || !callbacks_.packState)
    return OperationStatus::error(StatusCode::not_supported, "nekRS transactional callbacks are incomplete");
  auto enclosing = callbacks_.capture();
  if (!enclosing) return OperationStatus::error(StatusCode::solver_failure, "nekRS capture callback returned no snapshot");
  auto finish = [&](OperationStatus status) {
    const auto restoreStatus = callbacks_.restore(*enclosing);
    // Restoration is part of the operation's postcondition.  Never hide a
    // failed restore behind an earlier map/pack error.
    if (!restoreStatus) return restoreStatus;
    return status;
  };
  auto status = install(q, p); if (!status) return finish(status);
  status = callbacks_.advanceMap(spec); if (!status) return finish(status);
  status = callbacks_.packState(out); if (!status) return finish(status);
  return finish(OperationStatus::ok());
}

} // namespace bif::nekrs
