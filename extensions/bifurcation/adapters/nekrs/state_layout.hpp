#pragma once

#include "bif/state.hpp"
#include "bif/status.hpp"
#include <functional>

namespace bif::nekrs {

struct LayoutCallbacks {
  std::function<OperationStatus(StateLayout &)> describe;
  std::function<OperationStatus(ConstState, State &)> scatter;
  std::function<OperationStatus(State &)> gather;
};

inline OperationStatus validateLayout(const StateLayout &layout) {
  if (!layout.valid()) return OperationStatus::error(StatusCode::invalid_state, "nekRS layout is not a unique weighted free-DOF layout");
  for (const auto w : layout.metricWeights)
    if (!(w > 0.0)) return OperationStatus::error(StatusCode::invalid_state, "kinetic-energy metric weights must be positive");
  return OperationStatus::ok();
}

} // namespace bif::nekrs
