#pragma once

#include "bif/provider.hpp"
#include <functional>

namespace bif::nekrs {

struct SolverSnapshotCallbacks {
  std::function<std::unique_ptr<ProviderSnapshot>()> capture;
  std::function<OperationStatus(const ProviderSnapshot &)> restore;
};

// A snapshot must include clocks, BDF/EXT/pressure histories, coefficients,
// solver workspaces affecting convergence, and registered UDF state. The
// adapter refuses a callback set that cannot provide those semantics.
inline bool complete(const SolverSnapshotCallbacks &callbacks) {
  return static_cast<bool>(callbacks.capture) && static_cast<bool>(callbacks.restore);
}

} // namespace bif::nekrs
