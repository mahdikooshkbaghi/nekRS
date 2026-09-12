#pragma once

#include "bif/state.hpp"
#include "bif/status.hpp"
#include <memory>

namespace bif {

class ProviderSnapshot {
public:
  virtual ~ProviderSnapshot() = default;
};

class FlowMapProvider {
public:
  virtual ~FlowMapProvider() = default;
  virtual StateLayout layout() const = 0;
  virtual OperationStatus evaluate(ConstState q, const Parameters &p,
                                   const FlowMapSpec &spec, State &out) = 0;
  virtual std::unique_ptr<ProviderSnapshot> capture() = 0;
  virtual OperationStatus restore(const ProviderSnapshot &snapshot) = 0;

  // Accepted-state installation is deliberately separate from trial evaluation.
  virtual OperationStatus install(ConstState q, const Parameters &p) = 0;
  virtual OperationStatus validate(ConstState q) const;
};

class ParameterBinding {
public:
  virtual ~ParameterBinding() = default;
  virtual OperationStatus set(const Parameters &p) = 0;
};

} // namespace bif
