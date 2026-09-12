#pragma once

#include "bif/provider.hpp"
#include <functional>
#include <memory>

namespace bif::nekrs {

// The adapter owns no OCCA or nekRS object. The embedding layer supplies the
// narrow state-access operations; this keeps all nekRS concrete types here.
struct SessionState final : ProviderSnapshot {
  std::shared_ptr<void> payload;
};

struct SessionCallbacks {
  std::function<std::unique_ptr<ProviderSnapshot>()> capture;
  std::function<OperationStatus(const ProviderSnapshot &)> restore;
  std::function<OperationStatus(ConstState)> installState;
  std::function<OperationStatus(State &)> packState;
  std::function<OperationStatus(const Parameters &)> setParameters;
  std::function<OperationStatus(const FlowMapSpec &)> advanceMap;
};

class Session {
public:
  Session(StateLayout layout, SessionCallbacks callbacks);
  StateLayout layout() const { return layout_; }
  OperationStatus evaluate(ConstState q, const Parameters &p,
                           const FlowMapSpec &spec, State &out);
  std::unique_ptr<ProviderSnapshot> capture() { return callbacks_.capture(); }
  OperationStatus restore(const ProviderSnapshot &snapshot) { return callbacks_.restore(snapshot); }
  OperationStatus install(ConstState q, const Parameters &p);

private:
  StateLayout layout_;
  SessionCallbacks callbacks_;
};

} // namespace bif::nekrs
