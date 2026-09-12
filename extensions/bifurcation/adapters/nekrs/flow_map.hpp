#pragma once

#include "bif/provider.hpp"
#include "session.hpp"

namespace bif::nekrs {

class FlowMapProvider final : public bif::FlowMapProvider {
public:
  explicit FlowMapProvider(Session &session) : session_(session) {}
  StateLayout layout() const override { return session_.layout(); }
  OperationStatus evaluate(ConstState q, const Parameters &p,
                           const FlowMapSpec &spec, State &out) override {
    return session_.evaluate(q, p, spec, out);
  }
  std::unique_ptr<ProviderSnapshot> capture() override { return session_.capture(); }
  OperationStatus restore(const ProviderSnapshot &snapshot) override { return session_.restore(snapshot); }
  OperationStatus install(ConstState q, const Parameters &p) override { return session_.install(q, p); }

private:
  Session &session_;
};

} // namespace bif::nekrs
