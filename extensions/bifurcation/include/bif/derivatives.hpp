#pragma once

#include "bif/provider.hpp"

namespace bif {

class MapDerivative {
public:
  virtual ~MapDerivative() = default;
  virtual OperationStatus apply(ConstState q, const Parameters &p,
                                const FlowMapSpec &spec, ConstState direction,
                                State &out) = 0;
  virtual OperationStatus parameterAction(ConstState q, const Parameters &p,
                                          const FlowMapSpec &spec,
                                          const std::string &parameter,
                                          State &out) = 0;
};

class FiniteDifferenceMapDerivative final : public MapDerivative {
public:
  struct Options {
    double relativeStep{1.0e-5};
    double parameterStep{1.0e-5};
  };

  explicit FiniteDifferenceMapDerivative(FlowMapProvider &provider);
  FiniteDifferenceMapDerivative(FlowMapProvider &provider, Options options);
  OperationStatus apply(ConstState q, const Parameters &p,
                        const FlowMapSpec &spec, ConstState direction,
                        State &out) override;
  OperationStatus parameterAction(ConstState q, const Parameters &p,
                                 const FlowMapSpec &spec,
                                 const std::string &parameter,
                                 State &out) override;
  const Options &options() const { return options_; }

private:
  FlowMapProvider &provider_;
  Options options_;
};

class FlowMapResidual final : public MapDerivative {
public:
  FlowMapResidual(FlowMapProvider &provider, MapDerivative &mapDerivative)
      : provider_(provider), mapDerivative_(mapDerivative) {}

  OperationStatus residual(ConstState q, const Parameters &p,
                           const FlowMapSpec &spec, State &out);
  OperationStatus apply(ConstState q, const Parameters &p,
                        const FlowMapSpec &spec, ConstState direction,
                        State &out) override;
  OperationStatus parameterAction(ConstState q, const Parameters &p,
                                 const FlowMapSpec &spec,
                                 const std::string &parameter,
                                 State &out) override;

private:
  FlowMapProvider &provider_;
  MapDerivative &mapDerivative_;
};

} // namespace bif
