#pragma once

#include "bif/provider.hpp"
#include <functional>
#include <utility>

namespace bif::nekrs {

class ReynoldsBinding final : public ParameterBinding {
public:
  explicit ReynoldsBinding(std::function<OperationStatus(double)> setViscosity)
      : setViscosity_(std::move(setViscosity)) {}
  OperationStatus set(const Parameters &parameters) override {
    const double re = parameters.get("Re");
    if (!(re > 0.0)) return OperationStatus::error(StatusCode::invalid_parameter, "Re must be positive");
    if (!setViscosity_) return OperationStatus::error(StatusCode::not_supported, "viscosity callback is absent");
    return setViscosity_(1.0 / re);
  }
private:
  std::function<OperationStatus(double)> setViscosity_;
};

} // namespace bif::nekrs
