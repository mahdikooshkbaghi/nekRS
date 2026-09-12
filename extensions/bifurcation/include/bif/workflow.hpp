#pragma once

#include "bif/checkpoint.hpp"
#include <memory>

namespace bif {

struct NewtonOptions {
  int maxIterations{12};
  double tolerance{1.0e-7};
  double finiteDifferenceStep{1.0e-6};
};

struct CorrectionResult {
  OperationStatus status;
  State state;
  double residual{0.0};
  int iterations{0};
  int mapEvaluations{0};
};

CorrectionResult correctEquilibrium(FlowMapProvider &provider, ConstState initial,
                                    const Parameters &parameters,
                                    const FlowMapSpec &spec,
                                    const NewtonOptions &options = {});

struct LocalizationResult {
  OperationStatus status;
  double parameter{0.0};
  double growthRate{0.0};
  double bracketLo{0.0};
  double bracketHi{0.0};
  int evaluations{0};
  EventStatus event{EventStatus::none};
};

LocalizationResult localizeHopf(FlowMapProvider &provider, MapEigensolver &eigensolver,
                                MapDerivative &derivative, State loState,
                                Parameters loParameters, State hiState,
                                Parameters hiParameters, const FlowMapSpec &spec,
                                int maxIterations = 12, double bracketTolerance = 0.1);

} // namespace bif
