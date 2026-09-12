#pragma once

#include "bif/task.hpp"
#include "bif/workflow.hpp"
#include "nox_group.hpp"

namespace bif::trilinos {

struct LocaOptions {
  double scaledParameter{50.0};
  double initialStep{0.02};
  double maximumStep{0.04};
  int maximumSteps{128};
};

// The owner of LOCA objects is kept in the .cpp adapter. This result is a
// plain record so checkpoint files never serialize a LOCA object graph.
struct LocaPoint {
  State state;
  Parameters parameters;
  double arclength{0.0};
};

// Runs the real LOCA Stepper around the production matrix-free group. The
// returned status is false when LOCA does not finish the requested branch.
OperationStatus runLocaContinuation(FlowMapProvider& provider,
                                    FlowMapResidual& residual,
                                    MapDerivative& derivative,
                                    const FlowMapSpec& spec,
                                    ConstState initial,
                                    const LOCA::ParameterVector& params,
                                    const LocaOptions& options = {});

} // namespace bif::trilinos
