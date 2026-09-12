#include "loca_runner.hpp"

#include <LOCA_Abstract_Factory.H>
#include <LOCA_GlobalData.H>
#include <LOCA_Parameter_Vector.H>
#include <LOCA_Stepper.H>
#include <LOCA_Tpetra_Factory.hpp>
#include <NOX_StatusTest_Combo.H>
#include <NOX_StatusTest_MaxIters.H>
#include <NOX_StatusTest_NormF.H>
#include <Teuchos_ParameterList.hpp>

namespace bif::trilinos {

OperationStatus runLocaContinuation(FlowMapProvider& provider,
                                    FlowMapResidual& residual,
                                    MapDerivative& derivative,
                                    const FlowMapSpec& spec,
                                    ConstState initial,
                                    const LOCA::ParameterVector& params,
                                    const LocaOptions& options) {
  if (!spec.valid() || initial.empty() || !provider.layout().valid())
    return OperationStatus::error(StatusCode::invalid_state, "invalid LOCA continuation problem");
  if (!params.isParameter("Re"))
    return OperationStatus::error(StatusCode::invalid_parameter, "LOCA continuation requires parameter Re");

  try {
  auto parameterList = Teuchos::rcp(new Teuchos::ParameterList);
  auto& loca = parameterList->sublist("LOCA");
  auto& stepper = loca.sublist("Stepper");
  stepper.set("Continuation Parameter", "Re");
  const double initialValue = params.getValue("Re");
  const double span = std::max(options.maximumStep * options.maximumSteps, options.maximumStep);
  stepper.set("Initial Value", initialValue);
  stepper.set("Max Value", initialValue + span);
  stepper.set("Min Value", initialValue - span);
  stepper.set("Max Steps", options.maximumSteps);
  stepper.set("Max Nonlinear Iterations", 20);
  auto& stepSize = loca.sublist("Step Size");
  stepSize.set("Method", "Constant");
  stepSize.set("Initial Step Size", options.initialStep);
  stepSize.set("Max Step Size", options.maximumStep);
  stepSize.set("Min Step Size", options.initialStep);
  auto& nox = parameterList->sublist("NOX");
  nox.sublist("Printing").set("Output Information", 0);

  auto factory = Teuchos::rcp(new LOCA::Tpetra::Factory);
  auto globalData = LOCA::createGlobalData(parameterList, factory);
  NoxGroup concrete(residual, derivative, provider, spec, initial, params);
  Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup> group =
      Teuchos::rcp(new NoxGroup(concrete));
  group->setParams(params);

  auto normF = Teuchos::rcp(new NOX::StatusTest::NormF(1.0e-8));
  auto maxIters = Teuchos::rcp(new NOX::StatusTest::MaxIters(20));
  auto combo = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR,
                                                        normF, maxIters));
  LOCA::Stepper continuation(globalData, group, combo, parameterList);
  const auto result = continuation.run();
  LOCA::destroyGlobalData(globalData);
  if (result != LOCA::Abstract::Iterator::Finished)
    return OperationStatus::error(StatusCode::not_converged, "LOCA Stepper did not finish");
  return OperationStatus::ok();
  } catch (const std::exception& error) {
    return OperationStatus::error(StatusCode::solver_failure, error.what());
  }
}

int locaApiVersionAnchor() {
  return sizeof(LOCA::MultiContinuation::AbstractGroup) > 0 ? 0 : 1;
}
} // namespace bif::trilinos
