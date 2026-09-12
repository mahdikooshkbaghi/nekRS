#include "api.hpp"
#include "nekrs.hpp"
#include <mpi.h>
#include <memory>
#include <utility>

namespace bif::nekrs {
namespace {
class Snapshot final : public ProviderSnapshot {
public:
  explicit Snapshot(::nekrs::AnalysisState value) : state(std::move(value)) {}
  ::nekrs::AnalysisState state;
};

OperationStatus convert(const ::nekrs::AnalysisOperationStatus& status)
{
  return status ? OperationStatus::ok() : OperationStatus::error(StatusCode::solver_failure, status.message);
}
}

std::unique_ptr<Session> makeSession()
{
  ::nekrs::AnalysisLayout source;
  auto layoutStatus = ::nekrs::analysisLayout(source);
  if (!layoutStatus) return nullptr;
  StateLayout layout;
  layout.globalSize = source.globalSize; layout.localSize = source.localSize; layout.globalOffset = source.globalOffset;
  layout.globalIds = source.globalIds; layout.freeDofs = source.freeDofs; layout.metricWeights = source.metricWeights;
  layout.fieldDescription = "nekRS unique free velocity state";
  MPI_Comm_rank(MPI_COMM_WORLD, &layout.communicator.rank);
  MPI_Comm_size(MPI_COMM_WORLD, &layout.communicator.size);
  layout.communicator.globalSum = [](double local) {
    double global = 0.0;
    MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global;
  };

  SessionCallbacks callbacks;
  callbacks.capture = []() -> std::unique_ptr<ProviderSnapshot> {
    ::nekrs::AnalysisState state;
    if (!::nekrs::captureAnalysisState(state)) return nullptr;
    return std::make_unique<Snapshot>(std::move(state));
  };
  callbacks.restore = [](const ProviderSnapshot& base) -> OperationStatus {
    const auto* snapshot = dynamic_cast<const Snapshot*>(&base);
    if (!snapshot) return OperationStatus::error(StatusCode::invalid_state, "foreign nekRS snapshot");
    return convert(::nekrs::restoreAnalysisState(snapshot->state));
  };
  callbacks.installState = [](ConstState state) -> OperationStatus {
    return convert(::nekrs::unpackAnalysisState(state));
  };
  callbacks.packState = [](State& state) -> OperationStatus {
    return convert(::nekrs::packAnalysisState(state));
  };
  callbacks.setParameters = [](const Parameters& parameters) -> OperationStatus {
    return convert(::nekrs::setAnalysisReynolds(parameters.get("Re")));
  };
  callbacks.advanceMap = [](const FlowMapSpec& spec) -> OperationStatus {
    if (!spec.valid() || spec.steps <= 0) return OperationStatus::error(StatusCode::invalid_parameter, "invalid nekRS map specification");
    double time = ::nekrs::analysisTime();
    int step = ::nekrs::timeStep();
    if (spec.suppressOutput) ::nekrs::checkpointStep(0);
    for (int i = 0; i < spec.steps; ++i) {
      ++step;
      ::nekrs::initStep(time, spec.timestep, step);
      // runStep() performs one complete nekRS time-step solve.  Its integer
      // argument is the solver stage/corrector index, not a request to repeat
      // the whole step; repeating it here would advance the same timestep with
      // stale history and falsely turn the map into a nonlinear inner solve.
      if (!::nekrs::runStep(1) || !::nekrs::stepConverged())
        return OperationStatus::error(StatusCode::not_converged, "nekRS time-step solve did not converge");
      ::nekrs::finishStep();
      time += spec.timestep;
    }
    return OperationStatus::ok();
  };
  return std::make_unique<Session>(std::move(layout), std::move(callbacks));
}

std::unique_ptr<bif::FlowMapProvider> makeFlowMapProvider()
{
  auto session = makeSession();
  if (!session) return nullptr;
  // The provider must retain the session. The owning wrapper is intentionally
  // local to this adapter and never exposes a borrowed nekRS object.
  class OwnedProvider final : public bif::FlowMapProvider {
  public:
    explicit OwnedProvider(std::unique_ptr<Session> value) : session(std::move(value)) {}
    StateLayout layout() const override { return session->layout(); }
    OperationStatus evaluate(ConstState q, const Parameters& p, const FlowMapSpec& s, State& out) override { return session->evaluate(q, p, s, out); }
    std::unique_ptr<ProviderSnapshot> capture() override { return session->capture(); }
    OperationStatus restore(const ProviderSnapshot& snapshot) override { return session->restore(snapshot); }
    OperationStatus install(ConstState q, const Parameters& p) override { return session->install(q, p); }
  private:
    std::unique_ptr<Session> session;
  };
  return std::make_unique<OwnedProvider>(std::move(session));
}

} // namespace bif::nekrs
