#include "linear_solver.hpp"

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Operator.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <stdexcept>

namespace bif::trilinos {
namespace {
using Scalar = double;
using LO = int;
using GO = long long;
using Node = Tpetra::Map<LO, GO>::node_type;
using Map = Tpetra::Map<LO, GO, Node>;
using Vector = Tpetra::Vector<Scalar, LO, GO, Node>;
using MV = Tpetra::MultiVector<Scalar, LO, GO, Node>;
using OP = Tpetra::Operator<Scalar, LO, GO, Node>;

class CallbackOperator final : public OP {
public:
  CallbackOperator(Teuchos::RCP<const Map> map, Apply apply) : map_(std::move(map)), apply_(std::move(apply)) {}
  void apply(const MV &x, MV &y, Teuchos::ETransp trans = Teuchos::NO_TRANS,
             Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const override {
    if (trans != Teuchos::NO_TRANS) throw std::logic_error("transpose callback is not available");
    auto xv = x.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto yv = y.getLocalViewHost(Tpetra::Access::ReadWrite);
    State input(xv.extent(0)), output;
    for (std::size_t i = 0; i < input.size(); ++i) input[i] = xv(i, 0);
    auto status = apply_(input, output);
    if (!status || output.size() != input.size()) throw std::runtime_error(status.message.empty() ? "operator callback failed" : status.message);
    for (std::size_t i = 0; i < output.size(); ++i) yv(i, 0) = beta * yv(i, 0) + alpha * output[i];
  }
  bool hasTransposeApply() const override { return false; }
  Teuchos::RCP<const Map> getDomainMap() const override { return map_; }
  Teuchos::RCP<const Map> getRangeMap() const override { return map_; }
private:
  Teuchos::RCP<const Map> map_; Apply apply_;
};
}

LinearSolveResult solveGMRES(const StateLayout &layout, const Apply &apply,
                             ConstState rhs, State &solution,
                             const LinearSolverOptions &options) {
  LinearSolveResult result;
  if (rhs.size() != layout.localSize || !layout.valid()) { result.status = OperationStatus::error(StatusCode::invalid_state, "invalid linear-solver layout"); return result; }
  auto comm = Tpetra::getDefaultComm();
  if (comm->getSize() != layout.communicator.size || comm->getRank() != layout.communicator.rank) {
    result.status = OperationStatus::error(StatusCode::collective_failure, "Tpetra and nekRS communicators disagree");
    return result;
  }
  Teuchos::Array<GO> gids(layout.localSize); for (std::size_t i = 0; i < layout.localSize; ++i) gids[i] = static_cast<GO>(layout.globalIds[i]);
  auto map = Teuchos::rcp(new Map(static_cast<GO>(layout.globalSize), gids(), 0, comm));
  auto A = Teuchos::rcp(new CallbackOperator(map, apply));
  auto X = Teuchos::rcp(new MV(map, 1, true)), B = Teuchos::rcp(new MV(map, 1, false));
  auto xv = X->getLocalViewHost(Tpetra::Access::ReadWrite), bv = B->getLocalViewHost(Tpetra::Access::ReadWrite);
  for (std::size_t i = 0; i < layout.localSize; ++i) { xv(i, 0) = i < solution.size() ? solution[i] : 0.0; bv(i, 0) = rhs[i]; }
  auto problem = Teuchos::rcp(new Belos::LinearProblem<Scalar, MV, OP>(A, X, B));
  if (!problem->setProblem()) { result.status = OperationStatus::error(StatusCode::solver_failure, "Belos rejected linear problem"); return result; }
  Teuchos::ParameterList params; params.set("Maximum Iterations", options.maximumIterations); params.set("Num Blocks", options.restart); params.set("Convergence Tolerance", options.tolerance);
  params.set("Flexible Gmres", options.flexible);
  auto solver = Teuchos::rcp(new Belos::PseudoBlockGmresSolMgr<Scalar, MV, OP>(problem, Teuchos::rcpFromRef(params)));
  const auto solveStatus = solver->solve(); result.iterations = solver->getNumIters();
  auto xread = X->getLocalViewHost(Tpetra::Access::ReadOnly); solution.resize(layout.localSize); for (std::size_t i = 0; i < layout.localSize; ++i) solution[i] = xread(i, 0);
  result.status = solveStatus == Belos::Converged ? OperationStatus::ok() : OperationStatus::error(StatusCode::not_converged, "Belos GMRES did not converge");
  result.status.iterations = result.iterations; return result;
}

} // namespace bif::trilinos
