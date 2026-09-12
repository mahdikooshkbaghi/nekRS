#include "map_eigensolver.hpp"

#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziBlockKrylovSchurSolMgr.hpp>
#include <AnasaziTpetraAdapter.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_Operator.hpp>
#include <Tpetra_Vector.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <stdexcept>
#include <functional>

namespace bif::trilinos {
namespace {
using Scalar = double;
using LO = int;
using GO = long long;
using Map = Tpetra::Map<LO, GO>;
using MV = Tpetra::MultiVector<Scalar, LO, GO>;
using OP = Tpetra::Operator<Scalar, LO, GO>;

class CallbackOperator final : public OP {
  using ApplyFunction = std::function<OperationStatus(ConstState, State&)>;
public:
  CallbackOperator(Teuchos::RCP<const Map> map, ApplyFunction callback)
      : map_(std::move(map)), callback_(std::move(callback)) {}
  void apply(const MV& x, MV& y, Teuchos::ETransp trans = Teuchos::NO_TRANS,
             Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const override {
    if (trans != Teuchos::NO_TRANS) throw std::logic_error("map derivative has no transpose action");
    auto xView = x.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto yView = y.getLocalViewHost(Tpetra::Access::ReadWrite);
    for (std::size_t column = 0; column < x.getNumVectors(); ++column) {
      State input(xView.extent(0)), output;
      for (std::size_t row = 0; row < input.size(); ++row) input[row] = xView(row, column);
      auto status = callback_(input, output);
      if (!status || output.size() != input.size()) throw std::runtime_error(status.message.empty() ? "map derivative callback failed" : status.message);
      for (std::size_t row = 0; row < output.size(); ++row) yView(row, column) = alpha * output[row] + beta * yView(row, column);
    }
  }
  bool hasTransposeApply() const override { return false; }
  Teuchos::RCP<const Map> getDomainMap() const override { return map_; }
  Teuchos::RCP<const Map> getRangeMap() const override { return map_; }
private:
  Teuchos::RCP<const Map> map_;
  ApplyFunction callback_;
};
}

SpectrumResult MapEigensolver::solve(const StateLayout& layout, ConstState q,
                                     const Parameters& p, const FlowMapSpec& spec,
                                     bif::MapDerivative& derivative, int nev) {
  SpectrumResult result;
  if (q.empty() || nev <= 0 || !layout.valid() || layout.localSize != q.size()) { result.status = OperationStatus::error(StatusCode::invalid_state, "empty Anasazi map problem"); return result; }
  try {
    auto comm = Tpetra::getDefaultComm();
    if (comm->getSize() != layout.communicator.size || comm->getRank() != layout.communicator.rank) {
      result.status = OperationStatus::error(StatusCode::collective_failure, "Tpetra and nekRS communicators disagree");
      return result;
    }
    Teuchos::Array<GO> gids(layout.localSize); for (std::size_t i = 0; i < layout.localSize; ++i) gids[i] = static_cast<GO>(layout.globalIds[i]);
    auto map = Teuchos::rcp(new Map(static_cast<GO>(layout.globalSize), gids(), 0, comm));
    auto callback = [&](ConstState direction, State& output) { return derivative.apply(q, p, spec, direction, output); };
    auto operatorObject = Teuchos::rcp(new CallbackOperator(map, callback));
    auto initial = Teuchos::rcp(new MV(map, std::max(1, std::min(options_.eigenpairs, options_.subspace / 2)), true));
    initial->randomize();
    auto problem = Teuchos::rcp(new Anasazi::BasicEigenproblem<Scalar, MV, OP>(operatorObject, initial));
    problem->setNEV(std::min(nev, options_.eigenpairs));
    problem->setHermitian(false);
    if (!problem->setProblem()) { result.status = OperationStatus::error(StatusCode::solver_failure, "Anasazi rejected map eigenproblem"); return result; }
    Teuchos::ParameterList parameters;
    parameters.set("Which", "LM");
    // BlockKrylovSchur defaults to a scalar block. Keeping the default avoids
    // a 17.3 Anasazi/Teuchos type mismatch for the Block Size parameter.
    parameters.set("Num Blocks", std::max(options_.subspace, 3));
    parameters.set("Maximum Restarts", options_.maximumRestarts);
    parameters.set("Convergence Tolerance", options_.residualTolerance);
    parameters.set("Verbosity", 0);
    Anasazi::BlockKrylovSchurSolMgr<Scalar, MV, OP> solver(problem, parameters);
    const auto solveStatus = solver.solve();
    const auto solution = problem->getSolution();
    if (!solution.Evecs || solution.numVecs == 0) { result.status = OperationStatus::error(StatusCode::not_converged, "Anasazi returned no map multipliers"); return result; }
    for (int i = 0; i < solution.numVecs && i < static_cast<int>(solution.Evals.size()); ++i) {
      Eigenpair pair; pair.pairId = i; pair.value = std::complex<double>(solution.Evals[i].realpart, solution.Evals[i].imagpart);
      const int encoding = i < static_cast<int>(solution.index.size()) ? solution.index[i] : 0;
      int realIndex = i;
      int imagIndex = -1;
      double imagSign = 1.0;
      if (encoding > 0) imagIndex = i + 1;
      else if (encoding < 0) { realIndex = i - 1; imagIndex = i; imagSign = -1.0; }
      auto values = solution.Evecs->getLocalViewHost(Tpetra::Access::ReadOnly);
      if (realIndex < 0 || realIndex >= static_cast<int>(values.extent(1)) ||
          (imagIndex >= static_cast<int>(values.extent(1)))) {
        result.status = OperationStatus::error(StatusCode::solver_failure, "Anasazi returned an invalid eigenvector encoding");
        return result;
      }
      pair.vector.resize(q.size());
      State imaginary(imagIndex >= 0 ? q.size() : 0);
      for (std::size_t row = 0; row < q.size(); ++row) {
        pair.vector[row] = values(row, realIndex);
        if (imagIndex >= 0) imaginary[row] = imagSign * values(row, imagIndex);
      }
      State image; auto applyStatus = derivative.apply(q, p, spec, pair.vector, image);
      ++result.mapApplications;
      if (!applyStatus) { result.status = applyStatus; return result; }
      State residual(image.size());
      if (imagIndex < 0) {
        for (std::size_t row = 0; row < image.size(); ++row)
          residual[row] = image[row] - pair.value.real() * pair.vector[row];
      } else {
        State imageImaginary; applyStatus = derivative.apply(q, p, spec, imaginary, imageImaginary);
        ++result.mapApplications;
        if (!applyStatus) { result.status = applyStatus; return result; }
        const double a = pair.value.real(), b = pair.value.imag();
        for (std::size_t row = 0; row < image.size(); ++row)
          residual[row] = image[row] - a * pair.vector[row] + b * imaginary[row];
        const double realResidual = layout.norm(residual);
        for (std::size_t row = 0; row < imageImaginary.size(); ++row)
          residual[row] = imageImaginary[row] - a * imaginary[row] - b * pair.vector[row];
        pair.residual = std::hypot(realResidual, layout.norm(residual)) /
                        std::max(1.0, std::hypot(layout.norm(pair.vector), layout.norm(imaginary)));
      }
      if (imagIndex < 0)
        pair.residual = layout.norm(residual) / std::max(1.0, layout.norm(pair.vector));
      result.eigenpairs.push_back(std::move(pair));
    }
    result.status = solveStatus == Anasazi::Converged ? OperationStatus::ok() : OperationStatus::error(StatusCode::not_converged, "Anasazi map solve did not converge");
    return result;
  } catch (const std::exception& error) {
    result.status = OperationStatus::error(StatusCode::solver_failure, std::string("Anasazi map solve failed: ") + error.what()); return result;
  }
}

} // namespace bif::trilinos
