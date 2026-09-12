#include "nox_group.hpp"

#include <NOX_Abstract_MultiVector.H>
#include <Teuchos_ParameterList.hpp>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <random>
#include <stdexcept>

namespace bif::trilinos {
namespace {
NoxVector& vectorRef(NOX::Abstract::Vector& v) {
  auto* result = dynamic_cast<NoxVector*>(&v);
  if (!result) throw std::invalid_argument("NOX vector is not a bifurcation vector");
  return *result;
}
const NoxVector& vectorRef(const NOX::Abstract::Vector& v) {
  auto* result = dynamic_cast<const NoxVector*>(&v);
  if (!result) throw std::invalid_argument("NOX vector is not a bifurcation vector");
  return *result;
}
}

NOX::Abstract::Vector& NoxVector::init(double value) {
  std::fill(values_.begin(), values_.end(), value); return *this;
}
NOX::Abstract::Vector& NoxVector::random(bool useSeed, int seed) {
  std::mt19937 generator(useSeed ? static_cast<unsigned>(seed) : std::random_device{}());
  std::uniform_real_distribution<double> distribution(-1.0, 1.0);
  for (auto& value : values_) value = distribution(generator);
  return *this;
}
NOX::Abstract::Vector& NoxVector::abs(const NOX::Abstract::Vector& y) {
  const auto& source = vectorRef(y).values();
  if (source.size() != values_.size()) throw std::invalid_argument("NOX vector length mismatch");
  for (std::size_t i = 0; i < values_.size(); ++i) values_[i] = std::abs(source[i]);
  return *this;
}
NOX::Abstract::Vector& NoxVector::operator=(const NOX::Abstract::Vector& y) {
  values_ = vectorRef(y).values(); return *this;
}
NOX::Abstract::Vector& NoxVector::reciprocal(const NOX::Abstract::Vector& y) {
  const auto& source = vectorRef(y).values();
  if (source.size() != values_.size()) throw std::invalid_argument("NOX vector length mismatch");
  for (std::size_t i = 0; i < values_.size(); ++i) values_[i] = 1.0 / source[i];
  return *this;
}
NOX::Abstract::Vector& NoxVector::scale(double alpha) { for (auto& value : values_) value *= alpha; return *this; }
NOX::Abstract::Vector& NoxVector::scale(const NOX::Abstract::Vector& a) {
  const auto& source = vectorRef(a).values();
  if (source.size() != values_.size()) throw std::invalid_argument("NOX vector length mismatch");
  for (std::size_t i = 0; i < values_.size(); ++i) values_[i] *= source[i];
  return *this;
}
NOX::Abstract::Vector& NoxVector::update(double alpha, const NOX::Abstract::Vector& a, double gamma) {
  const auto& source = vectorRef(a).values();
  if (source.size() != values_.size()) throw std::invalid_argument("NOX vector length mismatch");
  for (std::size_t i = 0; i < values_.size(); ++i) values_[i] = alpha * source[i] + gamma * values_[i];
  return *this;
}
NOX::Abstract::Vector& NoxVector::update(double alpha, const NOX::Abstract::Vector& a, double beta,
                             const NOX::Abstract::Vector& b, double gamma) {
  const auto& av = vectorRef(a).values(); const auto& bv = vectorRef(b).values();
  if (av.size() != values_.size() || bv.size() != values_.size()) throw std::invalid_argument("NOX vector length mismatch");
  for (std::size_t i = 0; i < values_.size(); ++i) values_[i] = alpha * av[i] + beta * bv[i] + gamma * values_[i];
  return *this;
}
Teuchos::RCP<NOX::Abstract::Vector> NoxVector::clone(NOX::CopyType type) const {
  auto result = Teuchos::rcp(new NoxVector(values_.size()));
  if (type == NOX::DeepCopy) result->values_ = values_;
  return result;
}
double NoxVector::norm(NOX::Abstract::Vector::NormType type) const {
  double local = 0.0;
  if (type == MaxNorm) {
    for (double value : values_) local = std::max(local, std::abs(value));
    double global = local; MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); return global;
  }
  if (type == OneNorm) for (double value : values_) local += std::abs(value);
  else for (double value : values_) local += value * value;
  double global = 0.0; MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return type == OneNorm ? global : std::sqrt(global);
}
double NoxVector::norm(const NOX::Abstract::Vector& weights) const {
  const auto& w = vectorRef(weights).values();
  if (w.size() != values_.size()) throw std::invalid_argument("NOX vector length mismatch");
  double local = 0.0; for (std::size_t i = 0; i < values_.size(); ++i) local += w[i] * values_[i] * values_[i];
  double global = 0.0; MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); return std::sqrt(global);
}
double NoxVector::innerProduct(const NOX::Abstract::Vector& y) const {
  const auto& source = vectorRef(y).values();
  if (source.size() != values_.size()) throw std::invalid_argument("NOX vector length mismatch");
  double local = 0.0; for (std::size_t i = 0; i < values_.size(); ++i) local += values_[i] * source[i];
  double global = 0.0; MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); return global;
}
void NoxVector::print(std::ostream& stream) const { stream << "NoxVector(" << values_.size() << ")\n"; for (double v : values_) stream << std::setprecision(17) << v << '\n'; }

NoxGroup::NoxGroup(FlowMapResidual& residual, MapDerivative& derivative, FlowMapProvider& provider,
                   FlowMapSpec spec, ConstState initial, const LOCA::ParameterVector& params,
                   LinearSolverOptions linearOptions)
    : residual_(residual), derivative_(derivative), provider_(provider), spec_(std::move(spec)),
      params_(params), linearOptions_(linearOptions), x_(initial), f_(initial.size()), newton_(initial.size()) {
  if (!params_.isParameter("Re")) params_.addParameter("Re", 0.0);
}

NOX::Abstract::Group& NoxGroup::operator=(const NOX::Abstract::Group& source) { copy(source); return *this; }
void NoxGroup::copy(const NOX::Abstract::Group& source) {
  const auto* other = dynamic_cast<const NoxGroup*>(&source);
  if (!other) throw std::invalid_argument("cannot copy unrelated NOX group");
  x_.values() = other->x_.values(); f_.values() = other->f_.values(); newton_.values() = other->newton_.values(); params_ = other->params_;
  normF_ = other->normF_; fValid_ = other->fValid_; jacobianValid_ = other->jacobianValid_; newtonValid_ = other->newtonValid_;
}
void NoxGroup::setX(const NOX::Abstract::Vector& y) { x_.values() = vectorRef(y).values(); fValid_ = jacobianValid_ = newtonValid_ = false; }
void NoxGroup::computeX(const NOX::Abstract::Group& group, const NOX::Abstract::Vector& d, double step) {
  x_.values() = vectorRef(group.getX()).values(); x_.update(step, d, 1.0); fValid_ = jacobianValid_ = newtonValid_ = false;
}
NOX::Abstract::Group::ReturnType NoxGroup::status(const OperationStatus& result) const {
  return result ? NOX::Abstract::Group::Ok : NOX::Abstract::Group::Failed;
}
NOX::Abstract::Group::ReturnType NoxGroup::computeF() {
  Parameters parameters; parameters.reynolds = getParam("Re");
  const auto result = residual_.residual(x_.values(), parameters, spec_, f_.values());
  if (result) { normF_ = f_.norm(); fValid_ = true; }
  return status(result);
}
NOX::Abstract::Group::ReturnType NoxGroup::computeJacobian() { jacobianValid_ = true; return NOX::Abstract::Group::Ok; }
NOX::Abstract::Group::ReturnType NoxGroup::applyJacobian(const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const {
  const auto& in = vectorRef(input); auto& out = vectorRef(result);
  State mapped; Parameters parameters; parameters.reynolds = getParam("Re");
  const auto operation = residual_.apply(x_.values(), parameters, spec_, in.values(), mapped);
  if (!operation) return NOX::Abstract::Group::Failed; out.values() = std::move(mapped); return NOX::Abstract::Group::Ok;
}
NOX::Abstract::Group::ReturnType NoxGroup::applyJacobianInverse(Teuchos::ParameterList& params, const NOX::Abstract::Vector& input, NOX::Abstract::Vector& result) const {
  const auto& in = vectorRef(input); auto& out = vectorRef(result); State solution;
  auto options = linearOptions_; if (params.isParameter("Tolerance")) options.tolerance = params.get<double>("Tolerance");
  Parameters parameters; parameters.reynolds = getParam("Re");
  const auto operation = solveGMRES(provider_.layout(), [&](ConstState v, State& mapped) {
    return residual_.apply(x_.values(), parameters, spec_, v, mapped);
  }, in.values(), solution, options);
  if (!operation.status) return operation.status.code == StatusCode::not_converged ? NOX::Abstract::Group::NotConverged : NOX::Abstract::Group::Failed;
  out.values() = std::move(solution); return NOX::Abstract::Group::Ok;
}
NOX::Abstract::Group::ReturnType NoxGroup::computeNewton(Teuchos::ParameterList& params) {
  if (!fValid_ && computeF() != NOX::Abstract::Group::Ok) return NOX::Abstract::Group::Failed;
  NoxVector rhs = f_; rhs.scale(-1.0); const auto result = applyJacobianInverse(params, rhs, newton_); newtonValid_ = result == NOX::Abstract::Group::Ok; return result;
}

Teuchos::RCP<const NOX::Abstract::Vector> NoxGroup::getXPtr() const { return Teuchos::rcp(new NoxVector(x_)); }
Teuchos::RCP<const NOX::Abstract::Vector> NoxGroup::getFPtr() const { return Teuchos::rcp(new NoxVector(f_)); }
Teuchos::RCP<const NOX::Abstract::Vector> NoxGroup::getGradientPtr() const { return Teuchos::rcp(new NoxVector(f_)); }
Teuchos::RCP<const NOX::Abstract::Vector> NoxGroup::getNewtonPtr() const { return Teuchos::rcp(new NoxVector(newton_)); }
Teuchos::RCP<NOX::Abstract::Group> NoxGroup::clone(NOX::CopyType type) const {
  auto result = Teuchos::rcp(new NoxGroup(*this));
  if (type == NOX::ShapeCopy) { result->x_.init(0.0); result->f_.init(0.0); result->newton_.init(0.0); result->fValid_ = result->jacobianValid_ = result->newtonValid_ = false; }
  return result;
}

void NoxGroup::setParamsMulti(const std::vector<int>& ids, const NOX::Abstract::MultiVector::DenseMatrix& values) {
  for (std::size_t i = 0; i < ids.size(); ++i) setParam(ids[i], values(0, static_cast<int>(i)));
}
void NoxGroup::setParams(const LOCA::ParameterVector& params) { params_ = params; fValid_ = jacobianValid_ = newtonValid_ = false; }
void NoxGroup::setParam(int id, double value) { params_.setValue(static_cast<unsigned>(id), value); fValid_ = jacobianValid_ = newtonValid_ = false; }
void NoxGroup::setParam(std::string id, double value) { params_.setValue(std::move(id), value); fValid_ = jacobianValid_ = newtonValid_ = false; }
double NoxGroup::getParam(int id) const { return params_.getValue(static_cast<unsigned>(id)); }
double NoxGroup::getParam(std::string id) const { return params_.getValue(std::move(id)); }
NOX::Abstract::Group::ReturnType NoxGroup::computeDfDpMulti(const std::vector<int>& ids, NOX::Abstract::MultiVector& dfdp, bool isValidF) {
  if (!isValidF && computeF() != NOX::Abstract::Group::Ok) return NOX::Abstract::Group::Failed;
  vectorRef(dfdp[0]).values() = f_.values();
  for (std::size_t i = 0; i < ids.size(); ++i) {
    const auto label = params_.getLabel(static_cast<unsigned>(ids[i])); State action; Parameters parameters; parameters.reynolds = getParam("Re");
    const auto operation = derivative_.parameterAction(x_.values(), parameters, spec_, label, action);
    if (!operation) return NOX::Abstract::Group::Failed;
    auto& target = vectorRef(dfdp[static_cast<int>(i + 1)]); target.values() = std::move(action);
  }
  return NOX::Abstract::Group::Ok;
}

} // namespace bif::trilinos
