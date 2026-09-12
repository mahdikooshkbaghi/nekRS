#pragma once

#include "bif/derivatives.hpp"
#include "bif/provider.hpp"
#include "linear_solver.hpp"

#include <LOCA_MultiContinuation_AbstractGroup.H>
#include <LOCA_Parameter_Vector.H>
#include <NOX_Abstract_Vector.H>
#include <mpi.h>
#include <memory>
#include <vector>

namespace bif::trilinos {

// A small NOX vector implementation.  It deliberately stores only the free
// analysis state; nekRS remains behind FlowMapProvider and is never exposed to
// NOX/LOCA clients.
class NoxVector final : public NOX::Abstract::Vector {
public:
  explicit NoxVector(std::size_t n = 0) : values_(n, 0.0) {}
  explicit NoxVector(State values) : values_(std::move(values)) {}
  NoxVector& operator=(const NoxVector&) = default;
  const State& values() const { return values_; }
  State& values() { return values_; }

  NOX::Abstract::Vector& init(double) override;
  NOX::Abstract::Vector& random(bool useSeed = false, int seed = 1) override;
  NOX::Abstract::Vector& abs(const NOX::Abstract::Vector&) override;
  NOX::Abstract::Vector& operator=(const NOX::Abstract::Vector&) override;
  NOX::Abstract::Vector& reciprocal(const NOX::Abstract::Vector&) override;
  NOX::Abstract::Vector& scale(double) override;
  NOX::Abstract::Vector& scale(const NOX::Abstract::Vector&) override;
  NOX::Abstract::Vector& update(double, const NOX::Abstract::Vector&, double = 0.0) override;
  NOX::Abstract::Vector& update(double, const NOX::Abstract::Vector&, double,
                                const NOX::Abstract::Vector&, double = 0.0) override;
  Teuchos::RCP<NOX::Abstract::Vector> clone(NOX::CopyType = NOX::DeepCopy) const override;
  double norm(NOX::Abstract::Vector::NormType = TwoNorm) const override;
  double norm(const NOX::Abstract::Vector& weights) const override;
  double innerProduct(const NOX::Abstract::Vector&) const override;
  NOX::size_type length() const override { return values_.size(); }
  void print(std::ostream&) const override;

private:
  State values_;
};

// Concrete LOCA group for F(q,Re)=q-Phi_T(q,Re).  Jacobian actions and
// parameter actions are matrix-free and use the production derivative object.
class NoxGroup final : public LOCA::MultiContinuation::AbstractGroup {
public:
  NoxGroup(FlowMapResidual& residual, MapDerivative& derivative,
           FlowMapProvider& provider, FlowMapSpec spec, ConstState initial,
           const LOCA::ParameterVector& params,
           LinearSolverOptions linearOptions = {});
  NoxGroup(const NoxGroup&) = default;
  ~NoxGroup() override = default;

  NOX::Abstract::Group& operator=(const NOX::Abstract::Group&) override;
  void copy(const NOX::Abstract::Group&) override;
  void setX(const NOX::Abstract::Vector&) override;
  void computeX(const NOX::Abstract::Group&, const NOX::Abstract::Vector&, double) override;
  ReturnType computeF() override;
  ReturnType computeJacobian() override;
  ReturnType computeNewton(Teuchos::ParameterList&) override;
  ReturnType applyJacobian(const NOX::Abstract::Vector&, NOX::Abstract::Vector&) const override;
  ReturnType applyJacobianInverse(Teuchos::ParameterList&, const NOX::Abstract::Vector&, NOX::Abstract::Vector&) const override;
  bool isF() const override { return fValid_; }
  bool isJacobian() const override { return jacobianValid_; }
  bool isNewton() const override { return newtonValid_; }
  bool isGradient() const override { return false; }
  const NOX::Abstract::Vector& getX() const override { return x_; }
  const NOX::Abstract::Vector& getF() const override { return f_; }
  double getNormF() const override { return normF_; }
  const NOX::Abstract::Vector& getGradient() const override { return f_; }
  const NOX::Abstract::Vector& getNewton() const override { return newton_; }
  Teuchos::RCP<const NOX::Abstract::Vector> getXPtr() const override;
  Teuchos::RCP<const NOX::Abstract::Vector> getFPtr() const override;
  Teuchos::RCP<const NOX::Abstract::Vector> getGradientPtr() const override;
  Teuchos::RCP<const NOX::Abstract::Vector> getNewtonPtr() const override;
  Teuchos::RCP<NOX::Abstract::Group> clone(NOX::CopyType = NOX::DeepCopy) const override;

  void setParamsMulti(const std::vector<int>&, const NOX::Abstract::MultiVector::DenseMatrix&) override;
  void setParams(const LOCA::ParameterVector&) override;
  void setParam(int, double) override;
  void setParam(std::string, double) override;
  const LOCA::ParameterVector& getParams() const override { return params_; }
  double getParam(int) const override;
  double getParam(std::string) const override;
  ReturnType computeDfDpMulti(const std::vector<int>&, NOX::Abstract::MultiVector&, bool) override;

private:
  ReturnType status(const OperationStatus&) const;
  bool unpack(const NOX::Abstract::Vector&, ConstState&) const;
  FlowMapResidual& residual_;
  MapDerivative& derivative_;
  FlowMapProvider& provider_;
  FlowMapSpec spec_;
  LOCA::ParameterVector params_;
  LinearSolverOptions linearOptions_;
  NoxVector x_, f_, newton_;
  double normF_{0.0};
  bool fValid_{false}, jacobianValid_{false}, newtonValid_{false};
};

} // namespace bif::trilinos
