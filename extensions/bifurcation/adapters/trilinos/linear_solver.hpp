#pragma once

#include "bif/state.hpp"
#include "bif/status.hpp"
#include <functional>

namespace bif::trilinos {

using Apply = std::function<OperationStatus(ConstState, State &)>;

struct LinearSolverOptions {
  int maximumIterations{400};
  int restart{40};
  double tolerance{1.0e-8};
  bool flexible{false};
};

struct LinearSolveResult {
  OperationStatus status;
  int iterations{0};
  double residual{0.0};
};

// Belos-backed GMRES. The PDE operator remains an action callback: no global
// sparse matrix is assembled and no identity operation is used as a solver.
LinearSolveResult solveGMRES(const StateLayout &layout, const Apply &apply,
                             ConstState rhs, State &solution,
                             const LinearSolverOptions &options = {});

} // namespace bif::trilinos
