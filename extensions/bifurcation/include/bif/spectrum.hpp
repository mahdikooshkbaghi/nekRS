#pragma once

#include "bif/derivatives.hpp"
#include <complex>
#include <vector>

namespace bif {

struct Eigenpair {
  std::complex<double> value{0.0, 0.0};
  State vector;
  double residual{0.0};
  int pairId{-1};
};

struct SpectrumResult {
  OperationStatus status;
  std::vector<Eigenpair> eigenpairs;
  int mapApplications{0};
};

class MapEigensolver {
public:
  virtual ~MapEigensolver() = default;
  virtual SpectrumResult solve(const StateLayout& layout, ConstState q,
                               const Parameters &p, const FlowMapSpec &spec,
                               MapDerivative &derivative, int nev) = 0;
};

// A small dense reference implementation used by the solver-independent examples.
// Production PDE runs use the Anasazi adapter and never assemble this matrix.
class DenseMapEigensolver final : public MapEigensolver {
public:
  SpectrumResult solve(const StateLayout& layout, ConstState q,
                       const Parameters &p, const FlowMapSpec &spec,
                       MapDerivative &derivative, int nev) override;
};

struct PairMatch {
  int pairId{-1};
  std::complex<double> value{0.0, 0.0};
  double overlap{0.0};
  bool matched{false};
};

PairMatch matchPair(const SpectrumResult &current,
                    std::complex<double> previous,
                    ConstState previousVector,
                    const StateLayout &layout);

} // namespace bif
