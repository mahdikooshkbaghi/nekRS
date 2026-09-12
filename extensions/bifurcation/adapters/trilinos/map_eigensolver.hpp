#pragma once

#include "bif/spectrum.hpp"

namespace bif::trilinos {

struct AnasaziOptions {
  int eigenpairs{6};
  int subspace{50};
  int maximumRestarts{100};
  double residualTolerance{1.0e-5};
};

// Anasazi owns the distributed Krylov basis. The callback is a map derivative,
// so this adapter never diagonalizes J_F or assembles a PDE matrix.
class MapEigensolver final : public bif::MapEigensolver {
public:
  explicit MapEigensolver(AnasaziOptions options = {}) : options_(options) {}
  SpectrumResult solve(const StateLayout& layout, ConstState q,
                       const Parameters &p, const FlowMapSpec &spec,
                       bif::MapDerivative &derivative, int nev) override;
private:
  AnasaziOptions options_;
};

} // namespace bif::trilinos
