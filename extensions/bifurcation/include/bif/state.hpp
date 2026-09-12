#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace bif {

using Scalar = double;
using State = std::vector<Scalar>;
using ConstState = const State &;

struct Communicator {
  int rank{0};
  int size{1};
  // Optional collective reduction supplied by an adapter. The core contract
  // does not depend on MPI (or any other communication library).
  std::function<Scalar(Scalar)> globalSum;
};

struct StateLayout {
  std::size_t globalSize{0};
  std::size_t localSize{0};
  std::size_t globalOffset{0};
  std::vector<std::uint64_t> globalIds;
  std::vector<unsigned char> freeDofs;
  std::vector<Scalar> metricWeights;
  std::string fieldDescription;
  Communicator communicator;
  bool hostResident{true};

  bool valid() const {
    return globalSize != 0 && localSize == globalIds.size() &&
           localSize == freeDofs.size() && localSize == metricWeights.size();
  }

  Scalar norm(ConstState x) const;
  Scalar dot(ConstState x, ConstState y) const;
};

struct Parameters {
  double reynolds{0.0};
  std::vector<double> values;
  std::vector<std::string> names;

  double get(const std::string &name, double fallback = 0.0) const;
  void set(const std::string &name, double value);
};

struct FlowMapSpec {
  double horizon{2.0};
  double timestep{0.02};
  int steps{100};
  std::string integrator{"BDF1-EXT1"};
  std::string historyPolicy{"provider-defined"};
  bool suppressOutput{true};

  bool valid() const {
    if (!(horizon > 0.0) || !(timestep > 0.0) || steps <= 0) return false;
    const double integrated = timestep * static_cast<double>(steps);
    return std::abs(integrated - horizon) <=
           1.0e-10 * std::max(1.0, std::abs(horizon));
  }
};

} // namespace bif
