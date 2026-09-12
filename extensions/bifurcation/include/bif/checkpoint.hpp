#pragma once

#include "bif/state.hpp"
#include "bif/task.hpp"
#include <string>
#include <vector>

namespace bif {

struct Checkpoint {
  int schemaVersion{1};
  std::string branchId{"branch-0"};
  int acceptedStep{0};
  double arclength{0.0};
  double nextStep{0.0};
  double parameterScale{1.0};
  std::string sourceSha;
  std::string meshHash;
  State state;
  Parameters parameters;
  FlowMapSpec map;
  StateLayout layout;
  std::vector<Eigenpair> trackedSpectrum;
  EventStatus event{EventStatus::none};
  double bracketLo{0.0};
  double bracketHi{0.0};
  std::string metadata;
};

OperationStatus writeCheckpoint(const std::string &path, const Checkpoint &checkpoint);
OperationStatus readCheckpoint(const std::string &path, Checkpoint &checkpoint);

} // namespace bif
