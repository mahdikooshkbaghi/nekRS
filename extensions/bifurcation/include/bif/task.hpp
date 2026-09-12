#pragma once

#include "bif/spectrum.hpp"
#include <string>
#include <vector>

namespace bif {

enum class TaskKind { seed, continueBranch, resume, stability, dns };

enum class EventStatus { none, bracketed, hopfCandidate, hopf, ambiguous, rejected };

struct Task {
  TaskKind kind{TaskKind::seed};
  std::string caseFile;
  std::string configFile;
  std::string inputCheckpoint;
  std::string outputDirectory{"run"};
  std::string parameter{"Re"};
  double firstParameter{35.0};
  double lastParameter{60.0};
  double parameterStep{2.0};
  int stopAfterAccepted{-1};
  int eigenpairs{6};
};

struct BranchPoint {
  int branch{0};
  int step{0};
  double parameter{0.0};
  double arclength{0.0};
  double residual{0.0};
  double divergence{0.0};
  std::complex<double> multiplier{0.0, 0.0};
  double sigma{0.0};
  double omega{0.0};
  double Strouhal{0.0};
  double eigenResidual{0.0};
  EventStatus event{EventStatus::none};
  int newtonIterations{0};
  int linearIterations{0};
  int mapEvaluations{0};
  double elapsedSeconds{0.0};
  int mpiRanks{1};
};

const char *eventStatusName(EventStatus status);
const char *taskKindName(TaskKind kind);

} // namespace bif
