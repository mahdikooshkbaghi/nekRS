#if !defined(nekrs_nrs_hpp_)
#define nekrs_nrs_hpp_

#include <mpi.h>
#include <functional>
#include <map>
#include <string>
#include <vector>
#include <cstdint>

/*
  Basic high-level API
*/

namespace nekrs
{
void setup(MPI_Comm commg_in,
           MPI_Comm comm_in,
           int buildOnly,
           int commSizeTarget,
           int ciMode,
           const std::map<std::string, std::map<std::string, std::string>> &parKeyValuePairs,
           std::string casename,
           std::string _backend,
           std::string _deviceID,
           int nSessions,
           int sessionID,
           int debug);
void udfExecuteStep(double time, int tstep, int checkpointStep);
void writeCheckpoint(double time);
int checkpointStep(double time, int tStep);
void checkpointStep(int val);
int finalize();
int runTimeStatFreq();
int printStepInfoFreq();
int updateFileCheckFreq();
void printRuntimeStatistics(int step);
double writeInterval(void);
std::tuple<double, double> dt(int tStep);
double startTime(void);
double endTime(void);
int numSteps(void);
void lastStep(int val);
int lastStep(double time, int tstep, double elapsedTime);
int writeControlRunTime(void);
int exitValue(void);
bool stepConverged(void);
void processUpdFile();
void printStepInfo(double time, int tstep, bool printStepInfo, bool printVerboseInfo);
void updateTimer(const std::string &key, double time);
void resetTimer(const std::string &key);

void initStep(double time, double dt, int tstep);
bool runStep(std::function<bool(int)> convergenceCheck, int corrector);
bool runStep(int corrector);
void finishStep();
bool stepConverged();
int timeStep();

// Optional analysis-state boundary used by extensions. The payload is opaque
// and must only be restored with the same mesh, discretization and rank layout.
struct AnalysisState {
  std::vector<double> payload;
};
struct AnalysisLayout {
  std::uint64_t globalSize = 0;
  std::uint64_t localSize = 0;
  std::uint64_t globalOffset = 0;
  int components = 0;
  std::vector<std::uint64_t> globalIds;
  std::vector<unsigned char> freeDofs;
  std::vector<double> metricWeights;
};
struct AnalysisOperationStatus {
  bool success = true;
  std::string message;
  explicit operator bool() const { return success; }
};
AnalysisOperationStatus captureAnalysisState(AnalysisState &state);
AnalysisOperationStatus restoreAnalysisState(const AnalysisState &state);
AnalysisOperationStatus analysisLayout(AnalysisLayout &layout);
AnalysisOperationStatus packAnalysisVelocity(std::vector<double> &state);
AnalysisOperationStatus unpackAnalysisVelocity(const std::vector<double> &state);
// Pack the autonomous velocity/pressure phase state used by the optional
// flow-map adapter.  Histories and work arrays are reconstructed
// deterministically from this phase state before each trial map evaluation.
AnalysisOperationStatus packAnalysisState(std::vector<double> &state);
AnalysisOperationStatus unpackAnalysisState(const std::vector<double> &state);
AnalysisOperationStatus setAnalysisReynolds(double reynolds);
double analysisTime();
double finalTimeStepSize(double time);
} // namespace nekrs

#endif
