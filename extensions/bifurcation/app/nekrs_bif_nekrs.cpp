#include "bif/checkpoint.hpp"
#include "bif/derivatives.hpp"
#include "bif/workflow.hpp"
#include "map_eigensolver.hpp"
#include "loca_runner.hpp"
#include "api.hpp"
#include "nekrs.hpp"
#include "inipp.hpp"

#include <mpi.h>
#include <Tpetra_Core.hpp>
#include <unistd.h>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace fs = std::filesystem;

namespace {
std::string option(int argc, char** argv, const std::string& name, const std::string& fallback = {}) {
  for (int i = 1; i + 1 < argc; ++i) if (name == argv[i]) return argv[i + 1];
  return fallback;
}

std::map<std::string, std::string> readConfig(const std::string& path) {
  std::map<std::string, std::string> values; std::ifstream input(path); std::string line, section;
  while (std::getline(input, line)) {
    const auto comment = line.find_first_of("#;"); if (comment != std::string::npos) line.resize(comment);
    auto trim = [](std::string value) { const auto first = value.find_first_not_of(" \t\r\n"); const auto last = value.find_last_not_of(" \t\r\n"); return first == std::string::npos ? std::string{} : value.substr(first, last - first + 1); };
    line = trim(line); if (line.empty()) continue;
    if (line.front() == '[' && line.back() == ']') { section = line.substr(1, line.size() - 2); std::transform(section.begin(), section.end(), section.begin(), ::tolower); continue; }
    const auto equal = line.find('='); if (equal == std::string::npos) continue;
    auto key = trim(line.substr(0, equal)); auto value = trim(line.substr(equal + 1)); std::transform(key.begin(), key.end(), key.begin(), ::tolower); values[section + "." + key] = value;
  }
  return values;
}

template <typename T> T configValue(const std::map<std::string, std::string>& values, const std::string& key, T fallback) {
  const auto it = values.find(key); if (it == values.end()) return fallback; std::istringstream stream(it->second); T value = fallback; stream >> value; return value;
}

std::string rankCheckpoint(const std::string& path, int rank) {
  const fs::path value(path); const auto candidate = value.parent_path() / (value.stem().string() + ".rank" + std::to_string(rank) + value.extension().string());
  return fs::exists(candidate) ? candidate.string() : path;
}

bool choosePair(const bif::SpectrumResult& spectrum, bif::Eigenpair& selected) {
  bool found = false;
  for (const auto& pair : spectrum.eigenpairs) if (pair.value.imag() > 1.0e-8 && (!found || std::abs(pair.value) > std::abs(selected.value))) { selected = pair; found = true; }
  return found;
}

void usage() {
  std::cout << "nekrs-bif --case cylinder.par --config bifurcation.ini --task seed|continue|resume|stability \\\n --output RUN [--seed CHECKPOINT] [--restart CHECKPOINT] [--stop-after-accepted N] [--driver loca]\n";
}
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  Tpetra::ScopeGuard tpetraGuard(&argc, &argv, MPI_COMM_WORLD);
  int rank = 0, size = 1; MPI_Comm_rank(MPI_COMM_WORLD, &rank); MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (argc < 2 || option(argc, argv, "--help") == "1") { if (rank == 0) usage(); MPI_Finalize(); return argc < 2 ? 1 : 0; }
  const fs::path casePath = fs::absolute(option(argc, argv, "--case"));
  if (casePath.empty() || !fs::exists(casePath)) { if (rank == 0) std::cerr << "case file not found\n"; MPI_Finalize(); return 2; }
  const auto configPath = fs::absolute(option(argc, argv, "--config", (casePath.parent_path() / "bifurcation.ini").string()));
  const std::string task = option(argc, argv, "--task", "seed");
  const fs::path output = fs::absolute(option(argc, argv, "--output", "run-bifurcation"));
  fs::create_directories(output);
  const auto config = readConfig(configPath.string());

  const auto oldDirectory = fs::current_path(); fs::current_path(casePath.parent_path());
  std::ifstream parFile(casePath.filename()); std::stringstream parText; parText << parFile.rdbuf();
  inipp::Ini ini; ini.parse(parText); ini.interpolate();
  const std::string casename = casePath.stem().string();
  nekrs::setup(MPI_COMM_WORLD, MPI_COMM_WORLD, 0, size, 0, ini.sections, casename, "CPU", "LOCAL-RANK", 1, 0, 0);
  fs::current_path(oldDirectory);

  auto provider = bif::nekrs::makeFlowMapProvider();
  if (!provider) { if (rank == 0) std::cerr << "failed to construct nekRS analysis provider\n"; nekrs::finalize(); MPI_Finalize(); return 3; }
  const auto layout = provider->layout();
  bif::FlowMapSpec map; map.horizon = configValue(config, "map.horizon", 2.0); map.timestep = configValue(config, "map.dt", 0.02); map.steps = configValue(config, "map.steps", static_cast<int>(std::llround(map.horizon / map.timestep))); map.integrator = "BDF1-EXT1";
  if (!map.valid()) {
    if (rank == 0) std::cerr << "map.horizon must equal map.steps * map.dt within tolerance\n";
    nekrs::finalize(); MPI_Finalize(); return 4;
  }
  bif::Parameters parameters; parameters.reynolds = configValue(config, "continuation.first", 35.0);
  bif::State state; auto packed = nekrs::packAnalysisState(state); if (!packed) { if (rank == 0) std::cerr << packed.message << "\n"; nekrs::finalize(); MPI_Finalize(); return 4; }
  bif::Checkpoint restart; int accepted = 0;
  if (task == "resume" || task == "stability" || !option(argc, argv, "--seed").empty()) {
    const auto supplied = option(argc, argv, task == "resume" || task == "stability" ? "--restart" : "--seed");
    const auto path = rankCheckpoint(supplied, rank); auto status = bif::readCheckpoint(path, restart);
    if (!status) { if (rank == 0) std::cerr << status.message << "\n"; nekrs::finalize(); MPI_Finalize(); return 5; }
    state = restart.state; parameters = restart.parameters; accepted = restart.acceptedStep;
  }

  bif::FiniteDifferenceMapDerivative derivative(*provider, {configValue(config, "newton.finite_difference_step", 1.0e-5), 1.0e-5});
  if (option(argc, argv, "--driver") == "loca") {
    bif::FlowMapResidual residual(*provider, derivative);
    LOCA::ParameterVector locaParameters; locaParameters.addParameter("Re", parameters.reynolds);
    bif::trilinos::LocaOptions locaOptions;
    locaOptions.initialStep = configValue(config, "continuation.step", locaOptions.initialStep);
    locaOptions.maximumStep = configValue(config, "continuation.max_step", locaOptions.maximumStep);
    locaOptions.maximumSteps = configValue(config, "continuation.max_steps", locaOptions.maximumSteps);
    const auto locaStatus = bif::trilinos::runLocaContinuation(*provider, residual, derivative, map, state, locaParameters, locaOptions);
    if (rank == 0) {
      if (locaStatus) std::cout << "LOCA continuation completed\n";
      else std::cerr << "LOCA continuation failed: " << locaStatus.message << "\n";
    }
    nekrs::finalize(); MPI_Finalize(); return locaStatus ? 0 : 8;
  }
  bif::trilinos::MapEigensolver eigensolver({configValue(config, "eigensolver.nev", 6), configValue(config, "eigensolver.subspace", 50), configValue(config, "eigensolver.max_restarts", 100), configValue(config, "eigensolver.residual_tolerance", 1.0e-5)});
  const double first = task == "resume" ? parameters.reynolds + configValue(config, "continuation.step", 2.0) : parameters.reynolds;
  const double last = task == "stability" ? parameters.reynolds : configValue(config, "continuation.last", 60.0);
  const double step = configValue(config, "continuation.step", 2.0);
  const int stopAfter = std::stoi(option(argc, argv, "--stop-after-accepted", "-1"));
  std::ofstream csv; if (rank == 0) { csv.open(output / "branch.csv", task == "resume" ? std::ios::app : std::ios::out); if (task != "resume") csv << "branch,step,Re,arclength,residual,divergence,mu_real,mu_imag,sigma,omega,St,eigen_residual,event_status,newton_iterations,linear_iterations,map_evaluations,elapsed_seconds,mpi_ranks\n"; }
  bif::State previousState = state; bif::Parameters previousParameters = parameters; double previousSigma = 0.0;
  const bool onePoint = task == "seed" || task == "stability";
  for (double reynolds = first; onePoint || reynolds <= last + 1.0e-12; reynolds += step) {
    parameters.reynolds = reynolds; const auto started = std::chrono::steady_clock::now();
    const auto corrected = bif::correctEquilibrium(*provider, state, parameters, map, {configValue(config, "newton.max_iterations", 12), configValue(config, "newton.tolerance", 1.0e-7), configValue(config, "newton.finite_difference_step", 1.0e-5)});
    if (!corrected.status) { if (rank == 0) std::cerr << "correction failed: " << corrected.status.message << "\n"; nekrs::finalize(); MPI_Finalize(); return 6; }
    state = corrected.state; auto installed = provider->install(state, parameters); if (!installed) { if (rank == 0) std::cerr << installed.message << "\n"; nekrs::finalize(); MPI_Finalize(); return 7; }
    const auto spectrum = eigensolver.solve(layout, state, parameters, map, derivative, configValue(config, "eigensolver.nev", 6));
    if (!spectrum.status) { if (rank == 0) std::cerr << "spectrum failed: " << spectrum.status.message << "\n"; nekrs::finalize(); MPI_Finalize(); return 8; }
    bif::Eigenpair pair; const bool hasPair = choosePair(spectrum, pair); const double sigma = hasPair ? std::log(std::abs(pair.value)) / map.horizon : NAN; const double omega = hasPair ? std::atan2(pair.value.imag(), pair.value.real()) / map.horizon : NAN;
    bif::EventStatus event = (accepted > 0 && hasPair && previousSigma * sigma < 0.0) ? bif::EventStatus::hopfCandidate : bif::EventStatus::none;
    const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
    if (rank == 0) csv << "0," << accepted << "," << reynolds << "," << accepted * step << "," << corrected.residual << ",nan," << (hasPair ? pair.value.real() : NAN) << "," << (hasPair ? pair.value.imag() : NAN) << "," << sigma << "," << omega << "," << (hasPair ? std::abs(omega) / (2.0 * M_PI) : NAN) << "," << (hasPair ? pair.residual : NAN) << "," << bif::eventStatusName(event) << "," << corrected.iterations << ",0," << (corrected.mapEvaluations + spectrum.mapApplications) << "," << elapsed << "," << size << "\n";
    restart.state = state; restart.layout = layout; restart.parameters = parameters; restart.map = map; restart.acceptedStep = accepted; restart.nextStep = step; restart.event = event; restart.bracketLo = previousParameters.reynolds; restart.bracketHi = reynolds; restart.metadata = "time-discrete nekRS flow-map continuation";
    auto writeStatus = bif::writeCheckpoint((output / ("checkpoint.rank" + std::to_string(rank) + ".json")).string(), restart); if (!writeStatus) { if (rank == 0) std::cerr << writeStatus.message << "\n"; nekrs::finalize(); MPI_Finalize(); return 9; }
    if (rank == 0) bif::writeCheckpoint((output / "checkpoint.json").string(), restart);
    MPI_Barrier(MPI_COMM_WORLD);
    if (event == bif::EventStatus::hopfCandidate) { fs::create_directories(output / "events/HB1"); bif::writeCheckpoint((output / ("events/HB1/manifest.rank" + std::to_string(rank) + ".json")).string(), restart); if (rank == 0) bif::writeCheckpoint((output / "events/HB1/manifest.json").string(), restart); }
    previousState = state; previousParameters = parameters; previousSigma = sigma; ++accepted;
    if (onePoint || (stopAfter > 0 && accepted >= stopAfter)) break;
  }
  if (rank == 0) csv.close(); nekrs::finalize(); MPI_Finalize(); return 0;
}
