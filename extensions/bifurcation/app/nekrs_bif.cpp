#include "bif/checkpoint.hpp"
#include "bif/derivatives.hpp"
#include "bif/spectrum.hpp"
#include "bif/workflow.hpp"
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <string>

namespace {
class NormalFormSnapshot final : public bif::ProviderSnapshot {
public:
  explicit NormalFormSnapshot(bif::State value) : state(std::move(value)) {}
  bif::State state;
};

class HopfNormalForm final : public bif::FlowMapProvider {
public:
  bif::StateLayout layout() const override {
    bif::StateLayout l; l.globalSize = l.localSize = 2; l.globalIds = {0, 1};
    l.freeDofs = {1, 1}; l.metricWeights = {1.0, 1.0}; l.fieldDescription = "x,y"; return l;
  }
  bif::OperationStatus evaluate(bif::ConstState q, const bif::Parameters &p,
                                const bif::FlowMapSpec &spec, bif::State &out) override {
    if (q.size() != 2) return bif::OperationStatus::error(bif::StatusCode::invalid_state, "normal form requires two entries");
    const double a = p.get("Re"), w = 1.0, dt = spec.timestep;
    out = q;
    auto rhs = [a, w](const bif::State &z) {
      const double r2 = z[0] * z[0] + z[1] * z[1];
      return bif::State{a * z[0] - w * z[1] - r2 * z[0], w * z[0] + a * z[1] - r2 * z[1]};
    };
    for (int step = 0; step < spec.steps; ++step) {
      const auto k1 = rhs(out); bif::State z2{out[0] + 0.5 * dt * k1[0], out[1] + 0.5 * dt * k1[1]};
      const auto k2 = rhs(z2); bif::State z3{out[0] + 0.5 * dt * k2[0], out[1] + 0.5 * dt * k2[1]};
      const auto k3 = rhs(z3); bif::State z4{out[0] + dt * k3[0], out[1] + dt * k3[1]}; const auto k4 = rhs(z4);
      out[0] += dt * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) / 6.0;
      out[1] += dt * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6.0;
      if (!std::isfinite(out[0]) || !std::isfinite(out[1])) return bif::OperationStatus::error(bif::StatusCode::solver_failure, "normal-form integration overflow");
    }
    return bif::OperationStatus::ok();
  }
  std::unique_ptr<bif::ProviderSnapshot> capture() override { return std::make_unique<NormalFormSnapshot>(state_); }
  bif::OperationStatus restore(const bif::ProviderSnapshot &snapshot) override {
    const auto *saved = dynamic_cast<const NormalFormSnapshot *>(&snapshot); if (!saved) return bif::OperationStatus::error(bif::StatusCode::invalid_state, "foreign snapshot"); state_ = saved->state; return bif::OperationStatus::ok();
  }
  bif::OperationStatus install(bif::ConstState q, const bif::Parameters &) override { state_ = q; return bif::OperationStatus::ok(); }
private:
  bif::State state_{0.0, 0.0};
};

std::string option(int argc, char **argv, const std::string &name, const std::string &fallback = {}) {
  for (int i = 1; i + 1 < argc; ++i) if (argv[i] == name) return argv[i + 1]; return fallback;
}
void usage() {
  std::cout << "nekrs-bif --example hopf_normal_form [--task seed|continue|resume|stability] [--output DIR]\n"
               "         --case cylinder.par is accepted only by a nekRS-enabled build\n";
}
}

int main(int argc, char **argv) {
  if (argc == 1 || option(argc, argv, "--help") == "1") { usage(); return argc == 1 ? 1 : 0; }
  const auto example = option(argc, argv, "--example");
  if (example != "hopf_normal_form") {
    if (!option(argc, argv, "--case").empty()) {
      std::cerr << "cylinder flow requires the optional nekRS adapter; this executable was built without NEKRS support\n";
      return 2;
    }
    usage(); return 2;
  }
  const std::string task = option(argc, argv, "--task", "seed"), output = option(argc, argv, "--output", "run-normal-form");
  std::filesystem::create_directories(output);
  HopfNormalForm provider; bif::FlowMapSpec spec; spec.horizon = 2.0; spec.timestep = 0.01; spec.steps = 200;
  bif::Parameters p; p.reynolds = -0.1; bif::State seed{0.0, 0.0};
  bif::Checkpoint checkpoint;
  if (task == "resume" || task == "stability") {
    const auto input = option(argc, argv, "--restart", output + "/checkpoint.json"); auto status = bif::readCheckpoint(input, checkpoint);
    if (!status) { std::cerr << status.message << "\n"; return 3; } seed = checkpoint.state; p = checkpoint.parameters;
  }
  bif::FiniteDifferenceMapDerivative derivative(provider); bif::DenseMapEigensolver eigensolver;
  std::ofstream csv(output + "/branch.csv"); csv << "branch,step,Re,residual,mu_real,mu_imag,sigma,omega,St,event_status\n";
  double first = std::stod(option(argc, argv, "--first", "-0.1"));
  double last = std::stod(option(argc, argv, "--last", "0.1"));
  const double increment = std::stod(option(argc, argv, "--step", "0.02"));
  const int stopAfter = std::stoi(option(argc, argv, "--stop-after-accepted", "-1"));
  const bool branchTask = task == "continue" || task == "resume";
  if (task == "stability") first = last = p.reynolds;
  if (task == "resume") first = p.reynolds + increment;
  double previousSigma = 0.0; bif::State previousState = seed; bif::Parameters previousParameters = p;
  int accepted = task == "resume" ? checkpoint.acceptedStep : 0;
  const double direction = increment >= 0.0 ? 1.0 : -1.0;
  for (double parameter = first; direction * (parameter - last) <= 1.0e-12 || !branchTask; parameter += increment) {
    p.reynolds = parameter;
    auto corrected = bif::correctEquilibrium(provider, seed, p, spec);
    if (!corrected.status) { std::cerr << "correction failed: " << corrected.status.message << "\n"; return 4; }
    const auto spectrum = eigensolver.solve(provider.layout(), corrected.state, p, spec, derivative, 2);
    if (!spectrum.status) { std::cerr << "spectrum failed: " << spectrum.status.message << "\n"; return 5; }
    double sigma = 0.0, omega = 0.0, muReal = 0.0, muImag = 0.0; bif::EventStatus event = bif::EventStatus::none;
    for (const auto &pair : spectrum.eigenpairs) if (pair.value.imag() > 1.0e-8) { muReal = pair.value.real(); muImag = pair.value.imag(); sigma = std::log(std::abs(pair.value)) / spec.horizon; omega = std::atan2(pair.value.imag(), pair.value.real()) / spec.horizon; }
    if (accepted > 0 && sigma * previousSigma < 0.0) event = bif::EventStatus::hopfCandidate;
    csv << "0," << accepted << "," << parameter << "," << corrected.residual << "," << muReal << "," << muImag << "," << sigma << "," << omega << "," << std::abs(omega) / (2 * M_PI) << "," << bif::eventStatusName(event) << "\n";
    checkpoint.state = corrected.state; checkpoint.parameters = p; checkpoint.map = spec; checkpoint.layout = provider.layout(); checkpoint.acceptedStep = accepted; checkpoint.bracketLo = previousParameters.reynolds; checkpoint.bracketHi = p.reynolds; checkpoint.event = event; checkpoint.metadata = "hopf normal form; map eigenvalues are not Floquet multipliers";
    auto saved = bif::writeCheckpoint(output + "/checkpoint.json", checkpoint); if (!saved) { std::cerr << saved.message << "\n"; return 6; }
    if (event == bif::EventStatus::hopfCandidate) { std::filesystem::create_directories(output + "/events/HB1"); bif::writeCheckpoint(output + "/events/HB1/manifest.json", checkpoint); }
    previousSigma = sigma; previousState = seed = corrected.state; previousParameters = p; ++accepted;
    if (!branchTask || (stopAfter > 0 && accepted >= stopAfter)) break;
  }
  std::cout << "completed " << task << " for the solver-independent Hopf normal form; checkpoint=" << output << "/checkpoint.json\n";
  return 0;
}
