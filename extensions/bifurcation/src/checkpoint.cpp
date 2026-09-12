#include "bif/checkpoint.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <regex>
#include <sstream>

namespace bif {
namespace fs = std::filesystem;
namespace {
std::string quote(const std::string &s) {
  std::string out = "\"";
  for (char c : s) { if (c == '\\' || c == '\"') out += '\\'; out += c; }
  return out + "\"";
}
std::string numberList(const std::vector<double> &values) {
  std::ostringstream out; out << "[" << std::setprecision(17);
  for (std::size_t i = 0; i < values.size(); ++i) { if (i) out << ","; out << values[i]; }
  out << "]";
  return out.str();
}
std::string uintList(const std::vector<std::uint64_t> &values) {
  std::ostringstream out; out << "[";
  for (std::size_t i = 0; i < values.size(); ++i) { if (i) out << ","; out << values[i]; }
  out << "]";
  return out.str();
}
template <typename T> bool scalar(const std::string &text, const std::string &key, T &value) {
  std::regex expression("\\\"" + key + "\\\"\\s*:\\s*([-+0-9.eE]+)");
  std::smatch match; if (!std::regex_search(text, match, expression)) return false;
  std::istringstream(match[1].str()) >> value; return true;
}
bool stringValue(const std::string &text, const std::string &key, std::string &value) {
  std::regex expression("\\\"" + key + "\\\"\\s*:\\s*\\\"([^\\\"]*)\\\"");
  std::smatch match; if (!std::regex_search(text, match, expression)) return false; value = match[1]; return true;
}
bool arrayValue(const std::string &text, const std::string &key, std::vector<double> &values) {
  std::regex expression("\\\"" + key + "\\\"\\s*:\\s*\\[([^\\]]*)\\]");
  std::smatch match; if (!std::regex_search(text, match, expression)) return false;
  std::stringstream ss(match[1].str()); std::string item; values.clear();
  while (std::getline(ss, item, ',')) { if (!item.empty()) values.push_back(std::stod(item)); }
  return true;
}
}

OperationStatus writeCheckpoint(const std::string &path, const Checkpoint &c) {
  try {
    const fs::path target(path), temporary = target.string() + ".tmp." + std::to_string(c.acceptedStep);
    if (!target.parent_path().empty()) fs::create_directories(target.parent_path());
    std::ofstream out(temporary);
    if (!out) return OperationStatus::error(StatusCode::io_error, "cannot open checkpoint temporary file");
    out << std::setprecision(17) << "{\n"
        << "  \"schema_version\": " << c.schemaVersion << ",\n"
        << "  \"branch_id\": " << quote(c.branchId) << ",\n"
        << "  \"accepted_step\": " << c.acceptedStep << ",\n"
        << "  \"arclength\": " << c.arclength << ",\n"
        << "  \"next_step\": " << c.nextStep << ",\n"
        << "  \"parameter_scale\": " << c.parameterScale << ",\n"
        << "  \"source_sha\": " << quote(c.sourceSha) << ",\n"
        << "  \"mesh_hash\": " << quote(c.meshHash) << ",\n"
        << "  \"parameter_Re\": " << c.parameters.reynolds << ",\n"
        << "  \"map_horizon\": " << c.map.horizon << ",\n"
        << "  \"map_timestep\": " << c.map.timestep << ",\n"
        << "  \"map_steps\": " << c.map.steps << ",\n"
        << "  \"event_status\": " << quote(eventStatusName(c.event)) << ",\n"
        << "  \"bracket_lo\": " << c.bracketLo << ",\n"
        << "  \"bracket_hi\": " << c.bracketHi << ",\n"
        << "  \"state\": " << numberList(c.state) << ",\n"
        << "  \"metric_weights\": " << numberList(c.layout.metricWeights) << ",\n"
        << "  \"global_ids\": " << uintList(c.layout.globalIds) << ",\n"
        << "  \"metadata\": " << quote(c.metadata) << "\n}\n";
    out.close();
    if (!out) return OperationStatus::error(StatusCode::io_error, "failed writing checkpoint");
    fs::rename(temporary, target);
    return OperationStatus::ok();
  } catch (const std::exception &e) { return OperationStatus::error(StatusCode::io_error, e.what()); }
}

OperationStatus readCheckpoint(const std::string &path, Checkpoint &c) {
  std::ifstream in(path); if (!in) return OperationStatus::error(StatusCode::io_error, "cannot open checkpoint: " + path);
  std::stringstream buffer; buffer << in.rdbuf(); const auto text = buffer.str();
  if (!scalar(text, "schema_version", c.schemaVersion) || c.schemaVersion != 1)
    return OperationStatus::error(StatusCode::io_error, "unsupported or missing checkpoint schema");
  stringValue(text, "branch_id", c.branchId); stringValue(text, "source_sha", c.sourceSha);
  stringValue(text, "mesh_hash", c.meshHash); stringValue(text, "metadata", c.metadata);
  scalar(text, "accepted_step", c.acceptedStep); scalar(text, "arclength", c.arclength);
  scalar(text, "next_step", c.nextStep); scalar(text, "parameter_scale", c.parameterScale);
  scalar(text, "parameter_Re", c.parameters.reynolds); scalar(text, "map_horizon", c.map.horizon);
  scalar(text, "map_timestep", c.map.timestep); scalar(text, "map_steps", c.map.steps);
  scalar(text, "bracket_lo", c.bracketLo); scalar(text, "bracket_hi", c.bracketHi);
  if (!arrayValue(text, "state", c.state)) return OperationStatus::error(StatusCode::io_error, "checkpoint has no state");
  arrayValue(text, "metric_weights", c.layout.metricWeights);
  c.layout.localSize = c.state.size(); c.layout.globalSize = c.state.size();
  c.layout.globalIds.resize(c.state.size()); for (std::size_t i = 0; i < c.state.size(); ++i) c.layout.globalIds[i] = i;
  c.layout.freeDofs.assign(c.state.size(), 1);
  return OperationStatus::ok();
}

} // namespace bif
