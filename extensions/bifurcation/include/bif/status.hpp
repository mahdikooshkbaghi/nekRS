#pragma once

#include <string>
#include <utility>

namespace bif {

enum class StatusCode {
  success = 0,
  invalid_state,
  invalid_parameter,
  solver_failure,
  not_converged,
  not_supported,
  io_error,
  collective_failure,
  numerical_failure
};

struct OperationStatus {
  StatusCode code{StatusCode::success};
  double achieved_error{0.0};
  int iterations{0};
  std::string message;

  explicit operator bool() const { return code == StatusCode::success; }
  static OperationStatus ok() { return {}; }
  static OperationStatus error(StatusCode c, std::string text) {
    OperationStatus s; s.code = c; s.message = std::move(text); return s;
  }
};

const char *statusName(StatusCode code);

} // namespace bif
