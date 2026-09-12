#pragma once

#include "bif/provider.hpp"
#include "session.hpp"
#include <memory>

namespace bif::nekrs {

// Construct the in-process provider around the narrow analysis API exported by
// src/lib/nekrs.hpp. This translation unit is only linked into the optional
// nekRS-enabled target.
std::unique_ptr<Session> makeSession();
std::unique_ptr<bif::FlowMapProvider> makeFlowMapProvider();

} // namespace bif::nekrs
