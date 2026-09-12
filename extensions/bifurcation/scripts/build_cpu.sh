#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD=${BUILD_DIR:-"$ROOT/build-cpu"}
TRILINOS_ARGS=(-DBIF_ENABLE_TRILINOS=OFF)
if [[ -n "${TRILINOS_DIR:-}" ]]; then
  TRILINOS_ARGS=(-DBIF_ENABLE_TRILINOS=ON -DTrilinos_DIR="$TRILINOS_DIR")
fi
cmake -S "$ROOT" -B "$BUILD" -DCMAKE_BUILD_TYPE=Release -DBIF_BUILD_APP=ON "${TRILINOS_ARGS[@]}"
cmake --build "$BUILD" --parallel "${BUILD_JOBS:-2}"
printf 'built %s\n' "$BUILD/nekrs-bif"
