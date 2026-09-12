#!/usr/bin/env bash
set -euo pipefail

# This script is intentionally fail-fast. A cylinder run requires a nekRS build
# with -DNEKRS_BUILD_BIFURCATION=ON and an application-specific state callback.
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
NP=${NP:-1}
BUILD=${BUILD_DIR:-"$ROOT/build-cpu"}
OUT=${OUT_DIR:-"$ROOT/run-cylinder-np${NP}"}
"$ROOT/scripts/build_cpu.sh"
BIN="$BUILD/nekrs-bif"
if [[ ! -x "$BIN" ]]; then echo "missing $BIN" >&2; exit 1; fi
mkdir -p "$OUT"
"$BIN" --case "$ROOT/examples/cylinder_hopf/cylinder.par" \
  --config "$ROOT/examples/cylinder_hopf/bifurcation.ini" --task seed --output "$OUT/seed" 2>&1 | tee "$OUT/seed.log"
