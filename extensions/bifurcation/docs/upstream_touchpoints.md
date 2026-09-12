# nekRS touchpoints

The existing high-level API is `nekrs::setup`, `initStep`, one complete
`runStep(stage)`, `finishStep`, and `finalize` in `src/lib/nekrs.hpp`. The
adapter uses one `runStep(1)` per map timestep; the stage argument is not a
request to repeat the timestep.

The optional boundary adds an opaque immutable analysis snapshot and
unique-node layout/pack/scatter functions. The snapshot includes velocity,
pressure, explicit/advection histories, properties, relative velocity history,
BDF/EXT coefficients, clocks, and checkpoint flags. The phase-state adapter
packs free velocity DOFs plus pressure (with one pressure gauge DOF removed)
and reconstructs the other trial histories deterministically. UDF-owned state
still requires an embedding UDF that is pure with respect to the advertised
phase state.

The fluid change is limited to the analysis-state boundary and deterministic
phase-state history installation; ordinary solver algorithms are unchanged.
`CMakeLists.txt` adds only the disabled-by-default
`NEKRS_BUILD_BIFURCATION` option and subdirectory. Enabling it builds the
nekRS adapter without Trilinos; the CLI additionally requires compatible
NOX/LOCA/Belos/Anasazi packages.
