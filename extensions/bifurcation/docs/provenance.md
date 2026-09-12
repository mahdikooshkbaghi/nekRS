# Provenance

- nekRS working tree base observed by the implementation session: `master`
  (the handout records `995575be0ee6d7a88012c563f127c62673444eae`; verify with
  `git rev-parse HEAD` before publishing).
- Supplied Trilinos source tree: `Version.cmake` reports `17.3.0-dev`.
  Verify its full SHA and build prefix with `git -C Trilinos rev-parse HEAD`
  and the configured CMake cache; do not use an ABI-incompatible system MPI.
- The checked-in CMake project does not modify the supplied Trilinos tree.
  `scripts/build_cpu.sh` builds the dependency-free reference core. Set
  `BIF_ENABLE_TRILINOS=ON` only with a compatible installed/exported NOX
  package configuration.
- Validation performed in this tree: standalone Release CMake build with GCC
  16.2.1 and the Hopf normal-form smoke executable. No cylinder/MPI result is
  claimed by this repository state.
