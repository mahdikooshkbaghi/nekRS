# nekRS bifurcation extension

This is an optional, separately built continuation boundary. The core library
has no nekRS, OCCA, Trilinos, PETSc, or NOX headers. Build its verified
reference path with:

```bash
./scripts/build_cpu.sh
./build-cpu/nekrs-bif --example hopf_normal_form --task seed --output run
```

The normal-form run exercises the provider contract, map residual,
finite-difference derivative, small dense spectrum, and checkpoint writer. Its
multipliers are equilibrium-map multipliers, not periodic-orbit Floquet
multipliers.

For an installed compatible Trilinos build, configure separately with
`-DBIF_ENABLE_TRILINOS=ON -DTrilinos_DIR=/path/to/prefix` and inspect the
exported NOX package. `NOX::all_libs` supplies the configured Teuchos/Tpetra,
Belos, Anasazi and LOCA closure. The supplied uninstalled development build
must be configured/installed with paths valid in the current checkout; stale
absolute paths are rejected by CMake.

The nekRS adapter is deliberately callback-based. The upstream boundary now
provides an immutable snapshot plus a unique-node velocity/pressure phase-state
pack/scatter; multistep histories and work arrays are reconstructed
deterministically for each BDF1/EXT1 trial. The cylinder mesh and full
NOX/LOCA driver remain integration work. Therefore this change does **not**
claim a verified cylinder continuation, MPI result, DNS corroboration, or
literature-scale critical Reynolds number. See
`docs/implementation_status.md` and `docs/upstream_touchpoints.md`.
