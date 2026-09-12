# Implementation status

## Complete in this change

- Optional `extensions/bifurcation` library with provider/state contracts.
- Transaction-aware centered map finite differences and `q-Phi_T` actions.
- Weighted layout operations, pair matching, bounded dense reference spectrum
  path for small solver-independent examples.
- Atomic schema-1 checkpoint manifest containing state, map, layout, parameter,
  branch and bracket metadata.
- Generic correction and safeguarded Hopf localization workflow.
- nekRS adapter boundary with explicit callback-based transactional session.
- Trilinos boundary with source API probe and Belos GMRES callback operator;
  the adapter target is opt-in and must use a compatible exported build.
- Hopf normal-form executable smoke example and declarative cylinder case
  inputs. The ordinary nekRS build remains unchanged when the option is off.

## Not claimed / remaining integration work

The nekRS boundary now supplies an opaque transactional snapshot and an
explicit velocity/pressure phase-state pack/scatter. Histories and work arrays
are deterministically reconstructed for each advertised BDF1/EXT1 map. The
nekRS-backed CLI and Belos/Anasazi translation units compile when enabled, but
NOX/LOCA continuation is not yet wired into the production driver; correction
currently uses the provider-independent matrix-free path with communicator-
aware reductions. The cylinder mesh, serial/2-rank acceptance matrix, verified
Hopf number, and DNS comparison are not represented as results. This is
intentional: no fabricated cylinder evidence is recorded.
