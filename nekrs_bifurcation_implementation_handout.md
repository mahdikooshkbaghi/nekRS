# nekRS bifurcation implementation handout

Prepared 2026-09-11. This is an implementation specification, not a report of a
completed implementation or verified cylinder run.

## 1. Mission and non-negotiable boundary

Start from a clean checkout of the upstream nekRS default branch. Implement an
optional, separately built C++ continuation/bifurcation library using
NOX/LOCA, Belos, and Anasazi from the supplied Trilinos source. Finish with a
small, reproducible, one-parameter Hopf bifurcation of two-dimensional flow
past a stationary circular cylinder, running on both one and two MPI ranks.

The user's highest priority is a clear boundary and the minimum possible
changes to existing nekRS modules. Treat that as an acceptance requirement:

- Put all continuation, event detection, eigensolver, and checkpoint workflow
  code in a separate library under `extensions/bifurcation/`.
- Put every dependency on nekRS internals in its `adapters/nekrs/` directory.
- Keep Trilinos types in `adapters/trilinos/`; the mathematical provider
  contract must not require Trilinos, PETSc, OCCA, or a nekRS concrete type.
- Existing fluid, scalar, elliptic, mesh, and time-integration numerical
  algorithms must remain unchanged for this milestone.
- First use existing setup/step APIs. If state restoration cannot be correct
  through those APIs, add a narrow provider-state capture/restore facility.
  A small extraction of driver setup is also permissible if necessary to
  avoid duplicating it. Explain each existing-file edit individually.
- Do not add NOX/LOCA includes, bifurcation configuration, or event logic to
  fluid/scalar solvers, mesh code, UDF physics callbacks, or elliptic kernels.
- The ordinary nekRS build and executable must work without Trilinos. Prefer
  a standalone extension CMake project; a disabled-by-default `add_subdirectory`
  is acceptable only if upstream does not export usable library targets.
- No wholesale port of PR #4. Read it for lessons, not as a source of an
  intertwined implementation to copy.

The outcome is a working flow demonstration, not merely compilation, a scalar
example, natural parameter stepping with no stability calculation, or a DNS
movie. Do not declare completion without the acceptance evidence in section 13.

## 2. Scope decision: reuse the time integrator first

Use the existing incompressible nekRS time integrator as an in-process flow
map. This minimizes changes to existing numerical modules and avoids making
a new coupled steady momentum/pressure residual and preconditioner a
prerequisite for the first flow bifurcation.

For a deterministic numerical flow map `Phi_T(q, Re)`, solve

```text
F(q, Re) = q - Phi_T(q, Re) = 0
J_F v    = v - D_q Phi_T(q, Re) v
F_Re     = -D_Re Phi_T(q, Re)
```

Use LOCA to continue equilibria, including their unstable portion. Use
Anasazi on `D_q Phi_T`, NOT on `J_F`, for stability. At an equilibrium the
map multipliers give numerical growth rates

```text
sigma = log(abs(mu)) / T
omega = unwrapped_arg(mu) / T
St    = abs(omega) / (2*pi)       # D = U_infinity = 1
```

These approximate continuous-time growth rates; timestep convergence is
required. A complex conjugate pair crossing `abs(mu)=1`, away from `mu=+1`
and `mu=-1`, is the event sought. Do not call the equilibrium-map spectrum a
periodic-orbit Floquet calculation.

Separate future capabilities for direct steady operators and shooting. Do
not implement a direct spatial residual, analytic fluid tangent solver,
adjoint solver, normal forms, two-parameter continuation, or periodic-orbit
continuation in this milestone. They must not be prerequisites or claimed
features. A short perturbed DNS corroborates the Hopf result; it is not a
periodic-orbit continuation implementation.

Finite-difference map derivatives are explicitly allowed for this small CPU
demonstration. Hide them behind a derivative-action capability so an analytic
tangent implementation can replace them later. Do not inherit PR #4's
direct-residual-only or analytic-JVP-only plan: those were a different scope.

Do not confuse direct *steady formulation* with direct *linear factorization*.
The initial flow backend is explicitly iterative. If a direct algebra option
is not implemented, reject it clearly; do not silently select a strategy by
mesh size or MPI rank count. Dense operations on small Krylov/Hessenberg
matrices are permissible; global dense PDE matrices are not.

## 3. Inspected environment and source facts

These are observations from preparation, to be rechecked by the coordinator:

| Item | Observation |
|---|---|
| Machine | x86-64 Linux, Intel i7-8650U, 4 physical cores / 8 threads, approximately 16 GiB RAM |
| Existing nekRS | `/home/mahdik/workspace/nekRS`, branch `master`, SHA `995575be0ee6d7a88012c563f127c62673444eae` |
| Existing user data | That checkout has untracked `examples/bratu_gelfand/`; preserve it |
| Trilinos | `/home/mahdik/workspace/Trilinos`, SHA `5a4d907148920aef26172221bfd79c1de8295c27` |
| Trilinos version file | `Version.cmake` says `17.3.0-dev`; do not assume PR #4's 17.1.1 |
| nekRS high-level API | `src/lib/nekrs.hpp` declares `setup`, `initStep`, `runStep`, `finishStep`, `finalize` |
| Existing cylinder-named example | `examples/mv_cyl` is moving-geometry/low-Mach material, not a ready stationary external-cylinder Hopf benchmark |

The local Trilinos `git describe` can select a misleading old tag. Record the
full SHA and `Version.cmake`; use source headers as the API authority.

Use one consistent compiler/MPI/BLAS toolchain for both libraries and all
case code. Never link system MPI objects against a different Pixi/Conda MPI.
Follow applicable repository environment instructions. If using the PDECont
Pixi environment, build the entire C++ stack under its compatible toolchain;
do not assume its presence makes an independently built system Trilinos ABI
compatible. DOLFINx is not a dependency of this extension.

Initial resource limits: one build at a time, build parallelism 2, one
numerical run at a time, 1 or 2 MPI ranks, one BLAS/OpenMP thread per rank,
double precision, OCCA CPU backend. Target peak RSS below 8 GiB for the whole
job. Record actual timing and memory; no promised runtime is implied here.

## 4. Fresh source and Git procedure

The user intends to replace their fork. Repository deletion is not necessary
for implementation and is not an agent prerequisite. Do not delete a remote
fork, PR, or old local checkout as part of this handout. Preserve existing
work; use a fresh sibling directory such as `nekRS-bifurcation`.

Discover the upstream default branch instead of assuming it is literally
`main`. At preparation, upstream URLs use `master`. "From main" here means
the current upstream default branch, pinned to a recorded SHA.

Suggested discovery/checkout sequence, with the destination checked absent:

```bash
gh repo view Nek5000/nekRS --json defaultBranchRef,url
git clone https://github.com/Nek5000/nekRS.git /home/mahdik/workspace/nekRS-bifurcation
git -C /home/mahdik/workspace/nekRS-bifurcation status -sb
git -C /home/mahdik/workspace/nekRS-bifurcation rev-parse HEAD
git -C /home/mahdik/workspace/nekRS-bifurcation switch -c feature/external-bifurcation-cylinder
```

Read that checkout's `AGENTS.md`, build documentation, and release notes before
implementation. Record the base SHA in `docs/provenance.md`. Reconfigure
remotes for the user's replacement fork only when its actual URL exists;
publishing is not needed for local acceptance. Do not implement on the
default branch. Do not reset the old checkout.

## 5. Package layout and ownership

```text
extensions/bifurcation/
  CMakeLists.txt
  include/bif/
    status.hpp
    state.hpp
    provider.hpp
    derivatives.hpp
    spectrum.hpp
    task.hpp
    checkpoint.hpp
  src/
    finite_difference.cpp
    event_tracking.cpp
    hopf_localization.cpp
    workflow.cpp
    checkpoint.cpp
  adapters/trilinos/
    state_bridge.*
    nox_group.*
    linear_solver.*
    map_eigensolver.*
    loca_runner.*
  adapters/nekrs/
    session.*
    state_layout.*
    solver_snapshot.*
    flow_map.*
    parameter_binding.*
  app/nekrs_bif.cpp
  examples/brusselator_1d/
  examples/hopf_normal_form/
  examples/cylinder_hopf/
    mesh/                         # generator source + small generated mesh
    cylinder.par
    cylinder.udf
    cylinder.usr                  # only if needed for boundary/connectivity setup
    bifurcation.ini
    README.md
  scripts/
    build_cpu.sh
    run_cylinder_demo.sh
    plot_branch.py
  docs/
    provider_contract.md
    upstream_touchpoints.md
    provenance.md
    implementation_status.md
```

Names are proposed, not existing APIs. Keep the structure compact; combine
tiny implementation files where useful. The include dependency must point
from adapters toward contracts. Core compilation must not need nekRS headers.
A different PDE provider should link to the same algorithms without editing
them. Do not build Python bindings in this milestone; this boundary should
make a small C ABI binding possible later.

Keep case configuration declarative. `task` owns the source (initial/event/
restart), continuation parameter, branch type, and stopping conditions. The
cylinder UDF specifies fluid properties, geometry-independent boundary data,
and diagnostics; it must not instantiate NOX/LOCA, find snapshots, or run
continuation itself. Bind `Re -> nu = 1/Re` in the nekRS provider's parameter
binding, including coefficient/preconditioner invalidation.

### Provider contract to freeze before parallel coding

Use small C++ interfaces/composition compatible with the pinned toolchain.
Do not introduce a new compiler requirement solely for syntax. The following
is semantic pseudocode, not a promised existing header:

```cpp
struct OperationStatus {
  StatusCode code;             // success, invalid_state, solver_failure, ...
  double achieved_error;
  int iterations;
};

struct FlowMapProvider {
  StateLayout layout() const;
  OperationStatus evaluate(ConstState q, Parameters p,
                           FlowMapSpec spec, MutableState out);
  ProviderSnapshot capture();
  OperationStatus restore(const ProviderSnapshot&);
};

struct MapDerivative {
  OperationStatus apply(ConstState q, Parameters p, FlowMapSpec spec,
                        ConstState direction, MutableState out);
  OperationStatus parameter_action(ConstState q, Parameters p,
                                   FlowMapSpec spec, ParameterId,
                                   MutableState out);
};
```

Define communicator, scalar type, memory location, global IDs, free/constrained
degrees of freedom, and ghost ownership in `StateLayout`. Use opaque handles
or views plus operations; no `Tpetra::Vector` in this public contract. An
explicit CPU host-view capability is acceptable for this CPU-only milestone.
Do not pretend it is zero-copy GPU support.

Every vector action fills caller-owned output; input/output aliasing is
forbidden unless explicitly supported. Provider snapshots and workspaces
are owned objects; borrowed nekRS/OCCA objects are never destroyed by the
extension. All numerical calls are collective on the case communicator.
The provider is non-reentrant; evaluations of cloned NOX groups must serialize
through it and restore the corresponding state and parameters each time.
Caching keys include state generation, parameter generation, map specification,
layout, and discretization identity. Pointer equality is not a valid cache key.

## 6. The critical nekRS work: a correct numerical map

This is the highest-risk work package. A mathematical interface cannot hide
incorrect solver state restoration.

### 6.1 Audit before coding

Read the actual driver time loop, high-level API, nrs initialization, fluid
time-step routines, pressure solver, and UDF call sites at the pinned SHA.
Document the normal call order and return/error semantics. In particular,
`runStep` can participate in corrector iterations: do not infer that one call
is one complete timestep.

Inventory everything that can affect a repeated evaluation:

- velocity and any necessary pressure state;
- BDF velocity histories and EXT/nonlinear forcing histories;
- pressure extrapolation/correction state;
- time, timestep number, startup order, timestep history, coefficients;
- forcing/property callbacks, boundary data, and parameter-dependent caches;
- subcycling/interpolation buffers if enabled;
- projection initial guesses and solver caches that affect achieved tolerance;
- output timers/counters and UDF side effects.

The first baseline should use a supported first-order BDF/EXT configuration,
fixed timestep, stationary mesh, incompressible constant-density flow, no
scalars, no subcycling, and autonomous boundary conditions. Verify the actual
configuration spelling. A first-order setting reduces history requirements;
it does NOT prove that velocity is the complete state.

Choose and document one exact state model:

1. Prefer a one-step map whose independent input is the free velocity state,
   with all algebraic quantities reconstructed deterministically to controlled
   tolerance. Prove by the source audit and repeatability checks that no hidden
   history affects the result.
2. If the selected integrator needs independent history, include it in the
   continuation/map state and evolve the full numerical Markov state. Document
   history metrics and distinguish physical modes from computational modes.

Do not repeatedly fill arbitrary histories from `q` and present that startup
procedure as the exact flow map of the mature multistep integrator. A canonical
lifting would define a different numerical map and needs separate consistency
evidence; it is not the default shortcut.

Pressure is an algebraic variable unless the audited algorithm requires it as
numerical state. Do not give it a fictitious physical evolution equation.
Fix a pressure gauge only where a nullspace exists and preserve the actual
outlet pressure/traction convention; a fixed pressure outlet may already fix
the gauge. Homogeneous perturbation boundary conditions must follow from the
same prescribed base-flow boundary data.

### 6.2 Transactional evaluation

Implement `evaluate(q, Re, T, out)` as a transaction:

1. Capture/guard the enclosing solver session state.
2. Install `Re` and refresh affected viscosity/coefficient data.
3. Scatter free global state entries to conforming nekRS fields; impose fixed
   boundary values and reconstruct required algebraic state.
4. Restore the map's documented initial clock/history state.
5. Advance exactly `m` fixed steps, with `T=m*dt`, using the ordinary nekRS
   step sequence. No adaptive timestep, wall-time exit, or analysis output.
6. Pack the output, validate solver convergence and finiteness collectively.
7. Restore the enclosing session state on success and failure.

Keep accepted-state installation separate from trial map evaluation. A rejected
Newton trial must not change the next accepted output or checkpoint. Suppress
checkpoint writes and user diagnostics inside trial evaluations while still
executing callbacks needed for physical forcing and properties. Document the
allowed UDF behavior: hidden mutable user state must be registered with the
snapshot or unsupported explicitly.

A small `captureAnalysisState/restoreAnalysisState` hook is justified only if
existing public access cannot satisfy this contract. Its payload remains
opaque to continuation algorithms. Do not use state restoration as a reason
to move bifurcation code into `nrs_t`.

### 6.3 State layout, norm, and two-dimensional subspace

Use unique physical unknowns, not duplicated element-interface/halo values.
Reuse nekRS global numbering/gather-scatter facilities. Essential boundary
values are prescribed, not independent Newton unknowns. Handle periodic
identification and local singleton nodes explicitly.

Use a volume-normalized kinetic-energy inner product. For the initial CPU
Tpetra bridge, a diagonal transformation `z_i=sqrt(w_i)*q_i` is a useful way
to make standard Euclidean NOX/Belos/Anasazi operations equal the selected
physical metric. Here `w_i>0` are assembled unique-node quadrature weights,
including field scales. Apply the same transformation to state, residual,
derivatives, eigenvectors, and restart metadata. Do not apply mass weighting
twice. The weights must be independent of `Re` for this fixed-geometry case.

If nekRS needs 3D elements, extrude a 2D quadrilateral mesh by one periodic
element layer. Use a fixed modest span, e.g. `Lz=1`, since spanwise-independent
physics does not need a thin domain. Explicitly restrict base states, Krylov
seeds/directions, and outputs to spanwise-independent `u,v` with `w=0`.
One element layer alone does not eliminate spanwise polynomial modes.
Document whether the bridge reduces to 2D unknowns or supplies a consistent
projector. Verify the invariance and do not count excluded projector modes as
physical eigenvalues.

Use the FULL cross-stream domain. Do not impose centerline reflection on
perturbations: the cylinder shedding instability would be suppressed.

## 7. Trilinos integration: minimal and source-verified

Keep `/home/mahdik/workspace/Trilinos` read-only. Build out of source and record
the installation prefix, SHA, compiler versions, MPI library/version, BLAS,
scalar precision, and all CMake options. Reuse an existing compatible install
if its provenance can be established.

Start with the required stack: Teuchos, NOX/LOCA, Tpetra, Belos, Anasazi and
their actual dependency closure. Do not enable Thyra, Stratimikos, MueLu,
Ifpack2, Amesos2, or every optional package just because the old PR did.
Enable an extra dependency only when a selected API requires it.

Verified local paths to read before choosing flags or implementing wrappers:

```text
Version.cmake
packages/nox/cmake/Dependencies.cmake
packages/nox/CMakeLists.txt
packages/nox/src/NOX_Abstract_Group.H
packages/nox/src/NOX_Abstract_Vector.H
packages/nox/src-tpetra/NOX_Tpetra_Vector.hpp
packages/nox/src-tpetra/NOX_Tpetra_Vector_def.hpp
packages/nox/src-loca/src/LOCA_MultiContinuation_AbstractGroup.H
packages/nox/src-loca/src/LOCA_MultiPredictor_Restart.H
packages/nox/src-loca/src/LOCA_AnasaziOperator_AbstractStrategy.H
packages/nox/src-loca/src/Anasazi_LOCA_MultiVecTraits.H
packages/nox/src-loca/src/Anasazi_LOCA_OperatorTraits.H
packages/nox/examples/lapack/LOCA_Brusselator/BrusselatorContinuation.C
packages/nox/examples/lapack/LOCA_Brusselator/BrusselatorHopfContinuation.C
```

In this checkout `NOX_ENABLE_LOCA` is a package option, and the NOX Tpetra
abstract implementation is separately configurable. The optional Kokkos
solver-stack switch brings a much larger dependency set; it is not required
merely because a NOX Tpetra vector is used. Verify configuration output and
installed headers instead of guessing version-dependent switches.

Use CPU Serial initially and MPI distribution, double precision, an explicit
compatible local/global ordinal pair, Release builds, and no dependency test
suites. Produce `scripts/build_cpu.sh` with exact working CMake flags after
the build is verified. This handout intentionally does not invent a guaranteed

[493 more lines in file. Use offset=401 to continue.]
configure command for an unbuilt development SHA.

### Adapter behavior

- `NoxGroup`: state/parameter management, correct clone/copy/invalidation,
  residual, JVP, parameter derivative, and inverse application through the
  backend. No `nrs_t`, OCCA, mesh, or elliptic logic.
- `LinearSolver`: Belos factory-selected GMRES; enable flexible GMRES when
  the chosen preconditioner requires it. Report requested/achieved residual,
  convergence status, and iteration counts. Start with identity preconditioning
  as an explicit option; never implement `J^-1 v = v` as a solver.
- `LocaRunner`: use LOCA's predictor/corrector and arclength machinery.
  Correct the initial point before continuing. Cache each accepted state.
- `MapEigensolver`: nonsymmetric Anasazi Krylov-Schur through the available
  factory/traits interfaces, acting on the map derivative. Target largest
  magnitude map multipliers, retain complete complex pairs, and assess Ritz
  residuals using independent operator applications.

For `F=q-Phi`, do not use `J_F` eigenvalues as growth rates or a near-zero
steady-Jacobian spectral transform as a substitute for the map operator.
If an operator callback returns `void`, provide an outer checked failure
protocol. Arrange collective error agreement at safe boundaries; if a failure
is unrecoverable inside library collectives, terminate coherently with a clear
diagnostic. Never return zeros and permit apparent eigensolver success. A
post-hoc allreduce cannot repair ranks already taking different collectives.

Do not choose global sparse factorization as a fallback for a stalled flow
Krylov solve. First verify the map/JVP, timestep, horizon, norm, and tolerances.
Only add a provider-owned iterative preconditioning capability when measured
need warrants it; retain the same engine interface.

## 8. Derivatives, continuation, and event localization

### Finite differences

Use a centered map action initially:

```text
D Phi(q) v ~= [Phi(q+h*v) - Phi(q-h*v)] / (2*h)
h = eta * max(1, norm(q)) / norm(v)
```

Define this in the selected weighted state coordinates. Handle zero direction
without dividing by zero. For a normalized direction try `eta` around
`1e-4, 1e-5, 1e-6`; measure the error plateau and record the selected value.
Do not assume the smallest perturbation is best. Center parameter differences
similarly, respect positive `Re`, and restore coefficients between evaluations.

Pressure/velocity solve noise divided by `h` must remain below the derivative
accuracy budget. Use deterministic initialization and sufficiently tight inner
solves. Check directional consistency, approximate additivity, and repeatability
on a corrected nontrivial cylinder state. A constant perturbation alone is
not enough to exercise advection and incompressibility.

Initial numerical knobs, to be measured and finalized:

| Knob | Starting choice |
|---|---|
| Time horizon | `T=2` convective units, exactly an integer number of timesteps |
| Timestep | `dt=0.02` or smaller if actual CFL requires; compare with `dt/2` |
| Nonlinear correction | up to 12 iterations, normalized residual target `1e-7` |
| Krylov | restarted GMRES, initial dimension 40, explicit total-iteration limit |
| Newton inner accuracy | inexact tolerance initially `1e-2`, tighten with residual |
| Map eigensolve | 4–6 multipliers, subspace initially 40–60, complete pairs |
| Eigen residual target | `1e-5` normalized; tighten/check if crossing is uncertain |
| Continuation | `Re=35` to `60`, initial/max parameter increment around 2 |
| Event localization | Reynolds bracket width at most `0.1`, reliable opposite signs |

These are proposed defaults, not established results. A small map horizon can
make `I-DPhi` poorly conditioned; a long one can amplify unstable directions
and alias frequency. Measure these effects before tuning GMRES blindly. If
needed, adjust `T` while retaining enough phase margin and checking recovered
growth rates with a second horizon.

Continue with LOCA arclength and an explicitly scaled parameter, e.g.
`rho=Re/50`, while reporting physical `Re`. Ensure `F_rho=50*F_Re` if `rho` is
the internal parameter. Use a secant predictor after an explicitly configured
first-step predictor. Do not claim the cylinder branch proves fold traversal;
the requested cylinder event is Hopf.

At each accepted state, match the shedding conjugate pair using eigenvalue
proximity AND mass-weighted subspace overlap. Do not select an independently
sorted eigenvalue index at every step. Recycle converged vectors as a starting
subspace, not the previous operator's Hessenberg relation.

Bracket `g(Re)=log(abs(mu(Re)))/T`. Refine using safeguarded secant/bisection:
interpolate a nearby equilibrium, NOX-correct at the trial `Re`, recompute and
match the pair, then update the bracket. Do not interpolate eigenvalues and
call that a solved critical point. Local one-parameter refinement is legitimate
here; it does not require implementing LOCA's full augmented two-parameter
Hopf tracking interfaces.

An event is `HOPF_CANDIDATE` until the following are established:

- a simple conjugate pair, finite nonzero frequency, opposite growth-rate
  signs on corrected equilibria, and a nonzero crossing slope above numerical
  uncertainty;
- pair residuals and derivative sensitivity support those signs;
- absence of a competing unresolved neutral/unstable mode in the computed
  leading spectrum (otherwise report ambiguity and enlarge the subspace);
- timestep/horizon checks distinguish a physical instability from a map
  artifact; reject a `+1` or `-1` crossing as the requested Hopf event;
- independent perturbed DNS provides consistent decay/growth and frequency.

Use `HOPF` for the numerically validated crossing. Do not claim a computed
first Lyapunov coefficient or supercriticality from eigenvalues alone. The
literature identifies the canonical transition as supercritical; this milestone
does not calculate its normal form.

## 9. Canonical example: stationary circular-cylinder wake

### Physical problem and literature

Solve the nondimensional incompressible equations with `D=U_infinity=rho=1`:

```text
du/dt + (u.grad)u = -grad(p) + (1/Re)*Laplacian(u)
div(u) = 0
```

The first cylinder-wake instability is a two-dimensional oscillatory Hopf
bifurcation near `Re=46–47` for sufficiently unconfined/resolved flow. Jackson
locates onset using extended steady equations; Provansal, Mathis, and Boyer
study the near-threshold dynamics. A later cylinder-wake study uses `Re_c=46.6`.
These are reference scales, not an exact expected value for an arbitrary
coarse finite-domain discretization. See [R1]–[R3] in section 15.

Use a stationary unit-diameter cylinder centered at `(0,0)` in a domain
initially `[-15,30] x [-15,15]`. Boundary intent:

- inlet and distant top/bottom: prescribed uniform `(u,v)=(1,0)`;
- cylinder: stationary no-slip;
- outlet: one supported, documented open/outflow condition consistent with
  the pinned nekRS pressure treatment;
- extrusion planes: periodic, with an explicitly 2D state subspace.

Do not assume these boundaries reproduce a literature domain exactly. Record
them and assess domain sensitivity at the critical point. Do not substitute
the tightly confined Turek–Schäfer geometry and retain an unconfined `Re_c`
comparison. Do not reuse `mv_cyl` without removing its unrelated moving-mesh
and low-Mach physics through a NEW case, not edits to the existing example.

### Small reproducible mesh

Start with approximately 200–500 quadrilateral elements in the plane, one
periodic extrusion layer, polynomial degree 4 (`5^3` element-local nodes per
hex). At 500 hexes that is 62,500 element-local nodes before shared-node
identification, not millions. Count actual free/unique state unknowns in logs.
This is a starting budget, not a proven adequate mesh.

Use a body-fitted O-grid near the cylinder, at least 32–48 effective boundary
points around its circumference, smoothly graded cells, and a resolved near
wake out to approximately `x=10`. Far-field elements can be much larger.
Preserve curved geometry accurately; reducing polynomial order must not
turn the cylinder into a coarse polygon without documenting the effect.

Inspect existing read-only cylinder meshing sources in
`/home/mahdik/workspace/pdecont/nekStab/` and upstream Nek examples. The local
`nekStab/example/cylinder_re100/000_dns/fpcyl.box` is a mesh-construction
reference, not automatically the desired small ready-to-run 3D mesh. Choose
one reproducible supported generation route, e.g. structured quads plus
extrusion and an available converter. Do not write a new general mesher.

Check in the small generated `.re2` with generator inputs/script, tool versions,
boundary-ID table, periodic connectivity recipe, and provenance/license notice
when adapting another project. Do not require a GUI mesher, manually generated
untracked file, or downloaded flow solution to run the final example. Confirm
all boundary IDs, positive element Jacobians, and periodic pairing under MPI.

### Run sequence

1. Run the ordinary nekRS executable on the new case at `Re=35`, confirming
   the base solver works before using continuation. Use startup DNS only to
   obtain a seed, then NOX-correct it and save the result.
2. Continue corrected equilibria through `Re=35..60`. Near onset, the steady
   solution must continue after it loses stability; DNS alone cannot do this.
3. Locate the complex-pair crossing and persist its state/eigenvectors.
4. At a subcritical and supercritical point, restart ordinary time integration
   from the corrected equilibrium plus a small real shedding-mode perturbation.
   Use a relative perturbation around `1e-4` initially and verify an interval of
   linear behavior. Fit growth/decay from a mode projection or an oscillation
   envelope, not the instantaneous logarithm of a sign-changing signal.
5. Compare the frequency with the eigenvalue prediction using `St=omega/(2*pi)`.
   A few hundred convective time units may be necessary; measure and choose
   points several Reynolds units from onset to reduce critical slowing down.
   Do not require a fully saturated long-time limit cycle for this milestone.
6. Repeat the core continuation/localization on two MPI ranks and demonstrate
   interrupted continuation resume in a new process.

Mesh/time/domain sensitivity runs should correct and evaluate only near the
crossing where possible; do not repeat every far-from-onset continuation point.

## 10. Checkpoint and restart are part of the implementation

Implement restart early, before the full cylinder sweep. Distinguish the
provider's transient in-memory transaction snapshot from persistent accepted
continuation checkpoints. A saved CSV or one in-memory eigenvector is neither
a branch restart nor an event restart.

Manifest schema version stays 1. Persist:

- accepted state and all independent numerical-history components;
- physical parameters and internal parameter scaling;
- branch ID, accepted step number, arclength coordinate, tangent/secant history,
  orientation, next step size, and predictor/continuation method;
- map specification (`T`, `dt`, integrator, startup/history policy), state
  transformation/metric, constraints and pressure convention;
- bracket endpoints and their states, tracked complex eigenspace, pair identity,
  eigenvalues/residuals, and event status if localization is active;
- mesh hash, stable state IDs/layout, boundary definitions, provider metadata,
  source SHAs, build provenance, and numerical tolerances.

Use parallel field I/O already available in nekRS where suitable, or per-rank
binary shards plus an atomic manifest. Write new generation filenames, then
commit the manifest last after collective success. Preserve the previous
complete checkpoint if interrupted. Do not gather large fields to rank zero.

Minimum acceptance is a fresh-process resume at the same rank count for both
1-rank and 2-rank runs, plus seeding a separate stability task from the saved
Hopf event. If repartitioning is not supported, reject a different rank/layout
explicitly and document that limitation. Stable global IDs should leave a
clear path to repartitioning later; do not claim it without demonstrating it.

Changing Krylov limits and compatible solver tolerances on resume is allowed.
Changing the map timestep, horizon, mesh, or state definition changes the
numerical problem: require an explicit seed-and-recorrect workflow rather
than calling it an exact continuation resume. Use LOCA's restart predictor/
state mechanisms as supported by the pinned source; do not serialize raw
C++ pointers or an opaque LOCA object graph.

## 11. Swarm execution plan

### Model and coordination rule

The user requested `gpt-5.5-luna` workers. That exact identifier was NOT among
the models available in the preparation session; `gpt-5.5` and
`gpt-5.6-luna` were distinct available choices. The implementing coordinator
must check its own model catalog. If the requested identifier is unavailable,
ask the user which available model to use before launching those workers;
do not invent an alias or silently substitute a model. Planning and local
read-only preparation can continue while that choice is pending.

Use one coordinator plus at most three workers at a time. The coordinator
owns the provider contract, shared build entry points, integration branch,
acceptance decisions, and final report. Workers receive bounded tasks and
explicit file ownership. Keep every handoff in `docs/implementation_status.md`
with completed work, exact command, observed result, and next blocker.

Use isolated worktrees/branches if available. Freeze interface headers before
the implementation wave and have the coordinator merge/cherry-pick accepted
work. Workers must not concurrently edit the same tree/files, race builds,
change public interfaces unilaterally, or launch simultaneous MPI jobs.
One numerical runner owns execution on this small machine.

### Wave A — discovery and contract, three independent audits

| Worker | Bounded assignment | Deliverable and stop point |
|---|---|---|
| A: nekRS provider auditor | Read setup/step/state/pressure/UDF code; enumerate state, supported first-order mode, public access and minimum hooks | `provider_contract.md` draft and `upstream_touchpoints.md`; no numerical-kernel edits |
| B: Trilinos/build auditor | Inspect pinned APIs/options; establish minimal CPU-MPI dependency build and adapter signatures | Working build recipe, source/version manifest, small linkable interface probe; no CFD implementation |
| C: cylinder case designer | Select literature/domain/BCs and generate the small stationary-cylinder case | Reproducible mesh and plain-nekRS seed run with sizes, BC checks and timing; no continuation code |

Coordinator gate: resolve state representation, metric, parameter scaling,
map semantics, source base, target layout, and ownership. Approve a concrete
upstream edit list internally; this does not require routine user permission
when within scope. If the only proposed route requires rewriting the pressure
or momentum scheme, reject it and investigate a narrower adapter route.

### Wave B — implementation after contract freeze

| Worker | Owned files | Required result |
|---|---|---|
| D: provider implementation | `adapters/nekrs/` and only explicitly assigned minimal upstream hooks | Deterministic transactional map, correct layout/parameter updates, no library logic in solvers |
| E: numerical library adapter | `adapters/trilinos/`, `src/finite_difference.cpp` | NOX correction, real Belos inverse, LOCA arclength, Anasazi map spectrum using a solver-independent example |
| F: persistence/workflow | `src/checkpoint.cpp`, task/workflow implementation and app CLI | Atomic accepted-state/event persistence, restart/seed tasks, structured diagnostics; use mock provider until D is ready |

Coordinator integrates often and resolves only actual interface conflicts.
Do not let three agents invent three vector layouts or checkpoint schemas.

### Wave C — integrate and establish numerical correctness

Run the small 1D Brusselator first, using the shared adapter/engine rather
than merely launching a stock Trilinos example. A simple reaction-diffusion
provider is sufficient; record equations, domain, BCs, parameters and an
analytically checked threshold. The supplied LOCA Brusselator sources are
API references, not acceptance evidence for our provider.

Also provide a tiny solver-independent Hopf normal-form example
`x'=a*x-w*y-(x*x+y*y)*x`, `y'=w*x+a*y-(x*x+y*y)*y`, with `w=1`, threshold `a=0`.
Use it to verify map-vs-residual eigenvalue interpretation, frequency units,
pair matching, localization, and resume. It is a fast runnable numerical
example, not a new automated test suite.

Assign one worker to source/ownership review, one to map and derivative
diagnostics, and one to example/workflow documentation. Only the designated
runner executes builds/MPI. The coordinator ensures the stock solver still
runs with the extension disabled and the library builds without nekRS headers
when using only its independent provider.

For the Brusselator, one convenient analytically checkable choice is
`u_t=d*u_xx+A-(B+1)*u+u^2*v`,
`v_t=d*v_xx+B*u-u^2*v`, `A=1`, equal `d>0`, homogeneous Neumann BCs on a
short interval such as `[0,1]`. The equilibrium is `(1,B)`, homogeneous Hopf
is `B=2`, frequency 1, and equal diffusion shifts nonconstant modes left.
Document temporal discretization error if using a numerical flow map. This
choice avoids assuming the parameters in the stock Trilinos example have
the same threshold.

### Wave D — cylinder and MPI completion

The numerical runner performs section 9's sequence and section 13's evidence
matrix. A second worker reviews event matching, derivative/eigen residuals,
and mesh/time sensitivity without changing acceptance thresholds after seeing
results. A third worker prepares plots and README commands from real logs.

Coordinator must inspect actual data, not accept a worker's "PASS" alone.
Fix failures at their owning layer. Integrate, rerun only affected checks,
and finish the documented serial/MPI/restart demonstration. Do not stop at a
successful scalar example or leave the flow run as a recommendation.

### Standard worker prompt

```text
Read this handout and the applicable AGENTS.md. Implement only your assigned
work package and owned files on your assigned branch/worktree. Follow the
frozen provider contract. Existing nekRS numerical modules must remain
unchanged except the coordinator-listed narrow hooks. Do not modify reference
Trilinos/nekStab sources. Use numerical examples and diagnostics; do not add or
run an automated test suite without user instruction. No global PDE matrix
assembly, hidden operator failures, or fabricated acceptance data. Coordinate
build/MPI access with the single runner. Return your commit/diff, exact commands,
observed numerical evidence, limitations, and any interface request. Never
declare the overall task complete; that belongs to the coordinator.
```

## 12. Required executable workflow

Implement and document a compact `nekrs-bif` executable. The following is the
REQUIRED PROPOSED CLI to implement, not a claim that upstream provides it:

```bash
nekrs-bif --case cylinder.par --config bifurcation.ini --task seed --output run-seed
nekrs-bif --case cylinder.par --config bifurcation.ini --task continue --seed run-seed/checkpoint.json --output run-serial
nekrs-bif --case cylinder.par --config bifurcation.ini --task resume --restart run-serial/checkpoint.json --output run-resumed
nekrs-bif --case cylinder.par --config bifurcation.ini --task stability --seed run-serial/events/HB1/manifest.json --output run-event
mpirun -np 2 nekrs-bif --case cylinder.par --config bifurcation.ini --task seed --output run-seed-mpi2
mpirun -np 2 nekrs-bif --case cylinder.par --config bifurcation.ini --task continue --seed run-seed-mpi2/checkpoint.json --output run-mpi2
mpirun -np 2 nekrs-bif --case cylinder.par --config bifurcation.ini --task resume --restart run-mpi2/checkpoint.json --output run-mpi2-resumed
```

Provide a safe `--stop-after-accepted N` option for demonstrating interrupted
continuation. Save enough state that stopping before the Hopf crossing and
resuming discovers the same event. Preserve output history without duplicate
accepted step/event IDs. Use separate output paths for independent serial and
MPI runs. If exact CLI names change during implementation, update this recipe
and the README together and supply the actual executed commands.

`scripts/run_cylinder_demo.sh` should build case kernels once, run the required
tasks with configurable `NP=1`/`NP=2`, capture logs, and produce a final result
summary. It must use shell failure propagation so a plotting/logging pipeline
cannot hide a failed solver. Avoid unsupported global OpenMPI environment
workarounds unless a measured runtime problem requires one; document any such
workaround separately from the physics result.

Write one structured row per accepted point containing at least:

```text
branch, step, Re, arclength, residual, divergence,
mu_real, mu_imag, sigma, omega, St, eigen_residual,
event_status, newton_iterations, linear_iterations,
map_evaluations, elapsed_seconds, mpi_ranks
```

Include pair IDs/subspace information and all computed leading eigenvalues in
a spectrum file. Rank zero writes small metadata/CSV; all ranks participate
in numerical and field checkpoint operations. Generate a plot of growth rate
versus `Re` with bracket/localized point, frequency versus `Re`, and a base-flow
or perturbation vorticity image. Use ordinary plotting tools; these are static
numerical figures, not generated illustrations.

## 13. Acceptance matrix and numerical evidence

Tolerances below are initial acceptance targets, not measured claims. If a
target is unattainable, investigate and document the numerical error budget;
do not silently loosen it. A genuinely blocked run needs an exact reproducer
and explicit remaining work, not a "completed" label.

| Gate | Required evidence |
|---|---|
| Boundary | Reviewable diff of every existing nekRS file changed; no modifications to existing numerical algorithms; core headers independent of nekRS/Trilinos |
| Optional dependency | Ordinary nekRS build/run with extension off; extension build with pinned minimal Trilinos stack |
| Layout | Nonconstant field pack/scatter round-trip near roundoff; correct unique ownership, quadrature volume, periodic/essential BCs on 1 and 2 ranks |
| Map determinism | `Phi(q), Phi(q+v), Phi(q)` gives matching first/last outputs to a tolerance far below the FD signal; failures/rejected trials restore state |
| Derivatives | Step-size sweep and direction/additivity diagnostics on the nontrivial flow; selected FD action stable enough for `1e-5` Ritz residuals |
| Independent examples | Shared-engine Brusselator and Hopf normal form produce their documented equilibria, thresholds/frequencies, and restart behavior |
| Genuine steady flow | `||q-Phi_T(q)||/max(1,||q||)<=1e-7` AND normalized one-step change divided by `dt` around `1e-6` or better; no oscillation hidden by a period-length map |
| Incompressibility | Record normalized divergence and boundary mismatch; initial divergence target `1e-6`, consistent with discretization and pressure solve tolerance |
| Branch | Corrected nonzero cylinder equilibria on both sides of onset; continuation through the unstable segment |
| Hopf | Matched nonreal pair crossing the unit circle, opposite resolved signs, normalized Ritz residuals `<=1e-5`, `Re` bracket width `<=0.1` |
| Literature scale | Initial coarse screening target `Re_H` in `[43,50]`, `St_H` in `[0.10,0.14]`; this interval alone is NOT acceptance of accuracy |
| Timestep | Halve `dt` near the event and recorrect; target critical-Re change `<=0.5` and frequency change `<=3%`; reduce again if not met |
| Mesh | Raise degree from 4 to 5 or modestly refine near-cylinder/wake cells; target critical-Re change `<=1` and frequency change `<=5%` |
| Domain | Increase distant boundary placement on a comparable near-body mesh near the event; document change and target critical-Re shift `<=1` |
| DNS corroboration | A small perturbation decays below and grows above onset in a linear window; frequency agrees within approximately `10%` with the eigen prediction |
| MPI | Independent 1/2-rank results have the same crossing pair/event; target `|delta Re_H|<=0.1`, frequency difference `<=1%`, comparable residuals |
| Restart | Same-rank fresh-process continuation resumes on both 1 and 2 ranks, including a stop before onset; event state seeds a new stability task |
| Reproducibility | Small mesh plus generator, exact commands, configurations, SHAs, logs, numeric results, plots, timings and peak memory |

Do not compare eigenvectors entrywise across MPI layouts or complex phases;
compare physical observables and invariant subspaces after appropriate mapping.
At the event, a sign is only resolved if it exceeds the sensitivity to FD step,
inner tolerance, and eigen residual. A tiny computed `sigma` alone is not proof.

The `[43,50]` range is deliberately a coarse sanity screen around literature
onset, not a published universal bound. The sensitivity gates prevent using a
poorly resolved mesh just because a number falls inside it. If the final
accepted case needs slightly more than 500 elements, retain a small mesh and
state the measured size/runtime honestly; preserve near-body resolution before
shrinking the physical domain.

Keep the original coarse case as a cheap execution smoke case only if its
numerical inaccuracies are clearly labeled. Deliver at least one fully accepted
small flow configuration; a known-failing coarse case cannot be the final demo.

## 14. Failure handling and finishing rule

Investigate in this order when progress stalls:

1. Confirm a plain nekRS run of the same case, MPI, mesh and BCs.
2. Check state restoration, parameter invalidation, fixed timestep/startup,
   pressure constraints and 2D invariance.
3. Measure map repeatability and FD noise before altering Newton/eigen settings.
4. Verify metric/ownership and map-vs-residual operator sign/identity.
5. Check pair tracking, eigen residuals, spectral subspace size, horizon aliasing.
6. Adjust inner tolerances, horizon, or subspace based on evidence.
7. Refine the mesh/time/domain only after algebra and state semantics are sound.

Stop blind reruns. After repeated identical failures, produce a small numerical
reproducer and an owning-layer diagnosis, then fix that layer. Dependency
failures require the exact compiler/library versions, full command, failure
point, and expected result. Do not remove convergence checks to get a run to
exit successfully. Do not replace the requested cylinder flow with a scalar
normal form and call the objective complete.

The coordinator's final handoff must include:

- implementation branch/base SHA and concise architecture summary;
- exact inventory and rationale for every existing nekRS-module edit;
- a single working build recipe and serial/MPI demo commands;
- measured `Re_H`, `omega_H`, `St_H`, brackets and uncertainty/sensitivity data;
- fresh-process restart and event-seeding evidence;
- actual machine, rank count, state size, time and memory;
- README containing only current supported behavior and verified results;
- explicit limitations: CPU-only validation, finite-difference derivatives,
  time-discrete map, supported rank/layout restart, and no computed normal form
  or periodic-branch continuation unless separately implemented and verified.

No automated test suite is requested. Numerical examples and their recorded
acceptance diagnostics are the evidence. Do not publish a PR, close PR #4,
delete a fork, or modify the user's existing checkout merely to finish this
implementation. Those are separate repository-management actions.

## 15. References and how to use them

- [R1] C. P. Jackson (1987), *A finite-element study of the onset of vortex
  shedding in flow past variously shaped bodies*, JFM 182, 23–45.
  [Publisher](https://doi.org/10.1017/S0022112087002234).
  Supports the canonical flow and Hopf interpretation, not our coarse-mesh
  acceptance values or proposed numerical implementation.
- [R2] M. Provansal, C. Mathis, L. Boyer (1987), *Bénard–von Kármán instability:
  transient and forced regimes*, JFM 182, 1–22.
  [Publisher](https://doi.org/10.1017/S0022112087002222).
  Supports near-threshold transient/oscillation corroboration.
- [R3] S. Bagheri (2013), *Koopman-mode decomposition of the cylinder wake*, JFM.
  [Author-hosted paper](https://www.mech.kth.se/~shervin/pdfs/2013_jfm_koopman_cylinder.pdf).
  Gives the familiar `Re_c=46.6` reference; do not demand that exact number on
  a different finite domain and coarse numerical map.
- [R4] [NOX/LOCA architecture](https://trilinos.github.io/nox_and_loca.html).
  Abstract groups/vectors allow the separation specified here. Match actual
  signatures to the pinned local source, not an old online example.
- [R5] [nekRS upstream](https://github.com/Nek5000/nekRS) and
  [case documentation](https://nekrs.readthedocs.io/en/latest/problem_setup/case.html).
  Check supported CPU build, current case schema and mesh conventions against
  the pinned source. Documented `.re2` compatibility does not guarantee an old
  `.usr`/boundary map can be copied unchanged.
- [R6] [nekStab](https://nekstab.github.io/) and its supplied read-only local
  source. Use for CFD time-stepper/Krylov and mesh-design reference; do not
  port its Nek5000 internals into the generic continuation library.

The equations, contracts, task decomposition, proposed tolerances and mesh
budgets above are design recommendations. None is a claim of an already
performed cylinder computation.

## 16. Milestone handoff — 2026-09-11

This section is an execution log for the next agent. It supersedes neither the
acceptance matrix nor the no-fabrication rule above.

### Implemented in the working tree

- Added the optional `NEKRS_BUILD_BIFURCATION` build path. The default nekRS
  build remains independent of Trilinos.
- Added the provider/state/layout contracts, finite-difference map derivatives,
  transactional nekRS flow-map evaluation, BDF1/EXT1 history reconstruction,
  pressure gauge removal, MPI-aware layout norms, checkpointing, event records,
  and the standalone normal-form smoke workflow.
- Added a concrete `NoxVector` and `NoxGroup` in
  `extensions/bifurcation/adapters/trilinos/nox_group.{hpp,cpp}`. The group
  derives from `LOCA::MultiContinuation::AbstractGroup`, implements the NOX and
  LOCA parameter APIs, computes `q-Phi_T`, supplies matrix-free Jacobian and
  Reynolds actions, and delegates inverse actions to Belos GMRES.
- Added `runLocaContinuation()` in
  `extensions/bifurcation/adapters/trilinos/loca_runner.{hpp,cpp}`. It creates
  a LOCA Tpetra factory/global-data object, configures a real `LOCA::Stepper`
  and NOX status tests, installs the custom group, and calls `Stepper::run()`.
  The nekRS executable accepts `--driver loca` to select this path; the
  existing explicit map workflow remains available for diagnostics.
- Changed the Belos adapter to instantiate `PseudoBlockGmresSolMgr` directly.
  The supplied local Belos registration does not register the factory alias
  `GMRES` at runtime. The Anasazi block-size workaround uses its documented
  defaults and requires at least three Krylov blocks.
- Added the required `__okl__` UDF callbacks to
  `extensions/bifurcation/examples/cylinder_hopf/cylinder.udf`. This fixes the
  prior `Cannot find oudf or okl section in udf` compilation failure.

### Build evidence

On branch `feature/bifurcation-cylinder`, with build parallelism 2:

```text
cmake --build /tmp/nekrs-bif-full -j2
  passed: ordinary nekRS plus nekRS adapter, without Trilinos

cmake --build /tmp/nekrs-bif-trilinos --target nekrs-bif -j2
  passed: Trilinos adapter and nekrs-bif executable
```

The Trilinos link emits warnings about the local supplied tree mixing Kokkos
5.1 and 5.2 and stale/versioned LOCA-Tpetra library paths. Those warnings are
not resolved by this milestone.

### Runtime evidence and exact limits

The UDF now compiles and nekRS setup reaches fluid initialization when the
legacy mesh `/home/mahdik/workspace/pdecont/nekStab/example/cylinder_re1m/000_dns/1cyl.re2`
is copied into a temporary run directory. That mesh reports one boundary ID;
the repository case currently requests four boundary types and therefore aborts
with:

```text
Size of fluid velocity boundaryTypeMap (4) does not match number of boundary IDs in mesh (1)!
```

No timestep, continuation, Hopf localization, or MPI result was obtained from
that case. The repository still has no body-fitted stationary-cylinder `.re2`
mesh. `examples/mv_cyl/mv_cyl.re2` is not a substitute: it is a moving-
geometry/low-Mach material example and did not provide a usable boundary map
for this workflow. A KTH Nek5000 `ext_cyl.re2` was inspected/downloaded only as
an external reference; it is a 2-D binary mesh and was not inserted or claimed
as a nekRS 3-D validation mesh.

A temporary 1-rank legacy-mesh diagnostic with a single coordinate-classified
Dirichlet boundary and very low polynomial order reached the map eigensolver.
It exposed and fixed the Anasazi `Block Size` type mismatch and the minimum
Krylov-block requirement. It produced no complex-pair/Hopf result and must not
be used as cylinder evidence. A subsequent `--driver loca` diagnostic reached
Belos; the local Belos factory alias failed (`PSEUDOBLOCK GMRES` was not
registered), which is why direct solver-manager construction was added. The
run was stopped before a successful continuation result. The supplied
Trilinos/Kokkos ABI warnings and cleanup aborts remain an infrastructure risk.

### Required next actions

1. Provide or generate a proven body-fitted stationary external-cylinder mesh
   in `extensions/bifurcation/examples/cylinder_hopf/mesh/`, with documented
   boundary IDs, 2-D-like extrusion/periodicity, and provenance. Do not use
   `mv_cyl.re2` or the one-ID legacy mesh as the final benchmark.
2. Run plain nekRS on that mesh first, then run the adapter with one timestep,
   map repeatability, and derivative checks before enabling continuation.
3. Harden LOCA exception/cleanup handling and validate the custom group on a
   small provider-independent problem. Confirm `computeDfDpMulti` column
   semantics, accepted-state installation, predictor/corrector output, and
   restart behavior.
4. Repair the Trilinos installation/configuration so the final build has no
   stale absolute library paths or mixed Kokkos ABI. Re-run clean 1-rank and
   2-rank builds.
5. Only after those checks, execute fresh one- and two-rank cylinder runs and
   record the map residual, multiplier residual, Hopf Reynolds number,
   frequency/Strouhal number, timestep/mesh convergence, and DNS comparison.

There is currently no verified cylinder Hopf value, frequency, 1-rank result,
or 2-rank result in this repository.
