# Provider contract

`include/bif/provider.hpp` is the only dependency required by the numerical
core. A provider owns the discretization and is non-reentrant. Every
`evaluate` call is a transaction: parameters, state, histories, clocks, solver
workspaces, and user callbacks must be restored before returning. Trial calls
must not write output files.

`StateLayout` describes the local portion of a unique global state with stable
IDs, free-DOF flags, positive kinetic-energy metric weights, and communicator
metadata. The current reference examples are host-resident. GPU zero-copy and
rank repartitioning are deliberately not claimed.

The production map is `Phi_T`; `FlowMapResidual` exposes `q-Phi_T(q)` and its
Jacobian action `v-DPhi_T(v)`. Finite differences are centered and use the
selected weighted norm. Parameter derivatives are centered and reject a
non-positive Reynolds perturbation. Input/output aliasing is not supported.

The nekRS adapter accepts callbacks for pack/scatter and transactional state
capture instead of leaking `nrs_t` or OCCA into these headers. This is the
narrow integration point where a future upstream analysis-state hook belongs.
