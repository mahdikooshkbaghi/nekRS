# Brusselator contract example

The shared provider contract is intended for the one-dimensional check

`u_t=d u_xx + A -(B+1)u + u²v`

`v_t=d v_xx + Bu - u²v`

with homogeneous Neumann boundaries, `A=1`, equal diffusion, and equilibrium
`(u,v)=(1,B)`. The homogeneous Hopf threshold is `B=2`, frequency 1; the
spatial modes move left when the diffusion is equal. A provider implementation
can be added without changing the continuation library. No unverified numerical
result is included in this repository.
