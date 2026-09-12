# Hopf normal-form smoke example

This uses the shared provider, centered finite-difference map derivative, dense
reference spectrum path, checkpoint writer, and correction workflow. It is not
a cylinder result and its map spectrum is not a Floquet spectrum.

```bash
../../scripts/build_cpu.sh
../../build-cpu/nekrs-bif --example hopf_normal_form --task seed --output run
../../build-cpu/nekrs-bif --example hopf_normal_form --task stability \
  --restart run/checkpoint.json --output run-stability
```

The continuous system is `x'=a*x-y-(x²+y²)x`, `y'=x+a*y-(x²+y²)y`,
with `a=0` at onset and angular frequency 1. The example is deliberately
small; PDE states must use the Trilinos adapter rather than the dense path.
