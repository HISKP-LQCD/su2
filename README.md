# A MCMC sampling package for lattice gauge theories

This package implements Metropolis Monte Carlo and Hybrid Monte Carlo
for SU(2) and U(1) gauge theories in `d=2,3,4` dimensions. 

To compile, type

```
./init-cmake-build.sh
```

in the source directory. Then type make in `build/debug/` or
`build/release/` directories.

There are several executables being build. Invoking them with the `-h`
option will give instructions on command line parameters. 

There are also several observables implemented: the topological charge
in `d=2` and in `d=4`, Wilson loops and the gradient flow.
