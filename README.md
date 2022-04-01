# A MCMC sampling package for lattice gauge theories

This package implements Metropolis Monte Carlo and Hybrid Monte Carlo
for SU(2) and U(1) gauge theories in `d=2,3,4` dimensions. 

To compile, open the file

```
./init-cmake-build.sh
```

and set the value of the variable YAML_SRC_PATH to the path of your yaml-cpp source directory. Then run the above file with ```bash```. 

in the source directory. Then type make in `build/debug/` or
`build/release/` directories.

There are several executables being build. Invoking them with the `-h`
option will give instructions on command line parameters. 

There are also several observables implemented: the topological charge
in `d=2` and in `d=4`, Wilson loops and the gradient flow.
