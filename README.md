# A MCMC sampling package for lattice gauge theories

This package implements Metropolis Monte Carlo and Hybrid Monte Carlo
for SU(2) and U(1) gauge theories in `d=2,3,4` dimensions. 
The library supports also the measurements of several observables, e.g.: 

* topological charge
* Wilson loops
* gradient flow
* glueball correlators

## Software dependencies

* [Boost](https://www.boost.org/)
* [yaml-cpp](https://github.com/jbeder/yaml-cpp)
* [xtensor](https://github.com/xtensor-stack/xtensor)

## Installation instrucitons

This project is compiled using [cmake](https://cmake.org/). The installation process produces 2 directories, corresponding to the **debug** and **release** cmake build types.

In the following we assume that the source code is placed in `${HOME}/code/`.

1. Create a build directory in a custom location, e.g.: `${HOME}/build/su2/`
2. Copy the content of the directory `doc/installation_scripts/` in that directory.
3. Remove the `.example` extensions and change the fields of `do_cmake.sh` such that they match your choice.
4. Run 
   ``` bash 
   bash do_cmake.sh
   bash do_make.sh
   ```

Done. You'll find the **debug** and executables in:

-  **debug** : `${HOME}/build/su2/debug`
-  **release** : `${HOME}/build/su2/release`

