# test programs

This folder contains test programs for the routines implemented in this repository. 
If you want to test them all type in the shell:

```
./init_test.sh
cd build/
make -j8
```

## Conjugate Gradient

Conjugate Gradient (CG) methods implemented in the library: the standard CG and the BiCGStab. For that look at the files:

* ```CG.cpp```
* ```BiCGStab.cpp```

The library, ```../CG_solver.hpp```, is of general purpose. The user can define matrix and vector containers with the desired parallelizations and optimizations. There it is defined a base class for the CG and BiCGStab linear solvers, defined in ```../CG.hpp``` and ```BiCGStab.hpp```.

The file ```LA.hpp``` defines matrix and vector operations.

## yaml

Examples and tests of .yaml input files parsing.
