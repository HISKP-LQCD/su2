# Conjugate Gradient

In this folder is implemented an example of how to use the Conjugate Gradient (CG) methods implemented in the library: the standard CG and the BiCGStab. For that look at the files:

* ```CG.cpp```
* ```BiCGStab.cpp```

The library, ```../CG_solver.hpp```, is of general purpose. The user can define matrix and vector containers with the desired parallelizations and optimizations. There it is defined a base class for the CG and BiCGStab linear solvers, defined in ```../CG.hpp``` and ```BiCGStab.hpp```.

The file ```LA.hpp``` defines matrix and vector operations.
