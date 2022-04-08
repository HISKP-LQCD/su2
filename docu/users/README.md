# Users guide

This directory contains the users guide to this project.

## Theoretical background

* At the moment fermions are included as a degenerate doublet using the staggered discretization: see [./fermions/main.pdf](./fermions/main.pdf)


## Library usage

The library contains code for both $SU(2)$ and $U(1)$ in $d=2,3,4$ dimensions.

### $U(1)$ theory

#### Hybrid Monte Carlo

- Program: ```hmc-u1.cc```
- Parameters parsing: ```.yaml``` input file
- Input file example: see [./hmc-u1.input.example](hmc-u1.input.example)


#### Offline measurements

- Program: ```measure-u1.cc```
- Parameters parsing: ```.yaml``` input file
- Input file example: see [./measure-u1.input.example](measure-u1.input.example)
- Measurable observables: 
  * Planar Wilson Loop
  * Gradient Flow 
 
---