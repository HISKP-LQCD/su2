# Users guide

This directory contains the users guide to this project.
The latter is written in `Rmd` format.
In order to produce the `html` and `pdf` output run `make all` from the present directory.

## Theoretical background

* At the moment fermions are included as a degenerate doublet using the staggered discretization: see [./staggered_fermions.Rmd](./staggered_fermions.Rmd)


## Library usage

The library contains code for both $SU(2)$ and $U(1)$ in $d=2,3,4$ dimensions.

- Programs: `main-u1.cpp` and `main-su2.cpp`
- Parameters parsing: `yaml` input file
- Input file example: see [./main.input](main.input)
  Notes:
  - Each run produces/updates a file named `nconf_counter.txt` in the directory storing the gauge configurations. This file contains a header line and the values below. At the moment an example:
  ```
  heat i path_conf
  1 1140 ./confs/config_u1.4.4.1.16.b1.900000.1140
  ```
  - `restart` and `heat` cannot be passed simultaneously
  - If `restart: true` is passed, the last configuration id and path are read from the `nconf_counter.txt` file.
  - The program can do online and offline measurements

## tests

In the `test/` folder the user can find various simple programs implementing the routines of this repository. For further details please look at the `test/README.md` file.
