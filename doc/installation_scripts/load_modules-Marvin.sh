#!/bin/bash

# Modules to be loaded on the Marvin machine: https://www.hpc.uni-bonn.de/en/systems/marvin

module --force purge
module load GCCcore/13.3.0
module load Boost/1.85.0-GCC-13.3.0
module load CMake/3.29.3-GCCcore-13.3.0
module load Eigen/3.4.0-GCCcore-13.3.0

