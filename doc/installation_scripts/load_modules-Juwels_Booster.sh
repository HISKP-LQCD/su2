#!/bin/bash

# Modules to be loaded on the Juwels booster machine: https://apps.fz-juelich.de/jsc/hps/juwels/booster-overview.html
# NOTE: there are 2 options, see details below. Please select only the one you actually use

# Stages 2025 option (requires manual installation of xtensor)

module --force purge
module load Stages/2025
module load GCCcore/.13.3.0
module load Boost/1.86.0
module load CMake/3.30.3
module load Eigen/3.4.0

# Stages 2023 (no manual installation required, but Stage 2023 is deprecated)

module --force purge
module load Stages/2023
module load CMake/3.23.1
module load GCCcore/.11.3.0
module load xtensor/0.24.6
module load Boost/1.79.0
module load Eigen/3.4.0

