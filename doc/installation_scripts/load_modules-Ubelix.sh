#!/bin/bash

# Modules to be loaded on the Ubelix machine: https://hpc-unibe-ch.github.io/

module --force purge


module load GLib/2.77.1-GCCcore-12.3.0 
module load CMake/3.26.3-GCCcore-12.3.0


module load Boost/1.82.0-GCC-12.3.0 
module load Eigen/3.4.0-GCCcore-12.3.0 

module load yaml-cpp/0.7.0-GCCcore-12.3.0

module load xtensor/0.24.7-foss-2023a

